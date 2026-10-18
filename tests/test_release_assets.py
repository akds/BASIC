"""Release packaging preserves tested distributions and a portable skill."""
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import zipfile

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / 'packaging/release_assets.py'
spec = importlib.util.spec_from_file_location('basic_release_assets', SCRIPT)
assets = importlib.util.module_from_spec(spec)
spec.loader.exec_module(assets)


@pytest.fixture
def inputs(tmp_path):
    dist = tmp_path / 'built distributions'
    dist.mkdir()
    skill = tmp_path / 'source skill'
    (skill / 'references').mkdir(parents=True)
    (skill / 'SKILL.md').write_text('---\nname: basic-assemble\ndescription: Assemble receptors\n---\n')
    (skill / 'references/usage.md').write_text('Local workflow\n')
    _distributions(dist)
    return dist, skill


def _distributions(dist, wheel_version='2.0.0', sdist_version='2.0.0', project='basic-immune-assembler'):
    with zipfile.ZipFile(dist / 'basic_immune_assembler-2.0.0-py3-none-any.whl', 'w') as archive:
        archive.writestr('basic_immune_assembler-2.0.0.dist-info/METADATA',
                         f'Name: {project}\nVersion: {wheel_version}\n')
    with tarfile.open(dist / 'basic_immune_assembler-2.0.0.tar.gz', 'w:gz') as archive:
        data = f'Name: {project}\nVersion: {sdist_version}\n'.encode()
        entry = tarfile.TarInfo('basic_immune_assembler-2.0.0/PKG-INFO')
        entry.size = len(data)
        archive.addfile(entry, io.BytesIO(data))


def test_exact_distributions_portable_zip_and_checksums(inputs, tmp_path):
    dist, skill = inputs
    output = tmp_path / 'release assets'
    result = assets.prepare_assets(dist, skill, output)
    assert result['version'] == '2.0.0'
    assert set(path.name for path in output.iterdir()) == {
        'basic_immune_assembler-2.0.0-py3-none-any.whl',
        'basic_immune_assembler-2.0.0.tar.gz', 'basic-assemble-skill-2.0.0.zip', 'SHA256SUMS'}
    for source in dist.iterdir():
        assert (output / source.name).read_bytes() == source.read_bytes()
    with zipfile.ZipFile(output / 'basic-assemble-skill-2.0.0.zip') as archive:
        assert archive.namelist() == ['basic-assemble/SKILL.md', 'basic-assemble/references/usage.md']
        assert archive.read('basic-assemble/references/usage.md') == b'Local workflow\n'
        assert all(entry.date_time == (1980, 1, 1, 0, 0, 0) for entry in archive.infolist())
    lines = (output / 'SHA256SUMS').read_text().splitlines()
    assert len(lines) == 3
    for line in lines:
        digest, filename = line.split('  ')
        assert digest == hashlib.sha256((output / filename).read_bytes()).hexdigest()


def test_skill_zip_does_not_depend_on_source_mtime_or_permissions(inputs, tmp_path):
    dist, skill = inputs
    first, second = tmp_path / 'first', tmp_path / 'second'
    assets.prepare_assets(dist, skill, first)
    for source in skill.rglob('*.md'):
        os.utime(source, (1700000000, 1700000000))
        source.chmod(0o600)
    assets.prepare_assets(dist, skill, second)
    assert (first / 'basic-assemble-skill-2.0.0.zip').read_bytes() == (second / 'basic-assemble-skill-2.0.0.zip').read_bytes()


@pytest.mark.parametrize('pattern,extra', [('*.whl', False), ('*.tar.gz', False), ('*.whl', True), ('*.tar.gz', True)])
def test_wrong_distribution_count_does_not_create_output(inputs, tmp_path, pattern, extra):
    dist, skill = inputs
    source = next(dist.glob(pattern))
    if extra:
        (dist / ('extra-' + source.name)).write_bytes(source.read_bytes())
    else:
        source.unlink()
    output = tmp_path / 'assets'
    with pytest.raises(ValueError, match='exactly one'):
        assets.prepare_assets(dist, skill, output)
    assert not output.exists()


@pytest.mark.parametrize('kwargs', [{'sdist_version': '2.0.1'}, {'project': 'another-project'}, {'wheel_version': '../bad'}])
def test_mismatched_or_unsafe_metadata_rejected(inputs, tmp_path, kwargs):
    dist, skill = inputs
    _distributions(dist, **kwargs)
    output = tmp_path / 'assets'
    with pytest.raises(ValueError):
        assets.prepare_assets(dist, skill, output)
    assert not output.exists()


def test_existing_output_preserved_and_symlink_skill_rejected(inputs, tmp_path):
    dist, skill = inputs
    output = tmp_path / 'existing'
    output.mkdir()
    (output / 'keep').write_text('prior release')
    with pytest.raises(FileExistsError):
        assets.prepare_assets(dist, skill, output)
    assert (output / 'keep').read_text() == 'prior release'
    (skill / 'linked').symlink_to(output / 'keep')
    with pytest.raises(ValueError, match='symlink'):
        assets.prepare_assets(dist, skill, tmp_path / 'new')
    assert not (tmp_path / 'new').exists()


def test_cli_paths_with_spaces(inputs, tmp_path):
    dist, skill = inputs
    output = tmp_path / 'CLI assets'
    completed = subprocess.run([sys.executable, str(SCRIPT), '--dist', str(dist), '--skill', str(skill),
                                '--output', str(output)], capture_output=True, text=True, check=True)
    assert json.loads(completed.stdout)['output'] == str(output)


def test_copy_failure_removes_partial_assets_and_preserves_inputs(inputs, tmp_path, monkeypatch):
    dist, skill = inputs
    original_bytes = {path: path.read_bytes() for path in dist.iterdir()}
    output = tmp_path / 'failed assets'
    original_copy = assets.shutil.copyfile

    def fail_on_sdist(source, destination):
        if source.name.endswith('.tar.gz'):
            raise OSError('destination unavailable')
        return original_copy(source, destination)

    monkeypatch.setattr(assets.shutil, 'copyfile', fail_on_sdist)
    with pytest.raises(OSError, match='destination unavailable'):
        assets.prepare_assets(dist, skill, output)
    assert not output.exists()
    assert all(path.read_bytes() == data for path, data in original_bytes.items())
