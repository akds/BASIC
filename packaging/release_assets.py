"""Bundle exact built distributions and a standalone BASIC skill for a release."""
from __future__ import annotations

import argparse
from email.parser import BytesParser
import hashlib
import json
from pathlib import Path
import re
import shutil
import tarfile
import zipfile


def _metadata(path: Path) -> tuple[str, str]:
    if path.suffix == '.whl':
        with zipfile.ZipFile(path) as archive:
            entries = [name for name in archive.namelist()
                       if name.endswith('.dist-info/METADATA') and name.count('/') == 1]
            if len(entries) != 1:
                raise ValueError(f'Expected one wheel METADATA file: {path}')
            content = archive.read(entries[0])
    else:
        with tarfile.open(path, 'r:gz') as archive:
            entries = [member for member in archive.getmembers()
                       if member.name.endswith('/PKG-INFO')
                       and member.name.count('/') == 1 and member.isfile()]
            if len(entries) != 1:
                raise ValueError(f'Expected one top-level sdist PKG-INFO file: {path}')
            with archive.extractfile(entries[0]) as handle:
                content = handle.read()
    metadata = BytesParser().parsebytes(content)
    name, version = metadata.get('Name', ''), metadata.get('Version', '')
    if re.sub(r'[-_.]+', '-', name).lower() != 'basic-immune-assembler':
        raise ValueError(f'Unexpected distribution project in {path}: {name!r}')
    if not re.fullmatch(r'[0-9][A-Za-z0-9.!+_-]*', version):
        raise ValueError(f'Unsafe or missing distribution version in {path}: {version!r}')
    return name, version


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def prepare_assets(dist: Path, skill: Path, output: Path) -> dict:
    """Create a fresh asset directory; never rebuild or modify distributions."""
    dist, skill, output = (Path(path).expanduser().resolve() for path in (dist, skill, output))
    if output.exists():
        raise FileExistsError(f'Release asset directory already exists: {output}')
    if output.is_relative_to(skill):
        raise ValueError('Release assets must be outside the source skill directory')
    distributions = []
    for pattern in ('*.whl', '*.tar.gz'):
        matches = sorted(dist.glob(pattern))
        if len(matches) != 1 or not matches[0].is_file():
            raise ValueError(f'Expected exactly one {pattern} distribution in {dist}')
        distributions.append(matches[0])
    versions = {_metadata(path)[1] for path in distributions}
    if len(versions) != 1:
        raise ValueError('Wheel and source distribution versions do not match')
    version = versions.pop()
    if not (skill / 'SKILL.md').is_file():
        raise ValueError(f'Skill entrypoint not found: {skill / "SKILL.md"}')
    sources = []
    for source in sorted(skill.rglob('*')):
        if source.is_symlink():
            raise ValueError(f'Skill resource must not be a symlink: {source}')
        if source.is_file():
            sources.append(source)
    output.mkdir(parents=True)
    try:
        assets = []
        for source in distributions:
            destination = output / source.name
            shutil.copyfile(source, destination)
            assets.append(destination)
        skill_zip = output / f'basic-assemble-skill-{version}.zip'
        with zipfile.ZipFile(skill_zip, 'w', compression=zipfile.ZIP_STORED) as archive:
            for source in sources:
                name = 'basic-assemble/' + source.relative_to(skill).as_posix()
                entry = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                entry.create_system = 3
                entry.external_attr = 0o100644 << 16
                archive.writestr(entry, source.read_bytes())
        assets.append(skill_zip)
        checksums = {path.name: _sha256(path) for path in sorted(assets)}
        (output / 'SHA256SUMS').write_text(
            ''.join(f'{digest}  {name}\n' for name, digest in checksums.items()),
            encoding='utf-8')
    except BaseException:
        shutil.rmtree(output)
        raise
    return {'version': version, 'output': str(output), 'sha256': checksums}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dist', type=Path, default=Path('dist'))
    parser.add_argument('--skill', type=Path, default=Path('skills/basic-assemble'))
    parser.add_argument('--output', type=Path, default=Path('release-assets'))
    args = parser.parse_args(argv)
    try:
        result = prepare_assets(args.dist, args.skill, args.output)
    except (OSError, ValueError, tarfile.TarError, zipfile.BadZipFile) as error:
        parser.exit(2, f'release_assets: {error}\n')
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
