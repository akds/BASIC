from pathlib import Path

import pytest

from basic_assembler import skill_cli


def test_all_installs_two_identical_skills_and_is_idempotent(tmp_path):
    targets = skill_cli.install_skill(tmp_path, 'all', {'SKILL.md': 'test skill\n'})
    assert len(targets) == 2
    assert all(path.read_text() == 'test skill\n' for path in targets)
    assert skill_cli.install_skill(tmp_path, 'all', {'SKILL.md': 'test skill\n'}) == targets


@pytest.mark.parametrize('agent,folder', [('codex', '.agents'), ('agy', '.agents'), ('claude', '.claude')])
def test_agent_selects_discovery_directory(tmp_path, agent, folder):
    assert skill_cli.install_skill(tmp_path, agent, {'SKILL.md': 'skill'}) == [tmp_path / folder / 'skills/basic-assemble/SKILL.md']


def test_existing_different_skill_prevents_all_writes(tmp_path):
    existing = tmp_path / '.claude/skills/basic-assemble/SKILL.md'
    existing.parent.mkdir(parents=True)
    existing.write_text('user customized')
    with pytest.raises(FileExistsError):
        skill_cli.install_skill(tmp_path, 'all', {'SKILL.md': 'new'})
    assert existing.read_text() == 'user customized'
    assert not (tmp_path / '.agents').exists()


def test_symlink_cannot_redirect_install_outside_project(tmp_path):
    project, outside = tmp_path / 'project', tmp_path / 'outside'
    project.mkdir()
    outside.mkdir()
    (project / '.agents').symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError, match='inside'):
        skill_cli.install_skill(project, 'codex', {'SKILL.md': 'skill'})
    assert list(outside.iterdir()) == []


def test_second_write_failure_removes_only_new_skill(tmp_path, monkeypatch):
    original = skill_cli.os.open
    def fail_second(path, *args):
        if '.claude' in str(path):
            raise OSError('write unavailable')
        return original(path, *args)
    monkeypatch.setattr(skill_cli.os, 'open', fail_second)
    with pytest.raises(OSError):
        skill_cli.install_skill(tmp_path, 'all', {'SKILL.md': 'skill'})
    assert not list(tmp_path.rglob('SKILL.md'))


def test_print_and_missing_project(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(skill_cli, 'skill_text', lambda: 'printable skill\n')
    assert skill_cli.main(['--print']) == 0
    assert capsys.readouterr().out == 'printable skill\n'
    assert skill_cli.main(['--agent', 'all', '--project', str(tmp_path / 'missing')]) == 1
    assert not (tmp_path / 'missing').exists()


def test_bundled_skill_contains_standard_frontmatter():
    text = skill_cli.skill_text()
    assert text.startswith('---\nname: basic-assemble\n')
    assert '\ndescription:' in text


def test_referenced_guides_are_installed_with_entrypoint(tmp_path):
    content = skill_cli.skill_files()
    assert 'references/custom-and-annotation.md' in content
    targets = skill_cli.install_skill(tmp_path, 'claude', content)
    assert len(targets) == len(content)
    for name, text in content.items():
        assert (tmp_path / '.claude/skills/basic-assemble' / name).read_text() == text


def test_resource_path_cannot_escape_skill(tmp_path):
    with pytest.raises(ValueError, match='Unsafe'):
        skill_cli.install_skill(tmp_path, 'codex', {'../../outside': 'wrong'})
    assert not (tmp_path / '.agents').exists()
