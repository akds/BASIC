"""Export the bundled, provider-neutral BASIC skill into a chosen project."""
from __future__ import annotations

import argparse
from importlib.resources import files
import os
from pathlib import Path
import sys


def skill_root():
    checkout = Path(__file__).resolve().parent.parent / 'skills' / 'basic-assemble'
    return checkout if (checkout / 'SKILL.md').is_file() else files('basic_assembler.agent_skill')


def skill_files() -> dict[str, str]:
    root = skill_root()
    content = {'SKILL.md': root.joinpath('SKILL.md').read_text(encoding='utf-8')}
    references = root.joinpath('references')
    if references.is_dir():
        for item in sorted(references.iterdir(), key=lambda item: item.name):
            if item.is_file() and item.name.endswith('.md'):
                content[f'references/{item.name}'] = item.read_text(encoding='utf-8')
    return content


def skill_text() -> str:
    return skill_root().joinpath('SKILL.md').read_text(encoding='utf-8')


def install_skill(project: Path, agent: str, content: dict[str, str]) -> list[Path]:
    project = project.expanduser().resolve()
    if not project.is_dir():
        raise ValueError(f'Project directory does not exist: {project}')
    if agent not in ('codex', 'claude', 'agy', 'all'):
        raise ValueError(f'Unknown agent: {agent}')
    roots = ['.agents', '.claude'] if agent == 'all' else ['.claude' if agent == 'claude' else '.agents']
    for name in content:
        if Path(name).is_absolute() or '..' in Path(name).parts:
            raise ValueError(f'Unsafe skill resource name: {name}')
    targets = {project / root / 'skills' / 'basic-assemble' / name: text.encode('utf-8')
               for root in roots for name, text in content.items()}
    pending = []
    for target, payload in targets.items():
        if target.is_symlink() or not target.resolve().is_relative_to(project):
            raise ValueError(f'Skill destination must remain inside the project: {target}')
        if target.exists():
            if not target.is_file() or target.read_bytes() != payload:
                raise FileExistsError(f'Skill differs at {target}; back up and remove that file before reinstalling.')
        else:
            pending.append((target, payload))
    created = []
    try:
        for target, payload in pending:
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.resolve().is_relative_to(project):
                raise ValueError(f'Skill destination leaves the project: {target}')
            fd = os.open(target, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
            created.append(target)
            with os.fdopen(fd, 'wb') as handle:
                handle.write(payload)
    except BaseException:
        for target in reversed(created):
            target.unlink(missing_ok=True)
        raise
    return list(targets)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description='Install the bundled BASIC skill for Codex, Claude Code, or agy.')
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--agent', choices=['codex', 'claude', 'agy', 'all'])
    mode.add_argument('--print', dest='print_only', action='store_true', help='Print SKILL.md without creating files')
    parser.add_argument('--project', type=Path, default=Path('.'), help='Existing project directory (default: current directory)')
    args = parser.parse_args(argv)
    try:
        if args.print_only:
            print(skill_text(), end='')
        else:
            for path in install_skill(args.project, args.agent, skill_files()):
                print(path)
        return 0
    except (OSError, ValueError) as exc:
        print(f'basic-skill: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
