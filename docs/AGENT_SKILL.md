# Use BASIC 2.0 with an agent

BASIC includes a portable `basic-assemble` skill for running and interpreting
per-cell BCR/TCR assembly. Install it into an analysis project with `basic-skill`.
The skill contains instructions and examples; BASIC, Bowtie2, reference data,
and optional IgBLAST must be available separately.

## Install from your BASIC environment

If BASIC 2.0 is already installed, activate that environment and run:

```bash
basic-assemble --version
basic-skill --agent all --project '/path/to/analysis-project'
```

The project directory must already exist. Omit `--project` to use the current
directory. To install for a single agent, select it instead of `all`:

| Option | Destination relative to the project | Official discovery documentation |
| --- | --- | --- |
| `--agent codex` | `.agents/skills/basic-assemble/` | [OpenAI: local skills](https://learn.chatgpt.com/docs/build-skills#where-codex-loads-local-skills) |
| `--agent agy` | `.agents/skills/basic-assemble/` | [Antigravity: agent skills](https://antigravity.google/docs/skills) |
| `--agent claude` | `.claude/skills/basic-assemble/` | [Claude Code: skill locations](https://code.claude.com/docs/en/skills#choose-where-skills-load) |
| `--agent all` | Both directories above | Codex and agy share one copy. |

These are project installations. The command does not install a global skill,
change agent settings, or install an agent application.

To install a verified local BASIC 2.0 wheel in an isolated environment on a
POSIX system, use Python 3.10 or newer and replace the example paths below:

```bash
python3 -m venv '/path/to/basic-venv'
'/path/to/basic-venv/bin/python' -m pip install \
  '/path/to/basic_immune_assembler-2.0.1-py3-none-any.whl'
'/path/to/basic-venv/bin/basic-assemble' --version
'/path/to/basic-venv/bin/basic-skill' --agent all \
  --project '/path/to/analysis-project'
```

This example uses a wheel you already have. Installing it
supplies the `basic-skill` command. The skill installer copies only the bundled
instructions.

## Inspect and use the skill

Preview the main instructions without writing files:

```bash
basic-skill --print
```

`--print` displays only `SKILL.md`. A normal installation includes both:

```text
basic-assemble/
├── SKILL.md
└── references/
    └── custom-and-annotation.md
```

Keep the folder together: the main instructions link to the optional guide
relative to the installed skill folder. The guide covers custom references,
related-species seeding, and IgBLAST annotation.

Open the analysis project in your agent and ask, for example:

> Use the basic-assemble skill to assemble this human BCR cell from
> reads/cell_001.fastq.gz. Use the BASIC environment at /path/to/basic-venv,
> write to results/run_001, and check the inputs and dependencies first.

The official provider documentation above describes discovery and activation.
Check that `basic-assemble` appears in the agent's available skills before
relying on automatic selection. Make the BASIC environment accessible to the
agent, or provide its executable path. The skill does not change the agent's
permissions or grant access to files outside its environment.

## Source-checkout alternative

From a verified BASIC 2.0 source checkout, use its Python environment without
installing a wheel:

```bash
python -m basic_assembler.skill_cli --agent codex \
  --project '/path/to/analysis-project'
```

You can also copy `skills/basic-assemble/` directly into the desired discovery
directory. Copy the entire folder, not just `SKILL.md`. For example, this manual
copy refuses an existing destination instead of merging or replacing it:

```bash
python3 - <<'PY'
from pathlib import Path
import shutil

source = Path('/path/to/BASIC/skills/basic-assemble')
destination = Path('/path/to/analysis-project/.agents/skills/basic-assemble')
shutil.copytree(source, destination)
PY
```

Use `.claude/skills/basic-assemble` for Claude Code. The installer is preferred
when available because it also checks collisions and destination containment.

## Reinstall, customize, or update

Re-running `basic-skill` with identical installed files succeeds without
rewriting them. If any target file differs, or a file destination is occupied
by a directory, the command fails before writing any of the selected files.
Thus `--agent all` does not partially update the other agent when one copy has
local edits. Symbolic-link file destinations and paths escaping the project
are rejected. Unrelated files are not removed.

There is no overwrite flag. Before updating a customized skill, save your
changes outside its discovery directory, compare them with the new bundled
version, and remove the conflicting files or move the old skill
folder aside. Re-run the installer, then reapply relevant edits. If both
discovery directories are used, keep both copies up to date.

On an ordinary write failure, the installer attempts to remove files it
created during that attempt; existing files are preserved. This is not a
crash-atomic update mechanism, and newly created empty directories can remain.

See [the usage guide](USAGE.md) for command details and interpretation.
