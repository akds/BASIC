# BASIC 2.0 portable skill review

Review date: 2026-10-18.

The review covered [`SKILL.md`](../skills/basic-assemble/SKILL.md) and its
[custom references and annotation guide](../skills/basic-assemble/references/custom-and-annotation.md).
The aim was a self-contained skill folder that can be copied separately from
the repository and used with an independently installed BASIC 2.0 environment.

## External model reviews

Both requested reviewers received the main skill and optional guide through
the `agy` CLI. Each prompt requested a short, read-only critique of concrete
portability, command-use, installation, and scientific-interpretation issues.
The invocations used `--disable-slash-commands`, `--output-format stream-json`,
and a 120-second print timeout. Neither used planning mode.

| Reviewer model ID | Outcome |
| --- | --- |
| `claude-opus-5-5-high` | Substantive response on the first attempt, approximately 24 seconds. |
| `gemini-4-argon-high` | Substantive response on the first attempt, approximately 51 seconds. |

Both streams ended with a successful result containing review text. Neither
required a retry or showed tool execution. These were reviews of the supplied
instructions, not independent executions of BASIC or biological validation.

## Changes made after review

- **Installation:** explain how to locate an existing environment or obtain a
  verified BASIC 2.x checkout/wheel, use an isolated Python environment, install
  Bowtie2 separately, and recheck the installed commands. Do not guess a package
  name or assume the historical repository's default branch contains 2.x.
- **Reference paths:** state the exact `ROOT/<genome>/<receptor>/` resolution for
  `-d ROOT`, including the `--seed-species` override. Distinguish the skill's
  optional guide from the biological reference library.
- **Batch configuration:** put the shared annotation-configuration warning next
  to the mixed-species example, and explain explicit biological context and
  per-row sample names.
- **Stale locks:** require verification of the originating run, host/container,
  or scheduler job. A local PID check does not establish that a shared-filesystem
  lock is stale.
- **Ambiguous biology:** state that multiple assembled candidates do not resolve
  dual-chain expression, doublets, ambient RNA, or chimeras. Preserve the existing
  distinctions between assembly completion, productivity, and raw read support.

Codex checked suggestions against the implementation before adopting them.
Gemini's proposed blanket stop for every missing executable was not adopted:
IgBLAST and reference building are optional workflows. Its suggested prohibition
on `--batch` with `-n` was corrected to describe the actual behavior: manifest
`sample_id` values supply batch sample names. Claude's optional suggestion to
move all experimental settings out of the main skill was not adopted; the short
section keeps their experimental status visible.

## Validation and limits

The skill's YAML structure and local links were checked, including links after
copying the skill folder away from the repository. The example reference
manifest was parsed as JSON. Commands and caveats were checked against the
current BASIC 2.0 command-line implementation.

The skill does not bundle BASIC, Bowtie2, IgBLAST, or reference databases. These
reviews do not establish annotation accuracy, clinical suitability, or improved
recovery from experimental scoring or related-species seeding. Those require
appropriate data and separate validation.
