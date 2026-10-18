# README editorial review — 2026/10/18

The README was revised for someone encountering BASIC and immune-receptor
assembly for the first time. The aim was a clear explanation, a working first
run, and an accessible route to benchmarking.

## Independent reviewers

| Reviewer | Scope and outcome |
| --- | --- |
| Fresh Codex agent | Started without the earlier conversation; read the README, usage/algorithm/validation guides, CLI, pipeline, examples, and packaging. Supplied a proposed rewrite and checked the integrated version for factual, command, and link errors. |
| Claude Opus 5.5 High (`claude-opus-5-5-high`) via agy | Completed a read-only editorial critique in 46.28 seconds; successful result, no timeout or stderr. |
| Gemini 4 Argon High (`gemini-4-argon-high`) via agy | Completed a read-only editorial critique in 48.92 seconds; successful result, no timeout or stderr. |

The two external models received the original README, benchmark guide, example
paths, and factual context about per-cell inputs, dependencies, candidate
assemblies, the forthcoming patch, and deferred registry publication. They
were instructed not to run tools or edit files. These were editorial reviews,
not independent executions, full source reviews, or biological validation.

## Changes incorporated

- Explain reads, seeds, and contigs before listing release features. Describe
  Bowtie2 seed selection and extension through overlapping reads in three steps.
- State the per-cell input requirement and barcode/UMI scope near the opening.
- Put one main installation path and a bundled human B-cell example first.
  Include clone/cd steps, expected chain lengths, output paths, and rerun behavior.
- Describe the result files and distinguish assembly status from execution
  failure. Explain `-a` beside the own-data examples.
- Add an immediately runnable comparison with expected JSON metrics, plus a
  timed example in the benchmark guide. Explain identifier and orientation
  matching and the distinction between regression checks and biological truth.
- Move the changelog, experimental options, assistant skill, and development
  links below the first-use path.

## Suggestions corrected before adoption

External suggestions were checked against the implementation and data:

- Seeds are selected **reads** matching reference regions, not simply copied
  reference sequences.
- The full `cell_A1` example uses `-n result`, which makes the generated FASTA
  identifiers match the supplied result. Blanket advice to rename every
  example's IDs was unnecessary for this route.
- Proposed 50,000/100,000-read quick trials were not adopted as the default
  walkthrough: reduced coverage can produce incomplete or unsupported results.
- The timed example compares the output of the command it actually runs. It
  does not score a previously renamed FASTA that could have become stale.
- The run manifest marks output publication; its existence does not mean
  that the biological assembly status is `complete`.
- The bundled stored assembly is labeled a regression reference. Broader
  claims about biological accuracy or the provenance of every supplied
  sequence were not inferred from reviewers' wording.

## Verification

A fresh installation of the local 2.0.1 wheel outside the source checkout ran
the documented full single-end example. Its FASTA matched the supplied file
byte for byte: heavy chain 878 nt, light chain 829 nt. The benchmark reported
two exact matches, exact recall 1.0, and exact precision 1.0. These values are
specific to this regression example, not general accuracy estimates.

Commands and local links were checked against the current implementation.
The fresh Codex agent's final pass found no blocking issue; its clarification
of `partial` versus `unsupported` was incorporated. See the
[validation record](../benchmarks/results/validation.json) for software checks
and the [separate skill review](SKILL_REVIEW.md) for assistant instructions.
