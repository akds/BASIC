# Changelog

## 2.0.1 — 2026/10/18

- Preserve the original alignment/timeout error when process-group termination is denied.
- Fall back to terminating the direct child and bound failure-path process cleanup.
- Share cleanup across alignment, annotation, reference building, and benchmark subprocesses.
- Add deterministic permission-denial and process-exit regression tests.
- Include Bioconda's required downstream compatibility pin in generated recipes.
- Rewrite the introduction and first-run guide, with a runnable example and benchmark walkthrough reviewed by fresh Codex, Claude, and Gemini reviewers.

## 2.0 — 2026/10/18

This release modernizes BASIC 1.5.0. The indexed default assembler has
regression coverage; new biological scoring modes remain experimental.

### Reliability

- Checked, shell-free Bowtie2 execution with retained diagnostics.
- Streaming FASTQ validation, synchronized paired files, clear invalid-input errors.
- Unique temporary directories and per-sample locks; no uncontrolled queue waits.
- Atomic per-file publication with rollback on normal multi-file write failures.
- Existing results require `--overwrite`; failed reruns preserve previous results.
- Explicit complete, partial, unsupported, and failed outcomes; enforce `-a` after assembly.
- Validate all index shards, arguments, trimmed seeds, and optional annotation setup.

### Algorithms

- Shared read pool and sorted-prefix overlap index replace repeated full-pool scans.
- Stable tie handling with ambiguity diagnostics; palindromes are counted once.
- Bounded extension and correct handling of mixed read lengths and short seeds.
- Optional quality/mate evidence, bounded beam search, and multiple candidate seeds.
- Separate incompatible seed loci and label TCR contigs from actual reference loci.
- Experimental related-species and local-homology seeding for non-model organisms.

### Features

- FASTA, AIRR TSV, chain summary, run manifest, checksums, and alignment logs.
- Optional local IgBLAST annotation, preserving unknown biological fields as null.
- TSV batches, per-sample error handling, and dependency/path preflight.
- Reference import/download and checked index construction with source provenance.
- Reproducible overlap benchmark and synthetic/reference-based evaluation tools.
- Python package, installed commands, wheel/source distributions, and GitHub Actions.
- PyPI/TestPyPI publishing workflow and a checksum-verified Bioconda recipe generator.
- Portable Codex, Claude, and Antigravity skill with an explicit `basic-skill` installer.
- Standard MIT license, confirmed for the 2.0 release.

### Migration

- Python 3.10 or newer is required; Python 2 support is removed.
- CLI flags from 1.5 remain, but undocumented internal Python helpers are not an API.
- Paired `--subsample N` now means N synchronized pairs; formerly it limited the
  concatenated unpaired stream and could exclude mate 2 entirely.
- Low-quality trimming is optional; terminal Ns are trimmed before seed selection.
- Full FASTQ format and paired identifiers are validated. Arbitrary file formats
  that Bowtie2 tolerated may now fail validation.
- Results never silently append. Existing output bundles require `--overwrite`.
- TCR names follow observed loci; this corrects the historical mouse database ordering.
- `complete` means a contig connects both seeds. Productivity and correct V(D)J
  biology require annotation and validation; they are not implied by that status.
