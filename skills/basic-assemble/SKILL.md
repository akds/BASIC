---
name: basic-assemble
description: Runs and interprets BASIC 2.0 BCR/TCR assembly from per-cell RNA-seq FASTQ, including paired reads, batches, reference selection, custom species, and optional IgBLAST annotation. Use for BASIC immune-receptor assembly, troubleshooting, or result review; not for the BASIC programming language or generic RNA-seq alignment.
---

# BASIC immune-receptor assembly

Turn the user's per-cell reads into reproducible receptor-contig results. Use
ordinary local command execution and file inspection; no particular agent,
plugin, or external model is required. This folder can be copied independently
of the BASIC repository. BASIC, Bowtie2, reference data, and optional IgBLAST
must be available separately.

## Establish the run

1. Identify the input files, biological species, BCR or TCR, sample identifiers,
   and output destination from the request and available metadata. Ask only for
   consequential missing information; do not silently assume human BCR.
2. Check `basic-assemble --version` and `basic-assemble --help`. These workflows
   target **BASIC 2.0**. In a confirmed source checkout, `python BASIC.py` accepts
   the same arguments. Use that checkout's Python environment. For other
   versions, inspect their help before using these flags. Do not assume a
   package named `basic` installs this version or claim software is installed
   merely because the skill is present.
   If the command is missing, first locate an existing environment or checkout.
   Otherwise use a verified BASIC 2.x wheel/checkout in a Python 3.10+ virtual
   environment, install Bowtie2 separately with the environment's package
   manager, and recheck version/help. Do not assume the historical repository's
   default branch contains 2.x; request the intended source if unavailable.
3. Use one cell per run. Inputs are plain/gzipped four-line FASTQ with Phred+33
   quality scores. Paired files must have matching record IDs, counts, lane
   order, and mate direction. BASIC does not demultiplex cell barcodes, deduplicate
   UMIs, or assemble directly from a 10x count matrix. Do not pool distinct cells.
4. Select references deliberately. Bundled `human` and `mouse` libraries cover
   BCR and TCR; supported contigs depend on actual seeds and reads. TCR labels
   follow the reference locus: alpha/beta, and mouse gamma/delta where supported.
   With `-d ROOT`, BASIC uses `ROOT/<genome>/<receptor>/`: `<genome>` is `-g`
   (or `--seed-species` when set), and `<receptor>` is `-i`.
   For new species or annotation, read
   [references/custom-and-annotation.md](references/custom-and-annotation.md).
   Resolve this guide relative to the installed skill folder.
5. Use a distinct output sample name for independent runs. Preserve inputs and
   prior results. Do not add `--overwrite`, delete locks, or replace references
   merely to make a command pass. Honor explicit replacement authorization;
   otherwise use a new destination or ask if the requested destination is fixed.

## Execute the chosen workflow

The paths below are illustrative; substitute verified input paths. Quote
paths passed to a shell. Record the exact command and exit status.

Single-end BCR:

```bash
basic-assemble -SE 'reads/cell_001.fastq.gz' -g human -i BCR \
  -n cell_001 -o 'results/run_001' -p 4
```

Paired-end TCR:

```bash
basic-assemble -PE_1 'reads/cell_002_R1.fastq.gz' \
  -PE_2 'reads/cell_002_R2.fastq.gz' -g mouse -i TCR \
  -n cell_002 -o 'results/run_001' -p 4
```

Use comma-separated input paths for multiple lanes of the **same cell**;
paired lane lists must correspond. `--subsample N` takes the first N physical
single reads or N synchronized pairs across all lanes combined, before
filtering. It is not random sampling; later lanes may contribute nothing.
Only consumed records undergo FASTQ validation. Omit it for a full-data run.

For batches, create a tab-separated manifest with unique `sample_id` values
(also unique ignoring case). Use `fastq` for single-end rows or `left` and
`right` for paired rows. Paths resolve relative to the manifest. Optional
`genome` and `receptor` columns override `-g` and `-i`. Supply the biological
context explicitly through those columns or command-line flags. Each row's
`sample_id` supplies its output name, overriding `-n`:

```tsv
sample_id	left	right	genome	receptor
cell_001	reads/cell_001_R1.fastq.gz	reads/cell_001_R2.fastq.gz	human	BCR
cell_002	reads/cell_002_R1.fastq.gz	reads/cell_002_R2.fastq.gz	mouse	TCR
```

This mixed-species/receptor example has **no IgBLAST annotation**. With
`--igblast`, every row must use the same species/receptor-matched databases and
auxiliary data; split batches when those configurations differ.

```bash
basic-assemble --batch cells.tsv -o 'results/run_001' -p 4
```

Cells run sequentially; `-p` controls each aligner's threads. Add
`--continue-on-error` when the user wants other cells processed after a failure;
the batch still exits nonzero if any cell fails. Do not combine `--batch` with
command-line read inputs.

Append `--check` to the intended command for a cheap path/index/dependency
preflight, then run without it. This does not inspect every FASTQ record or
establish biological suitability of references. Use `-b /path/to/bowtie2` if
the intended aligner is not on PATH. `-t /scratch` chooses the parent for a
unique temporary directory. Use `--keep-temp` for troubleshooting and consider
its disk cost. Set `--alignment-timeout SECONDS` only when an explicit limit is
needed; mapping has no default timeout.

## Interpret and recover

Inspect `<sample>.summary.tsv`, `<sample>.airr.tsv`, `<sample>.run.json`, and
`<sample>.fasta`; successful runs retain aligner logs in `<sample>.logs/`.
The run manifest records parameters, input/reference provenance, read-filter
statistics, and output checksums. It is published last; verify checksums when
consuming an existing bundle. A prior FASTA alone does not prove the latest run
succeeded. Keep the manifest with exported sequences.

| Observation | Interpretation and next action |
| --- | --- |
| `complete` | Assembly connects usable seeds. It does **not** establish a full coding sequence, productive receptor, correct CDR3, or absence of a chimera. |
| Partial or missing chain | Inspect seed coverage, filtering, references, and logs. Use `-a` only when partial/empty biological results are acceptable; it does not suppress execution errors. |
| Empty annotation fields | Unknown/unavailable. Do not convert an empty `productive` field to false or infer productivity from contig length. |
| Ambiguity or multiple candidates | Report alternatives and their support; do not silently pick biological truth from a score. Read support is not UMI count, molecule count, or calibrated confidence. |
| Nonzero exit or `<sample>.failed.json` | Diagnose the recorded failure; retain prior successful results and any recovery directory. Do not relabel the attempt successful because old outputs exist. |
| Existing lock | A local PID check is meaningful only on the originating host/container. Confirm that run or scheduler job has ended before treating the lock as stale; if its origin is uncertain, preserve it and use a fresh name. |

Multiple recovered chains alone cannot distinguish true dual-receptor
expression from doublets, ambient RNA, or chimeric assembly.

Report sample-level success/failure, chain/locus status, output locations,
software/reference versions, and material caveats. Distinguish an actual run
from a preflight or a proposed command. Do not claim improved biological
accuracy from a successful execution, synthetic test, or assembly score.

## Optional exploration

Keep default scoring/search settings unless the user requests a comparison or
an experiment. `--scoring evidence`, `--quality-weighted`, `--mate-aware`, wider
`--beam-width`, and multiple `--max-anchors`/`--max-candidates` are experimental.
Quality/mate weighting requires evidence scoring; mate awareness also requires
paired reads. Broader searches cost time/memory. Compare with a separately
named default run and independent truth where available. `--trim-quality 20`
is optional terminal trimming; inspect retained-read statistics rather than
assuming it improves recovery. Related-species/local seeding is also
experimental; follow the reference guide and report target and seed species.
