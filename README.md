# BASIC 2.0

**Assemble B-cell and T-cell receptor sequences from single-cell RNA-seq.**

Version 2.0 (2026/10/18) · patch **2.0.1**

BASIC reconstructs candidate immune-receptor transcripts from the short RNA
sequencing reads of individual cells. Give it reads from one cell, choose a
species and receptor type, and get assembled sequences with supporting reports.
It supports B-cell receptor (**BCR**) heavy and light chains and T-cell receptor
(**TCR**) chains. Human and mouse reference sequences are included.

Reads must already be separated by cell. BASIC does not separate cell barcodes
or deduplicate UMIs.

## The idea behind BASIC

Sequencing reads cover short pieces of a transcript. **Assembly** joins reads
that overlap into a longer sequence, called a **contig**. BASIC focuses this
search on immune receptors:

1. **Find starting reads.** Bowtie2 identifies reads matching known variable
   and constant regions of a receptor. These reads become assembly *seeds*.
2. **Extend through overlaps.** BASIC builds outward from the seeds using exact
   sequence overlaps with other reads from the same cell.
3. **Report candidate chains.** It writes assembled sequences and records
   whether the variable and constant seeds connect, along with read support
   and ambiguity information.

A **complete** assembly connects its seeds. It does not by itself establish a
full coding sequence, a productive receptor, or a correct CDR3. Optional
IgBLAST annotation can add gene calls, junction/CDR3 information, and predicted
productivity.

## Install and try an example

BASIC needs **Python 3.10+** and **Bowtie2**. On Linux or macOS, start from a
checkout with the included examples. With Git and Conda installed:

```bash
git clone https://github.com/akds/BASIC.git
cd BASIC
conda env create -f environment.yml
conda activate basic-2
python -m pip install .
basic-assemble --version
```

Already have Python and Bowtie2? Run `python -m pip install .` in your Python
environment from the checkout. The [installation guide](docs/INSTALLATION.md)
also covers [release wheels and checksums](https://github.com/akds/BASIC/releases/tag/v2.0.1),
pipx, uv, and platform details. PyPI publication is deferred; the Bioconda
`basic` update requires separate review and publication.

From the repository root, run the included human B-cell example:

```bash
basic-assemble -SE examples/cell_A1/ERR1421621_1M_reads.fastq.gz \
  -g human -i BCR -n result -o results/cell_A1 -p 2
```

This example has one million single-end reads. Expect two contigs: a heavy
chain of **878 nucleotides** and a light chain of **829 nucleotides**, saved in
`results/cell_A1/result.fasta`. The reports are beside it. The comparison below
checks both sequences against the supplied result.

Add `--check` to check paths, indexes, and dependencies before assembly; it
does not read every FASTQ record. To rerun, choose a new output directory or
add `--overwrite` to replace existing results explicitly.

## Use your own reads

Use reads with coverage across the receptor transcript. BASIC accepts
single-end or paired-end FASTQ, plain or gzipped,
with Phred+33 quality scores. Paired files must have matching identifiers and
counts.

Replace the filenames below with reads from one cell:

```bash
# Human BCR, paired-end reads
basic-assemble -PE_1 cell_R1.fastq.gz -PE_2 cell_R2.fastq.gz \
  -g human -i BCR -n cell_001 -o results

# Human TCR, single-end reads
basic-assemble -SE t_cell.fastq.gz \
  -g human -i TCR -n t_cell_001 -o results
```

Use `-g mouse` for mouse reads. TCR outputs carry locus labels such as `TRA`
and `TRB`; matching mouse reads and references also support gamma/delta
(`TRG`/`TRD`) candidates. Recovery depends on the reads and available seeds.

Incomplete assembly exits with an error by default. Add **`-a`** to retain
partial or empty biological results; invalid input and execution errors still
fail. Use `-p 4` for four Bowtie2 threads or `-b /path/to/bowtie2` to select the
executable. Run `basic-assemble --help` for all options.

## Read the results

For the included example, look in `results/cell_A1/`:

| File | What it tells you |
| --- | --- |
| `result.fasta` | Assembled nucleotide sequences |
| `result.summary.tsv` | Status and ambiguity by chain or reference group |
| `result.airr.tsv` | Sequences and assembly metadata in an AIRR-compatible table |
| `result.run.json` | Parameters, versions, diagnostics, and file checksums |
| `result.logs/` | Bowtie2 diagnostics |

A chain is **complete** when its seeds connect, **partial** when a seed yields
a contig but the seeds do not connect, or **unsupported** when neither seed is
usable. Execution failures are recorded separately when possible. Gene calls,
CDR3, and productivity remain empty
unless annotation supplies them. Read support is a count of observations, not
a molecule count or a confidence score.

For automated workflows, use the run manifest as the completion marker and
verify its output checksums. See [output and algorithm details](docs/ALGORITHMS.md).

## Check an example, then benchmark

Compare the example you just assembled with the supplied BASIC result:

```bash
python benchmarks/benchmark.py \
  --expected-fasta examples/cell_A1/result.fasta \
  --actual-fasta results/cell_A1/result.fasta \
  --label included-human-bcr-regression \
  --output results/cell_A1/comparison.json
```

Expect `exact_match_count: 2`, `exact_recall: 1.0`, and `exact_precision: 1.0`
inside the JSON report's `sequences` object. Keep `-n result` in the example
assembly command so the FASTA identifiers match. This checks reproducibility
against an existing assembly; biological accuracy needs independent truth.

For your own benchmark, replace `--expected-fasta` with independently validated
sequences and `--actual-fasta` with your assemblies. Match FASTA identifiers
and sequence orientation first: the script does not align or automatically
pair contigs. Optional AIRR tables add gene-call, junction/CDR3, and productivity
comparisons. The [benchmark guide](benchmarks/README.md) includes a timed example,
metric definitions, and advice for fair comparisons. See the
[validation report](docs/VALIDATION.md) for measured results and their limits.

## Go further

| Task | Where to start |
| --- | --- |
| Process many cells | [Batch TSV format and commands](docs/USAGE.md#multiple-cells) |
| Add gene and CDR3 annotations | [IgBLAST setup and usage](docs/USAGE.md#optional-igblast-annotation) |
| Work with another organism | [Build references from FASTA, HTTPS sources, or NCBI accessions](docs/NON_MODEL_ORGANISMS.md) |
| Use an AI coding assistant | [Portable skill for Codex, Claude, and Antigravity](docs/AGENT_SKILL.md) |
| Explore assembly settings | [Algorithms and experimental options](docs/ALGORITHMS.md) |

Custom references must identify the correct receptor loci and orientation.
Related-species seeding is experimental; it does not infer a validated germline
catalogue or justify definitive gene calls across species.

## What's new in 2.0

- **BCR and TCR assembly:** indexed overlap searches, shared read pools,
  deterministic ties, and bounded searches; the default retains BASIC's
  original entropy-based extension rule.
- **TCR reporting:** locus-aware chain labels, corrected historical mouse
  ordering, and support for matching mouse gamma/delta seeds.
- **Everyday workflows:** stricter input checks, protected outputs, structured
  reports, per-cell batches, optional annotation, and custom reference tools.
- **Exploratory options:** quality and mate evidence, alternative assembly
  paths, and multiple candidates. These remain experimental; no general
  biological accuracy improvement is established.

Patch 2.0.1 fixes error handling when the operating system denies process-group
cleanup. See the [changelog and migration notes](CHANGELOG.md),
[review record](docs/REVIEW.md), and [README editorial reviews](docs/README_REVIEW.md).
For development, see [testing and release instructions](docs/RELEASING.md).

## Citation and license

Canzar S, Neu KE, Tang Q, Wilson PC, Khan AA. **BASIC: BCR assembly from single
cells.** *Bioinformatics* 33(3), 1 February 2017.

BASIC uses the [MIT license](LICENSE). Bundled references are preserved from
BASIC 1.5.0; their original download dates are unknown. See
[reference provenance](db/provenance.json) and the [historical license notice](LICENSE-NOTICE.md).
Reference data and external tools retain their respective terms.
