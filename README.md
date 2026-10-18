# BASIC 2.0

**BCR and TCR assembly from single-cell RNA-seq**

**Version 2.0 (2026/10/18)**

BASIC reconstructs candidate immune-receptor transcripts from per-cell RNA-seq
reads. It assembles B-cell receptor heavy and light chains and T-cell receptor
chains using variable- and constant-region seeds. Human and mouse references
are included; custom references support additional organisms.

## What's new in 2.0

| Update | What it provides |
| --- | --- |
| Indexed assembly | A shared read pool and indexed overlap search replace repeated full-pool scans. Ties are deterministic and searches have explicit bounds. |
| Locus-aware TCR reporting | Chain labels follow the mapped reference loci, correcting the historical mouse TCR ordering. Matching mouse references also support gamma/delta candidates. |
| Reliable execution | Validated FASTQ and paired-read inputs, checked Bowtie2 failures, isolated temporary files, and protection against accidental output overwrites. |
| Structured results | FASTA sequences, AIRR-compatible tables, chain summaries, logs, and manifests with parameters, versions, and checksums. |
| More workflows | Per-cell batch processing, optional IgBLAST annotation, and tools to find, import, and build custom references. |
| Experimental search options | Quality and mate evidence, alternative assembly paths, multiple candidates, and related-species seeding for exploratory work. |

The default retains BASIC's original entropy-based extension rule. The new
experimental modes do not establish an improvement in biological accuracy;
see the [validation results and limitations](docs/VALIDATION.md).

## Install

Requires **Python 3.10+** and **Bowtie2**. From this checkout:

```bash
conda env create -f environment.yml
conda activate basic-2
python -m pip install .
```

Alternatively, install Bowtie2 separately and run `python -m pip install .`.
The Python package has no third-party runtime dependencies; Bowtie2 is a
separate executable. Installed commands are `basic-assemble`, `BASIC.py`,
`basic-references`, and `basic-skill`. `python BASIC.py` also works directly
from the checkout.

You can also install the published 2.0 wheel into a Python environment with
Bowtie2 available:

```bash
python -m pip install https://github.com/akds/BASIC/releases/download/v2.0.0/basic_immune_assembler-2.0.0-py3-none-any.whl
```

[Release downloads and checksums](https://github.com/akds/BASIC/releases/tag/v2.0.0)
include the wheel, source, and standalone agent skill. PyPI publication as
`basic-immune-assembler` is deferred. The Bioconda update keeps the package name
`basic` and requires a separate review and build; check its available version
before installing. See the [installation guide](docs/INSTALLATION.md) for pip,
pipx, uv, Conda/Mamba, and platform details.

## Quick start

```bash
# Human BCR: paired-end reads from one cell
basic-assemble -PE_1 cell_R1.fastq.gz -PE_2 cell_R2.fastq.gz \
  -g human -i BCR -n cell_001 -o results

# Human TCR: single-end reads
basic-assemble -SE t_cell.fastq.gz \
  -g human -i TCR -n t_cell_001 -o results

# Mouse TCR: paired-end reads
basic-assemble -PE_1 mouse_R1.fastq.gz -PE_2 mouse_R2.fastq.gz \
  -g mouse -i TCR -n mouse_cell_001 -o results
```

For TCRs, matching reads and references can produce alpha/beta (`TRA`/`TRB`)
contigs, or mouse gamma/delta (`TRG`/`TRD`) contigs where those seeds are
supported. The command does not assume every cell expresses every chain.

Input files must already be separated by cell. BASIC accepts plain or gzipped
four-line FASTQ with Phred+33 quality, including comma-separated lanes. Paired
files must have matching identifiers and counts. BASIC does not demultiplex
barcodes or deduplicate UMIs.

Useful options:

- `-p 4` sets Bowtie2 threads; `-b /path/to/bowtie2` overrides PATH lookup.
- `--check` checks paths, indexes, and dependencies without assembling; it does
  not parse the complete FASTQ files.
- `-a` retains partial or empty biological results. Execution errors still fail.
- `--subsample N` uses the first N reads or N synchronized pairs across lanes.
- `--overwrite` explicitly replaces existing results; failed reruns preserve
  the previous successful output bundle.

Run `basic-assemble --help` for all options, including quality trimming,
alignment timeouts, temporary files, and experimental search settings.

## Results

Each sample produces:

| Output | Contents |
| --- | --- |
| `cell_001.fasta` | Assembled sequences with unique contig identifiers |
| `cell_001.airr.tsv` | AIRR columns and `basic_*` assembly metadata |
| `cell_001.summary.tsv` | Chain or reference-group status and ambiguity counts |
| `cell_001.run.json` | Parameters, versions, diagnostics, and input/reference/output checksums |
| `cell_001.logs/` | Bowtie2 diagnostics |

Results distinguish **complete**, **partial**, **unsupported**, and **failed**
outcomes. Complete means that a contig connects its two seeds; it does not
establish a productive receptor, a correct CDR3, or an error-free sequence.
Read support is not a molecule count or a calibrated confidence score.
Biological annotation fields remain empty unless annotation supplies them.

Use the run manifest as the completion marker and verify its output checksums
when consuming results automatically. [Algorithm details](docs/ALGORITHMS.md)
explain assembly status, ambiguity, scoring, and search limits.

## Multiple cells and annotation

Batch input is a tab-separated file with `sample_id` and either `fastq` or
`left`/`right` columns. Optional `genome` and `receptor` columns override the
command defaults. File paths are relative to the batch file.

```bash
basic-assemble --batch cells.tsv -o results --continue-on-error
```

Each cell has separate results; a batch summary records failures. Cells are
processed sequentially, with `-p` controlling each alignment's threads.

Optional local **IgBLAST** annotation adds available gene calls, junction/CDR3
information, and predicted productivity. It requires a separate IgBLAST
installation and appropriate germline databases; BASIC's Bowtie2 indexes are
not IgBLAST databases. See the [usage and annotation guide](docs/USAGE.md).

## Other organisms

Build references from selected local FASTA files, HTTPS sources, or NCBI
nucleotide accessions, then use the new species directory:

```bash
basic-references build --manifest species-references.json --output refs \
  --species my_species --receptor BCR
basic-assemble -SE cell.fastq.gz -g my_species -i BCR \
  -d refs -n custom_cell -o results
```

Reference sources must declare the appropriate receptor loci and orientation.
For organisms lacking suitable references, `--seed-species` with
`--local-seeding` can explore related-species homology. This is experimental:
it does not discover or validate a germline catalogue or justify transferring
definitive gene calls between species. The [non-model organism guide](docs/NON_MODEL_ORGANISMS.md)
covers searching, reference manifests, supported loci, and these limitations.

## Documentation and development

- [Installation and upgrades](docs/INSTALLATION.md)
- [Codex, Claude, and Antigravity skill](docs/AGENT_SKILL.md), with [external review feedback](docs/SKILL_REVIEW.md)
- [Extended usage, batches, and annotation](docs/USAGE.md)
- [Assembly algorithms and experimental options](docs/ALGORITHMS.md)
- [Validation results](docs/VALIDATION.md), [review record](docs/REVIEW.md), and [benchmark tools](benchmarks/README.md)
- [Changes and migration from BASIC 1.5](CHANGELOG.md)
- [Development checks and release instructions](docs/RELEASING.md), including [Bioconda](docs/BIOCONDA.md)

Run tests with `python -m pytest` after installing `requirements-test.lock`.
Real-aligner tests skip when Bowtie2 is unavailable; set
`BASIC_BOWTIE2=/path/to/bowtie2` to require a particular executable.

## Citation and provenance

Canzar S, Neu KE, Tang Q, Wilson PC, Khan AA. **BASIC: BCR assembly from single
cells.** *Bioinformatics* 33(3), 1 February 2017.

The software is distributed under the [MIT license](LICENSE). Bundled references
are preserved from BASIC 1.5.0; their original download dates are unknown. See
[reference provenance](db/provenance.json) and the [historical license notice](LICENSE-NOTICE.md).
Reference data and external tools retain their respective terms.
