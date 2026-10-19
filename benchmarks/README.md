# Assembly evaluation

## Try the included example

After installing BASIC and Bowtie2, run this from the repository root. It
assembles the bundled human B-cell reads, measures the command's elapsed time,
and compares the result with the supplied BASIC assembly:

```bash
python benchmarks/benchmark.py \
  --expected-fasta examples/cell_A1/result.fasta \
  --actual-fasta results/cell_A1-timed/result.fasta \
  --label included-human-bcr-regression \
  --output results/cell_A1-timed/comparison.json \
  --run basic-assemble -SE examples/cell_A1/ERR1421621_1M_reads.fastq.gz \
    -g human -i BCR -n result -o results/cell_A1-timed -p 2
```

Use a fresh output directory for each run. Keep `-n result` so the generated
FASTA identifiers match the supplied file. In `comparison.json`, expect
`sequences.exact_match_count` to be `2` and exact recall and precision to be
`1.0`. `execution.wall_seconds` records this run's elapsed time, including
input preparation, alignment, assembly, and report writing. It varies with
the machine and workload.

This is a regression check against an existing assembly. To measure biological
accuracy, use independently validated receptor sequences as described below.
For comparisons between tools or settings, use the same reads and reference
versions, repeat timing runs under comparable conditions, and record the
hardware and thread counts. The run manifest records BASIC's configuration.

## Compare your own assemblies

`benchmark.py` compares assemblies with independently curated truth and writes
JSON. It uses only Python's standard library. FASTA identifiers (the first token
after `>`) must match between files; AIRR rows use the same `sequence_id`.
Normalize identifiers in a separate copy when comparing tools. Nucleotide
comparisons are case-insensitive and orientation-sensitive; `N` is a literal
base, not a wildcard. FASTA input may be gzipped.

```sh
python benchmarks/benchmark.py \
  --expected-fasta truth.fasta --actual-fasta assembled.fasta \
  --expected-airr truth.tsv --actual-airr assembled.airr.tsv \
  --label human-bcr-low-coverage --output metrics.json
```

The report includes exact sequence recall/precision, expected sequences fully
contained in assembled contigs, strict partial reconstructions, mismatches,
and missing/extra IDs. Containment includes exact matches. It does not establish
that additional flanking sequence is correct or that a contig is free of
chimeric joins. Strict partial sequences match a proper substring of the truth.
Missing truth and empty denominators yield JSON `null`, not artificial scores.

Annotation metrics separately assess nucleotide/amino-acid junction and CDR3,
productivity, gene calls, and locus. Only nonempty truth values are scored;
missing predictions count against recall. Reported-value accuracy excludes
missing predictions and must be read alongside recall. Gene call comparisons
are exact strings, including alleles and ordering of ambiguous calls. According
to the [AIRR definitions](https://docs.airr-community.org/en/stable/datarep/rearrangements.html),
junction includes the conserved flanking residues whereas CDR3 excludes them;
the script never substitutes one for the other.

Optionally append `--run COMMAND ARGUMENTS...` to execute a command with a timeout
before comparison and record wall-clock time. The command is invoked directly
without a shell; quote each argument as required by your shell. All benchmark
options must precede `--run`. Failed commands stop evaluation, so stale FASTA
files are not accidentally scored. Use an isolated output directory for each
run. Peak memory is not measured by this wrapper; measure the actual assembler
with a platform-appropriate resource profiler (for example, GNU `/usr/bin/time
-v` on Linux), not the Python comparison process.

## Dataset matrix and interpretation

Before accepting changes to assembly defaults, collect independently validated
human and mouse BCR and TCR sequences with their raw reads and provenance. Run
fixed-seed coverage subsamples and include low coverage, sequencing errors,
mixed read lengths, multiple expressed chains, and ambiguous branches. Compare
default and experimental modes on the same reads and reference versions.
Compare results separately for each combination of species, receptor, and
coverage so abundant easy cases cannot hide failures in difficult cases.
Curated junction/CDR3 truth and targeted inspection of joins are needed to
evaluate biological correctness beyond containment.

The small fixtures in `tests/test_benchmarks.py` verify metric arithmetic. They
are synthetic software tests, not biological performance evidence. No expanded
biological benchmark panel or accuracy gain is claimed by this release.

## Compare assembly settings

Use the same per-cell reads and sample names, with a fresh output directory
for every configuration. Keep `-a` enabled so partial and empty results remain
available for scoring. Start with the default and combined runs below; use
the separate beam and anchor runs to understand which change helps. These
are suggested trial settings, not validated accuracy presets.

| Run | Additional assembly options | Question |
| --- | --- | --- |
| Default | None | What does the standard search recover? |
| Beam | `--beam-width 4` | Does keeping alternative paths help? |
| Anchors | `--max-anchors 2` | Do other starting reads help? |
| Combined: broader search | `--beam-width 4 --max-anchors 2` | Do the two changes help together? |
| Evidence | Combined options plus `--scoring evidence --quality-weighted --mate-aware` | Does quality and mate support improve paired-read results? |
| Local seeds | `--local-seeding` | Does local reference alignment recover useful seeds missed by the default alignment? |

For single-end data, omit `--mate-aware`. Test `--max-candidates 2` separately
if secondary predictions matter. It retains additional seed-pair results;
it does not expand the search, and it is not a two-contig cap. Freeze a rule
for selecting the primary prediction before looking at truth. Report primary
predictions and all-candidate results separately rather than selecting whichever
contig happens to match the expected sequence.

For example, from the downloaded source, time the combined settings on the
included example without installing BASIC as a package:

```bash
python3 benchmarks/benchmark.py \
  --expected-fasta examples/cell_A1/result.fasta \
  --actual-fasta results/cell_A1-beam/result.fasta \
  --label included-human-bcr-beam \
  --output results/cell_A1-beam/comparison.json \
  --run python3 BASIC.py -SE examples/cell_A1/ERR1421621_1M_reads.fastq.gz \
    -g human -i BCR -n result -o results/cell_A1-beam -p 2 -a \
    --beam-width 4 --max-anchors 2
```

This is still a regression comparison with the stored BASIC result, not an
independent test of accuracy. Alternative paths can introduce additional
FASTA identifiers; inspect missing and extra IDs instead of silently removing
them to improve a score.

Choose settings on a development subset and report performance on held-out
cells. Include failures and missing chains. Measure correct recovery,
junction/CDR3 accuracy, extra predictions, sequence errors, runtime, and peak
memory together. Longer contigs or more `complete` calls alone do not establish
an improvement. Use independent truth and an external alignment workflow if
you need sequence identity or constant-region coverage; this script does not
align sequences or automatically match chains across tools.

Use the same IgBLAST version and annotation databases across configurations.
For coverage experiments, subsample with fixed random seeds while preserving
pairs: BASIC's `--subsample` takes a file prefix, not a random sample. More
search costs time and memory, so retain the smallest configuration that
improves correct recovery without unacceptable errors or resource use.

## Optional IgBLAST annotation

The pipeline can run a local `igblastn` with AIRR output (`-outfmt 19`, supported
[since IgBLAST 1.9](https://ncbiinsights.ncbi.nlm.nih.gov/2018/05/24/igblast-1-9-airr-rearrangement-reporting/)).
Supply V and J BLAST database prefixes, a D database for
VDJ loci, and a matched J auxiliary file. An optional C database enables C gene
calls that may identify isotype when sufficient constant-region sequence is
present. Annotation never infers productivity or CDR3 from an assembly status.

Follow the [official installation instructions](https://ncbi.github.io/igblast/cook/How-to-set-up.html)
to obtain/build your databases. BASIC does not download them. Set `--igdata`
to the directory containing `internal_data` when IgBLAST cannot find its
installation data. The wrapper sets `IGDATA` only for its child process.
Without a matching J auxiliary file, frame/CDR3 annotations can be missing or
incorrect; do not interpret empty annotations as negative biological findings.
The [official examples](https://ncbi.github.io/igblast/cook/examples.html)
describe `-ig_seqtype TCR`, database selection, and `-c_region_db`.
