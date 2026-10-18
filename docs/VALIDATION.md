# BASIC 2.0 validation — 2026/10/18

This release was reviewed and exercised on macOS arm64, Python 3.12.14, pytest
9.1.1, and Bowtie2 2.5.5. The results below separate software checks from
biological evidence. Passing these checks does not guarantee that every input
is handled correctly or establish a general improvement in receptor accuracy.

## Software checks

The test suite covers:

- Exact indexed overlap behavior, entropy ties, ambiguous branches, mixed read
  lengths, IUPAC symbols, short seeds, support thresholds, search bounds,
  quality/mate scoring, beam paths, and legacy helper compatibility.
- Plain/gzipped FASTQ, Phred+33 validation, synchronized mates, terminal
  trimming, lane subsampling, input aliases, and sparse restoration of IUPAC
  bases that Bowtie2 represents as `N`.
- Bowtie2 command construction, SAM parsing, failed subprocesses, timeout
  cleanup, optional annotation, and empty biological results.
- Sample and batch CLI execution, reference/target species separation, locus
  matching, output IDs, overwrite protection, publication rollback, and
  retention of recovery files when rollback itself fails.
- AIRR serialization, unknown annotations, local IgBLAST command and error
  handling, database fingerprints, and benchmark metric calculations.
- Reference manifests, orientation/locus validation, bounded HTTPS downloads,
  truncated responses, checksums, NCBI search parsing, index-build failures,
  and preservation of existing references.

The complete suite is rerun after review fixes. Tests of subprocess error
paths use executable fixtures; those checks are distinct from the real
Bowtie2 tests described below. Coverage includes statements and branches;
it is a diagnostic, not proof that all possible behaviors have been tested.
The final test totals and coverage are recorded in [the review record](REVIEW.md).

Reproduce the suite from the checkout:

```bash
python -m pip install -r requirements-test.lock
python -m pip install --no-build-isolation -e .
BASIC_BOWTIE2=/path/to/bowtie2 python -m pytest \
  --cov=basic_assembler --cov-branch --cov-report=term-missing
```

Without Bowtie2, real-aligner tests skip. The explicit environment variable
requires the named executable. GitHub Actions defines additional Linux/macOS
and Python 3.10/3.12/3.14 checks; those remote jobs have not been run from this
source snapshot.

## Assembly regressions with real Bowtie2

| Case | Observed result |
| --- | --- |
| Supplied human BCR single-end example | Both chains exactly match the supplied reference result and the BASIC 1.5 baseline: heavy 878 nt, light 829 nt. |
| Supplied human BCR paired-end example | Both contigs contain their supplied Sanger sequences; heavy chain matches the 1.5 baseline exactly. |
| Synthetic human/mouse BCR and TCR | Exact recovery of constructed variable–junction–constant targets in the tested cases. |
| Synthetic mouse gamma/delta TCR | Correct `TRG`/`TRD` sequence recovery and locus labeling with the matching bundled reference entries. |
| Mixed 50/75 nt paired reads and low support | Tested complete/partial behavior and search diagnostics. |
| New custom reference database | All four indexes built with native `bowtie2-build`; real reads mapped to the intended loci, including the required constant-region orientation. |

The paired example's new light contig is **1,631 nt**, versus **1,640 nt** in
the original run. They share the first 1,584 nt and then differ in the extra
3′ tail. Both contain the available Sanger truth. The additional tail has not
been independently validated, so exact full-contig parity and improved tail
accuracy are not claimed. See the [comparison data](../benchmarks/results/biological-regression.json).

The synthetic cases validate implementation on constructed sequences. They
do not substitute for independently validated biological TCR datasets, diverse
organisms, somatic hypermutation, low-coverage samples, or sequencing errors.

## Measured efficiency

The [overlap benchmark](../benchmarks/results/overlap.json) uses 102,401 unique
100 nt reads and reconstructs the same 2,500 nt sequence with the original-style
scan and indexed implementations. In one local run:

| Operation | Seconds |
| --- | ---: |
| Repeated full-pool scan | 0.7076 |
| Build overlap index | 0.0916 |
| Indexed extension | 0.0109 |
| Index build plus extension | 0.1025 |

That is **6.90× faster including index construction** for this unambiguous
synthetic workload. It is one run, with the scan first; it is not an end-to-end
speedup claim. The reproducible command and comparison method are documented
in [the algorithm guide](ALGORITHMS.md).

Profiling the million-read single-end example identified FASTQ validation as
a major remaining cost. Equivalent compiled character-range validation now
replaces a Python-level check of every quality character, with boundary tests
for accepted and rejected qualities. The [profile summary](../benchmarks/results/profile-summary.json)
records the pre-optimization profile and a separate microbenchmark. No peak
memory number is reported because the local resource measurement was blocked.
Original/new example runtimes were collected under different concurrency
conditions and should not be used as a controlled speed comparison.

## External services and remaining limits

- A live NCBI nucleotide search returned structured accession candidates. The
  search requires biological curation before reference building; download
  success does not validate a sequence's locus or completeness.
- IgBLAST integration is tested with executable fixtures, including successful
  AIRR output and failures. A native IgBLAST run with curated biological
  databases has not been validated in this environment.
- Evidence scoring, beam search, multiple anchors/candidates, and related-species
  seeding remain experimental. No general sensitivity/specificity or CDR3
  accuracy gain is established for them.
- Exact overlaps require literal sequence agreement. IUPAC ambiguity is
  preserved, not treated as a wildcard; this is not an indel-tolerant assembler.
- `complete` means that seeds are connected. Receptor productivity, clonotype,
  alleles, and true dual-chain expression need additional evidence.
- Output publication supports recovery from ordinary write failures. Multiple
  file renames do not form a crash-atomic filesystem transaction.

The [benchmark guide](../benchmarks/README.md) describes the broader biological
evaluation needed before changing assembly defaults or claiming higher accuracy.
