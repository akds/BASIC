# BASIC 2.0 review record — 2026/10/18

## Scope

Codex reviewed the assembler, read preparation, Bowtie2 integration, CLI and
batch execution, output publication, optional annotation, reference builder,
benchmarks, packaging, and documentation. Separate review passes examined
the algorithm/reference layer, I/O/integration layer, and annotation/metrics
layer. Confirmed findings received focused regression tests.

This is a software review, not certification of biological performance.
See [validation evidence and limits](VALIDATION.md) before interpreting the
new experimental assembly options.

## Confirmed findings addressed

| Area | Changes and regression coverage |
| --- | --- |
| Assembly | Indexed exact overlap lookup; deterministic entropy ties; short/mixed-length reads; bounded state expansion; supported alternatives; distinct complete/partial/unsupported outcomes. |
| Seed selection | Aggregate identical seeds across reference alleles; prevent joining incompatible loci; preserve unmatched one-sided candidates; retain all four possible TCR loci. |
| TCR labels | Infer loci from mapped references, including mouse ordering and gamma/delta candidates; summary labels follow the retained assembly. |
| FASTQ preparation | Validate paired identifiers and directions; synchronize subsampling; trim terminal Ns and low qualities together; preserve IUPAC bases through mapping; protect input aliases. |
| External programs | Check exit status, retain stderr, validate SAM/AIRR output, and terminate process groups on timeout/error where supported. |
| Batch safety | Check all proposed outputs against all batch inputs and the batch manifest before running any sample; include directory, symlink, hardlink, and case-insensitive aliases. |
| Result publication | Require explicit overwrite, retain prior successful results after ordinary failures, and preserve recoverable originals if rollback itself fails. |
| Annotation | Preserve unknown values, validate sequence identity and output completeness, distinguish seed locus from IgBLAST locus, and fingerprint discoverable annotation references. |
| Reference construction | Validate sources, loci and orientation; bound downloads; detect truncated responses; record hashes; check each index build; preserve prior indexes on failure. |
| Efficiency | Reuse one prepared read pool and overlap index; retain optional evidence only when requested; optimize the profiled quality-character validation bottleneck. |

## Validation record

**307 pytest cases passed, with no failures or skips**: the 295-case assembly
and skill suite in 45.00 seconds, followed by 12 new release-asset tests in
0.06 seconds. Nine additional Bioconda-generator unit tests passed. The main
suite includes real Bowtie2 2.5.5 assembly regressions. Combined statement/branch
coverage of `basic_assembler` is **93.38%**; statement coverage is 94.56%, and
branch coverage is 90.63%. The assembler has every statement covered and 99% combined coverage.
Module summaries and hashes of the tested implementation are saved in
[the 2.0.0 validation record](../benchmarks/results/validation-2.0.0.json).

Both the wheel and source distribution build successfully. A fresh environment
outside the checkout passed all four installed-command checks, discovered all
four human/mouse BCR/TCR index sets, and reproduced the single-end example
FASTA exactly after normalizing only the sample-name header field. The source archive includes tests, examples,
documentation, the agent skill, release helpers, and GitHub Actions workflows. Native IgBLAST biological
validation is still pending.

Coverage is measured in the pytest process. Some module entrypoints and
subprocess execution are tested externally and therefore are not credited by
that measurement. Untested paths include platform-specific process cleanup
and unusual filesystem races. The goal is meaningful behavioral coverage,
not a claim that every possible failure has been tested.

## Independent Gemini review

**Independent source review remains incomplete.** Gemini 4 Argon (High) was
requested through the `agy` CLI. Full-source requests and a final set of three
component requests returned empty responses after their time limits; each
component was allowed five minutes. The CLI returned exit code 0 with a
`SUCCESS` envelope while stderr explicitly reported a print timeout. Those
responses are not treated as completed reviews or approvals.

A shorter, high-level request did return three priorities: branching/chimeric
joins, repeat-driven cycles, and strand/stitch/annotation consistency. These
were compared with the tests and documented limitations, but that response did
not inspect the implementation and does not replace a source review. Branching,
search bounds, reverse complements, disconnected seeds, and AIRR serialization
have test coverage; biological chimera rates and native IgBLAST annotation
accuracy remain unestablished. See the [Gemini request record](GEMINI_REVIEW.md).

## Release scope

The source is prepared for `akds/BASIC`, with a clean README identifying
**Version 2.0 (2026/10/18)** and distribution version `2.0.0`. The maintainer
confirmed the standard MIT license. PyPI/TestPyPI publication and GitHub release
asset uploads have a prepared workflow; Bioconda has a checksum-verified recipe
generator. Local tests are not evidence of a completed remote CI run, registry
upload, or conda-build. See [release instructions](RELEASING.md).

The new portable skill received substantive, successful critiques from both
Claude Opus 5.5 High and Gemini 4 Argon High via agy. Their concrete suggestions
were incorporated and the standalone bundle validated. This **skill review**
is separate from the incomplete earlier Gemini full-source review above; see
[the skill review record](SKILL_REVIEW.md).

## 2.0.1 correction after release

A later macOS/Python 3.10 CI run exposed an intermittent process-group signal
permission error during malformed-SAM cleanup. The cleanup error masked the
original `AlignmentError`, despite earlier green CI runs. Version 2.0.1 adds
a shared best-effort cleanup helper: it attempts group termination, falls back
to the direct child, bounds failure-path reaping, and preserves the original
error. Alignment, annotation, reference building, and benchmarking use it.

Deterministic tests force denied signals, vanished processes, failed fallback
termination, and failed pipe closure. Real subprocess tests cover malformed
SAM and alignment/version/annotation/builder/benchmark timeouts. Existing tests
for descendants that outlive the wrapper remain in place.

The operating system still controls whether a process can be signalled. This
change cannot guarantee termination of inaccessible descendants; if such a
descendant holds the streaming aligner's stdout open, the timeout can still
be exceeded. It does not change assembly scoring or establish new biological
accuracy. Latest validation is recorded in
[validation.json](../benchmarks/results/validation.json).

The full local 2.0.1 suite passed **319 tests** with real Bowtie2 2.5.5 in
45.63 seconds. Statement coverage was **95.08%**, branch coverage **91.20%**,
and combined coverage **93.92%**. The shared cleanup helper had 100% statement
and branch coverage in these tests; this does not imply all operating-system
conditions are covered. The nine Bioconda recipe-generator tests also passed.

A fresh local wheel installation passed all four entrypoints, all sixteen
reference index prefixes, both assistant-skill exports, MIT metadata/license
checks, and the installed-package Bioconda smoke harness. The full single-end
example was byte-identical to the supplied result, and the benchmark reported
two exact matches. This is wheel validation, not a completed conda build.

The README and benchmark introduction were separately copyedited using a
fresh Codex review plus successful Claude and Gemini critiques through agy.
The [editorial review record](README_REVIEW.md) documents the changes and the
suggestions corrected before adoption.
