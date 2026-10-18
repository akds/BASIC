# Assembly algorithms and limits

This document describes `basic_assembler/assembly.py`. Defaults preserve the
historical BASIC entropy rule while making ties deterministic and execution
bounded. Optional evidence scores and beam search are experimental assembly
heuristics. They are not calibrated measures of biological correctness.

## Sequence and read representation

Sequences are converted to uppercase, whitespace is removed, and `U` becomes
`T`. The accepted alphabet is `ACGTRYSWKMBDHVN`; complements cover all these
IUPAC symbols. Overlap matching is exact and literal: `R` matches `R`, but does
not act as a wildcard for `A` or `G`. Terminal `N` bases are stripped from
assembly seeds and after extensions. A seed shorter than `min_overlap` after
trimming is unsupported.

The index receives counts for already-oriented reads. It does not add reverse
complements. The FASTQ layer adds a physical read to its sequence and reverse
complement entries, counting a self-reverse-complement read once. Counts are
read observations, not UMI-deduplicated molecule counts. Normalized duplicate
keys are combined; zero-count or empty reads are omitted.

## Exact overlap lookup

`OverlapIndex` holds a snapshot of the counts and a sorted tuple of distinct
read strings. For each suffix of the current contig, two binary searches find
the contiguous range of reads with that prefix. The index supports every
eligible overlap length, rather than looking up only the shortest suffix as
a read prefix. For example, suffix `CATCG` can match read `CATCGAA` even though
the contig's last three bases, `TCG`, are internal to that read.

With maximum read length `L`, contig length `C`, and minimum overlap `k`,
eligible overlaps are:

```text
k <= overlap <= min(C, L - 1)
```

This permits an entire short seed to overlap the prefix of a longer read. For
seeds at least `L` bases long, the maximum overlap remains `L - 1`, as in the
original algorithm. All matching reads contribute to the overlap's entropy
statistic, including reads below the extension support threshold.

For `N` distinct reads, index construction validates the input bases and sorts
`N` strings: `O(N log N)` string comparisons, plus validation. Auxiliary storage
is `O(N)` count entries and references; already-normalized read strings are
reused. It does not materialize every read prefix.

One extension examines at most `L - k` suffix lengths. Ignoring string
comparison and copying costs, lookup costs `O(L log N + M)`, where `M` counts
all matches across those lengths. Ranking `K` supported candidates costs
`O(K log K)`. Actual string comparisons can inspect up to `L` bases, and
constructing all suffix strings can copy `O(L^2)` characters. Repeat checks
also search the current contig. Mate scoring adds a scan of the path's source
reads for each candidate. Repetitive pools can have many matches and remain
expensive; this index does not guarantee a fixed speedup. It removes the
historical full-pool substring scan from every extension step.

## Default legacy scoring

For the reads matching an overlap, let `c_i` be a read's count. The historical
statistic is:

```text
p_i = c_i / (sum(c_j) + 0.5)
H = -sum(p_i * log2(p_i))
```

The probabilities do not sum to one because the extra `0.5` appears only in
the denominator. This is the original entropy-like rule, not a newly
calibrated Shannon-entropy model.

A candidate needs a raw count at least `min_support` (default `2`). Candidates
are ordered by:

1. Lowest `H` for their overlap.
2. Shortest overlap when `H` ties, preserving the original `<=` update rule.
3. Highest raw read count within that overlap.
4. Lexicographically smallest read sequence for a reproducible count tie.

With `beam_width=1`, the first candidate is followed. Its nonoverlapping suffix
is appended. If its sequence with terminal Ns removed already occurs anywhere
in the contig, the branch stops. This protects against repetition but can also
stop extension through a genuine biological repeat. Floating-point entropy
calculations are deterministic for the same normalized pool and interpreter;
they are not an exact rational comparison.

## Optional quality and mate evidence

With `scoring="evidence"`, the raw-count threshold still applies. For a
candidate read `r` on the current directional path `P`:

```text
E(r, P) = Q(r) + sum(mate_links[source][r] for source in P)
```

`Q(r)` is the supplied summed reliability weight. Missing weights default to
the raw count. Supplied weights must be finite, nonnegative, and no larger
than their read counts. A zero-score candidate is excluded.

When enabled in the FASTQ layer, each observation's weight is the product of
`1 - 10^(-Phred/10)` over its retained bases, summed by oriented sequence. This
uses an independence assumption and favors reads with higher base quality;
it is not the probability that an assembled receptor is correct.

Mate links count observed physical pairs. The FASTQ layer links all retained
orientation combinations, with each directed edge counted at most once per
physical pair. At assembly time, sources are the initial seed and reads
selected earlier in the **current directional search**. This source list
resets when the contig is reversed for extension from the other end. Mate
support does not enforce insert size, library orientation, or molecule-level
phasing. Correlated observations can contribute to several extension steps.

Evidence candidates are ordered by highest `E`, longest overlap, highest raw
count, lowest historical `H`, then lexicographic read sequence. The reported
score is uncalibrated. In particular, the sum of quality weight and mate count
is a ranking heuristic, not a likelihood ratio or confidence percentage.

## Optional bounded beam search

`beam_width > 1` retains alternative paths. For each expanded state, candidates
use the selected scoring rule above; identical successor sequences are
deduplicated. At most `beam_width` distinct successors are considered from
that state. A candidate's zero-based rank contributes `log2(rank + 1)` to a
cumulative path cost. Skipped duplicate successors retain their original
candidate ranks.

Across states, paths are ordered by:

1. Presence of the opposite seed in either orientation.
2. Lowest cumulative rank cost.
3. Longest sequence.
4. Lexicographic sequence and termination reason for remaining ties.

Active and finished path sets each retain at most `beam_width` distinct
sequences. Candidate generation can temporarily hold more states before
pruning. Beam search can rescue a lower-ranked branch that reaches the other
seed, but it is not exhaustive and does not enumerate every receptor chain.

## Seed selection and reported chains

For each historical reference group (`hv`/`hc` and `lv`/`lc`), the pipeline
retains up to `--max-anchors` seeds per locus on each side. Known, different
loci are never paired. An unknown locus can pair with either side, and seeds
without a compatible partner remain one-sided attempts. All selected loci
are retained, including custom TCR groups containing all four TCR loci.
For `V` selected variable seeds and `C` selected constant seeds, the number
of attempts is at most `V*C + V + C`, or one empty attempt if neither side has
seeds. The per-locus limit bounds this work; no fixed three-locus cap applies.

The pipeline ranks completed attempts first, followed by combined anchor
support and sequence order, and retains up to `--max-candidates` distinct
assemblies per reference group. The group's reported status and chain label
come from its best retained attempt. The run manifest keeps `region_group`
so a reported TCR locus can be traced to the historical reference group.

## Chain status, ambiguity, and stopping

`assemble_chain()` starts from the reverse-oriented constant seed, extends,
reverses, and extends again. If no resulting path includes the variable seed,
it tries the variable seed first and repeats the two-direction procedure.
The supplied index is reused in each search. The original constant-first
preference remains: a complete constant-first result does not trigger an
additional variable-first search.

Results have one of three statuses:

| Status | Meaning |
| --- | --- |
| `complete` | A returned contig contains both usable seeds, allowing either seed orientation in the containment check. |
| `partial` | At least one usable seed produced a contig, but no returned path connects both seeds. |
| `unsupported` | Neither seed is present and long enough to start assembly. |

“Complete” means seed-to-seed containment. It does **not** establish a complete
coding sequence, productive V(D)J rearrangement, correct CDR3, or absence of a
chimera. Annotation and external biological truth are separate checks.

A contig has `ambiguous=true` if an encountered state had supported, distinct
continuations whose appended sequences were not merely prefixes of one
another. The flag uses candidates admitted by the selected scoring mode; a
zero-evidence candidate is excluded in evidence mode. It considers supported
alternatives even if the beam later prunes them. Multiple distinct complete
contigs retained from a search are also marked ambiguous. `false` does not
prove uniqueness or correctness. `branch_points` counts branching states on
the chosen path, not alternative receptors in the sample.

Default bounds and stop conditions are:

| Setting or reason | Behavior |
| --- | --- |
| `min_overlap=15` | Shorter seed or overlap is unsupported. |
| `min_support=2` | A candidate needs at least two raw observations. |
| `max_steps=10000` | Shared cap on expanded states for the entire chain, including all beam paths, both directions, and fallback. A dead-end lookup also uses a state expansion. |
| `max_contig_length=30000` | An extension that would exceed the cap is not appended. An initial seed above the cap raises an error. |
| `no_supported_overlap` | No eligible supported candidate remains. |
| `repeated_read` | The selected trimmed read already occurs in the contig. |
| `short_seed` | A directional search's seed is too short. |

Contig `termination` combines the reasons from its two directional searches.
`diagnostics.expanded_states` records the shared budget use. A returned contig
that hits a step or length bound causes `diagnostics.bounded=true`. If the
bridge between both seeds was already found, a bounded contig can still have
status `complete`; the bound limits further extension and exploration.
Otherwise its status remains partial. Sums and minima of support describe
extension steps and must not be interpreted as distinct molecule counts.

## Limitations and reproducible performance check

There is no UMI/barcode-aware pooling, cell demultiplexing, insertion/deletion
tolerant overlap matching, or base-quality-aware consensus sequence. Exact
overlaps can fail across sequencing errors or true variants. Sequence counts
can reflect PCR duplication. Short shared overlaps can connect unrelated
transcripts. Evidence and beam modes require evaluation on independent truth
before making claims about accuracy, sensitivity, or chimera rates.

Run the self-contained overlap benchmark from the repository root:

```bash
python benchmarks/overlap_benchmark.py \
  --random-seed 7281 --background-reads 100000 \
  --read-length 100 --contig-length 2500 --min-overlap 15 \
  --output benchmarks/results/overlap.json
```

The script includes the original scan-based extension rule with safe bounds.
Its synthetic target has unique minimum-length k-mers. Every target read has
count three; distinct background reads have count one, so only the target is
supported for extension. It asserts that both algorithms return exactly the
same sequence and that the sequence equals the generated target.

JSON records exact parameters, sequence and pool checksums, the assembled
sequence, Python/platform information, termination, and separate timings for
the scan, index construction, and indexed extension. Generation and hashing
are outside the timed sections. The scan runs first, once; timings vary by
machine and workload and do not provide uncertainty intervals. This measures
neither peak memory nor FASTQ parsing, alignment, paired-read evidence, beam
search, or end-to-end performance. `benchmarks/benchmark.py` separately
compares output sequences and annotations to supplied biological truth.
