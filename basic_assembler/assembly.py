"""Deterministic, bounded overlap assembly without external dependencies.

IUPAC symbols are preserved and compared literally during overlap matching.
The default scoring retains BASIC's entropy rule; evidence scores and beam
search are optional heuristics, not calibrated probabilities of correctness.
"""

from __future__ import annotations

from bisect import bisect_left
from dataclasses import dataclass, replace
import math
from numbers import Integral, Real
from types import MappingProxyType
from typing import Mapping


_COMPLEMENTS = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
_ALPHABET = frozenset("ACGTRYSWKMBDHVN")


def normalize_sequence(sequence: str) -> str:
    """Normalize case, whitespace, and RNA U; reject non-IUPAC symbols."""
    if not isinstance(sequence, str):
        raise TypeError("A nucleotide sequence must be a string")
    if set(sequence) <= _ALPHABET:
        # Retain already normalized read objects when indexing a large pool.
        return sequence
    sequence = "".join(sequence.split()).upper().replace("U", "T")
    invalid = set(sequence) - _ALPHABET
    if invalid:
        raise ValueError("Invalid nucleotide symbol(s): " + ", ".join(sorted(invalid)))
    return sequence


def complement(sequence: str) -> str:
    return normalize_sequence(sequence).translate(_COMPLEMENTS)


def revcom(sequence: str) -> str:
    return complement(sequence)[::-1]


@dataclass(frozen=True)
class AssemblyOptions:
    min_overlap: int = 15
    min_support: int = 2
    max_steps: int = 10000
    max_contig_length: int = 30000
    scoring: str = "legacy"
    beam_width: int = 1

    def __post_init__(self) -> None:
        for name in ("min_overlap", "min_support", "max_steps", "max_contig_length", "beam_width"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Integral) or value < 1:
                raise ValueError(f"{name} must be a positive integer")
        if self.max_contig_length < self.min_overlap:
            raise ValueError("max_contig_length must be at least min_overlap")
        if self.scoring not in {"legacy", "evidence"}:
            raise ValueError("scoring must be 'legacy' or 'evidence'")


class OverlapIndex:
    """One sorted read pool supporting arbitrary exact prefix lookups.

    Two binary searches bound each lookup. Unlike indexing only the final
    minimum-length suffix, this finds every read starting with a longer
    contig suffix, including matches in which its last k-mer is internal.
    Reads must already have the desired strand counts; this class never adds
    reverse complements or doubles palindromic reads.
    """

    def __init__(self, counts: Mapping[str, int]):
        normalized: dict[str, int] = {}
        for raw_sequence, count in counts.items():
            if isinstance(count, bool) or not isinstance(count, Integral) or count < 0:
                raise ValueError("Read counts must be nonnegative integers")
            sequence = normalize_sequence(raw_sequence)
            if sequence and count:
                normalized[sequence] = normalized.get(sequence, 0) + int(count)
        self.counts: Mapping[str, int] = MappingProxyType(normalized)
        self.sequences: tuple[str, ...] = tuple(sorted(normalized))
        self.max_read_length = max(map(len, self.sequences), default=0)

    def lookup(self, prefix: str) -> tuple[str, ...]:
        """Return read sequences beginning with an already normalized prefix."""
        lower = bisect_left(self.sequences, prefix)
        upper = bisect_left(self.sequences, prefix + "Z")
        return self.sequences[lower:upper]


@dataclass(frozen=True)
class _Candidate:
    sequence: str
    overlap: int
    support: int
    entropy: float
    quality_support: float
    mate_support: int
    evidence_score: float


@dataclass(frozen=True)
class _Path:
    sequence: str
    steps: int = 0
    ambiguous: bool = False
    branch_points: int = 0
    support_sum: int = 0
    support_min: int | None = None
    evidence_sum: float = 0.0
    cost: float = 0.0
    termination: str = ""
    reads: tuple[str, ...] = ()


@dataclass
class _ExpansionBudget:
    limit: int
    used: int = 0


def _prepare_evidence(index, quality_weights, mate_links):
    quality: dict[str, float] = {}
    for raw_sequence, weight in (quality_weights or {}).items():
        sequence = normalize_sequence(raw_sequence)
        if isinstance(weight, bool) or not isinstance(weight, Real) or not math.isfinite(weight) or weight < 0:
            raise ValueError("Quality weights must be finite nonnegative numbers")
        quality[sequence] = quality.get(sequence, 0.0) + float(weight)
    for sequence, weight in quality.items():
        if weight > index.counts.get(sequence, 0) + 1e-9:
            raise ValueError("A quality weight cannot exceed its read count")
    links: dict[str, dict[str, int]] = {}
    for raw_source, neighbors in (mate_links or {}).items():
        source = normalize_sequence(raw_source)
        target_counts = links.setdefault(source, {})
        for raw_target, support in neighbors.items():
            if isinstance(support, bool) or not isinstance(support, Integral) or support < 0:
                raise ValueError("Mate link counts must be nonnegative integers")
            target = normalize_sequence(raw_target)
            if support:
                target_counts[target] = target_counts.get(target, 0) + int(support)
    return quality, links


def _candidates(path, index, options, quality, links, read_length_cap=None):
    sequence = path.sequence
    maximum = index.max_read_length
    if read_length_cap is not None:
        maximum = min(maximum, read_length_cap)
    # Full-seed overlaps are valid when a seed is shorter than an input read.
    upper = min(len(sequence), maximum - 1)
    candidates: list[_Candidate] = []
    for overlap in range(upper, options.min_overlap - 1, -1):
        matches = index.lookup(sequence[-overlap:])
        if not matches:
            continue
        denominator = float(sum(index.counts[read] for read in matches)) + 0.5
        probabilities = [index.counts[read] / denominator for read in matches]
        entropy = -sum(p * math.log(p, 2) for p in probabilities)
        for read in matches:
            support = index.counts[read]
            if support < options.min_support:
                continue
            quality_support = quality.get(read, float(support))
            mate_support = sum(links.get(source, {}).get(read, 0) for source in path.reads)
            evidence_score = quality_support + mate_support
            if options.scoring == "evidence" and evidence_score <= 0:
                continue
            candidates.append(_Candidate(read, overlap, support, entropy,
                                         quality_support, mate_support, evidence_score))
    if options.scoring == "legacy":
        # Original <= entropy comparison prefers the shortest overlap on a
        # tie. Lexical ordering makes ties within that overlap reproducible.
        candidates.sort(key=lambda c: (c.entropy, c.overlap, -c.support, c.sequence))
    else:
        candidates.sort(key=lambda c: (-c.evidence_score, -c.overlap,
                                       -c.support, c.entropy, c.sequence))
    return candidates


def _has_branch(candidates, sequence):
    additions = sorted({c.sequence[c.overlap:] for c in candidates
                        if c.sequence.strip("N") not in sequence and len(c.sequence) > c.overlap})
    # Nested prefixes merely represent different lengths of the same path.
    return any(not right.startswith(left) for left, right in zip(additions, additions[1:]))


def _has_target(sequence: str, targets: tuple[str, ...]) -> bool:
    return any(target and target in sequence for target in targets)


def _path_key(path, targets):
    return (not _has_target(path.sequence, targets), path.cost,
            -len(path.sequence), path.sequence, path.termination)


def _best_paths(paths, width, targets):
    distinct = {}
    for path in sorted(paths, key=lambda item: _path_key(item, targets)):
        distinct.setdefault(path.sequence, path)
    return list(distinct.values())[:width]


def _extend_paths(seed, index, options, quality, links, *, targets=(), read_length_cap=None,
                  budget=None):
    """Return up to beam_width terminal paths within a shared expansion budget."""
    seed = normalize_sequence(seed).strip("N")
    if len(seed) > options.max_contig_length:
        raise ValueError("Seed exceeds max_contig_length")
    initial = _Path(seed, reads=(seed,) if links else ())
    if len(seed) < options.min_overlap:
        return [replace(initial, termination="short_seed")]
    budget = budget if budget is not None else _ExpansionBudget(options.max_steps)
    active = [initial]
    finished: list[_Path] = []
    while active:
        next_paths: list[_Path] = []
        for position, path in enumerate(active):
            if budget.used >= budget.limit:
                finished.extend(replace(item, termination="max_steps")
                                for item in active[position:] + next_paths)
                return _best_paths(finished, options.beam_width, targets)
            budget.used += 1
            candidates = _candidates(path, index, options, quality, links, read_length_cap)
            if not candidates:
                finished.append(replace(path, termination="no_supported_overlap"))
                continue
            branching = _has_branch(candidates, path.sequence)
            # Keep at most one occurrence of an identical next sequence;
            # different overlap lengths must not consume the whole beam.
            successors = set()
            accepted = 0
            for rank, candidate in enumerate(candidates):
                repeated = candidate.sequence.strip("N") in path.sequence
                sequence = path.sequence if repeated else (
                    path.sequence + candidate.sequence[candidate.overlap:]).strip("N")
                if sequence in successors:
                    continue
                successors.add(sequence)
                cost = path.cost + math.log2(rank + 1)
                common = dict(ambiguous=path.ambiguous or branching,
                              branch_points=path.branch_points + int(branching), cost=cost)
                if repeated or sequence == path.sequence:
                    finished.append(replace(path, termination="repeated_read", **common))
                elif len(sequence) > options.max_contig_length:
                    finished.append(replace(path, termination="max_contig_length", **common))
                else:
                    next_paths.append(_Path(
                        sequence=sequence, steps=path.steps + 1,
                        support_sum=path.support_sum + candidate.support,
                        support_min=candidate.support if path.support_min is None else min(path.support_min, candidate.support),
                        evidence_sum=path.evidence_sum + candidate.evidence_score,
                        reads=path.reads + (candidate.sequence,) if links else (),
                        **common))
                accepted += 1
                if accepted >= options.beam_width:
                    break
        finished = _best_paths(finished, options.beam_width, targets)
        active = _best_paths(next_paths, options.beam_width, targets)
    return _best_paths(finished, options.beam_width, targets)


def _join_metadata(first, second, sequence):
    minimums = [p.support_min for p in (first, second) if p.support_min is not None]
    return replace(second, sequence=sequence, steps=first.steps + second.steps,
                   ambiguous=first.ambiguous or second.ambiguous,
                   branch_points=first.branch_points + second.branch_points,
                   support_sum=first.support_sum + second.support_sum,
                   support_min=min(minimums) if minimums else None,
                   evidence_sum=first.evidence_sum + second.evidence_sum,
                   cost=first.cost + second.cost,
                   termination=(first.termination + ";" + second.termination))


def _bidirectional(seed, reverse_first, target, index, options, quality, links, budget):
    targets = (target, revcom(target)) if target else ()
    first = _extend_paths(seed, index, options, quality, links, targets=targets, budget=budget)
    results = []
    for path in first:
        second = _extend_paths(revcom(path.sequence), index, options, quality, links,
                               targets=targets, budget=budget)
        for continuation in second:
            sequence = continuation.sequence if reverse_first else revcom(continuation.sequence)
            results.append(_join_metadata(path, continuation, sequence))
    return _best_paths(results, options.beam_width, targets)


def _contig(path, region):
    return {"sequence": path.sequence, "region": region,
            "ambiguous": path.ambiguous, "termination": path.termination,
            "steps": path.steps, "branch_points": path.branch_points,
            "read_support_min": path.support_min, "read_support_sum": path.support_sum,
            "evidence_score_sum": path.evidence_sum}


def assemble_chain(
    chain_start: str,
    chain_end: str,
    counts: Mapping[str, int],
    options: AssemblyOptions | None = None,
    *,
    index: OverlapIndex | None = None,
    quality_weights: Mapping[str, float] | None = None,
    mate_links: Mapping[str, Mapping[str, int]] | None = None,
) -> dict:
    """Assemble one chain, returning sequence records and explicit status.

    ``chain_start`` is forward oriented; ``chain_end`` is reverse complemented,
    matching the historical BASIC anchor convention. ``counts`` already holds
    any required strand orientations, counting each palindrome only once.
    A supplied index is authoritative and must have been built from this pool.
    Complete means seed-to-seed containment, not validated receptor productivity.
    ``max_steps`` caps the total number of expanded states across all beam
    paths, both directions, and the variable-seed fallback for this chain.
    """
    options = options or AssemblyOptions()
    index = index if index is not None else OverlapIndex(counts)
    quality, links = _prepare_evidence(index, quality_weights, mate_links)
    budget = _ExpansionBudget(options.max_steps)
    diagnostics = {
        "scoring": options.scoring, "beam_width": options.beam_width,
        "unique_oriented_reads": len(index.counts),
        "quality_weighted": options.scoring == "evidence" and bool(quality),
        "mate_weighted": options.scoring == "evidence" and bool(links),
        "score_interpretation": "heuristic, not calibrated confidence",
        "support_interpretation": "sum/minimum over extension steps, not unique molecules",
        "completeness_definition": "seed-to-seed containment",
        "unsupported_seeds": [],
    }
    seeds = []
    for label, raw_sequence in (("variable", chain_start), ("constant", chain_end)):
        seed = normalize_sequence(raw_sequence).strip("N")
        if len(seed) > options.max_contig_length:
            raise ValueError(f"{label} seed exceeds max_contig_length")
        if not seed or len(seed) < options.min_overlap:
            diagnostics["unsupported_seeds"].append(label)
            seed = ""
        seeds.append(seed)
    start, end = seeds
    constant_paths = _bidirectional(end, True, start, index, options, quality, links, budget) if end else []
    full_paths = [p for p in constant_paths if start and _has_target(p.sequence, (start, revcom(start)))]
    variable_paths = []
    if not full_paths and start:
        variable_paths = _bidirectional(start, False, end, index, options, quality, links, budget)
        full_paths = [p for p in variable_paths if end and _has_target(p.sequence, (end, revcom(end)))]
    if full_paths:
        selected = _best_paths(full_paths, options.beam_width, ())
        if len(selected) > 1:
            selected = [replace(path, ambiguous=True) for path in selected]
        contigs = [_contig(path, "full") for path in selected]
        status = "complete"
    else:
        contigs = [_contig(path, "variable") for path in variable_paths]
        contigs += [_contig(path, "constant") for path in constant_paths]
        status = "partial" if contigs else "unsupported"
    diagnostics["bounded"] = any(
        "max_steps" in c["termination"] or "max_contig_length" in c["termination"] for c in contigs)
    diagnostics["contig_count"] = len(contigs)
    diagnostics["expanded_states"] = budget.used
    diagnostics["max_steps"] = options.max_steps
    return {"status": status, "contigs": contigs, "diagnostics": diagnostics}


def splitread(s: str, rl: int) -> list[str]:
    """Historical suffix helper, safe for short or empty input."""
    sequence = normalize_sequence(s)
    return [sequence[-length:] for length in range(min(rl, len(sequence)) - 1, 14, -1)]


def extend(s, d, verb=False, rl=None, reverse_comp=False):
    """Compatibility wrapper for BASIC.extend with deterministic legacy scoring."""
    index = OverlapIndex(d)
    options = AssemblyOptions()
    paths = _extend_paths(s, index, options, {}, {}, read_length_cap=rl)
    result = paths[0].sequence
    if verb:
        print(revcom(result) if reverse_comp else result)
    return result


def stitch(chain_end, chain_start, const, variable, verb, read_length, cellid,
           chain_name, result_queue):
    """Historical queue/FASTA interface using the bounded extension wrapper."""
    const_token = extend(revcom(extend(chain_end, const, verb, read_length, True)),
                         variable, verb, read_length) if chain_end else ""
    if chain_start and _has_target(const_token, (chain_start, revcom(chain_start))):
        records = [("", const_token)]
    else:
        var_token = revcom(extend(revcom(extend(chain_start, variable, verb, read_length)),
                                 const, verb, read_length, True)) if chain_start else ""
        if chain_end and _has_target(var_token, (chain_end, revcom(chain_end))):
            records = [("", var_token)]
        else:
            records = [("[variable_region_contig]", var_token),
                       ("[constant_region_contig]", const_token)]
    result = "".join(f">cell_id={cellid};{chain_name}_chain{region};BASIC\n{sequence}\n"
                     for region, sequence in records if sequence)
    result_queue.put({chain_name: result})
