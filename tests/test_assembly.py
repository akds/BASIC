"""Synthetic assembly cases with known paths, branches, and stop conditions."""

from __future__ import annotations

from itertools import permutations
from queue import SimpleQueue

import pytest

from basic_assembler.assembly import (
    AssemblyOptions,
    OverlapIndex,
    _extend_paths,
    _prepare_evidence,
    assemble_chain,
    complement,
    extend,
    normalize_sequence,
    revcom,
    splitread,
    stitch,
)


def extend_test_path(seed, counts, *, quality=None, mates=None, **option_values):
    options = AssemblyOptions(min_overlap=3, **option_values)
    index = OverlapIndex(counts)
    weights, links = _prepare_evidence(index, quality, mates)
    return _extend_paths(seed, index, options, weights, links)[0]


def reference_pool(reference, read_length=10, support=3):
    counts = {}
    for offset in range(len(reference) - read_length + 1):
        read = reference[offset:offset + read_length]
        counts[read] = support
        counts[revcom(read)] = support
    return counts


def test_iupac_complement_and_normalization():
    assert normalize_sequence(" aCg u\nry ") == "ACGTRY"
    assert complement("ACGTRYSWKMBDHVN") == "TGCAYRSWMKVHDBN"
    sequence = "ACGTRYSWKMBDHVN"
    assert revcom(revcom(sequence)) == sequence
    assert revcom("") == ""
    with pytest.raises(ValueError, match="Invalid nucleotide"):
        revcom("ACGT-")


@pytest.mark.parametrize("values", [
    {"min_overlap": 0}, {"min_support": -1}, {"max_steps": 0},
    {"max_contig_length": 1}, {"beam_width": 0}, {"beam_width": True},
    {"scoring": "probability"},
])
def test_invalid_options_fail_early(values):
    with pytest.raises(ValueError):
        AssemblyOptions(**values)


def test_index_normalizes_without_adding_strands_or_support():
    index = OverlapIndex({" acgt ": 1, "ACGT": 2, "CGTA": 4, "TTTT": 0})
    assert dict(index.counts) == {"ACGT": 3, "CGTA": 4}
    assert index.lookup("AC") == ("ACGT",)
    assert index.lookup("TA") == ()
    assert index.max_read_length == 4
    with pytest.raises(TypeError):
        index.counts["ACGT"] = 10
    with pytest.raises(ValueError):
        OverlapIndex({"ACGT": -1})


def test_long_overlap_is_found_when_minimum_suffix_is_internal_to_read():
    # CATCG is the real overlap. Indexing only read prefixes of length 3
    # against the contig's TCG suffix would incorrectly miss this extension.
    path = extend_test_path("TGCATCG", {"CATCGAA": 3})
    assert path.sequence == "TGCATCGAA"
    assert path.steps == 1


def test_shorter_seed_can_overlap_an_entire_read_prefix():
    path = extend_test_path("ACG", {"ACGTT": 2})
    assert path.sequence == "ACGTT"


def test_known_reference_reconstructs_exactly_in_forward_orientation():
    reference = "ATGACCTTAGCGTACAGTCCGATCTGACGTTACGGA"
    counts = reference_pool(reference)
    result = assemble_chain(reference[:10], revcom(reference[-10:]), counts,
                            AssemblyOptions(min_overlap=7))
    assert result["status"] == "complete"
    assert [c["sequence"] for c in result["contigs"]] == [reference]
    assert result["contigs"][0]["region"] == "full"
    assert not result["contigs"][0]["ambiguous"]


def test_compatible_prebuilt_index_can_be_shared_between_chains():
    counts = {"ACGTT": 2, "CGTTA": 2}
    index = OverlapIndex(counts)
    options = AssemblyOptions(min_overlap=3)
    expected = assemble_chain("ACGT", "TAAC", counts, options)
    assert assemble_chain("ACGT", "TAAC", counts, options, index=index) == expected


def test_equal_support_branch_is_deterministic_under_input_permutations():
    seed = "ACGTTGCACTGATCGT"
    candidates = [("CGTTGCACTGATCGTA", 2), ("CGTTGCACTGATCGTC", 2)]
    results = {extend(seed, dict(order), False, 16) for order in permutations(candidates)}
    assert results == {seed + "A"}


def test_branch_ambiguity_is_reported_without_probability_claim():
    result = assemble_chain("AACG", "", {"ACGA": 3, "ACGT": 2},
                            AssemblyOptions(min_overlap=3))
    assert result["status"] == "partial"
    assert any(c["ambiguous"] for c in result["contigs"])
    assert result["diagnostics"]["score_interpretation"] == "heuristic, not calibrated confidence"


def test_palindromic_read_does_not_gain_artificial_support():
    palindrome = "AGCTAGCTAGCTAGCT"
    seed = "T" + palindrome[:-1]
    assert revcom(palindrome) == palindrome
    assert extend(seed, {palindrome: 1}, False, 16) == seed
    assert extend(seed, {palindrome: 2}, False, 16) == seed + "T"


def test_terminal_ns_and_short_seeds_are_safe():
    assert extend("NNNACGTACGTACGTN", {}, False, 16) == "ACGTACGTACGT"
    assert splitread("ACG", 16) == []
    result = assemble_chain("NNNACGTN", "NNNN", {}, AssemblyOptions())
    assert result["status"] == "unsupported"
    assert result["contigs"] == []
    assert result["diagnostics"]["unsupported_seeds"] == ["variable", "constant"]


def test_disconnected_anchors_have_explicit_partial_status():
    result = assemble_chain("AACG", "GGTA", {}, AssemblyOptions(min_overlap=3))
    assert result["status"] == "partial"
    assert {(c["sequence"], c["region"]) for c in result["contigs"]} == {
        ("AACG", "variable"), ("TACC", "constant")}


def test_existing_seed_to_seed_containment_is_complete_without_extension():
    result = assemble_chain("AACG", revcom("TTAACGCC"), {},
                            AssemblyOptions(min_overlap=3))
    assert result["status"] == "complete"
    assert result["contigs"][0]["sequence"] == "TTAACGCC"


def test_search_stops_at_step_bound():
    path = extend_test_path("AACG", {"ACGT": 2, "CGTT": 2, "GTTC": 2}, max_steps=1)
    assert path.sequence == "AACGT"
    assert path.steps == 1
    assert path.termination == "max_steps"


def test_beam_and_bidirectional_search_share_one_expansion_budget():
    counts = {"CGTA": 3, "CGTC": 2, "GTCC": 2, "TCCT": 2, "CCTA": 2}
    result = assemble_chain("ACGT", revcom("CCTA"), counts,
                            AssemblyOptions(min_overlap=3, beam_width=3, max_steps=4))
    assert result["diagnostics"]["expanded_states"] == 4
    assert result["diagnostics"]["bounded"]


def test_budget_exhaustion_before_seed_bridge_remains_partial():
    reference = "ATGACCTTAGCGTACAGTCCGATCTGACGTTACGGA"
    result = assemble_chain(reference[:10], revcom(reference[-10:]), reference_pool(reference),
                            AssemblyOptions(min_overlap=7, max_steps=1))
    assert result["status"] == "partial"
    assert result["diagnostics"]["bounded"]
    assert all(c["region"] != "full" for c in result["contigs"])


def test_bridge_reached_at_budget_limit_is_complete_but_marked_bounded():
    result = assemble_chain("AAGC", "TCAA", {"TCAAGCTT": 2},
                            AssemblyOptions(min_overlap=3, max_steps=1))
    assert result["status"] == "complete"
    assert result["contigs"][0]["sequence"] == "AAGCTTGA"
    assert result["diagnostics"]["bounded"]
    assert result["diagnostics"]["expanded_states"] == 1


def test_search_stops_before_exceeding_length_bound():
    path = extend_test_path("AACG", {"ACGT": 2, "CGTT": 2}, max_contig_length=5)
    assert path.sequence == "AACGT"
    assert path.termination == "max_contig_length"
    with pytest.raises(ValueError, match="exceeds max_contig_length"):
        assemble_chain("AACGTT", "", {}, AssemblyOptions(min_overlap=3, max_contig_length=5))


def test_repeat_stops_without_looping():
    path = extend_test_path("AACGACG", {"ACGA": 2})
    assert path.sequence == "AACGACG"
    assert path.steps == 0
    assert path.termination == "repeated_read"


def test_beam_rescues_seed_to_seed_path_lost_by_greedy_search():
    # CGTA has higher support but dead-ends. The lower-support CGTC branch
    # reaches CCTA through GTCC and TCCT.
    counts = {"CGTA": 3, "CGTC": 2, "GTCC": 2, "TCCT": 2, "CCTA": 2}
    greedy = assemble_chain("ACGT", revcom("CCTA"), counts,
                            AssemblyOptions(min_overlap=3))
    beam = assemble_chain("ACGT", revcom("CCTA"), counts,
                          AssemblyOptions(min_overlap=3, beam_width=2))
    assert greedy["status"] == "partial"
    assert beam["status"] == "complete"
    assert any("ACGTCCTA" in c["sequence"] for c in beam["contigs"])
    assert all(c["ambiguous"] for c in beam["contigs"])
    reordered = assemble_chain("ACGT", revcom("CCTA"), dict(reversed(list(counts.items()))),
                               AssemblyOptions(min_overlap=3, beam_width=2))
    assert beam == reordered


def test_evidence_scoring_uses_quality_instead_of_raw_duplicate_count():
    counts = {"ACGA": 3, "ACGT": 2}
    quality = {"ACGA": 0.3, "ACGT": 1.9}
    legacy = extend_test_path("AACG", counts, quality=quality)
    evidence = extend_test_path("AACG", counts, quality=quality, scoring="evidence")
    assert legacy.sequence == "AACGA"
    assert evidence.sequence == "AACGT"


def test_evidence_scoring_uses_mates_from_the_current_path():
    counts = {"ACGA": 3, "ACGT": 2}
    mate_links = {"AACG": {"ACGT": 3}}
    legacy = extend_test_path("AACG", counts, mates=mate_links)
    evidence = extend_test_path("AACG", counts, mates=mate_links, scoring="evidence")
    unrelated = extend_test_path("AACG", counts, mates={"CCCC": {"ACGT": 3}}, scoring="evidence")
    assert legacy.sequence == unrelated.sequence == "AACGA"
    assert evidence.sequence == "AACGT"


@pytest.mark.parametrize("weight", [-1, 3, float("nan"), float("inf")])
def test_invalid_quality_weights_are_rejected(weight):
    with pytest.raises(ValueError):
        assemble_chain("AACG", "", {"ACGT": 2}, AssemblyOptions(min_overlap=3),
                       quality_weights={"ACGT": weight})


def test_zero_quality_evidence_does_not_extend():
    path = extend_test_path("AACG", {"ACGT": 2}, quality={"ACGT": 0}, scoring="evidence")
    assert path.sequence == "AACG"


def test_invalid_mate_count_is_rejected():
    with pytest.raises(ValueError, match="Mate link"):
        assemble_chain("AACG", "", {"ACGT": 2}, AssemblyOptions(min_overlap=3),
                       mate_links={"AACG": {"ACGT": -1}})


def test_sequence_normalization_rejects_nonstring_input():
    with pytest.raises(TypeError, match="string"):
        normalize_sequence(None)


def test_verbose_legacy_extend_displays_requested_orientation(capsys):
    seed = "AACGTAGCTACGATCG"
    assert extend(seed, {}, True, 16, reverse_comp=True) == seed
    assert capsys.readouterr().out.strip() == revcom(seed)
    assert extend(seed, {}, True, 16) == seed
    assert capsys.readouterr().out.strip() == seed


def test_legacy_extend_obeys_default_length_bound():
    with pytest.raises(ValueError, match="exceeds max_contig_length"):
        extend("ACGT" * 8000, {}, False, 100)


def test_beam_deduplicates_successors_from_different_overlap_lengths():
    result = assemble_chain("AACG", "", {"ACGT": 2, "CGT": 2},
                            AssemblyOptions(min_overlap=2, beam_width=2))
    assert [contig["sequence"] for contig in result["contigs"]] == ["AACGT"]
    assert not result["contigs"][0]["ambiguous"]


def test_legacy_stitch_finds_constant_first_path():
    sequence = "AACGTAGCTACGATCGTTAGCGTACGATGCTAG"
    output = SimpleQueue()
    stitch(revcom(sequence), sequence[:16], {}, {}, False, len(sequence),
           "cell1", "heavy", output)
    assert output.get()["heavy"] == ">cell_id=cell1;heavy_chain;BASIC\n" + sequence + "\n"


def test_legacy_stitch_falls_back_to_variable_path():
    sequence = "AACGTAGCTACGATCGTTAGCGTACGATGCTAG"
    output = SimpleQueue()
    stitch(revcom(sequence[-16:]), sequence[:16], {}, {sequence: 2}, False,
           len(sequence), "cell1", "heavy", output)
    assert output.get()["heavy"] == ">cell_id=cell1;heavy_chain;BASIC\n" + sequence + "\n"


def test_legacy_stitch_returns_labeled_partial_or_empty_results():
    output = SimpleQueue()
    start, end = "AACGTAGCTACGATCG", "TTAGCGTACGATGCTA"
    stitch(revcom(end), start, {}, {}, False, 16, "cell1", "heavy", output)
    result = output.get()["heavy"]
    assert "[variable_region_contig]" in result and start in result
    assert "[constant_region_contig]" in result and end in result
    stitch("", "", {}, {}, False, 16, "cell1", "heavy", output)
    assert output.get() == {"heavy": ""}
