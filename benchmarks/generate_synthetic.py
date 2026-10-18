#!/usr/bin/env python3
"""Generate explicitly synthetic V + random junction + C assembly test reads.

These constructs do not model V(D)J recombination, productivity, expression,
library biases, or real error distributions. They are software test fixtures,
not empirical evidence of biological accuracy. Bundled constant sequences are
reverse-complemented to match the assembler's configured reference orientation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import random
import re


REGIONS = ("hv", "hc", "lv", "lc")
DEFAULT_REFERENCE_INDICES = {region: 0 for region in REGIONS}


def _revcom(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _references(path):
    records, header, parts = [], None, []
    with Path(path).open(encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(parts).upper()))
                header, parts = line[1:], []
            else:
                if header is None:
                    raise ValueError(f"Sequence before FASTA header: {path}")
                parts.append(line)
    if header is not None:
        records.append((header, "".join(parts).upper()))
    return records


def _mutate(sequence, probability, rng):
    bases, qualities = [], []
    for base in sequence:
        changed = rng.random() < probability
        bases.append(rng.choice([other for other in "ACGT" if other != base]) if changed else base)
        qualities.append("+" if changed else "I")
    return "".join(bases), "".join(qualities)


def generate_cell(
    output_dir,
    *,
    database_path=None,
    species="human",
    receptor="BCR",
    name=None,
    read_lengths=(75,),
    step=20,
    copies=3,
    substitution_rate=0.0,
    noise_reads=0,
    seed=1,
    paired=False,
    junction_length=30,
    candidates_per_chain=1,
    reference_indices=None,
):
    """Write a deterministic per-cell fixture and return its manifest.

    References use zero-based FASTA record indexes (default: record 0 in each
    region); every chosen header, index, and source checksum is recorded.
    ``step`` controls tiled start positions and ``copies`` repeats each tile.
    Paired fragments are twice the read length with mate 2 reverse-complemented.
    ``noise_reads`` counts reads in single-end mode and pairs in paired mode.
    Multiple candidates share their V/C references but have different random
    synthetic junctions; their truth order need not equal assembler rank order.
    """
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", species):
        raise ValueError("species must be a directory basename")
    if receptor not in ("BCR", "TCR"):
        raise ValueError("receptor must be BCR or TCR")
    name = name or f"synthetic_{species}_{receptor}"
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name):
        raise ValueError("name must be a safe sample basename")
    read_lengths = tuple(read_lengths)
    if not read_lengths or any(not isinstance(value, int) or value < 16 for value in read_lengths):
        raise ValueError("read_lengths must contain integers of at least 16")
    for key, value in (("step", step), ("copies", copies),
                       ("junction_length", junction_length), ("candidates_per_chain", candidates_per_chain)):
        if not isinstance(value, int) or value < 1:
            raise ValueError(f"{key} must be a positive integer")
    if not isinstance(noise_reads, int) or noise_reads < 0:
        raise ValueError("noise_reads must be a nonnegative integer")
    if not math.isfinite(substitution_rate) or not 0 <= substitution_rate <= 1:
        raise ValueError("substitution_rate must be between 0 and 1")
    selected_indices = dict(DEFAULT_REFERENCE_INDICES)
    selected_indices.update(reference_indices or {})
    if set(selected_indices) != set(REGIONS):
        raise ValueError("reference_indices keys must be hv, hc, lv, lc")
    database = Path(database_path or Path(__file__).resolve().parents[1] / "db")
    references = {}
    provenance = {}
    for region in REGIONS:
        path = database / species / receptor / region
        index = selected_indices[region]
        records = _references(path)
        if not isinstance(index, int) or index < 0 or index >= len(records):
            raise ValueError(f"Invalid zero-based reference index {index!r} for {path}")
        header, raw_sequence = records[index]
        sequence = raw_sequence.strip("N")
        if not sequence or set(sequence) - set("ACGT"):
            raise ValueError(f"Reference {region}[{index}] is empty or contains internal ambiguous bases")
        match = re.search(r"(?:^|[|:_])(IGH|IGK|IGL|TRA|TRB|TRG|TRD)", header)
        if not match:
            raise ValueError(f"Cannot derive locus from reference header: {header}")
        references[region] = sequence
        provenance[region] = {"record_index": index, "header": header, "locus": match.group(1),
                              "source_file": str(path.resolve()),
                              "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                              "terminal_ns_removed": len(raw_sequence) - len(sequence),
                              "construct_orientation": "forward" if region.endswith("v") else "reverse_complement"}
    rng = random.Random(seed)
    constructs = []
    for position, (variable, constant) in enumerate((("hv", "hc"), ("lv", "lc"))):
        locus = provenance[variable]["locus"]
        if provenance[constant]["locus"] != locus:
            raise ValueError(f"Selected {variable}/{constant} references have incompatible loci")
        chain = ("heavy", "light")[position] if receptor == "BCR" else locus
        for candidate in range(1, candidates_per_chain + 1):
            junction = "".join(rng.choice("ACGT") for _ in range(junction_length))
            sequence = references[variable] + junction + _revcom(references[constant])
            if max(read_lengths) * (2 if paired else 1) > len(sequence):
                raise ValueError("Read/fragment length exceeds synthetic construct length")
            identifier = f"cell_id={name};{chain}_chain;BASIC"
            if candidate > 1:
                identifier += f";candidate={candidate}"
            constructs.append({"sequence_id": identifier, "sequence": sequence, "chain": chain,
                               "locus": locus, "candidate": candidate,
                               "variable_region": variable, "constant_region": constant,
                               "synthetic_junction": junction,
                               "junction_start": len(references[variable]),
                               "junction_end": len(references[variable]) + len(junction)})
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    read_names = ["R1.fastq", "R2.fastq"] if paired else ["reads.fastq"]
    files = ["truth.fasta", "manifest.json", "samples.tsv", *read_names]
    if any((output / filename).exists() for filename in files):
        raise FileExistsError("Synthetic fixture outputs already exist; use a fresh output directory")
    physical_count = 0
    with (output / "truth.fasta").open("w", encoding="ascii") as truth:
        for construct in constructs:
            truth.write(f">{construct['sequence_id']}\n{construct['sequence']}\n")
    from contextlib import ExitStack
    with ExitStack() as stack:
        handles = [stack.enter_context((output / filename).open("w", encoding="ascii")) for filename in read_names]
        number = 0
        for construct in constructs:
            template = construct["sequence"]
            for length in read_lengths:
                fragment_length = length * (2 if paired else 1)
                last_start = len(template) - fragment_length
                starts = sorted(set(range(0, last_start + 1, step)) | {last_start})
                for start in starts:
                    for copy in range(copies):
                        number += 1
                        fragments = [template[start:start + length]]
                        if paired:
                            fragments.append(_revcom(template[start + fragment_length - length:start + fragment_length]))
                        elif rng.random() < 0.5:
                            fragments[0] = _revcom(fragments[0])
                        for mate, fragment in enumerate(fragments, 1):
                            sequence, quality = _mutate(fragment, substitution_rate, rng)
                            identifier = f"synthetic_{number}" + (f"/{mate}" if paired else "")
                            handles[mate - 1].write(f"@{identifier}\n{sequence}\n+\n{quality}\n")
                            physical_count += 1
        for noise in range(noise_reads):
            number += 1
            for mate, handle in enumerate(handles, 1):
                length = rng.choice(read_lengths)
                sequence = "".join(rng.choice("ACGT") for _ in range(length))
                identifier = f"synthetic_noise_{number}" + (f"/{mate}" if paired else "")
                handle.write(f"@{identifier}\n{sequence}\n+\n{'I' * length}\n")
                physical_count += 1
    if paired:
        table = f"sample_id\tleft\tright\tgenome\treceptor\n{name}\tR1.fastq\tR2.fastq\t{species}\t{receptor}\n"
    else:
        table = f"sample_id\tfastq\tgenome\treceptor\n{name}\treads.fastq\t{species}\t{receptor}\n"
    (output / "samples.tsv").write_text(table, encoding="ascii")
    manifest = {"dataset_kind": "synthetic_software_validation", "sample_id": name,
                "limitations": "V + random junction + reverse-complemented bundled C. No V(D)J, productivity, expression, or empirical accuracy model.",
                "species": species, "receptor": receptor, "references": provenance,
                "options": {"read_lengths": list(read_lengths), "step": step, "copies": copies,
                            "substitution_rate": substitution_rate, "noise_reads": noise_reads,
                            "noise_unit": "pairs" if paired else "reads", "seed": seed,
                            "paired": paired, "junction_length": junction_length,
                            "candidates_per_chain": candidates_per_chain},
                "physical_reads": physical_count, "constructs": constructs,
                "files": {"truth": "truth.fasta", "reads": read_names, "batch": "samples.tsv"}}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--db", dest="database_path")
    parser.add_argument("--species", default="human")
    parser.add_argument("--receptor", choices=("BCR", "TCR"), default="BCR")
    parser.add_argument("--name")
    parser.add_argument("--read-lengths", default="75", help="Comma-separated lengths, each >=16")
    parser.add_argument("--step", type=int, default=20)
    parser.add_argument("--copies", type=int, default=3)
    parser.add_argument("--substitution-rate", type=float, default=0.0)
    parser.add_argument("--noise-reads", type=int, default=0)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--paired", action="store_true")
    parser.add_argument("--junction-length", type=int, default=30)
    parser.add_argument("--candidates-per-chain", type=int, default=1)
    parser.add_argument("--reference-index", action="append", default=[], metavar="REGION=INDEX")
    values = vars(parser.parse_args(argv))
    try:
        values["read_lengths"] = tuple(int(value) for value in values["read_lengths"].split(","))
        indexes = {}
        for selection in values.pop("reference_index"):
            region, index = selection.split("=", 1)
            indexes[region] = int(index)
        values["reference_indices"] = indexes
        result = generate_cell(**values)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"synthetic fixture: {exc}\n")
    print(json.dumps({"dataset_kind": result["dataset_kind"], "physical_reads": result["physical_reads"],
                      "sample_id": result["sample_id"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
