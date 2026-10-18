# Custom references and optional annotation

Use this guide only when the selected workflow needs new references,
related-species seeds, or IgBLAST. Check the installed commands' `--help` before
using these BASIC 2.0 examples. No online documentation is required to use this
skill. The historical [BASIC repository](https://github.com/akds/BASIC) may serve
an older release; match documentation to the installed version.

## Custom species

Prefer curated references for the actual species. If references are missing,
`basic-references search --organism 'Danio rerio' --query immunoglobulin
--limit 10 --output candidates.json` searches public nucleotide metadata.
Inspect candidates before selecting them. Hits can be unrelated superfamily
genes, partial records, rearranged transcripts, or whole chromosomes; the
search does not identify a validated germline catalog. It sends the organism
and query to NCBI, not the sample reads.

Create a manifest using selected regional nucleotide FASTA sources. All four
groups are required. This BCR example uses illustrative local paths, relative
to the manifest:

```json
{
  "schema_version": 1,
  "regions": {
    "hv": [{"path": "igh_v.fa", "orientation": "forward", "locus": "IGH"}],
    "hc": [{"path": "igh_c.fa", "orientation": "forward", "locus": "IGH"}],
    "lv": [{"path": "igk_v.fa", "orientation": "forward", "locus": "IGK"}],
    "lc": [{"path": "igk_c.fa", "orientation": "forward", "locus": "IGK"}]
  }
}
```

Each source needs exactly one `path`, HTTPS `url`, or selected NCBI `accession`
(prefer accession.version), plus explicit orientation and locus. Source FASTA
must be uncompressed. `forward` means biological 5′-to-3′ sequence: the builder
reverse-complements C references into BASIC's required orientation. `as_is`
means already BASIC-oriented, not unknown orientation. BASIC does not extract
coordinates or splice exons from whole genomic records.

BCR uses IGH in `hv`/`hc`, IGK or IGL in `lv`/`lc`. TCR permits TRA, TRB, TRG,
or TRD; the `hv`/`hc` and `lv`/`lc` locus sets must match within each pair.
Declare one locus per source; split mixed-locus FASTA files. These declarations
are not independent biological validation and cannot represent every species'
receptor nomenclature.

```bash
basic-references build --manifest sources.json --output references \
  --species target_species --receptor BCR --bowtie-build /path/to/bowtie2-build
basic-assemble -SE cell.fastq.gz -g target_species -i BCR \
  -d references -n target_cell -o results
```

The builder checks all four Bowtie2 indexes and records sources, orientation,
loci, checksums, and versions in `reference_manifest.json`. Preserve it. Existing
libraries require explicit `--overwrite`; avoid replacing one while assemblies
use it. For a new comparison, prefer a separate reference root.

If the user elects exploratory related-species seeding, keep the biological
target in `-g target_species` and choose the library with
`--seed-species related_species -d references --local-seeding`. Report both
species. BASIC does not select the nearest species automatically. Local
seeding can admit less specific seeds; read-to-read overlaps remain exact.
Recovered expressed contigs do not define novel validated germline alleles.

## IgBLAST annotation

Annotation is optional and separate from assembly. Use local `igblastn` and
species/receptor-matched **nucleotide BLAST database prefixes** for V/J, and D
for VDJ loci. BASIC's Bowtie2 seed indexes are not IgBLAST databases. Obtain
the appropriate J auxiliary file and IgBLAST internal data. Use the
[NCBI setup guide](https://ncbi.github.io/igblast/cook/How-to-set-up.html) when
preparing databases rather than inventing identifiers or files.

Add these options to the assembly command, replacing paths with verified data:

```bash
--igblast /path/to/igblastn \
--germline-v /data/human_ig_V --germline-d /data/human_ig_D \
--germline-j /data/human_ig_J --auxiliary-data /data/human_gl.aux \
--igdata /path/to/igblast --annotation-timeout 300
```

These example databases are for human BCR. Use TR databases and matching
auxiliary data for TCR, and actual target-species data for other organisms.
`--igdata` names the parent of `internal_data`. Optional `--germline-c PREFIX`
requests C assignment when supported by the installed IgBLAST version.

Retain database versions and the run's annotation provenance. BASIC preserves
available gene calls, junction/CDR3, and predicted productivity; missing fields
remain unknown. Junction includes conserved flanking codons, whereas CDR3
excludes them. Missing/mismatched J auxiliary data can prevent or distort
frame/CDR3 calls. An annotation failure is a failed run, not permission to
invent calls. Homology-based annotation for an unsupported species is not
definitive germline assignment.
