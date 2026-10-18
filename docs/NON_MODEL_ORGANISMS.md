# References and assembly for non-model organisms

BASIC can search public nucleotide records, download selected reference
sequences, and build a custom reference directory. It can also use another
species' references as experimental seeds for assembling the target animal's
reads. These are different operations: neither creates a validated germline
catalog automatically.

## Find candidate sequences by organism

After installing BASIC, search NCBI nucleotide records:

```bash
basic-references search \
  --organism 'Danio rerio' --query 'immunoglobulin' --limit 10 \
  --output candidates.json
```

The command combines an organism field restriction with your query, uses
ESearch to retrieve record identifiers, then ESummary to retrieve their
metadata. Both interfaces are documented by
[NCBI's E-utilities reference](https://www.ncbi.nlm.nih.gov/books/NBK25499/).
The result includes accession/version when available, title, nucleotide
length, taxonomic ID, record URL, the query, total match count, and response
checksums. `--limit` is between 1 and 100. Searches are bounded by a default
30-second timeout and a 5 MB response limit; `--timeout` and `--max-bytes`
adjust these limits.

Search hits are **unvalidated candidates**. A broad query can return an
immunoglobulin-superfamily gene unrelated to an antibody, a whole genomic
region, an expressed rearrangement, or a partial record. BASIC does not assign
V/C regions, receptor loci, or germline status from titles. Inspect the linked
records and select appropriate regional nucleotide sequences before adding
them to a build manifest. An empty search does not establish biological
absence. Search results and annotations can change as NCBI is updated.

## Build from files, URLs, or selected accessions

A manifest declares sources for all four reference groups. For BCR assembly,
`hv` and `hc` contain heavy variable and constant references; `lv` and `lc`
contain the corresponding light references. Each source must declare its
locus and orientation. For example, this manifest imports four local FASTA
files relative to the manifest's directory:

```json
{
  "schema_version": 1,
  "metadata": {
    "scientific_name": "Oryctolagus cuniculus",
    "notes": "References selected and checked by the analyst"
  },
  "download": {
    "timeout": 30,
    "max_bytes": 50000000
  },
  "regions": {
    "hv": [{"path": "rabbit_igh_v.fa", "orientation": "forward", "locus": "IGH"}],
    "hc": [{"path": "rabbit_igh_c.fa", "orientation": "forward", "locus": "IGH"}],
    "lv": [{"path": "rabbit_igk_v.fa", "orientation": "forward", "locus": "IGK"}],
    "lc": [{"path": "rabbit_igk_c.fa", "orientation": "forward", "locus": "IGK"}]
  }
}
```

The example's paths are placeholders for analyst-selected sequences. Save it
as `rabbit_sources.json`, supply those files, then run:

```bash
basic-references build --manifest rabbit_sources.json \
  --output ./references --species rabbit --receptor BCR \
  --bowtie-build /path/to/bowtie2-build
```

The result is `./references/rabbit/BCR/`, containing the four oriented FASTA
files, all four Bowtie2 indexes, build logs, and `reference_manifest.json`.
Use that directory root for assembly:

```bash
basic-assemble -SE cell.fastq.gz -g rabbit -i BCR \
  -d ./references -n rabbit_cell
```

To avoid downloading or uploading files manually, replace a source's `path`
with either an HTTPS `url` or a selected NCBI `accession`:

```json
{
  "accession": "SELECTED_ACCESSION.VERSION",
  "orientation": "forward",
  "locus": "IGH",
  "region": "V"
}
```

Replace `SELECTED_ACCESSION.VERSION` with an actual selected nucleotide
accession, preferably including its version. Do not copy an arbitrary search
hit into this declaration merely because its title contains a receptor term.
The builder retrieves the complete record as FASTA using NCBI EFetch with
`db=nuccore`, `rettype=fasta`, and `retmode=text`.
[NCBI documents accession.version input and EFetch output formats](https://www.ncbi.nlm.nih.gov/books/NBK25499/).
There is no automatic coordinate extraction or exon splicing: a whole
chromosome accession is not a substitute for a regional reference.

An HTTPS source has the same fields, with `url` replacing `accession`:

```json
{
  "url": "https://your-reference-host.example/selected_igh_v.fa",
  "orientation": "forward",
  "locus": "IGH",
  "sha256": "REPLACE_WITH_64_HEXADECIMAL_CHARACTERS"
}
```

The URL and checksum above are placeholders. `sha256` is optional; when
provided, it must be the SHA256 of the exact source bytes and must match.
The builder always records the actual source checksum, even without a pinned
value. HTTPS is required for URLs and redirects. The downloader limits each
source to 50 MB and 30 seconds by default, checks declared response lengths,
and does not extract archives. Inputs must be uncompressed ASCII FASTA.
Requests to NCBI are spaced below three per second within this process,
consistent with [NCBI's guidance for requests without an API key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/).

## Orientation and locus declarations

`orientation="forward"` means each input sequence is written in its
biological 5′-to-3′ direction. Variable references remain in that direction.
**Constant references are reverse-complemented** to match BASIC's historical
`hc`/`lc` reference orientation. `orientation="as_is"` means the input is
already in BASIC's reference orientation; no reverse complementation occurs.
Use `as_is` only for previously oriented references, not as a way to avoid
checking strand direction.

The containing group explicitly declares whether a source is V or C. An
optional `region` value must be `V` inside `hv`/`lv`, or `C` inside `hc`/`lc`.
This is a declaration by the analyst; BASIC does not infer functional region
boundaries from the downloaded sequence.

Supported declarations are:

| Receptor | Reference groups | Allowed loci |
| --- | --- | --- |
| BCR | `hv`, `hc` | `IGH` |
| BCR | `lv`, `lc` | `IGK`, `IGL` |
| TCR | Any of the four groups | `TRA`, `TRB`, `TRG`, `TRD` |

Within each channel pair, the locus sets must match between `hv` and `hc`,
and between `lv` and `lc`. A light-chain library can contain separate IGK and
IGL sources if both are supplied for V and C. TCR groups are not hardcoded to
a human alpha/beta layout: a matched V/C pair can contain, for example, TRB
and TRD sources, with TRA and TRG in the other pair. Actual reported locus
depends on the mapped seeds.

Each source has one explicit locus. Split a mixed-locus FASTA into curated
sources rather than assigning one locus to the whole file. FASTA identifiers
are prefixed with that locus, such as `IGH|original_identifier`, so downstream
seed mapping retains the declaration. Duplicate identifiers within a source
or merged group, conflicting locus-bearing identifiers, empty records, and
non-IUPAC sequence characters are rejected. These checks cannot establish
the biological correctness of an analyst's locus assignment. The fixed locus
vocabulary does not represent every non-model organism's receptor nomenclature.

Optional `metadata` accepts `scientific_name`, a positive integer `taxon_id`,
and `notes`. These are recorded declarations, not independently verified
taxonomy. The `--species` argument is a local lowercase directory slug using
letters, numbers, underscores, and hyphens.

## Installation, failure handling, and provenance

All source retrieval and validation happens in a sibling staging directory.
The builder then runs Bowtie2-build separately for all four groups and checks
that each has six readable, nonempty `.bt2` or `.bt2l` shards. The default
build timeout is 600 seconds per group; change it with `--build-timeout`.
Timed-out or interrupted builds terminate their process groups on POSIX.

An existing reference directory is preserved unless `--overwrite` is given.
Source errors, checksum failures, builder failures, missing shards, and
timeouts leave an existing directory intact. Installation occurs only after
all four indexes pass. For overwrite, the previous directory is temporarily
renamed to a rollback backup, and the staged directory is renamed into place;
a failed installation restores the previous directory. This uses atomic
individual renames but has a brief directory-swap window, so do not run
assemblies against the target while replacing it. A lock prevents another
BASIC reference builder from replacing the same target concurrently.

`reference_manifest.json` records source URLs, accession IDs or resolved local
paths, source and oriented-FASTA checksums, locus and orientation declarations,
record counts, BASIC version, and Bowtie2-build version. Keep this file with
the library. Build stdout reports the installed path; if an old backup could
not be removed after successful installation, it also reports that backup.

## Experimental assembly using related-species seeds

When suitable target-species references are missing, choose an available
related-species reference library and explicitly request its use:

```bash
basic-assemble -SE target_cell.fastq.gz -g target_species -i BCR \
  --seed-species related_species --local-seeding \
  -d ./references -n target_cell
```

`--seed-species` chooses the reference directory; it does not replace the
sample's biological species. Select the related species yourself; BASIC does
not estimate phylogenetic distance or automatically choose the closest taxon.
`--local-seeding` allows local read-to-reference alignments during seed
selection. It is useful to try when a target's reads diverge from the supplied
reference, but it can also admit less specific seeds. Read-to-read assembly
still uses exact overlaps; local seeding does not make that assembly
indel-tolerant.

The output is a set of candidate **expressed receptor contigs from the input
reads**, not newly validated germline alleles. Receptor rearrangement,
somatic mutation, sequencing errors, paralogous transcripts, and incomplete
coverage complicate interpretation. Seed-to-seed completion alone does not
resolve these issues. Inspect ambiguity/partial-status diagnostics and use
independent sequence or annotation evidence appropriate to the organism.

Automatic taxon-wide germline discovery is not supported: public records may
lack reliable locus/segment annotation, and expressed rearrangements do not
by themselves define a species' germline repertoire. The supported workflow
is organism-name candidate search, explicit selection and declarations,
automatic retrieval/index building, and optional related-reference seeding.
