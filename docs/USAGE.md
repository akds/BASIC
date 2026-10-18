# Extended usage

These examples use the installed `basic-assemble` command. From a checkout,
`python BASIC.py` accepts the same arguments. See the [README](../README.md)
for installation and a first run.

## Input files, lanes, and subsampling

Prepare one cell's reads per invocation. Input must be plain or gzipped,
four-line FASTQ with Phred+33 quality scores. Choose `-SE` or both `-PE_1`
and `-PE_2`. Paired records must have matching identifiers, with optional
`/1` and `/2` or standard Illumina mate labels.

Comma-separated files are processed in the order supplied. For paired data,
the first left file pairs with the first right file, and so on:

```bash
basic-assemble \
  -PE_1 cell_L001_R1.fastq.gz,cell_L002_R1.fastq.gz \
  -PE_2 cell_L001_R2.fastq.gz,cell_L002_R2.fastq.gz \
  -g human -i BCR -n cell_001 -o results --subsample 100000
```

`--subsample 100000` consumes at most the first 100,000 single-end records or
100,000 synchronized pairs **across all lanes combined**, before trimming and
filtering. It does not take 100,000 from each lane, and it is not random
sampling. If lane 1 supplies the entire prefix, later lanes contribute no
reads. A filtered read or pair still counts toward this limit, so fewer reads
may be retained. Only consumed records undergo FASTQ validation; input files
are also read for provenance checksums. Omit `--subsample` to process and
validate all records.

## Multiple cells

Save the following as `cells.tsv`. The separators in this example are actual
tabs. Paths are relative to the TSV file; absolute paths also work.

```tsv
sample_id	left	right	genome	receptor
b_cell_001	reads/b001_R1.fastq.gz	reads/b001_R2.fastq.gz	human	BCR
t_cell_001	reads/t001_R1.fastq.gz	reads/t001_R2.fastq.gz	mouse	TCR
```

For single-end data, use a `fastq` column instead of `left` and `right`:

```tsv
sample_id	fastq	genome	receptor
b_cell_002	reads/b002.fastq.gz	human	BCR
```

```bash
basic-assemble --batch cells.tsv -o results -p 4 --continue-on-error
```

Sample IDs must be unique even when case is ignored. Optional `genome` and
`receptor` columns override `-g` and `-i`; omit those columns to use command
defaults. Input fields can contain comma-separated lane lists. Each row must
specify either single-end or paired input, and batch mode cannot be combined with
command-line `-SE`/`-PE_1`/`-PE_2` inputs.

Cells run sequentially; `-p` controls each Bowtie2 alignment's threads.
`--continue-on-error` permits later cells to run after a sample failure, while
the final exit status remains nonzero if any sample failed. A uniquely named
`batch-*.json` records the processed samples and their outcomes. Structural
TSV errors, such as duplicate sample IDs or incorrect column counts, stop
the batch before any sample runs.

## Optional IgBLAST annotation

Install `igblastn` separately and obtain suitable germline databases following
the [NCBI installation guide](https://ncbi.github.io/igblast/cook/How-to-set-up.html).
BASIC does not download annotation databases. Its bundled Bowtie2 seed indexes
cannot serve as IgBLAST databases.

Supply nucleotide BLAST database **prefixes** for V and J, a D database for
VDJ loci, and a J auxiliary file matched to the sequences and identifiers in
your J database. If starting with germline FASTA, follow NCBI's identifier
preparation and `makeblastdb` instructions. Set `--igdata` to the directory
containing IgBLAST's `internal_data` folder when installation data cannot be
found automatically.

Replace the illustrative paths below with your installation and databases:

```bash
basic-assemble -PE_1 cell_R1.fastq.gz -PE_2 cell_R2.fastq.gz \
  -g human -i BCR -n cell_annotated -o results \
  --igblast /opt/igblast/bin/igblastn \
  --germline-v /data/igblast/human_ig_V \
  --germline-d /data/igblast/human_ig_D \
  --germline-j /data/igblast/human_ig_J \
  --auxiliary-data /opt/igblast/optional_file/human_gl.aux \
  --igdata /opt/igblast --annotation-timeout 300
```

For TCRs, select `-i TCR` and use **TCR germline databases**, not the IG databases
in the preceding example:

```bash
basic-assemble -SE t_cell.fastq.gz -g human -i TCR \
  -n t_cell_annotated -o results \
  --igblast /opt/igblast/bin/igblastn \
  --germline-v /data/igblast/human_tr_V \
  --germline-d /data/igblast/human_tr_D \
  --germline-j /data/igblast/human_tr_J \
  --auxiliary-data /data/igblast/human_tr.aux --igdata /opt/igblast
```

Add `--germline-c /data/igblast/human_ig_C` to request C-gene assignment when
your IgBLAST version supports it. C calls may inform isotype if the assembled
sequence contains sufficient constant-region evidence. See the
[official examples](https://ncbi.github.io/igblast/cook/examples.html).

The wrapper requests AIRR output and retains available gene calls,
junction/CDR3 annotations, and predicted productivity. Unavailable fields
stay empty; an empty `productive` field does not mean false. Junction includes
the conserved flanking codons, whereas CDR3 excludes them. Omitting or
mismatching the J auxiliary file can leave frame/CDR3 annotation unavailable
or incorrect. If annotation fails, the run stops without publishing a partial
annotation table.

Annotation settings are shared across batch rows. Use separate invocations
for batches requiring different species, receptor databases, or auxiliary
files. Discoverable annotation reference files receive checksums; unresolved
database/configuration paths are reported in the run manifest. Retain the
database versions and installation data used for reproducible annotation.

## Trimming and experimental assembly options

Terminal `N` bases are trimmed by default. `--trim-quality 20` also trims bases
below Phred 20 from both ends. The default minimum retained length is 16 bases
(`--min-overlap` plus one); reads with at least 50% of one A/C/G/T base are
filtered. These filters can reduce coverage, so inspect `read_statistics` in
the run manifest.

```bash
basic-assemble -SE cell.fastq.gz -g human -i BCR \
  -n cell_trimmed -o results --trim-quality 20
```

To try quality and mate evidence while limiting alternative paths:

```bash
basic-assemble -PE_1 cell_R1.fastq.gz -PE_2 cell_R2.fastq.gz \
  -g human -i BCR -n cell_evidence -o results -a \
  --scoring evidence --quality-weighted --mate-aware \
  --beam-width 3 --max-anchors 2 --max-candidates 2
```

`--quality-weighted` and `--mate-aware` require `--scoring evidence`;
`--mate-aware` also requires paired input. More paths/seeds can consume more
time and memory. These modes are experimental, and their scores are not
confidence probabilities. Compare them with defaults using independent
truth before changing an analysis workflow. See [algorithm details](ALGORITHMS.md)
and [benchmark tools](../benchmarks/README.md).

## Checks, outputs, and recovery

Use `--check` with the intended arguments to check input paths, indexes, and
dependencies without assembly. It does not validate every FASTQ record or
fully validate IgBLAST database contents. `--alignment-timeout 600` sets a
600-second limit for each regional Bowtie2 run; no mapping timeout is imposed
by default. `--annotation-timeout` defaults to 300 seconds.

By default, an incomplete reconstruction stops publication. Add `-a` to retain
partial or empty biological results; execution errors still fail. A
`complete` assembly connects its seeds but is not necessarily a complete or
productive receptor. Inspect contig status and ambiguity, not only whether
a FASTA file exists.

Existing sample outputs require `--overwrite`. Results are staged together,
and `<sample>.run.json` is published last as the completion marker. Automated
readers should verify the manifest's output checksums. Ordinary publication
errors trigger rollback to the previous bundle. Individual renames are
atomic, but the multi-file bundle is **not crash-atomic**: a power loss or
forced termination can interrupt publication. If rollback itself fails,
BASIC reports a `.basic-recovery-*` directory containing recoverable prior
files. Preserve that directory until recovery is complete.

Failures during processing write `<sample>.failed.json` with diagnostics.
`-t scratch --keep-temp` preserves the unique temporary run directory, whose
path is recorded in the success or failure manifest. Without `--keep-temp`,
temporary reads are removed; successful runs retain Bowtie2 logs in
`<sample>.logs/`.

A `.<sample>.basic.lock` prevents overlapping runs with the same output name.
After a forcibly terminated process it may remain. Inspect the PID recorded
inside and confirm that the original run is no longer active before removing
a stale lock and retrying. Use a new sample name for independent parameter
comparisons. For custom or related-species seeds, follow the
[non-model organism guide](NON_MODEL_ORGANISMS.md).
