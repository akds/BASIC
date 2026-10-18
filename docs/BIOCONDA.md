# Preparing a Bioconda update

This repository contains a **staged recipe generator**, not a published BASIC
2.0.1 conda package. It preserves the existing conda package name `basic`.
The Python distribution is named `basic-immune-assembler`.

## Existing package and upstream comparison

The official Bioconda recipe currently packages `basic` **1.5.1, build 1** as
`noarch: python`, with runtime dependencies on Python and Bowtie2. Its build
script copies the references to `$PREFIX/db/basic` and patches the old
`BASIC.py` reference path. This update instead installs the Python package and
its bundled reference resources. See the official [recipe](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/basic/meta.yaml),
[build script](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/basic/build.sh),
and [package history](https://bioconda.github.io/recipes/basic/README.html).

Upstream's latest release is [v1.5.1](https://github.com/akds/BASIC/releases/tag/v1.5.1),
published July 12, 2019. Compared with the starting 1.5.0 snapshot, the sole
change captures Bowtie2's version output instead of printing it; there are no
algorithm or reference changes in that comparison. The 1.5.1 README still
reports 1.5.0. See the [upstream comparison](https://github.com/akds/BASIC/compare/v1.5.0...v1.5.1).

The old recipe labels the license `MIT`, but the [upstream README](https://github.com/akds/BASIC/blob/v1.5.1/README.md)
only says “Software provided to academic users under MIT License.” That
historical release omitted the complete license text. BASIC 2 supplies the
standard [MIT license](../LICENSE), confirmed by its maintainer, and the new
recipe uses `license: MIT` with `license_file: LICENSE`. Bundled reference
provenance and third-party tool terms remain separate from the software
license. Bioconda requires a [license allowing redistribution](https://bioconda.github.io/contributor/guidelines.html#bioconda-recipe-checklist).

## Build a recipe from the actual sdist

Build the final source distribution using the release instructions, then pass
that exact file to the generator. For the current version:

```bash
python packaging/bioconda/prepare.py \
  --sdist dist/basic_immune_assembler-2.0.1.tar.gz \
  --output /tmp/basic-conda-recipe
```

The output directory must not already exist. The generator reads the package
name and version from `PKG-INFO`, verifies required source files, hashes the
actual archive, and writes `meta.yaml`, `run_test.py`, and `preparation.json`.
The default source is a local `file:` URI, suitable for local testing. There
is no invented release URL or placeholder checksum. Regenerate after every
sdist change; the checksum binds the recipe to exact bytes.

The recipe requires Python >=3.10 and Bowtie2, installs `basic-assemble`,
`BASIC.py`, `basic-references`, and `basic-skill`, and packages human/mouse
BCR/TCR references and the AI-agent skill.
IgBLAST and its germline databases remain separate, optional installations.
It builds with pip without dependency downloads or build isolation; conda
supplies setuptools, wheel, and pip. The [conda metadata documentation](https://docs.conda.io/projects/conda-build/en/stable/resources/define-metadata.html)
describes source checksums, entry points, and build dependencies.

With conda-build available in a dedicated environment, test the local recipe:

```bash
conda build /tmp/basic-conda-recipe --override-channels \
  -c conda-forge -c bioconda
```

The recipe checks installed entry points and skill content, all 16 packaged reference indexes,
their discovery from outside a checkout, and a small assembly using the
installed Bowtie2 executable. It does not require pytest, IgBLAST, downloads,
or manually supplied biological data. For the complete Bioconda lint/build
workflow, follow [Testing Recipes Locally](https://bioconda.github.io/contributor/building-locally.html).

## Stage the public update after release

After the final sdist is actually hosted, rerun the generator with
`--source-url` set to its **exact existing HTTPS URL** and a new output
directory. The generator downloads it and rejects any size or SHA256
difference from the local archive. The default license metadata is `MIT`
and `LICENSE`; the generator requires that file to exist inside the sdist.

Review the generated recipe and tests, remove any obsolete old `build.sh`
when preparing the update to Bioconda's existing `recipes/basic/` directory,
and run the Bioconda checks. Publication and opening a pull request are
separate maintainer actions. Version 2.0.1 is not available through
`conda install basic` merely because these preparation files exist.

## Local validation

The generator's tests require only Python:

```bash
python -m unittest discover -s packaging/bioconda -p 'test_*.py'
```

No conda, mamba, or micromamba executable was available in the preparation
environment. A successful generator or installed-package smoke test is not
a completed conda-build or Bioconda CI run.
