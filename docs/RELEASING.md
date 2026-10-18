# Release BASIC 2.0.1

The canonical repository is [akds/BASIC](https://github.com/akds/BASIC).
The README presents **2.0 (2026/10/18)**; Python and Conda use **2.0.1**, and the
Git tag is **v2.0.1**. The software uses the standard [MIT license](../LICENSE).
Reference data and external tools retain their terms.

PyPI publication can be deferred independently of the GitHub release.
Automatic PyPI uploads are disabled unless the repository Actions variable
`PYPI_PUBLISH_ENABLED` is set to `true`. GitHub release assets still publish.

Publication has three steps: merge the source update on GitHub, publish Python
distributions to PyPI, and update the existing Bioconda `basic` recipe.
Preparing these files does not publish any registry package.

## Preserve GitHub history

The starting folder is a source snapshot without `.git`. Apply the prepared
source to a branch in a real clone, preserving existing history:

```bash
git clone https://github.com/akds/BASIC.git BASIC-release
cd BASIC-release
git switch -c basic-2.0
rsync -av \
  --exclude='.git/' --exclude='.venv/' --exclude='dist/' --exclude='build/' \
  --exclude='release-assets/' \
  --exclude='*.egg-info/' --exclude='__pycache__/' --exclude='.pytest_cache/' \
  --exclude='.coverage*' --exclude='htmlcov/' --exclude='/results/' \
  --exclude='.agents/' --exclude='.codex/' --exclude='.aws/' --exclude='.cache/' \
  /path/to/updated/BASIC-master/ ./
git rm --ignore-unmatch .travis.yml
git diff --check
git diff --stat
```

Review changes against upstream before committing. The public v1.5.1 change
only adjusts the Bowtie2 version probe relative to v1.5.0; the new checked
subprocess wrapper supersedes that probe. Do not force-push or replace history.

## Run release checks

Install a working Bowtie2, then run:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-release.lock
python -m pip install --no-build-isolation -e .
BASIC_BOWTIE2=/absolute/path/to/bowtie2 python -m pytest
python -m unittest discover -s packaging/bioconda -p 'test_*.py'
python -m build --no-isolation
python -m twine check --strict dist/*.whl dist/*.tar.gz
python packaging/release_assets.py --dist dist --output release-assets
```

Use a clean `dist/` when building; keep only the intended version's wheel and
source archive. Install the wheel in another fresh environment outside the
checkout. Check all four commands (`basic-assemble`, `BASIC.py`,
`basic-references`, `basic-skill`), all four reference index sets, and a real
assembly. See [validation results and limits](VALIDATION.md).
Native IgBLAST biological validation remains pending.

Stage the reviewed source, documentation, tests, skills, packaging helpers,
license, dependency files, and workflows. Exclude environments, credentials,
results, and build outputs. Commit and push `basic-2.0`, then open a pull request
against `master`. GitHub Actions must pass on the submitted commit.

## Configure PyPI once

The Python distribution name is **basic-immune-assembler**. A registry lookup
returned no project during preparation; this does not reserve the name or
guarantee that PyPI will permit it. The existing Conda package remains `basic`.

In PyPI, register a [pending Trusted Publisher](https://docs.pypi.org/trusted-publishers/creating-a-project-through-oidc/)
for a new project, or add a publisher to the existing project if you own it:

| Field | Value |
| --- | --- |
| PyPI project | `basic-immune-assembler` |
| GitHub owner | `akds` |
| GitHub repository | `BASIC` |
| Workflow filename | `publish.yml` |
| Environment | `pypi` |

Create the `pypi` environment in the GitHub repository. Configure TestPyPI
separately with the same owner/repository/workflow and environment `testpypi`;
create that GitHub environment too. Apply the repository's desired reviewer
and tag restrictions. The workflow uses short-lived OIDC credentials rather
than a stored API token, as described by the
[PyPI publishing guide](https://docs.pypi.org/trusted-publishers/using-a-publisher/).

After configuring and testing the publisher, set the repository Actions
variable `PYPI_PUBLISH_ENABLED` to `true` if future published releases should
upload automatically. To publish this existing release later, manually run
the workflow with destination `pypi` and ref `v2.0.1`; that explicit run does
not require the variable. Until then, leave the variable unset or `false`.

## Publish the checked version

1. Merge the reviewed branch. In Actions, run **Build and publish Python
   distributions** with destination **build-only**. It tests with real Bowtie2,
   builds wheel/source archives, checks metadata, and tests wheel installation.
2. Run with destination **testpypi**, then verify installation in a fresh
   environment using only the TestPyPI index:

   ```bash
   python -m pip install --index-url https://test.pypi.org/simple/ basic-immune-assembler==2.0.1
   basic-assemble --version
   ```

3. Tag the exact reviewed commit `v2.0.1`, push the tag, and publish its GitHub
   release. The `release: published` event runs checks and uploads GitHub assets.
   It also publishes to PyPI when `PYPI_PUBLISH_ENABLED` is `true`.
   A manual **pypi** run is also supported, but it must select the matching
   version tag. Production rejects a branch or mismatched tag.
4. Verify `python -m pip install basic-immune-assembler==2.0.1` in a new
   environment and confirm the PyPI file checksums. The workflow also attaches
   the exact tested wheel, source archive, standalone skill ZIP, and checksum
   list to the GitHub release. It refuses to overwrite existing assets.
5. Generate and submit the [Bioconda recipe update](BIOCONDA.md) using the
   actual published source URL and SHA256. Verify installation from Bioconda
   after its maintainers merge and publish the recipe.

A PyPI release version's uploaded filenames cannot be replaced. If an uploaded
release needs changed code or artifacts, choose a new version and matching tag;
do not mask discrepancies with `skip-existing`. TestPyPI has the same constraint,
so repeated pre-release experiments should use distinct pre-release versions.

The prepared workflow does not open a Bioconda pull request or publish Conda
packages automatically. That separate build supplies Bowtie2 as an installation
dependency and updates the Anaconda.org channel listing.
