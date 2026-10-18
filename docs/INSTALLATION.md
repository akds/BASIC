# Install or upgrade BASIC 2.0

BASIC needs Python 3.10+ and the **Bowtie2 executable**. IgBLAST is optional
and requires its own germline databases. The wheel contains BASIC, human/mouse
reference indexes, and the agent skill; it does not contain Bowtie2 or IgBLAST.

## Install 2.0 now

From the source checkout:

```bash
conda env create -f environment.yml
conda activate basic-2
python -m pip install .
basic-assemble --version
bowtie2 --version
```

The environment file installs Python and Bowtie2; pip installs BASIC.
Alternatively, install Bowtie2 separately, then install the source or wheel
in a Python environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install https://github.com/akds/BASIC/releases/download/v2.0.0/basic_immune_assembler-2.0.0-py3-none-any.whl
basic-assemble --version
```

The command downloads the [published GitHub wheel](https://github.com/akds/BASIC/releases/tag/v2.0.0);
you can also pass the path to a previously downloaded wheel. If Bowtie2 is
outside `PATH`, supply `-b /path/to/bowtie2` when assembling.
Reference construction similarly needs `bowtie2-build`.

## After the public registries publish 2.0.0

These are the intended registry commands. PyPI publication is deferred, and
the Bioconda update requires its own review and build. Use the GitHub wheel or
source until your chosen registry lists 2.0.0; do not assume these commands
provide it before then.

| Installer | Install or upgrade command | Source |
| --- | --- | --- |
| pip | `python -m pip install --upgrade basic-immune-assembler==2.0.0` | PyPI |
| pipx, first install | `pipx install basic-immune-assembler==2.0.0` | PyPI, isolated application environment |
| pipx, later updates | `pipx upgrade basic-immune-assembler` | PyPI, latest compatible release |
| uv | `uv tool install --upgrade basic-immune-assembler==2.0.0` | PyPI, isolated application environment |
| Conda | `conda create -n basic-2 -c conda-forge -c bioconda --strict-channel-priority basic=2.0.0` | Bioconda, including Bowtie2 dependency |
| Mamba | `mamba create -n basic-2 -c conda-forge -c bioconda --strict-channel-priority basic=2.0.0` | Same Bioconda package |

For an existing Conda environment, activate it and use
`conda install -c conda-forge -c bioconda --strict-channel-priority basic=2.0.0`.
Let the solver report whether the environment's dependencies permit the upgrade.
A new environment avoids mixing BASIC 1.x and 2.x installations.

pip, pipx, and uv use the **same PyPI publication**. Python itself is not a
separate registry. Conda/Mamba use a separate recipe and build published by
Bioconda; its [Anaconda.org listing](https://anaconda.org/bioconda/basic) updates
after that process completes. See [release instructions](RELEASING.md) and the
[Bioconda guide](BIOCONDA.md).

## Platforms and verification

Use Linux or macOS with a compatible Bowtie2 build. Native Windows is not a
validated target; use WSL2 with the Linux route. A platform-neutral Python wheel
does not guarantee external aligners are available on every architecture.
Conda dependency availability also varies by platform.

After upgrading, verify the executable actually being used:

```bash
basic-assemble --version
bowtie2 --version
basic-assemble -SE cell.fastq.gz -g human -i BCR -n check -o results --check
```

The first command should print `2.0.0`. `--check` validates setup and paths;
it does not assemble reads or parse the complete FASTQ. Use a representative
cell to validate a real workflow. BASIC 2 has stricter input checks and different
output handling; read the [migration notes](../CHANGELOG.md#migration).

To add the bundled skill to an existing project, follow the
[Codex, Claude, and Antigravity guide](AGENT_SKILL.md).
