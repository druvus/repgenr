# Installing RepGenR

RepGenR is a Python package that calls external bioinformatics tools. The
package does not install those tools. This page covers both parts and three
ways to provide the tools.

## What gets installed

- **The package.** `pip install .` installs the `repgenr` command. It needs
  Python 3.12 or later. Two extras exist: `sparse` adds the
  `sourmash_plugin_branchwater` back-end for sparse sourmash dereplication
  (without it the adapter falls back to dense `sourmash compare`), and `dev`
  adds the test, lint and type-check tools.
- **The tools.** Dereplicators, aligners, SNP typers, the masker, tree
  builders, assemblers and polishers are separate programs. A stage that needs
  one that is missing stops with exit code 4 and names it, except under `assemble --assembler auto`, which excuses runs whose assembler is missing (reason `assembler_not_installed`, with a warning) and exits 4 only when no run can be assembled. A missing `auto` polisher is skipped (see [usage.md](usage.md#starting-from-sequencing-reads)). See
  [the per-tool table](#per-tool-table) for what each adapter needs.

```bash
pip install .            # the package only
pip install ".[sparse]"  # with the sparse sourmash back-end
pip install -e ".[dev]"  # for development
```

## Which method to use

| Situation | Recommended method | Why |
|---|---|---|
| Linux x86_64 workstation | One conda environment from `environment.yml` | The single-environment route. `environment.yml` covers every tool except Cactus and the databases; it is not expected to solve on macOS (see section 2). |
| HPC cluster | `--container singularity` with `--container-cache` on shared storage | Tools run from pinned images without site installs. Images are pulled once and reused. |
| macOS on Apple Silicon | Per-tool conda environments on `PATH` for the core tools; containers for progressiveMauve, Cactus, snippy, dRep, skesa, shovill, flye, medaka, racon and CheckM2 | Except racon and CheckM2, which have offline tests only, these were run inside containers on the audit machine (`verification.md`). skder, galah, sourmash, SibeliaZ, `simple`, parsnp, ska2, Gubbins and the tree builders ran natively. |
| Nextflow on a cluster | A site image or conda per profile; `-profile slurm,singularity` | The `slurm` profile sets only the executor. The container profiles set `--container` for every stage. |
| Nextflow on a cloud executor | A site config with the executor, queue and an image that provides `repgenr` and the tools | No cloud profile ships, because the region, queue and image are site-specific (see [usage.md](usage.md#profiles)). |

## Three ways to provide the tools

### 1. One conda environment

`environment.yml` lists the package and every tool except Cactus and the
databases. It is the single-environment route and is not expected to solve on macOS (see section 2). progressiveMauve from bioconda needs the adapter's `boost-cpp=1.74.0` pin, or the pinned image, so check that tool first with `repgenr list-tools --check`.

```bash
mamba env create -f environment.yml
mamba activate repgenr
```

The versions in `environment.yml` are lower bounds. They match the minimum
versions RepGenR checks before a run and guard against known-incompatible old
tools. They do not make an environment reproducible. For a reproducible
environment, create a pinned lock file from a solved environment, one per
platform, for example:

```bash
conda env export --no-builds > environment.lock.yml
```

(`conda-lock` serves the same purpose.)

### 2. Several conda environments and PATH

On macOS the tools are split over several environments. The known reasons are:

- `mashtree` depends on `perl-bio-samtools`, which pins samtools 0.1.x. The
  `simple` SNP typer needs samtools and bcftools 1.10 or later, so the two
  live in different environments.
- On exFAT or NTFS volumes, macOS writes `._*` AppleDouble files. The adapters
  ignore them, but skDER and dRep must run on a local filesystem, so keep their
  workdir off an exFAT disk.
- Some tools have no osx-arm64 build. parsnp and harvesttools run from an
  osx-64 (Rosetta) environment.
- progressiveMauve is not packaged for macOS at all. Use a container (below)
  or a Linux host.

Create one environment per group and put the environments' `bin` directories
on `PATH`. The first match wins, so list the environment that holds the newer
samtools and bcftools before the one that holds mashtree:

```bash
export PATH=$HOME/miniforge3/envs/repgenr_snp/bin:$HOME/miniforge3/envs/repgenr_core/bin:$PATH
export PATH=$HOME/miniforge3/envs/repgenr_parsnp_x64/bin:$PATH   # osx-64 env for parsnp

repgenr list-tools --check
```

The environment names above are examples. A single command can also be run
inside a named environment without activating it:

```bash
conda run -n repgenr_core --no-capture-output repgenr phylo -wd $WD --treebuilder mashtree
```

The live test suite uses the same idea: its config has a `[bin_dirs]` table
that puts tools from other environments on `PATH` (see `tests/live/README.md`).

### 3. Containers

RepGenR can run each tool inside a pinned container instead of on `PATH`. The
package itself runs on the host and only the tool subprocess is containerized,
so this works in the plain CLI and inside Nextflow. The default is `none`
(native execution). The quick start and the option table are in
[usage.md](usage.md#running-tools-in-containers); the rest is here.

```bash
repgenr --container docker dereplicate -wd $WD --tool drep
repgenr --container singularity --container-cache /data/repgenr_sif phylo -wd $WD
```

The options are global and go before the subcommand. Each has an environment
variable: `REPGENR_CONTAINER`, `REPGENR_CONTAINER_ENGINE`,
`REPGENR_CONTAINER_CACHE`, `REPGENR_CONTAINER_PLATFORM` and `REPGENR_WAVE`.

- `--container docker|singularity` selects the backend.
  `--container-engine` replaces the binary (apptainer, podman).
- `--container-cache <dir>` is where Singularity images and the Wave cache
  live. Put it on a large disk.
- `--platform linux/amd64` runs amd64 images on an arm64 host. Docker Desktop
  on Apple Silicon needs the Apple Virtualization framework with "Use Rosetta
  for x86/amd64 emulation" enabled. Docker's QEMU emulation cannot run some
  SIMD-heavy binaries; the bundled `vg` in Cactus hangs even on `vg version`.
  On a native amd64 Linux host no emulation is involved.
- `--wave` mints an image from the adapter's conda spec with the Seqera Wave
  CLI. It is the only image route for the four adapters that carry no pinned
  image: `simple`, parsnp, skder and SibeliaZ. Without `--wave` those four run
  on the host and the log says so.

#### Image sources

Each adapter declares its container metadata in `ToolCapabilities`:

- `container` is a pinned image URI, used as is when Wave is off. The
  single-package adapters pin a BioContainer (galah, sourmash, dRep, snippy,
  ska2, Gubbins, IQ-TREE, FastTree, RAxML-NG, mashtree, and the assemblers and
  polishers). progressiveMauve pins a BioContainer with a compatible boost, and
  Cactus pins its project image. Four adapters have no pin and run in an image
  only under `--wave`: the `simple` SNP typer (minimap2, samtools, bcftools)
  and parsnp (parsnp, harvesttools) span several packages, and the skder and
  SibeliaZ BioContainers are BusyBox-based. The GNU-only calls in their shell
  wrappers (`sort --parallel` in skder's greedy mode, `mktemp --suffix` in
  SibeliaZ) fail there and the run silently yields nothing.
- `conda` is a conda spec (for example `bioconda::skder`). With `--wave`,
  RepGenR mints an image for it. This is arm64-native and the only route for
  the multi-package adapters. The pin stays the default without Wave.

BioContainers are `linux/amd64`. On Apple Silicon pass
`--platform linux/amd64` with Rosetta, or use `--wave`.

#### Storage location

- **Singularity or Apptainer.** `--container-cache` sets `APPTAINER_CACHEDIR`
  and `SINGULARITY_CACHEDIR` (and the matching `*_TMPDIR`). `docker://` images
  are pulled once to `<cache>/<name>.sif` and reused. Put this on a large or
  external disk.
- **Docker.** The daemon manages image storage (Docker Desktop's disk image).
  It is set there, not per run.
- Large run-time data (the Cactus job store, scratch, downloads) lives under
  `--workdir` and `TMPDIR`.

#### Notes

- Before a tool runs in an image, RepGenR checks that the engine binary is
  on `PATH` and, for Docker, that `docker info` reaches the daemon. A stopped
  Docker Desktop therefore exits 4 before the stage starts, and
  `repgenr --container docker list-tools --check` reports it per tool. The
  image itself is pulled on first use.
- Docker runs as the host UID and GID so outputs are owned by you. The workdir
  and `TMPDIR` are bind-mounted at identical paths.
- Symlinked inputs (a `genomes/` directory staged by `repgenr ingest`) are
  followed. The directory each link points to is bound as well, so the tool
  sees the same paths inside the container.
- On Apple Silicon most BioContainers are `linux/amd64`. They run under Docker
  emulation, or use `--wave` for arm64-native images.
- The macOS SibeliaZ BSD-wrapper workaround is skipped automatically when
  running in a (Linux) container.
- dRep's CheckM step needs its reference data at run time. Set
  `CHECKM_DATA_PATH` on the host: RepGenR passes it into the container and
  binds that directory at the same path. Otherwise run dRep with
  `--ignoreGenomeQuality`.
- The bioconda `mauve` (progressiveMauve) build is broken upstream (boost ABI,
  `undefined symbol`). Use the pinned image, or run that tool natively on
  Linux.

## Per-tool table

"Native on macOS arm64" records what was run on the audit machine (Apple
Silicon) and is taken from the adapter status table in
[verification.md](verification.md#adapter-status). "No" means not run natively
there, not that the tool cannot work on another platform. The container pin is
the image in the adapter's declared capabilities; `none` means the adapter
carries no pin and needs `--wave` to run in an image. Minimum versions are in
`environment.yml`.

| Family | Tool | conda package | Native on macOS arm64 | Container pin | Databases or extras |
|---|---|---|---|---|---|
| Dereplicator | drep | `drep`, `checkm-genome` | no | `drep:3.7.1` | CheckM reference data through `CHECKM_DATA_PATH`, or `--ignoreGenomeQuality` |
| Dereplicator | skder | `skder`, `skani` | yes | none (Wave) | none |
| Dereplicator | galah | `galah` | yes | `galah:0.4.2` | none |
| Dereplicator | sourmash | `sourmash` | yes | `sourmash:4.9.4` | `sourmash_plugin_branchwater` for the sparse back-end (extra `sparse`) |
| Aligner | progressivemauve | `progressivemauve` | no (unpackaged on macOS) | `mauve:2.4.0.snapshot_2015_02_13` | none |
| Aligner | sibeliaz | `sibeliaz` | yes | none (Wave) | none |
| Aligner | cactus | not on conda | no | `cactus:v2.9.3` | none |
| SNP typer | simple | `minimap2`, `samtools`, `bcftools` | yes | none (Wave) | samtools and bcftools 1.10 or later |
| SNP typer | snippy | `snippy` | no | `snippy:4.6.0` | none |
| SNP typer | parsnp | `parsnp`, `harvesttools` | yes, from an osx-64 environment | none (Wave) | none |
| SNP typer | ska2 | `ska2` | yes | `ska2:0.5.1` | none |
| Masker | gubbins | `gubbins` | yes | `gubbins:3.4.3` | a multi-threaded RAxML build for more than one thread (see [usage.md](usage.md#snp-typing-and-masking)) |
| Tree builder | iqtree | `iqtree` | yes | `iqtree:3.1.3` | none |
| Tree builder | fasttree | `fasttree` | yes | `fasttree:2.2.0` | none |
| Tree builder | raxmlng | `raxml-ng` | yes | `raxml-ng:2.0.3` | none |
| Tree builder | mashtree | `mashtree` | yes | `mashtree:1.4.6` | none |
| Tree builder | sourmash | `sourmash` | yes | `sourmash:4.9.4` | none |
| Assembler | skesa | `skesa` | no | `skesa:2.5.1` | none |
| Assembler | shovill | `shovill` | no | `shovill:1.4.2` | none |
| Assembler | flye | `flye` | no | `flye:2.9.6` | none |
| Classifier | sourmash | `sourmash` | yes | `sourmash:4.9.4` | a GTDB sourmash sketch and its lineages CSV |
| Polisher | medaka | `medaka` | no | `medaka:2.2.2` | medaka models. Reads from SRA have no basecaller in their headers; the adapter then assumes ONT's bacterial R10.4.1 model and records it. Give another with `--tool-arg model=...` |
| Polisher | racon | `racon`, `minimap2` | no | `racon:1.5.0` | none |
| Quality (reads chain) | checkm2 | `checkm2` | no | `checkm2:1.1.0` | CheckM2 database: `checkm2 database --download`, then `--checkm2-db` or `CHECKM2DB` |

The genome download steps also need `ncbi-datasets-cli`.

## Cactus

Cactus has no conda package in `environment.yml` and no native install path
in these docs. The adapter pins the project image
(`quay.io/comparative-genomics-toolkit/cactus:v2.9.3`), so use
`--container docker|singularity`. On Apple Silicon this needs the Rosetta
setting described above.

## Nextflow

The Nextflow layer needs Nextflow 26.04 or later (`nextflowVersion =
'!>=26.04.0'` in `nextflow/nextflow.config`). The nf-test suite uses nf-test
0.9.5, the version CI installs. The `docker`, `singularity` and `wave`
profiles set `params.repgenr_opts`, so every stage calls `repgenr` with
`--container docker`, `--container singularity` or `--container singularity
--wave`. The `repgenr` command itself must be available to each process. See
[usage.md](usage.md#nextflow) for parameters and profiles.

## Verify the installation

```bash
repgenr --version
repgenr list-tools --check
```

`list-tools` prints each family with its adapters and, where one is declared,
the genome limit. Its last line names the dereplicators that `glance` can run
(`drep`, `sourmash`). `--check` runs each adapter's preflight and prints `ok` with
the versions found, or `missing` or `error` with the reason. It reports and
exits with status 0. A stage run with a missing or outdated tool exits with
status 4, except under `assemble --assembler auto`, which excuses runs whose assembler is missing (reason `assembler_not_installed`, with a warning) and exits 4 only when no run can be assembled. A missing `auto` polisher is skipped (see [usage.md](usage.md#starting-from-sequencing-reads)). The full exit-code table is in
[usage.md](usage.md#troubleshooting).
