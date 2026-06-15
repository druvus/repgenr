# Running tools in containers

RepGenR can run each external tool inside a pinned container instead of relying
on tools installed on `PATH`. This unblocks Linux-only tools on any host (e.g.
macOS) and makes tool versions reproducible. RepGenR itself runs on the host;
only the tool subprocess is containerized, so it works in the plain CLI and
inside Nextflow.

## Quick start

```bash
# Run every tool in a container via Docker:
repgenr --container docker dereplicate -wd $WD --tool drep

# Singularity/Apptainer (HPC), with .sif images on a large disk:
repgenr --container singularity --container-cache /Volumes/LaCie/repgenr_sif \
        phylo -wd $WD --aligner progressivemauve --treebuilder iqtree
```

`--container` is a top-level option (place it before the subcommand). Default is
`none` (native execution, unchanged).

## Options (top-level, or env var)

| Option | Env var | Meaning |
|--------|---------|---------|
| `--container {none,docker,singularity}` | `REPGENR_CONTAINER` | execution backend |
| `--container-engine <bin>` | `REPGENR_CONTAINER_ENGINE` | engine override (apptainer, podman) |
| `--container-cache <dir>` | `REPGENR_CONTAINER_CACHE` | Singularity `.sif` / Wave cache (can be external) |
| `--platform <plat>` | `REPGENR_CONTAINER_PLATFORM` | e.g. `linux/amd64` to emulate BioContainers on Apple Silicon |
| `--wave / --no-wave` | `REPGENR_WAVE` | resolve multi-tool/arm64 images via the Seqera Wave CLI |

## Image sources

Each adapter declares its container metadata in `ToolCapabilities`:
- `container` — a pinned image URI (e.g. the Cactus image), used as-is.
- `conda` — a conda spec (e.g. `bioconda::skder`). With `--wave`, RepGenR mints an
  image for it via the Wave CLI (arm64-native, and handles multi-tool adapters
  such as the `simple` SNP typer = minimap2+samtools+bcftools); without Wave, set
  an explicit `container` to use a specific BioContainer.

## Storage location

- **Singularity/Apptainer:** `--container-cache` sets `APPTAINER_CACHEDIR` /
  `SINGULARITY_CACHEDIR` (+ `*_TMPDIR`); `docker://` images are pulled once to
  `<cache>/<name>.sif` and reused. Put this on a large/external disk.
- **Docker:** image storage is managed by the Docker daemon (Docker Desktop's
  disk-image location) and is set there, not per-run.
- Large run-time data (Cactus jobstore, scratch, downloads) lives under
  `--workdir` / `TMPDIR`.

## Notes

- Docker runs as the host UID/GID so outputs are owned by you; the workdir and
  `TMPDIR` are bind-mounted at identical paths.
- On Apple Silicon most BioContainers are `linux/amd64` (run via Docker
  emulation, or use `--wave` for arm64-native images).
- The macOS SibeliaZ BSD-wrapper workaround is skipped automatically when running
  in a (Linux) container.
- dRep's CheckM needs its reference DB at run time — mount it via
  `CHECKM_DATA_PATH`, or run dRep with `--ignoreGenomeQuality`.
- The bioconda `mauve` (progressiveMauve) build is broken upstream (boost ABI
  `undefined symbol`); pin a known-good image or run that tool natively on Linux.

## Nextflow

`nextflow run nextflow/main.nf -profile docker` (or `singularity`, `wave`) sets
`params.repgenr_opts` so every stage calls `repgenr --container …`. On HPC,
Singularity is the clean target (RepGenR dispatches each tool to Singularity);
stacking with Nextflow's own Docker engine implies docker-in-docker.
