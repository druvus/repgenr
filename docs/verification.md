# Adapter verification status

Status of each pluggable tool adapter against its real external binary. "Unit"
means covered by offline tests (mock/parse); "Live" means run end-to-end against
the installed tool on real or closely-related data.

Live verification used a small GTDB/NCBI **Francisella** set (real download) and
a closely-related **synthetic** set (for tools that need within-species
divergence). Tools were installed via conda/mamba on macOS (Apple Silicon).

## Dereplicators

| Tool | Unit | Live | Notes |
|------|------|------|-------|
| skder | yes | yes | native + container (Wave); local scratch dir (exFAT-safe); membership from skani edges |
| sourmash | yes | yes | native; k-mer compare + greedy clustering |
| galah | yes | yes | native; cluster-definition parsed to representatives + members |
| drep | — | yes | container (Wave, amd64); a full run needs CheckM + its DB, so verified with `--virus` (sets `--ignoreGenomeQuality`, skipping CheckM) -> 8 reps |

## Aligners

| Tool | Unit | Live | Notes |
|------|------|------|-------|
| progressivemauve | yes (converter) | yes | container; full XMFA -> MSA -> tree on the synthetic set. progressiveMauve (libMems) needs boost-cpp 1.74; a naive `bioconda::mauve` Wave build solves against current conda-forge boost and fails at runtime (`undefined symbol _ZNK5boost...path8filenameEv`). **Both image paths are fixed and verified:** the adapter pins a **BioContainer** (`mauve:2.4.0.snapshot_2015_02_13--hdfd78af_4`, ships boost-cpp 1.74) as the default, and the `conda` spec pins `conda-forge::boost-cpp=1.74.0` so the **Wave** build also works. `container` wins over `conda`, so the BioContainer is used unless that pin is removed |
| sibeliaz | yes (converter) | yes | native end-to-end on closely-related genomes → MAF → MSA → tree. Required a macOS fix: SibeliaZ's wrapper uses Linux-only `free`/`find -printf`/`stat -c`/`mktemp --suffix`; the adapter passes `-f` and runs a BSD-patched wrapper, and `maf_to_fasta` takes a seqid→genome name-map |
| cactus | — | partial | container (`cactus:v2.9.3`, amd64). The backend now runs the full Minigraph-Cactus pipeline after two fixes: a writable `HOME` for Toil's config dir, and per-call bind mounts for genome paths listed inside `seqfile.txt`. On Apple Silicon it then fails inside `vg` with a QEMU amd64-emulation segfault; expected to run on a native amd64/Linux host |

## SNP typers

| Tool | Unit | Live | Notes |
|------|------|------|-------|
| simple | yes | yes | native (samtools/bcftools 1.23 from a dedicated env ahead on PATH) + container (Wave multi-tool image); minimap2 + samtools/bcftools → core-SNP alignment + distance matrix |
| parsnp | — | yes | native (parsnp 2.1.5 + harvesttools env) → core-SNP FASTA |
| snippy | — | yes | container (Wave, amd64); per-genome calling + snippy-core → 2406 core SNP sites |
| gubbins (mask) | — | yes | native; `--mask gubbins` in its own Python 3.10 env, converged and filtered |

## Tree builders

| Tool | Unit | Live | Notes |
|------|------|------|-------|
| iqtree | — | yes | ML tree from the SNP MSA |
| fasttree | — | yes | approximate-ML tree from the SNP MSA |
| raxmlng | — | yes | `--threads auto{N}` + `--redo` (fixes oversubscription / re-run) |
| mashtree | — | yes | alignment-free; validated on the real Francisella set |
| sourmash | yes | yes | alignment-free; k-mer distance + neighbor-joining |

## Front-end stages

| Stage | Live | Notes |
|-------|------|-------|
| metadata (GTDB) | yes | downloaded + parsed r207 bac120 (62,291 rep accessions) |
| genome (NCBI datasets) | yes | downloaded 8 Francisella genomes + outgroup |
| tree2tax | yes | FlexTaxD relations + genome map |
| vmetadata / vgenome (viral) | yes | live Hepeviridae run: BV-BRC FTPS download + NCBI Entrez -> 1256 genomes -> skder 799 reps -> mashtree -> tree2tax. Required an FTPS/TLS-session-reuse fix (BV-BRC dropped plain FTP) |

## Containers

RepGenR can run any tool in a pinned container (`--container docker|singularity`;
see `docs/containers.md`), pinning versions and unblocking tools that don't
install on the host. The backend has unit tests (argv construction, mounts, UID,
native vs wrapped) and was validated **live on macOS + Docker + Wave** across
every tool family:
- `dereplicate --tool skder` in a Wave-built single-tool image → 8 reps.
- `dereplicate --tool drep` in a Wave image (amd64) → 8 reps (`--virus` skips CheckM).
- `snptype --tool simple` in a Wave-built **multi-tool** image (minimap2 +
  samtools + bcftools) → 2413 core SNP sites.
- `snptype --tool snippy` in a Wave image (amd64) → 2406 core SNP sites.
- `phylo --treebuilder mashtree` in a Wave image → tree with all leaves.

Two backend fixes came out of this sweep:
- **Writable HOME.** Containers run as the host UID with no passwd entry, so HOME
  defaults to `/` and is not writable. The Docker wrapper now sets
  `-e HOME=<workdir>` (the mounted, writable working dir), which unblocks tools
  that touch HOME — e.g. Toil/Cactus creating its config dir.
- **Per-call extra mounts.** `run_tool(..., extra_mounts=[...])` lets an adapter
  declare input directories that are referenced indirectly (paths listed inside a
  manifest file rather than passed as argv tokens). Cactus uses this for the
  genome paths in `seqfile.txt`.

Singularity/Apptainer is Linux-only (no macOS build), so it can't run natively on
the macOS dev box. The exact command forms the backend emits were validated
against real **apptainer 1.4.4** in a Linux container: `apptainer pull <sif>
docker://<image>` (the `--container-cache` `.sif` behavior) and `apptainer exec
--bind <dir> --pwd <wd> <sif> <argv>`. On HPC/Linux this is the production engine.

Notes: macOS firmlinked temp/home paths must be bind-mounted un-resolved (handled
in the backend). The bioconda `mauve` (progressiveMauve) image is broken upstream
(boost ABI `undefined symbol`), so that tool stays unvalidated regardless of
container; the container mechanism itself ran it.

## Platform notes (macOS / Apple Silicon)

- Several tools lack osx-arm64 builds; some run via an osx-64 (Rosetta) conda env
  (e.g. parsnp/harvesttools). `mauve` is unpackaged on macOS entirely.
- `mashtree` pulls `perl-bio-samtools`, which pins an ancient samtools 0.1.x;
  modern samtools/bcftools for the `simple` SNP typer live in a separate env and
  are placed ahead on PATH.
- On exFAT/NTFS volumes, macOS `._*` AppleDouble files must be ignored (handled
  in the adapters/stages) and skDER must run on a local filesystem.
