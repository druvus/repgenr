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
| sibeliaz | yes (converter) | yes (native + container) | native end-to-end → MAF → MSA → tree (closely-related genomes; 26 genus-level Francisella reps -> 168k-col MSA) and **container** verified on real divergent Francisella genomes (4 seqs -> 80-col MSA -> tree). Required a macOS fix: SibeliaZ's wrapper uses Linux-only `free`/`find -printf`/`stat -c`/`mktemp --suffix`; the adapter passes `-f` and runs a BSD-patched wrapper, and `maf_to_fasta` takes a seqid→genome name-map. **Memory note:** SibeliaZ's per-block alignment calls `spoa` with global (O(n^2)-memory) alignment and suppresses its stderr; on a single very large collinear block (e.g. the ~50 kb synthetic genomes, one block) spoa can be OOM-killed in a memory-limited container VM, leaving an empty MAF. The adapter now detects an empty MAF and raises a clear error (raise the container VM RAM or run natively); real genomes fragment into short LCBs and align fine |
| cactus | — | yes | container (`cactus:v2.9.3`, amd64); full Minigraph-Cactus run -> HAL -> MAF -> MSA -> tree. Verified on the synthetic set and live on 5 real **F. tularensis** strains (intraspecific is the right granularity for a pangenome; genus-level inputs are too divergent and align only to the reference). Fixes: writable `HOME` for Toil; per-call bind mounts for `seqfile.txt` genome paths; `_find_hal` selects the combined `*.full.hal` (not a per-chromosome HAL under `chrom-alignments/`, which omits genomes); and the `_MINIGRAPH_` backbone pseudo-genome is excluded in `maf_to_fasta` so it is not a taxon (genomes Minigraph-Cactus drops from the graph are tolerated, not errored). Plus Rosetta emulation: the bundled `vg` cannot run under Docker's QEMU (even `vg version` hangs), so on Apple Silicon Docker Desktop must use the **Apple Virtualization framework with "Use Rosetta for x86/amd64"** enabled. On native amd64/Linux (HPC, the Singularity target) no emulation is involved |

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
native vs wrapped) and **every tool was validated live in a container** on macOS
+ Docker + Wave/BioContainers (all 16 tools; SibeliaZ has a memory note for
single very large blocks, see its row above):
- Dereplicators (4/4): `skder` -> 8 reps; `drep` (amd64, `--virus` skips CheckM)
  -> 8 reps; `sourmash` -> 8 reps; `galah` -> 8 reps.
- SNP typers (4/4): `simple` (Wave multi-tool image minimap2+samtools+bcftools)
  -> 2413 sites; `snippy` (amd64) -> 2406 sites; `parsnp` -> 80 sites;
  `--mask gubbins` -> 2413 sites.
- Tree builders (5/5): `mashtree`, `sourmash`, `iqtree`, `fasttree`, `raxmlng`
  all -> trees.
- Aligners (3/3): `progressivemauve` (BioContainer); `cactus` (amd64 + Rosetta);
  `sibeliaz` on real genomes (the synthetic single-large-block case OOM-kills
  spoa in the VM -- see its row).

Three backend fixes came out of these sweeps:
- **Writable HOME.** Containers run as the host UID with no passwd entry, so HOME
  defaults to `/` and is not writable. The Docker wrapper now sets
  `-e HOME=<workdir>` (the mounted, writable working dir), which unblocks tools
  that touch HOME — e.g. Toil/Cactus creating its config dir.
- **Per-call extra mounts.** `run_tool(..., extra_mounts=[...])` lets an adapter
  declare input directories that are referenced indirectly (paths listed inside a
  manifest file rather than passed as argv tokens). Cactus uses this for the
  genome paths in `seqfile.txt`; the sourmash dereplicator/tree builder use it for
  the genomes listed in their `--from-file` fofn.
- **Un-resolved fofn paths.** `write_fofn` emits `os.path.abspath` (not
  `Path.resolve()`) paths so a fofn read inside a container matches the backend's
  un-resolved bind mounts (macOS firmlinks resolve `/Users` to a path outside
  Docker's shared dirs).

Singularity/Apptainer is Linux-only (no macOS build), so it can't run natively on
the macOS dev box. The exact command forms the backend emits were validated
against real **apptainer 1.4.4** in a Linux container: `apptainer pull <sif>
docker://<image>` (the `--container-cache` `.sif` behavior) and `apptainer exec
--bind <dir> --pwd <wd> <sif> <argv>`. On HPC/Linux this is the production engine.

Notes: macOS firmlinked temp/home paths must be bind-mounted un-resolved (handled
in the backend). Containers run as the host UID with no passwd entry, so the
backend sets a writable `HOME` (the mounted workdir); adapters whose inputs are
listed inside a manifest file declare those directories via `extra_mounts`. The
freshly Wave-built `mauve` image is broken (boost ABI mismatch); the adapter pins
a working BioContainer and a boost-pinned conda spec instead (see the aligners
table).

## Platform notes (macOS / Apple Silicon)

- Several tools lack osx-arm64 builds; some run via an osx-64 (Rosetta) conda env
  (e.g. parsnp/harvesttools). `mauve` is unpackaged on macOS entirely.
- `mashtree` pulls `perl-bio-samtools`, which pins an ancient samtools 0.1.x;
  modern samtools/bcftools for the `simple` SNP typer live in a separate env and
  are placed ahead on PATH.
- On exFAT/NTFS volumes, macOS `._*` AppleDouble files must be ignored (handled
  in the adapters/stages) and skDER must run on a local filesystem.
- amd64-only container tools run under emulation. Docker's QEMU emulation cannot
  run some SIMD-heavy binaries (e.g. Cactus's bundled `vg` hangs even on
  `vg version`). Set Docker Desktop to the **Apple Virtualization framework** with
  **"Use Rosetta for x86/amd64 emulation"** enabled; Rosetta runs these binaries
  correctly. On native amd64/Linux hosts no emulation is involved.
