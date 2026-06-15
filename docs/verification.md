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
| skder | yes | yes | run in a local scratch dir (exFAT-safe); membership from skani edges |
| sourmash | yes | yes | k-mer compare + greedy clustering |
| galah | yes | yes | cluster-definition parsed to representatives + members |
| drep | — | no | needs CheckM + its reference DB (heavy); not installed |

## Aligners

| Tool | Unit | Live | Notes |
|------|------|------|-------|
| progressivemauve | yes (converter) | no | `mauve`/`mauvealigner` not packaged for macOS on bioconda; `xmfa_to_fasta` unit-tested |
| sibeliaz | yes (converter) | yes | end-to-end on closely-related genomes → MAF → MSA → tree. Required a macOS fix: SibeliaZ's wrapper uses Linux-only `free`/`find -printf`/`stat -c`/`mktemp --suffix`; the adapter passes `-f` and runs a BSD-patched wrapper, and `maf_to_fasta` takes a seqid→genome name-map |
| cactus | — | no | Minigraph-Cactus is heavy (Toil); distributed separately |

## SNP typers

| Tool | Unit | Live | Notes |
|------|------|------|-------|
| simple | yes | yes | minimap2 + samtools/bcftools 1.23 → core-SNP alignment + SNP distance matrix |
| parsnp | — | yes | parsnp 2.1.5 + harvesttools → core-SNP FASTA |
| snippy | — | no | heavy dependency stack; not installed |
| ksnp | — | no | not installed |
| gubbins (mask) | — | yes | `--mask gubbins`; runs in its own Python 3.10 env, converged and filtered |

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
| vmetadata / vgenome (viral) | no | network-bound (BV-BRC FTP + NCBI Entrez); offline-parsing unit-tested |

## Platform notes (macOS / Apple Silicon)

- Several tools lack osx-arm64 builds; some run via an osx-64 (Rosetta) conda env
  (e.g. parsnp/harvesttools). `mauve` is unpackaged on macOS entirely.
- `mashtree` pulls `perl-bio-samtools`, which pins an ancient samtools 0.1.x;
  modern samtools/bcftools for the `simple` SNP typer live in a separate env and
  are placed ahead on PATH.
- On exFAT/NTFS volumes, macOS `._*` AppleDouble files must be ignored (handled
  in the adapters/stages) and skDER must run on a local filesystem.
