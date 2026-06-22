# RepGenR

RepGenR (Representative-Genome Repositories) builds dereplicated, representative
genome repositories for large-scale genomic studies. It selects taxa, downloads
genomes, clusters them by average nucleotide identity (ANI), computes
phylogenetic trees, and emits taxonomy files for downstream tools such as
FlexTaxD.

Version 2 is a modular, importable Python package (Python 3.12+) with four
pluggable tool families and an optional Nextflow pipeline for scaling to
thousands of genomes.

## Pipeline

```
metadata -> genome -> dereplicate -> (align | snptype) -> phylo -> tree2tax
```

Each stage is a `repgenr` subcommand and an importable function. Stages
communicate through a single working directory whose state is recorded in
`repgenr.yaml` (provenance) and a SQLite genome manifest.

## Pluggable tools

Tools are discovered through Python entry points; adding one needs no change to
the core (see `docs/adding-tools.md`).

| Family | Built-in adapters |
|--------|-------------------|
| Dereplicators | `drep`, `skder`, `galah`, `sourmash` |
| Aligners | `progressivemauve`, `sibeliaz`, `cactus` |
| SNP typers | `simple` (samtools/bcftools), `snippy`, `parsnp` (optional Gubbins masking) |
| Tree builders | `iqtree`, `fasttree`, `raxmlng` (MSA), `mashtree`, `sourmash` (alignment-free) |

`repgenr list-tools` prints what is available in your environment.

## Installation

```bash
# Tools + Python environment (conda/mamba)
mamba env create -f environment.yml
mamba activate repgenr

# Or just the Python package (tools must be on PATH separately)
pip install .
```

Cactus is distributed separately (containers/binaries); see its documentation.

## Usage

```bash
WD=./francisella

# Full GTDB table (release-pinned):
repgenr metadata -wd $WD -r 232.0 --gtdb-version bac120 -d rep -l genus -tg francisella
# Or query just the target taxon via the GTDB API (no full-table download):
repgenr metadata -wd $WD --source api -d rep -l genus -tg francisella

repgenr genome -wd $WD
repgenr dereplicate -wd $WD --tool skder -t 16
repgenr phylo -wd $WD --aligner progressivemauve --treebuilder iqtree
repgenr tree2tax -wd $WD --include-dereplicated
```

Or run the whole chain in one command (bacterial by default; `--viral` for the
NCBI Virus path), then check progress at any time:

```bash
repgenr run -wd $WD -d rep -l genus -tg francisella --tool skder --treebuilder iqtree
repgenr status -wd $WD     # which stages are done, and what to run next
```

### Resume and `--force`

Each stage records its parameters in `repgenr.yaml`; re-running a stage that
already completed with the same parameters is a safe no-op (it logs that it
skipped). Change a parameter, or pass `--force`, to re-run it. A stage that
crashed mid-run has no completion stamp and so always re-runs. If you re-run an
upstream stage (e.g. `dereplicate --force`) and then a downstream stage whose
parameters are unchanged, the downstream skip is flagged as potentially stale —
pass `--force` there too to rebuild against the new inputs.

Alternatives:

```bash
# Scalable dereplication then an alignment-free tree
repgenr dereplicate -wd $WD --tool skder
repgenr phylo -wd $WD --treebuilder mashtree

# SNP-based phylogeny (core-SNP alignment as the MSA source)
repgenr snptype -wd $WD --tool simple
repgenr phylo -wd $WD --msa-source snptype --treebuilder iqtree --mask gubbins
```

### Viruses

The viral path selects from NCBI Virus by default (via the `datasets` CLI);
`vmetadata --source bvbrc` uses the legacy BV-BRC FTP path instead.

```bash
WD=./hav
repgenr vmetadata -wd $WD --target hepatovirus            # NCBI Virus (default)
repgenr vgenome   -wd $WD --target-genus Hepatovirus      # add --group-segments for segmented viruses
repgenr dereplicate -wd $WD --tool skder --virus
repgenr phylo -wd $WD --treebuilder mashtree
repgenr tree2tax -wd $WD --include-dereplicated
# or: repgenr run -wd $WD --viral -t hepatovirus -tg Hepatovirus --treebuilder mashtree
```

## Troubleshooting

- **`MissingBinaryError` / a tool is not found.** The Python package does not
  install the bioinformatics tools. Use the conda environment
  (`mamba env create -f environment.yml`) or put the tool on `PATH`. Run
  `repgenr list-tools` to see the adapters and `--container docker` (or
  `singularity`) to run tools in pinned images instead.
- **Apple Silicon / arm64.** BioContainers are amd64; pass
  `--platform linux/amd64` (and enable Rosetta) so emulated images run.
- **A stage failed; where are the details?** Errors print a concise message; the
  full traceback is in `<workdir>/repgenr.log`. Re-run with `--verbose` to see it
  on the console. `repgenr status -wd <WD>` shows what completed and what is next.
- **GTDB download fails.** Check `--release` (e.g. `232.0`) and `--gtdb-version`
  (`bac120`/`ar53`); transient HTTP errors are retried automatically. The
  `--source api` mode fetches only the target taxon (no full-table download).
- **NCBI Entrez throttling (viral BV-BRC path).** Set `NCBI_API_KEY` (and
  optionally `NCBI_EMAIL`) to raise the request-rate limit.
- **A tool hangs.** Set `REPGENR_SUBPROCESS_TIMEOUT=<seconds>` to cap every
  external tool; on expiry the process group is killed with a clear error.

## Nextflow

The Nextflow layer runs the pipeline as typed data channels (no shared working
directory); results are published under `--outdir`.

```bash
nextflow run nextflow/main.nf -profile standard --outdir results \
    --metadata_args "-r 232.0 --gtdb-version bac120 -d rep -l genus -tg francisella" \
    --derep_tool sourmash --phylo_args "--treebuilder mashtree"
```

`--mode viral` runs the viral path instead (NCBI Virus by default; BV-BRC is
available via `--vmetadata_args "--source bvbrc"`). Profiles: `standard` (local),
`slurm`, `cloud`, `test` (add a container profile such as `singularity` to run
the tools in pinned images). Resource labels (`process_low/medium/high`) are
tuned per profile; heavy aligners such as Cactus use `process_high`. Set
`--derep_process_size` to scatter dereplication across tasks for large inputs.

## Development

```bash
pip install -e ".[dev]"
ruff check src/ tests/
mypy src/repgenr
pytest -q
```

See `docs/architecture.md` for the design and `docs/adding-tools.md` for writing
a new adapter.

## License

MIT.
