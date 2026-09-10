# Live verification suite

Tests under this directory run the installed `repgenr` console script as a
real process against real tools. They are marked `live` and are deselected by
default (`addopts = "-m 'not live'"` in `pyproject.toml`), so neither plain
`pytest` nor CI collects them.

```bash
cp tests/live/live.example.toml tests/live/live.local.toml   # edit paths
conda run -n repgenr_dev --no-capture-output \
    pytest tests/live -m live --live-config tests/live/live.local.toml -ra \
    --basetemp /private/tmp/repgenr-live-tmp
```

Pass a fixed `--basetemp`: pytest clears it at every start, so only one
run's scratch exists at a time (the species-set and container tests leave a
few hundred MB per test, and the default temp directory keeps the last three
runs). Keep it on the boot disk: external volumes with non-native
filesystems lost files under samtools/bcftools, and Docker cannot bind them.

Subsets: `-m "live and not network and not container"` (offline, synthetic
data), `-m "live and network"` (GTDB/NCBI, cached after the first run),
`-m "live and container"` (Docker). A test whose tool is absent from PATH is
skipped by its `requires_binary` marker, so a partial environment still gives
a meaningful report.

## Network tests and the cache

`-m "live and network"` reaches GTDB (API and the release table, about 290
MB), NCBI Datasets, NCBI Virus and BV-BRC. Reference sets are built once
under `cache_dir` from the live config (default `~/.cache/repgenr-live`) and
reused on later runs; every test works on a copy under pytest's tmp_path.
Delete a subdirectory of the cache to rebuild that set. The reference sets
are genus Francisella (GTDB representatives), Francisella tularensis limited
to ten genomes, the r232 table, hepatovirus from NCBI Virus and Hepatitis E
from BV-BRC. First run about ten minutes, later runs about five.

## Tool environments for the species set

`tests/live/test_species_set.py` needs ska (`ska`), parsnp (`parsnp`,
`harvesttools`), gubbins (`run_gubbins.py`), minimap2, samtools, bcftools,
IQ-TREE, FastTree and RAxML-NG. Tools that live in other conda environments
are put on PATH through the `[bin_dirs]` table of the live config; a test
whose tool is still missing is skipped, not failed.

## Containers

`tests/live/test_container_runs.py` (`-m "live and container"`) runs adapters
through `--container docker --platform linux/amd64`. The module skips itself
unless `docker info` succeeds, `docker run --platform linux/amd64 alpine
uname -m` prints `x86_64` (Docker Desktop on Apple Silicon: Virtualization
framework with Rosetta) and the Seqera `wave` CLI is on PATH. The pinned
progressiveMauve (475 MB) and Cactus (1 GB) images are pulled on first use;
Wave mints and caches images for the conda-spec adapters (skder, sourmash,
dRep, snippy, sibeliaz, the simple SNP typer).
