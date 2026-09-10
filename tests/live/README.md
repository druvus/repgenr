# Live verification suite

Tests under this directory run the installed `repgenr` console script as a
real process against real tools. They are marked `live` and are deselected by
default (`addopts = "-m 'not live'"` in `pyproject.toml`), so neither plain
`pytest` nor CI collects them.

```bash
cp tests/live/live.example.toml tests/live/live.local.toml   # edit paths
conda run -n repgenr_dev --no-capture-output \
    pytest tests/live -m live --live-config tests/live/live.local.toml -ra
```

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
