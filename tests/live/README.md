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
