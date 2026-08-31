# Scaling and bias benchmarks

Audit tooling for the scaling/bias study (`docs/scaling-audit.md`). Lives outside
`src/repgenr` on purpose: it measures the pipeline, it is not part of it.

## Layout

- `genomegen.py` — seeded synthetic genome-set generator (pure stdlib). Scenarios
  `balanced`, `clonal` (over-represented clone block + diverse background), and
  `mixed`; `--order {clustered,interleaved,random}` controls where clone
  accessions land in the alphabetical sort, which is itself a mechanism under
  test. Each set gets a `truth.json` with the true cluster membership.
- `cells.py` — the benchmark matrix. A cell is one measured run: dereplication
  steps ({skder, galah, sourmash-sparse, sourmash-dense} x {100, 1000, 5000} x
  {balanced, clonal}), chunked-vs-single-pass at n=5000 with both accession
  orders, and tree builders (mashtree to 5000; the pure-Python sourmash NJ only
  to 1000; ML builders on small subsets via sibeliaz).
- `run_bench.py` — resumable runner. Wraps each cell in `/usr/bin/time -l`
  (macOS: max RSS in bytes), writes one JSON per completed cell to
  `benchmarks/results/cells/`, and aggregates with `--collect`.
- `dense_driver.py` — forces the sourmash dense (N x N) path in-process; there
  is no CLI switch that disables the branchwater sparse path.
- `metrics.py` — stdlib scoring: adjusted Rand index vs `truth.json`,
  representative counts, and which genome won the clone block.

## Data locations

Genome sets and per-cell workdirs live on external storage
(`/Volumes/sekvens2/repgenr/{sets,work,logs}`), not in the repo. Only the small
result JSONs and `summary.tsv` are committed under `benchmarks/results/`.

## Running

```bash
PY=$(conda run -n repgenr_dev python -c "import sys; print(sys.executable)")

# Generate sets (one-off, ~2 min per 5000-genome set)
$PY -m benchmarks.genomegen --scenario clonal --n 1000 \
    --out /Volumes/sekvens2/repgenr/sets/clonal_1000_clustered --seed 1

# See what would run, then run resumably (tiers: smoke -> mid -> heavy)
$PY -m benchmarks.run_bench --tier smoke --dry-run
$PY -m benchmarks.run_bench --tier smoke --resume
caffeinate -i $PY -m benchmarks.run_bench --tier heavy --resume   # overnight

# Aggregate committed cell JSONs into results/summary.tsv
$PY -m benchmarks.run_bench --collect
```

A cell is considered done when `results/cells/<id>.json` exists with
`"status": "ok"`; `--resume` skips those, so an interrupted tier can simply be
restarted. `--only '<glob>'` selects cells by id.

## Caveats

Synthetic genomes are mutated copies of a random ancestor: k-mer and
alignment-identity tools (the in-scope chain) see exactly the designed ANI
structure, but gene-model tools would not find real genes, and whole-genome
aligner timings (sibeliaz) are a lower bound. skani-validated: clone blocks
measure 99.98-99.99 percent ANI, background clusters about 95 percent.
