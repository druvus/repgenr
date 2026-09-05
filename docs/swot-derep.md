# SWOT and gap analysis: acquisition and dereplication

## Acquisition (stages/metadata.py, stages/genome.py)

### SWOT

**Strengths**
- Dual metadata sources: release-pinned GTDB TSV with MD5 verification and cleanup of corrupt transfers (metadata.py:155-178), or a low-transfer API path (metadata.py:371-398). Selection is replace-not-upsert so re-selection removes stale rows (metadata.py:298-300), and a portable `selection.tsv` hand-off supports the Nextflow data channel (metadata.py:307-319).
- Download robustness: fixed 5000-accession sub-batches so a failure loses one batch and re-runs resume by recomputing what is missing (genome.py:32, genome.py:69-75, genome.py:171-195); FASTA-shape validation rejects HTML-error-page bodies before they land at a genome path (genome.py:69-72, genome.py:226-234); undelivered accessions are recorded for downstream completeness checks (genome.py:91-93); disk-space floor and estimate guard (genome.py:127-144); stale-genome pruning is restricted to FASTA files (genome.py:117-124).
- `dataset=rep` restricts to GTDB species representatives (metadata.py:215-217), the single strongest built-in mitigation against clone blocks at genus/family level.

**Weaknesses**
- `--limit` takes the first N genomes in GTDB file order (metadata.py:268-275) and the API path slices `rows[:limit]` (metadata.py:381-382): a biased, order-dependent subsample, not a random or stratified one.
- No assembly-quality filtering at acquisition: parsing reads only accession, taxonomy, and representative flag (metadata.py:210-225), although the GTDB table carries CheckM completeness/contamination. All downstream stages inherit whatever quality mix NCBI holds.
- Flat 3600 s timeout per `datasets` invocation regardless of batch size (genome.py:37): a legitimate 5000-genome rehydrate on a slow link can be killed and retried indefinitely.
- Automatic outgroup selection returns the first matching representative in dict iteration order (metadata.py:290-294; API: first qualifying row, metadata.py:422-430) -- deterministic per release file, but sensitive to GTDB row order between releases.
- No content checksums for downloaded genomes; only the "starts with `>`" heuristic (genome.py:147-160), so a truncated FASTA passes.

**Opportunities**
- Stratified or random `--limit` sampling (per-species quota) would remove the ordering bias cheaply.
- Reading GTDB quality columns during parsing would enable quality-aware selection and quality-informed dereplication downstream at no extra download cost.
- Batches are downloaded serially (genome.py:189-195); bounded parallel batches would help at 5000+.

**Threats**
- GTDB layout drift is already real (tsv.gz vs tar.gz fallback, metadata.py:152-178); the API path tracks only the current release, so `api`-based runs are not reproducible across GTDB releases (metadata.py:64-69).
- NCBI throttling/outages during a 5000-genome batch surface as one opaque tool failure; per-accession retry is not possible inside `datasets`.

## Dereplication (stages/dereplicate.py + adapters)

### SWOT

**Strengths**
- Clean plugin contract (adapter class + entry point, docs/adding-tools.md:10-129) with auto-selection and scale warnings (dereplicate.py:77-88), unconsumed-extra warnings (dereplicate.py:95-100), and a normalized `DerepResult` contract with zero-length-representative guards (dereplicate.py:460-478).
- Recursive chunking bounds any single tool call at `--process-size` genomes, with parallel chunk workers and a measured ~2-3x speedup rationale (dereplicate.py:214-226, dereplicate.py:287-354); membership is composed in O(N) (dereplicate.py:357-389).
- Genuine scale engineering: fofn inputs avoid ARG_MAX (galah.py:54-56, drep.py:92-94, sourmash.py:179-180), sourmash refuses the dense N x N path above 5000 and requires the sparse branchwater backend (sourmash.py:47-49, sourmash.py:128-134), and the target-reps search reuses a digest-keyed sketch cache (dereplicate.py:244-261, sourmash.py:296-300).

**Weaknesses (scientific defensibility)**
- Representative choice is quality-blind in three of four tools: sourmash picks the most-connected genome with alphabetical-filename tie-break (dereplicators/sourmash.py:446-486, :459), skder membership is best-ANI first-seen-wins (skder.py:142-152), and the galah adapter passes no quality file so galah's quality-aware scoring runs on defaults (galah.py:58-64). Clone blocks (over-sequenced outbreak genotypes) therefore dominate connectivity and capture representative choice.
- The `secondary_ani` knob means different things per tool: skder/dRep treat it as ANI, sourmash treats it as Jaccard similarity (sourmash.py:75-76), which is numerically far stricter (99% ANI is nowhere near 0.99 Jaccard at k=31). `aligned_fraction` is silently ignored by galah and sourmash. Cross-tool results are not comparable at the same CLI settings.
- Stage-1 chunk collapse is irreversible and chunks are sliced alphabetically (dereplicate.py:310 over the sorted list from dereplicate.py:159-160): filename order correlates with taxonomy, so chunk composition is systematically non-random and chunk-local greedy choices can differ from a single-pass run.
- `--target-reps` bisects secondary ANI in [0.80, 0.9999] for at most 12 iterations and silently returns the closest count (dereplicate.py:229-284); taxonomy `--reduce` keeps the largest cluster per taxon (dereplicate.py:422) -- again a most-sequenced-genotype criterion.
- dRep `--virus` weights size only, with single-linkage chaining (drep.py:104-112), and appends a second `--S_algorithm` flag relying on argparse last-wins (drep.py:101, drep.py:106-112).
- skder genomes absent from the edge table are marked contained with no representative (skder.py:154-156); galah returns an empty result if `clusters.tsv` has no rows, with no completeness check (galah.py:72-87).

**Opportunities** -- quality-aware keeper scoring at the contract layer; wiring CheckM2/GTDB quality into galah and dRep; randomized or minhash-ordered chunking.

**Threats** -- all `recommended_max_genomes` values are unbenchmarked, so `--tool auto` and the scale warnings rest on guesses; skder still passes genome paths on argv with only a warning above 5000 (skder.py:46, skder.py:78-84).

### Gap list

1. **High** -- Per-tool semantic drift of `secondary_ani` (ANI vs Jaccard) and silently ignored `aligned_fraction` (sourmash.py:75-76, galah.py:51-64). Remediation: convert Jaccard to ANI-equivalent (or document per-tool semantics) and warn when a threshold parameter is unused.
2. **High -- CLOSED** -- Quality-blind, clone-block-sensitive representative selection (sourmash.py:446-486; skder.py:142-152; galah.py:58-64). Remediation: optional quality table (GTDB/CheckM2) consumed by a shared keeper-scoring step after adapter output. Closed by `stages/derep_keeper.py`: after any adapter runs, `rescore_representatives` re-picks each cluster's representative by `completeness - 5 x contamination` from `Manifest.quality()` (GTDB CheckM values), gated by `--keeper quality` (default) / `--keeper tool` to keep the adapter's own pick.
3. **High** -- skder argv-based invocation breaks past ARG_MAX at 5000+ genomes (skder.py:78-96). Remediation: chunk automatically or use skder's fofn-style input if available.
4. **Medium** -- Alphabetical chunk slicing biases stage-1 collapse (dereplicate.py:310). Remediation: shuffle with a recorded seed, or interleave by species.
5. **Medium -- CLOSED** -- `--limit` first-N sampling bias (metadata.py:268-275, :381-382). Closed 2026-09-05: both paths now cut by species-stratified round-robin with CheckM quality ranking within each species (`_stratified_limit`); file order plays no part.
6. **Medium** -- Flat 3600 s download timeout independent of batch size (genome.py:37). Remediation: scale timeout with batch size.
7. **Medium** -- `--target-reps` silently accepts the closest count and taxonomy reduce keeps the largest cluster (dereplicate.py:229-284, :422). Remediation: report the miss explicitly; make the reduce keeper criterion pluggable.
8. **Medium** -- skder/galah membership completeness holes (skder.py:154-156; galah.py:72-87). Remediation: assert every input genome receives a status and a representative, else fail.
9. **Low** -- Unbenchmarked `recommended_max_genomes` everywhere (e.g. drep.py:44). Remediation: benchmark on synthetic 1k/5k/20k sets and record measured limits.
10. **Low** -- Outgroup choice depends on GTDB row order (metadata.py:290-294). Remediation: sort candidates before picking.

### Dereplication alternatives matrix

| Tool | Algorithmic approach | Expected scaling | Representative criterion (quality-aware? deterministic?) | Clone-block behavior | Adapter status / feasibility |
|---|---|---|---|---|---|
| skder (in-tree) | skani sketch ANI + dynamic/greedy dereplication | Near-linear in practice; native to 10k+ (skder.py:2-4) but argv-limited ~5000 (skder.py:46) | skDER internal score (length/N50-based); not externally quality-aware; deterministic; membership first-seen ties (skder.py:142-152) | One rep per block; tie assignment arbitrary within block | Shipped (skder.py:49) |
| galah (in-tree) | dashing/skani prefilter + greedy clustering | Built for large sets; native scaling (galah.py:38-39) | Quality-aware only with a quality file, which the adapter does not pass (galah.py:58-64); deterministic | Rep chosen by galah defaults, not block size | Shipped; quality-file wiring is a small extension |
| sourmash (in-tree) | MinHash Jaccard; greedy most-connected | Sparse path ~linear in close pairs; dense capped at 5000 (sourmash.py:49) | Most-connected, alphabetical tie-break (sourmash.py:446-486, :459); not quality-aware; deterministic | Biased toward the densest (most-sequenced) genotype | Shipped |
| dRep (in-tree) | Mash primary + fastANI/ANImf secondary; weighted scoring | O(N^2) within primary clusters; recommended max 2000 (drep.py:44), chunk-wrapped | Score of completeness/contamination/N50/size -- quality-aware, deterministic; virus branch size-only (drep.py:104-112) | Best-scored genome per cluster, robust to block size | Shipped |
| dashing2 | SetSketch/ProbMinHash + LSH all-pairs | Sublinear candidate generation; demonstrated at millions of sequences | None built in -- adapter must implement greedy over the edge list (choice ours, can be quality-aware) | Depends on our greedy; same risks as sourmash unless scored | Feasible (edge-list output fits the sparse-cluster path); repo active, but release recency unverified |
| fastANI + greedy | Exact-ish alignment-free ANI, all pairs | O(N^2) compute; impractical past ~2-5k without prefilter | Ours to define; can be quality-aware and deterministic | Depends on our clustering | Feasible but low value: last release 1.34 (Jul 2023), repo semi-maintained |
| Assembly-Dereplicator | Iterative Mash distance pruning, keeps best assembly | Mash-based, handles thousands; fofn input since v0.3.1 | Assembly-quality heuristic (N50-based) -- partially quality-aware, deterministic | Keeps one good assembly per block | Feasible; v0.3.1 is the latest release, low activity; emits no cluster membership, so `clusters.tsv` needs reconstruction |
| CheckM2-informed dRep | dRep with `--genomeInfo` from CheckM2 | Same as dRep (chunk-wrapped) | Completeness/contamination-weighted -- fully quality-aware, deterministic | Best-quality genome per cluster; strongest clone-block answer | Straightforward extension of drep.py; CheckM2 actively maintained (v1.1.x, Feb 2025) |

External-tool status sources (checked 2026-08-31): [dashing2](https://github.com/dnbaker/dashing2), [Assembly-Dereplicator](https://github.com/rrwick/Assembly-Dereplicator) ([v0.3.1](https://zenodo.org/records/7894123)), [FastANI releases](https://github.com/ParBLiSS/FastANI/releases), [CheckM2 releases](https://github.com/chklovski/CheckM2/releases).
