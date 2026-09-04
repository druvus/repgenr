# Scaling and bias audit: alignment/SNP typing, tree building, tree2tax

Section of the RepGenR scaling/bias audit. All behavioral claims cite file:line.
External tool status verified by web search on 2026-08-31.

## Alignment + SNP typing

### SWOT

**Strengths.** The progressiveMauve path is deliberately linear in n: one pairwise
reference alignment per query, parallelized across the thread budget
(aligners/progressivemauve.py:74-88). SibeliaZ and Cactus are multi-genome and
reference-free internally; the reference enters only at projection
(sibeliaz.py:117, cactus.py:84-87). Adapters are uniform behind the plugin ABC
with declared scale limits and containers (docs/adding-tools.md). The simple
typer needs only minimap2/samtools/bcftools and emits a SNP distance matrix
(snptypers/simple.py:71-72). Divergence warnings exist for the too-diverse
direction (phylo.py:359-384).

**Weaknesses.** Every MSA path is reference-projected, and the reference defaults
to genomes[0] -- the alphabetically first file from list_fasta -- in five places
(progressivemauve.py:52-53, sibeliaz.py:62-63, cactus.py:46-47,
phylo.py:329-341, stages/snptype.py:84,186): an arbitrary, order-dependent bias
with no notice to the user. The simple typer is strictly serial with ~8
subprocesses per genome, holds all consensuses in RAM plus an O(n^2 S) distance
matrix (simple.py:61-68,120-147); at 5000 genomes that is ~40k process launches
and multi-GB dicts. Its core-SNP step also silently truncates all rows to the
shortest consensus (simple.py:123). Snippy runs per-genome serially (each run
gets the full thread budget, snippy.py:42-55) and snippy-core receives every
sample dir on argv -- ARG_MAX risk at scale, as do mashtree and sibeliaz genome
lists (snippy.py:59, mashtree.py:57, sibeliaz.py:86). ParSNP copies every input
genome (parsnp.py:45-48), doubling disk. The zero-SNP guard exists only in the
simple typer (simple.py:76-77); parsnp and snippy will hand an empty or
degenerate core alignment silently downstream. _warn_divergence never warns on
the opposite failure mode -- a clonal set where the reference-based core is
nearly invariant -- and infers taxonomy from filenames only (phylo.py:344-384).
Gubbins is fed the SNP-only alignment (stages/snptype.py:107-109,
maskers/gubbins.py:31-37); Gubbins expects a whole-genome alignment, so
recombination detection on concatenated SNP sites is methodologically shaky for
simple/parsnp/snippy (snippy's core.aln is SNP-only, snippy.py:63-67).

**Opportunities.** The SnpTyper ABC (call(genomes, reference, out_dir, params,
logger) -> SnpResult, docs/adding-tools.md) fits reference-free k-mer callers
directly; ska2 would remove reference bias for clonal sets. Parallelizing
simple's per-genome loop via the parallel_map already used by progressivemauve
is a small change. FOFN via write_fofn is the sanctioned pattern for long input
lists (adding-tools.md:73-76).

**Threats.** At 1000-5000 taxa the serial typers become multi-day runs before
any tree is built; argv overflows abort late; a wrong (alphabetical) reference
skews every downstream SNP distance and tree topology, worst on clonal data
where a handful of reference-private errors dominate the signal.

### Gaps

1. **High** -- Alphabetical genomes[0] reference default, unlogged (5 sites
   above). Log the choice loudly; pick a central genome (e.g. mash medoid) or
   require --reference above a size threshold.
2. **High** -- simple typer serial + O(n^2 S) in-RAM matrix
   (simple.py:61-68,120-147). parallel_map the per-genome calls; stream the
   matrix or make it optional above ~1000 genomes.
3. **High** -- argv path lists in snippy-core/mashtree/sibeliaz (snippy.py:59,
   mashtree.py:57, sibeliaz.py:86). Preflight length check + FOFN where the
   tool supports one.
4. **Medium** -- No zero/low-SNP guard in parsnp/snippy (guard only at
   simple.py:76-77). Count variable sites on the returned core FASTA in
   stages/snptype.py and fail/warn centrally.
5. **Medium** -- Gubbins gets SNP-only input (stages/snptype.py:107-109). Feed
   a whole-genome alignment (snippy core.full.aln) or document the caveat.
6. **Medium** -- simple truncates to min consensus length silently
   (simple.py:123). Assert equal lengths; error otherwise.
7. **Low** -- ParSNP input copying doubles disk (parsnp.py:45-48). Symlink
   instead.
8. **Low** -- No low-diversity warning anywhere (phylo.py:359-384 covers only
   divergence). Warn when the core SNP count is near zero per genome.

## Tree building

### SWOT

**Strengths.** Clean two-kind builder API (MSA vs GENOMES, treebuilders/base via
phylo.py:118-131) with auto-selection and scale warnings (phylo.py:96-105). The
RAxML-NG adapter caps the runaway autoMRE{1000} default at autoMRE{200}
(raxmlng.py:50-58) -- the only builder with support values out of the box.
IQ-TREE threads are capped sanely (-T auto --threads-max, iqtree.py:41).

**Weaknesses.** bootstrap defaults to 0 everywhere (phylo.py:50); IQ-TREE then
emits no support values (iqtree.py:46-47) and FastTree ignores both
params.bootstrap and params.outgroup entirely (fasttree.py:44-49) -- users get
support-free trees by default from 4 of 5 builders, and IQ-TREE rejects
-B < 1000 if a small value is passed. mashtree uses --mindepth 0 and honors the
genomesize extra in distance_matrix but not in build (mashtree.py:39 vs :57) --
the tree and the matrix are computed under different settings. The sourmash
builder converts Jaccard similarity to distance as 1 - similarity
(sourmash.py:82), which is not an evolutionary distance and compresses exactly
in the clonal regime, then feeds a pure-Python O(n^3) neighbor-joining with an
O(n^2)-per-iteration Q-scan in interpreted loops (tree/newick.py:30-67) --
behind a declared 10000-genome limit (sourmash.py:30); at 5000 taxa this is on
the order of 10^11 Python operations. Declared limits (iqtree 500 / raxmlng
1000 / fasttree 5000 / mashtree 10000 / sourmash 10000; iqtree.py:23,
raxmlng.py:26, fasttree.py:22, mashtree.py:22, sourmash.py:30) are
unbenchmarked. Label sanitization can collide distinct names (newick.py:74-79).

**Opportunities.** The FastTree adapter already anticipates VeryFastTree ("both
accepted", fasttree.py:17); preferring VeryFastTree when present is trivial.
Replacing the Python NJ with decenttree or rapidNJ behind the same
distance-matrix path removes the O(n^3) wall. attotree is a drop-in faster
mashtree.

**Threats.** For 1000-5000 clonal genomes the default path (progressivemauve
MSA -> iqtree, phylo.py:43-45) is the least scalable combination in the
codebase; near-zero-length branches plus no supports yield polytomous,
unstable topologies that tree2tax then freezes into taxonomy.

### Gaps

1. **High** -- Pure-Python NJ O(n^3) behind a 10000 limit (newick.py:30-67,
   sourmash.py:30). Swap in decenttree/rapidNJ or lower the declared limit to
   ~500.
2. **High** -- 1-minus-Jaccard used as a tree distance (sourmash.py:82). Use a
   proper distance (e.g. sourmash's ANI estimate) or document as topology-only.
3. **Medium** -- No support values by default; FastTree silently drops
   bootstrap and outgroup (phylo.py:50, iqtree.py:46-47, fasttree.py:44-49).
   Warn when bootstrap=0 with an MSA builder; log ignored params.
4. **Medium -- CLOSED** -- mashtree build/distance_matrix parameter mismatch
   (mashtree.py:39 vs :57). Both entry points now share one argument builder
   and a single mashtree call writes the tree and the matrix together.
5. **Medium** -- Declared limits unbenchmarked (all treebuilder capabilities).
   Benchmark once at 500/1000/5000 and correct the numbers.
6. **Low** -- Newick label sanitization can merge distinct labels
   (newick.py:74-79). Deduplicate after sanitizing.

## tree2tax

### SWOT

**Strengths.** Stateless dendropy core shared by the workdir and data-channel
paths (tree2tax.py:63-136); outgroup rooting resolved by exact accession before
substring (tree2tax.py:191-198); internal-node names are content-derived hashes,
stable across reruns (tree2tax.py:237-238); dereplicated genomes are re-attached
under their representative from clusters.tsv (tree2tax.py:268-290).

**Weaknesses.** Node naming hashes the sorted descendant leaf list of every
internal node (tree2tax.py:237-238) and _build_paths walks every leaf's full
ancestor chain (tree2tax.py:242-253): ~O(n^2) time and O(n^2) edge tuples on
ladder-shaped trees -- precisely the shape clonal data produces -- so 5000
leaves means ~12.5M hashed label joins and a large duplicated edge list
(tree2tax.py:255-261). On clonal input the upstream tree has near-zero branch
lengths and no supports, but tree2tax freezes every internal node into a taxon
regardless -- no collapse threshold for zero-length or unsupported branches. A
missing outgroup degrades silently to an unrooted tree with only a warning
(tree2tax.py:176-200), and rooting then depends on the builder's arbitrary
output rooting.

**Opportunities.** A single post-order pass computing child hashes
incrementally makes naming O(n); collapsing branches below a length/support
threshold before emission would make the taxonomy robust on clonal sets and is
purely local to this stage.

**Threats.** A taxonomy built from an unsupported, reference-biased tree
inherits every upstream artifact; at 5000 leaves the quadratic pass is minutes
of CPU and a swollen tsv, not a crash -- easy to miss until FlexTaxD import.

### Gaps

1. **Medium** -- O(n^2) hash naming + path building on ladder trees
   (tree2tax.py:237-238,242-261). Compute hashes bottom-up in one pass; emit
   each edge once.
2. **Medium** -- No collapse of zero-length/unsupported internal nodes; every
   node becomes a taxon (tree2tax.py:215-239). Add a collapse threshold option.
3. **Low** -- Unrooted fallback only warns (tree2tax.py:176-200). Consider
   midpoint rooting as an explicit fallback.

## Alternatives matrix

| Tool | Approach | Scaling shape | Clonal vs diverse | Adapter feasibility (plugin API) |
|---|---|---|---|---|
| iqtree (in-tree) | ML from MSA | poor >500 taxa (declared, iqtree.py:23) | good both; needs -B>=1000 for supports | n/a (present) |
| fasttree (in-tree) | approximate ML | good to ~5000 | weak resolution on clonal (approximate) | n/a (present) |
| raxmlng (in-tree) | ML + autoMRE{200} | poor >1000 | best supports; slow | n/a (present) |
| mashtree (in-tree) | MinHash distance + NJ | good ~10k (external NJ) | poor on clonal (Mash saturates near identity) | n/a (present) |
| sourmash+NJ (in-tree) | k-mer Jaccard + Python NJ | breaks well before declared 10k (newick.py:30-67) | poor on clonal (1-Jaccard compression) | n/a (present) |
| VeryFastTree | vectorized FastTree-2 reimplementation | very good; authors report 1M taxa; v4.0.5 (Apr 2026), bioconda | same regime as FastTree, much faster | trivial: same CLI family; fasttree.py:17 already anticipates it |
| decenttree | vectorized NJ/BIONJ from matrix | 64k SARS-CoV-2 genomes demonstrated; active (iqtree org) | distance-quality-bound; fine both | easy: replace neighbor_joining() call; matrix already exists |
| rapidNJ | heuristic NJ from matrix | good to tens of thousands; stable, older | as decenttree | easy, same slot |
| attotree | mashtree reimplementation (mash + NJ) | fraction of mashtree runtime; v0.1.6, bioconda | as mashtree | trivial: InputKind.GENOMES builder, ~60 lines |
| ska2 | split-k-mer reference-free SNP calling | very fast, Rust; incremental add-to-callset; Genome Research 2024, active 2026 | designed for clonal/outbreak sets; no reference bias | easy: SnpTyper with requires_reference=False; emits alignment directly |
| snp-sites | extract SNP sites from any MSA | linear, trivial | regime-neutral utility | easy: replaces simple.py's Python column scan; long-stable (not re-verified) |
| parsnp-as-aligner | core-genome MSA (not just SNPs) | ~2000 (its own guidance) | intraspecific/clonal only | moderate: already in-tree as SnpTyper; an Aligner adapter would reuse parsnp.xmfa via harvesttools -M |

External verification sources: VeryFastTree
(https://github.com/citiususc/veryfasttree/releases), SKA2
(https://genome.cshlp.org/content/34/10/1661), attotree
(https://github.com/karel-brinda/attotree), decenttree
(https://github.com/iqtree/decenttree).
