# Viral chain and the scale-declaration surface

## 1. SWOT: the viral chain

**Strengths**

- **`median_of_medians` is a real over-representation defense.** The default length window
  is built from one median per species, then the median of those
  (`viral/selection.py:151-160`), so 900 outbreak isolates of one species contribute a single
  vote. The BV-BRC path does the same over per-taxid medians (`viral/bvbrc.py:216-234`).
  This is the single most important bias control in the chain.
- **Lazy sequence access.** `_seq_map` returns `SeqIO.index` rather than a dict, so only kept
  and outgroup accessions are materialized (`viral/selection.py:165-172`); BV-BRC does the same
  and scans the group FASTA once for metadata instead of four times
  (`viral/bvbrc.py:55-58,150-164,167-181`).
- **Tolerance divergence is deliberate and recorded.** `RECORDS_LENGTH_TOLERANCE = 0.15` vs
  `BVBRC_LENGTH_TOLERANCE = 0.10` differ because they are measured against different centres
  (selection-range midpoint vs candidate-taxid median), and the module docstring says so and
  says not to unify them (`viral/_outgroup.py:10-21,33-34`).
- **Fail-soft outgroup on the default path.** Missing candidates, no length-compatible
  candidates, or an unassignable label all warn and proceed
  (`viral/selection.py:252,279,284`), so a thin genus still produces a repository.
- **Single-source metadata.** NCBI Virus `datasets` yields sequences and per-record taxonomy in
  one call (`viral/ncbi_virus.py:112-165`), removing the BV-BRC FTPS + Entrez batching path
  from the default flow.

**Weaknesses**

- **The defense is undocumented.** `--length-method` has no help text
  (`cli/cmd_viral.py:66`) and appears nowhere in `README.md` or `docs/`. A user who sets
  `--length-method mean` silently disables it, and on the records path `mean` is taken over
  *every record* (`mean(all_lens)`, `viral/selection.py:154,158`), so the most-sequenced
  species sets the window. The BV-BRC `mean` path averages per-taxid means
  (`viral/bvbrc.py:216-234`) and is therefore *not* equally vulnerable -- the two back-ends
  disagree about what `mean` means.
- **The `S_` panel is order-dependent, not representative.** `_emit("S", list(kept), 12)`
  (`viral/selection.py:273`) takes the first 12 kept records in download order. For an outbreak
  set that is 12 near-clones, so the mashtree distance matrix compares candidates against one
  clade rather than against the selection's diversity.
- **`group_segments` bypasses length filtering and the outgroup entirely.** Kept records are
  everything with a sequence (`viral/selection.py:58-62`) and outgroup selection is skipped
  (`viral/selection.py:92`). Segmented viruses get no length QC at all: a truncated segment
  enters a concatenated genome unchecked.
- **Segment concatenation order is length-ranked, not segment-ranked.**
  `sorted(recs, key=lambda r: -r.length)` (`viral/selection.py:214`). Two isolates whose
  segments are close in length can concatenate in different orders, producing non-homologous
  columns downstream.
- **Serotype matching is substring-based.** `_norm(v) in _norm(rec.organism)`
  (`viral/selection.py:131`) -- `h1` matches `h11n9`, `h13`, etc. Over-broad selection is
  indistinguishable from a correct one in the logs.
- **Species granularity depends on NCBI's `organismName`.** `_classify` sets species to the
  organism leaf (`viral/ncbi_virus.py:59-75`). Where NCBI reports strain-level organism names,
  `by_species` fragments into singletons and median-of-medians degrades toward a plain median
  over records -- the defense weakens exactly on the datasets that need it.
- **Off-by-one in the BV-BRC candidate cap.** `if written.get(rec.taxid, 0) > max_per_taxid`
  with `max_per_taxid = 3` (`viral/bvbrc.py:300-303`) admits **four** per taxid; the records
  path's `written >= cap` (`viral/selection.py:263`) admits exactly three.
- **No global `S_` cap on the BV-BRC path.** Records-path `S_` is capped at 12; BV-BRC caps per
  taxid only, so a 100-taxid selection stages ~400 `S_` files into mashtree
  (`viral/bvbrc.py:302-314`).
- **`SeqIO.index` handles are never closed** (`viral/selection.py:172`,
  `viral/bvbrc.py:175`) -- a file descriptor per stage invocation in an embedding caller.

**Opportunities**

- Log the chosen length method, the number of species voting, and the resulting window at
  INFO; today only the window is logged (`viral/selection.py:161`).
- Replace `_emit("S", ...)` first-12 with a diversity-spread pick (e.g. length-stratified or a
  cheap sketch-distance spread) -- the sequences are already indexed, so this is local.
- Length-filter segments *per segment name* under `group_segments` (the `segment` field is
  already parsed, `viral/ncbi_virus.py:105`) and order concatenation by segment name, falling
  back to length only for `ANONYMOUS`.
- Make serotype matching token-boundary-aware, or add `--target-serotype-exact`.
- Expose an outgroup panel size (currently the literals `12` and `3`,
  `viral/selection.py:273-275`) so large genera can widen the comparison.

**Threats**

- **Undocumented default that users will change.** Because `--length-method` is
  undiscoverable, the most likely way it gets exercised is copy-paste from a forum, at which
  point the over-representation guard is off and nothing in the output records it --
  `record_stage` for vgenome writes `source`, `selected`, `group_segments`, `no_outgroup`
  only (`viral/selection.py:104-112`), not the length method or window.
- **Outgroup candidacy is popularity-gated.** Candidate species need
  `>= outgroup_candidates_taxid_min_genomes` (default 5) records
  (`viral/selection.py:247-250`), so a well-sampled but phylogenetically closer species beats
  a correct but rare sister taxon.
- **Scale is unproven above ~1.3k.** The largest recorded viral run is 1256 hepatitis E
  genomes (`docs/verification.md:54`). The records path holds all `VirusRecord` objects in
  memory plus the `SeqIO.index` offset table; `datasets download virus genome taxon` for a
  large family (e.g. all `Flaviviridae`) is 10^5-10^6 sequences, and the whole
  `virus_records.json` is loaded eagerly (`viral/ncbi_virus.py:172-173`).
- **Nothing validates concatenated-genome sanity.** No warning fires when isolates end up with
  differing segment counts, which is the normal case for incomplete GenBank submissions.

## 2. SWOT: the scale-declaration system

**Strengths**

- One declarative field drives three behaviours (auto-selection, over-scale warning,
  alternative suggestions) with no per-stage tables: `ToolCapabilities.recommended_max_genomes`
  (`core/plugins.py:51`), `auto_select` (`core/plugins.py:207-241`), `scale_warning`
  (`core/plugins.py:244-266`). Third-party adapters participate for free.
- `auto_select` is environment-aware: availability is the primary sort key and is
  container-aware, so a declared image counts as available under `--container`
  (`core/plugins.py:180-192,232-238`).
- Tightest-fit ranking is correct and tested. Among fitting tools the *smallest* limit wins
  (`core/plugins.py:229-237`), so the accurate-but-limited tool beats the scalable-but-coarse
  one on small inputs; unbounded tools sort last among fitters
  (`tests/unit/test_auto_select.py:90`).
- Graceful degradation: a broken plugin is skipped with a warning rather than aborting
  auto-selection (`core/plugins.py:220-226`), and never fits by accident.

**Weaknesses**

- **All 15 declared values are unbenchmarked round numbers** -- 500/1000/2000/5000/10000 across
  `aligners/cactus.py:34`, `sibeliaz.py:50`, `progressivemauve.py:40`, `treebuilders/fasttree.py:22`,
  `mashtree.py:22`, `raxmlng.py:26`, `sourmash.py:30`, `iqtree.py:23`, `snptypers/snippy.py:24`,
  `simple.py:34`, `parsnp.py:25`, `dereplicators/drep.py:44`, plus `None` for `galah.py:38`,
  `sourmash.py:61`, `skder.py:57`. No measurement, no hardware assumption, no citation is
  recorded anywhere. `iqtree=500` and `progressivemauve=500` being equal is coincidence, not
  a finding.
- **Warnings only log.** Every call site emits `logger.warning` and proceeds
  (`stages/dereplicate.py:81-88`, `stages/phylo.py:99-105,288-294`, `stages/snptype.py:72-78`).
  There is no refuse threshold, no `--allow-over-scale` gate, and nothing lands in
  `repgenr.yaml`, so a run that blew past a limit is indistinguishable afterwards.
- **The `auto` path never warns.** In `dereplicate.py:77-88` and `phylo.py:95-105` the warning
  lives in the `else` branch. When nothing fits, `auto_select` silently returns the
  largest-capacity tool (`core/plugins.py:214-215`) and the user is told only which tool was
  picked, not that it is over scale.
- **`scale_warning` suggests unavailable tools.** Unlike `auto_select`, it never calls
  `_tool_available` (`core/plugins.py:257-265`), so "consider: cactus, sibeliaz" can name tools
  that are not installed and have no image.
- **Limits are invisible until violated.** `repgenr list-tools` prints names and a
  `(broken)` marker only (`cli/cmd_misc.py:174-189`); nothing surfaces the declared scale.
- **Genome count is the only dimension.** No memory, no genome *size*, no thread count. 2000
  bacterial genomes and 2000 viral genomes are treated identically by `drep=2000`, and
  `cactus=2000` carries the comment "Toil manages its own parallelism"
  (`aligners/cactus.py:34`) -- an assumption about the executor, encoded as a genome count.
- **`_PREFERRED_ORDER` hardcodes in-tree names in core** (`core/plugins.py:197`), so a
  third-party adapter can never win a tie. This is the one plugin-neutrality leak left in
  `plugins.py`.
- **Maskers have no scale declaration and no call site.** `scale_warning` is never invoked for
  the masker family.
- **The README claim outruns the evidence.** "scaling to thousands of genomes"
  (`README.md:11`) against a largest recorded run of 1256 (`docs/verification.md:54`).

**Opportunities**

- Record the tool's declared limit and the actual `n_items` in `record_stage` provenance, so
  over-scale runs are auditable after the fact and the limits become empirically refutable.
- Add `recommended_max_genomes` (and whether it fits the current input) to `list-tools`.
- Add a soft/hard split: `recommended_max_genomes` keeps warning; a new
  `max_genomes` refuses without an explicit override flag.
- Derive one benchmark per family from the existing 1256-genome hepatitis E set and replace
  the roundest values (`iqtree=500`, `drep=2000`) with measured ones plus a comment naming the
  hardware.
- Make `scale_warning` filter alternatives through `_tool_available`, reusing the function
  `auto_select` already calls.

**Threats**

- Values that look authoritative but were never measured are worse than absent values: they
  make `auto` produce confident wrong choices, and users will cite them.
- Tightest-fit ranking amplifies this. Because the *smallest* declared limit wins,
  the least-benchmarked small number (`iqtree=500`) becomes the default tree builder for the
  common case, purely on the strength of an unverified constant.
- A warning that never blocks trains users to ignore it, so the eventual real OOM at 20k
  genomes arrives with no distinguishing signal.

## 3. Gap list

| # | Sev | Gap | One-line remediation |
|---|-----|-----|----------------------|
| 1 | High | `--length-method` undocumented and unhelped; `mean` disables the over-representation guard (`cli/cmd_viral.py:66`, `viral/selection.py:158`) | Add help text plus a `docs/usage.md` subsection stating that `median_of_medians` is one-vote-per-species and `mean` is not |
| 2 | High | Records-path `mean` averages per record, BV-BRC `mean` averages per taxid (`viral/selection.py:154-158` vs `viral/bvbrc.py:216-234`) | Make the records path average per-species means so `mean` means the same thing on both back-ends |
| 3 | High | All 15 `recommended_max_genomes` values unbenchmarked (`core/plugins.py:51` + 15 adapters) | Benchmark one input size per family; annotate each constant with the measurement and hardware, or set it to `None` |
| 4 | High | `--tool auto` never emits a scale warning (`stages/dereplicate.py:77-88`, `stages/phylo.py:95-105`) | Call `scale_warning` on the resolved tool after auto-selection, not only in the manual branch |
| 5 | High | `group_segments` skips length filtering entirely (`viral/selection.py:58-62`) | Filter per segment name using the already-parsed `segment` field |
| 6 | Med | `S_` outgroup panel is the first 12 kept records (`viral/selection.py:273`) | Select the 12 by length-stratified spread rather than input order |
| 7 | Med | Segment concatenation ordered by descending length (`viral/selection.py:214`) | Order by segment name, falling back to length only for `ANONYMOUS` |
| 8 | Med | `scale_warning` suggests tools that are not installed (`core/plugins.py:257-265`) | Filter alternatives through `_tool_available` |
| 9 | Med | Serotype match is substring-based (`viral/selection.py:131`) | Match on token boundaries within `organism`, or add an exact-match flag |
| 10 | Med | Over-scale runs leave no trace in provenance (all four call sites) | Record `{tool, limit, n_items}` in `record_stage` whenever a warning fires |
| 11 | Med | Outgroup candidacy gated on record count >= 5 (`viral/selection.py:247-250`) | Document the bias; consider ranking by taxonomic distance before applying the count gate |
| 12 | Low | BV-BRC per-taxid cap admits 4 not 3 (`viral/bvbrc.py:300-303`) | Change `>` to `>=` to match `viral/selection.py:263` |
| 13 | Low | No global `S_` cap on the BV-BRC path (`viral/bvbrc.py:302-314`) | Apply the records-path global cap of 12 |
| 14 | Low | `_PREFERRED_ORDER` hardcodes in-tree names in core (`core/plugins.py:197`) | Move tie-break preference onto `ToolCapabilities` as an integer field |
| 15 | Low | Declared limits invisible in `list-tools` (`cli/cmd_misc.py:174-189`) | Print `recommended_max_genomes` per tool |
| 16 | Low | `SeqIO.index` handles never closed (`viral/selection.py:172`, `viral/bvbrc.py:175`) | Close the index in a `finally`, or return it from a context manager |
| 17 | Low | README claims "thousands of genomes" against a 1256-genome largest run (`README.md:11`, `docs/verification.md:54`) | State the verified scale and the date it was measured |

## 4. Adapter-candidate assessment

Every candidate needs the same three things (`docs/adding-tools.md`, sections "1. Subclass the
family ABC", "ToolCapabilities", "2. Register the entry point"): a subclass implementing the
family's single method, a `capabilities = ToolCapabilities(...)` with `required_binaries` and a
`container`/`conda` spec, and one entry-point line in `pyproject.toml`. Note the doc's hard
requirement: *"Set `container` and/or `conda` or the tool cannot run under
`--container`/`--wave`"*. All six were verified by web search on 2026-08-31.

| Tool | Verified | Adapter shape | Blocker |
|------|----------|---------------|---------|
| **VeryFastTree** v4.0.5 (Apr 2025), GPL-3.0, bioconda | yes | `TreeBuilder`, default `InputKind.MSA_FASTA`, `build(msa, out_dir, params, logger) -> Newick Path`; `conda=("bioconda::veryfasttree",)` | None. Cleanest candidate -- FastTree-2 CLI, Newick out, BioContainer available. |
| **decenttree** 1.0.0 (bioconda recipe Sep 2025), GPL-2.0 | yes | `TreeBuilder`, `build` from an MSA or a PHYLIP distance matrix; Newick out | None functional. No upstream GitHub releases, so pin the bioconda build, not a tag. |
| **rapidNJ** 2.3.3 via the `johnlees` fork, GPL-2.0, bioconda | yes | `TreeBuilder`, MSA in, Newick out | **Format blocker**: requires single-line (unwrapped) FASTA or Stockholm, so the adapter must rewrite the MSA before invoking. Upstream (`somme89/rapidNJ`) is dormant; bioconda tracks the fork. |
| **ska2** v0.5.1 (Jan 2025), Apache-2.0, bioconda | yes | Best fit is `Aligner`: `ska build -f <list>` then `ska align`, returning `AlignResult`. `ska distance` could additionally back `TreeBuilder.distance_matrix` for the viral outgroup step | **fofn shape mismatch**: `-f` takes a two-column `name<TAB>path` file, while `core.process.write_fofn` writes bare paths (used by `dereplicators/drep.py:60`, `galah.py:56`, `sourmash.py:149`). Needs its own writer. |
| **dashing2** 2.1.20 (bioconda Apr 2025), MIT | yes | `Dereplicator` with the optional `compare()` hook (`docs/adding-tools.md`, "Optional capability hooks"), and/or a `TreeBuilder.distance_matrix` provider. Native fofn via `-F paths.txt` | **Not a dereplicator by itself** -- it emits distances, so the clustering must live in the adapter, as `dereplicators/sourmash.py` already does. Output is PHYLIP or TSV and must be converted to the `genome1`/`genome2`/`similarity` CSV that `CompareResult` expects, and to the mashtree-style labelled TSV that `viral/_common.py:56-101` parses. |
| **Assembly-Dereplicator** v0.3.2 (Aug 2023), GPL-3.0 | yes | `Dereplicator`, `dereplicate(genomes, out_dir, params, logger) -> DerepResult`; representatives are copied into an output directory, which maps onto `DerepResult` directly | **Two blockers.** (a) Input is a *directory*, not argv paths and not a fofn -- the adapter must stage symlinks, for which `dereplicators/drep.py:25,94` (`link_or_copy`) is the precedent. (b) **Not on bioconda and no BioContainer**, so neither `container` nor `conda` can be set and the tool is unusable under `--container`/`--wave`; it also shells out to Mash, which would need its own `BinarySpec`. Low upstream churn (last release 2023). |
