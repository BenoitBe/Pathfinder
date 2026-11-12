# Changelog

## Unreleased

### Added
- **Day / night button** in the report header (2026-11-09). Night stays the default; the day
  theme puts the page and every figure on white. Page colours are CSS variables; figure
  neutrals (background, panels, grid, text, white lines) are swapped by the browser from one
  map (`LIGHT_SWAP`), light hues are darkened on text only (`LIGHT_TEXT`), marks keep their
  colours. SVG / PNG exports follow the theme on screen (same maps written into the copy).
  The choice is stored in `localStorage` (ignored when the browser blocks it).
- Test: every grey used by a figure must have a day equivalent.
- **Variance shares per protein** (2026-11-10, figure audit, batch 1): split violins under
  *What drives the response*, one pair per design term. Right half = selected proteins,
  left half = all quantified proteins without missing value (background). The selected
  proteins were picked on between-condition differences, which inflates their factor
  shares: a factor is read by its excess over the background. Densities are bounded to
  [0, 1] (reflection); each half is scaled to its own maximum. New sheets
  `Variance_quantiles` (quartiles per set and term) and `Variance_shares` (shares of every
  quantified protein).
- **Shared and divergent responses** (2026-11-10, figure audit, batch 1): log2FC against
  log2FC. A matrix of every pair of contrasts (3 or more) coloured by the standardised
  major-axis slope on the shared DEP, then a scatter for each of the 3 pairs sharing the
  most DEP (status of each protein, quadrant counts, sign concordance, fit drawn over the
  span of the shared DEP only). The slope is fitted on the shared DEP alone: on the union,
  proteins significant in one contrast only form a cross along the axes and pull the slope
  towards it. Same rule as the FC vs FC page of the Deimos dashboard. New sheet `FC_vs_FC`
  (one row per pair).
- **Enriched terms and their redundancy** (2026-11-10, figure audit, batch 1), before the
  functional groups: a squarified treemap of the enriched terms (area = −log10 q, colour =
  Meta-GO group, each group outlined and named) and, under it, the redundancy tree itself —
  terms joined by 1 − Jaccard on their proteins, average linkage, with the `--jaccard` cut
  drawn. Cutting that tree is what defines the groups, so the pair shows how much of an
  enrichment is one redundant block rather than independent evidence. No new dependency:
  the layout is `_squarify` (Bruls 2000) in numpy and the tree reuses
  `axis_meta_go._jaccard_condensed`.
- **Order and end-to-end change** (2026-11-10, figure audit, batch 1), for factorial
  designs, before *What drives the response*: a slope graph of the proteins moving most
  between the end levels of the ordered factor (mean log2 abundance over the conditions
  sharing the level, labels spread by isotonic regression with leader lines), and a bump
  chart of the programmes ranked by median z at each level, one panel per level of the
  lines factor. Crossing lines mean the programmes swap order, which an amplitude plot
  hides.
- **Dynamic range** (2026-11-10, figure audit, batch 1), before the enriched terms: the
  proteins ranked by mean log2 abundance with the DEP marked, and a strip giving the share
  of DEP in each decile of abundance. It answers whether the differential proteins are
  simply the abundant ones; the reading bullet says so when the top decile holds more than
  2.5× the share of the bottom one.
- **Programmes and functional groups** (2026-11-10, figure audit, batch 1), after the
  volcanos: bipartite graph of the proteins shared by each response programme and each
  Meta-GO group, line width = proteins in common, solid when more than expected
  (hypergeometric, BH over the whole table). A group tied to several programmes is one the
  clustering does not split. New sheet `Programme_groups` and a new public
  `programmes.theme_overlaps()` (the long table behind the existing `_theme_enrichment`,
  which is unchanged).
- **Shape of each programme** (2026-11-10, figure audit, batch 1): radar, one spoke per
  condition in design order, one polygon per programme, radii = median z on one common
  scale with the z = 0 ring drawn. Only for a design with 3 to 8 conditions; beyond that
  the programme profiles read better and the figure is skipped.

### Fixed
- **The whole *Biological reading* tab vanished on Python 3.11** (`SyntaxError:
  unterminated string literal`, `reading/svg.py`): the dynamic-range figure used an
  f-string whose expression spanned two source lines, which only became legal in 3.12
  (PEP 701). Written on a 3.13 machine, the module would not parse on 3.11, so the
  import of `reading.run` failed and the tab was dropped — with the reason only in
  `pathfinder_log.txt`. `tests/test_reading.py` now compiles every module with the
  oldest interpreter available on the machine, so a 3.12-only construct fails the
  suite instead of reaching a user.
- **A failed biological reading now says so in the report.** The tab used to vanish
  entirely, leaving the reason in `pathfinder_log.txt` only. It is now kept and carries
  a red banner naming the exception, stating that the rest of the report is unaffected,
  and pointing at `PATHFINDER_DEBUG=1` for the traceback; an `ImportError` or
  `SyntaxError` additionally says to check that the whole `reading/` folder matches
  `pathfinder.py`. New `red` theme variable so the banner reads in both themes.
- Biological reading skipped without a clear reason: when the `reading/` package is not
  next to `pathfinder.py` (an install updated file by file), the import failed and the
  whole *Biological reading* tab disappeared with only a generic ImportError. The
  message now names the missing folder and says to copy all of `reading/`.
- Day theme: the edges of the functional-group network (`#57606a`) had no day equivalent and
  stayed dark on white. Added to `LIGHT_SWAP` (the invariant test only caught it on a run
  whose groups are linked).

### Changed
- Factor profiles: the reference level is called the *neutral line* (white at night, black
  by day) instead of "in white".

## v1.3.0 — 2026-11-07

### Added
- **Biological reading tab** — InSight merged into Pathfinder: one report. Read from what
  Pathfinder already loads (same selection rule as the Proteins tab, same Meta-GO groups)
  plus Deimos' `matrices/` folder (DPC or imputed abundances, observed flags, samples.tsv),
  `PCA_Scores` and `last_config.yaml`. Interactive SVG in the report style (hover = name,
  log2FC, p, programme, group; clickable volcano legends; SVG / PNG export):
  how much changes; condition map (trajectories for ordered designs); response programmes;
  functional groups (enrichment map + group × contrast matrix, 16 most significant groups);
  annotated volcanos; spotlight proteins. Every figure is guarded: one failing figure never
  stops the others or the report.
- **Functional groups = the Meta-GO redundancy groups** (average-linkage Jaccard, `--jaccard`
  cut), named by the Meta-GO representative term: the same group carries the same name in
  every tab. InSight's own themes (Jaccard ≥ 0.35 + modularity) are gone.
- **GSEA when the ORA is empty** (H. glaber case): the `GSEA_*` sheets of Deimos (manifest
  section `gsea_sheet`, schema 1.3; leading edge read from `Lead_proteins`, else matched
  through the Genes column) become the functional groups — pathways grouped by Jaccard on
  their leading-edge proteins, NES as direction. `--groups-from auto|ora|gsea`.
- **Stability of the response programmes**: replicates resampled within each condition,
  k-means re-run, clusters matched to the programmes (Hungarian on Jaccard); mean Jaccard
  per programme (Hennig's clusterboot thresholds: ≥ 0.75 stable, < 0.6 dissolved) and, per
  protein, the fraction of resamplings where it stays. `--prog-boot` (50; 0 = off).
- **Factorial designs** (Hordeum: 9 genotypes × T0 / T3 / T7 × light / dark = 45 groups).
  - Factors from `--factor-regex` (named groups; absent group = `—`), else `samples.tsv`
    columns, else the design file, else the condition names: time-like tokens form `time`
    wherever they sit, shorter names (`Tron_T0` beside `Tron_T3_Light`) are aligned on the
    vocabulary of each position, an unknown token stops the guess; aliased factors dropped.
    A bad regex is reported and the automatic reading takes over. If the separators give
    fewer than two factors, letters glued to a number are split (`Grp_Male4` -> sex x age,
    `embAviso24` -> cultivar x time); bare numbers are never called `time`. Contrasts match
    their conditions whatever the separator (`Grp_Male4` in the contrast, `Grp-Male4` in
    samples.tsv). Checked on a real H. glaber run (2 x 2, 6 contrasts).
  - Roles: x = the ordered factor (time), lines = the factor with most levels (genotype),
    panels = the rest (light); a level of x seen in one panel only (T0) is the shared start.
  - *Effects of each design factor*: comparisons grouped by the factor they vary, rows and
    columns = the context levels, cells = down | up counts on a shared scale.
  - Condition map coded by factor (colour = genotype, shape = light, size = time); programme
    heatmap ordered genotype > light > time; *Programmes by design factor* = interaction
    plots per programme (reference genotype in white, clickable legend).
  - *What drives the response*: share of each protein's variance carried by each factor and
    two-factor interaction (type-II sums of squares, sum-to-zero coding), higher-order
    condition effects and replicates; median per programme. Sheets `Variance_partition`,
    `Design_factors`, `var_*` columns; provenance `reading_factors`, `reading_factor_roles`.
- Synthesis workbook sheets `Reading_facts`, `Programmes` (with `stability` and `var_*`),
  `Programme_summary`, `Variance_partition`, `Design_factors`, `Groups`, `Group_terms`,
  `Group_x_contrast`, `Reading_columns`; `pathfinder_programmes.tsv`; summary section
  *Response programmes*; provenance keys `reading_*`.
- CLI: `--no-reading`, `--k`, `--volcano`, `--max-volcanos`, `--categorical`,
  `--order-regex`, `--reference-fasta`, `--groups-from`, `--prog-boot`, `--factor-regex`.
- No new dependency (no matplotlib, networkx or adjustText: layouts and label placement
  are done in numpy). Repository invariants (Python 3.9 syntax, English strings) extended
  to `reading/`.

## v1.2.0 — 2026-11-06

### Changed
- **Proteins selected = the Deimos ones.** The protein signature, the core / specific /
  reversed patterns, the WGCNA core set and the Publications axis now apply the volcano rule
  declared in the Deimos `_manifest`: statistic (`p.adj` or raw `p.val`), threshold and ratio,
  `p < thr & |log2FC| >= log2(ratio)` — the comparisons of Deimos' `volcano_status`. Counts
  now match the Deimos volcanos, UpSet and dashboard (previously `p.adj < 0.05` without a
  fold-change criterion, so the two counts differed). Structured manifest keys (Deimos schema
  1.2: `volcano_p_column`, `volcano_p_thresh`, `volcano_ratio_min`) are read first, else the
  `volcano_rule` text of schemas 1.0/1.1. Without a manifest (Phobos, older workbooks) the
  historical rule is kept.
- **GO / pathway terms at q < 0.05**, the threshold of the upstream gProfiler / STRING
  enrichment, decoupled from the protein rule (`--pthr` used to drive both).
- `--pthr` now defaults to the upstream rule; new `--min-ratio`. An explicit override is
  printed, written in the provenance (`significance_rule_matches_upstream: no`) and in the
  summary.

### Added
- `io_layer.SelectionRule` / `selection_rule()`; tests on the rule, on overrides, and a
  protein-by-protein check against the Deimos contract workbook.

## v1.1.1 — 2026-10

CI: the v1.1.0 `io_layer.py` / `species.py` (manifest reader, replicate index) were
missing from the pushed tree; actions bumped to checkout@v5 / setup-python@v6 (Node 24).

Corrections found on a real time-course × cultivar run (16 contrasts).

### Fixed
- **"Opposite regulation" was mostly an artefact of contrast orientation.** Any protein
  (or GO term) up in one contrast and down in another was called opposite. Between
  contrasts sharing a condition this only expresses an ordering (`A24_vs_A0` up and
  `A24_vs_B24` down means A0 < A24 < B24); between unrelated contrasts it compares
  different things. A reversal is now reported only between **comparable** contrasts:
  no shared condition and log2FC profiles correlated over all proteins (|r| ≥ 0.3,
  sign-aligned when one contrast is written the other way round). The others are
  labelled `mixed`. On the run that revealed it: 279 "opposite" proteins → 1 reversal
  and 278 mixed. New column `opposite_in` (the contrasts and fold-changes involved);
  rule stated in the provenance. Without fold-changes the previous rule is kept.
- **"Core" meant three different sets.** The headline card counts proteins significant
  in at least half the contrasts, the module enrichment tested the top 200 ranked
  proteins, and the verdict mixed the two ("the 1 recurrent proteins concentrate in
  module blue, 92/676"). The chain and verdict now name the set actually tested
  ("top 200 proteins ranked by recurrence then Pi-score").
- **KPI "Core proteins" showed the number of significant proteins**: renamed
  "Significant proteins".
- **Contrast names cut at 22 characters** in the summary chain: the font now shrinks to
  fit (cut only beyond 40 characters).

### Added
- **Protein annotation** (`axis_proteins.annotate`): for each significant protein, the
  significant GO/pathway terms it belongs to (`go_terms`, `n_go_terms`) and its WGCNA
  module (`module`), matched on accession or gene symbol; columns in Proteins_core and
  the Proteins tab.
- **Contrast families**: above 8 contrasts, the per-contrast detail groups contrasts
  whose GO signatures correlate at r ≥ 0.5 into families instead of one line each.
- **Hub descriptions**: the WGCNA hub table carries the protein description from the
  statistics workbook when available.
- **Quality: heterogeneous conditions.** The correlation check runs on all proteins and
  misses a condition whose replicates sit at different stages of the response
  (asynchronous development): r stayed at 0.95. The spread of each condition along the
  PCs carrying the design is now compared with the median condition (flag at ≥ 2.5×,
  about 1–2 % false alarms on homogeneous simulated designs), naming the replicates that
  lie closer to another condition.

## v1.1.0 — 2026-10

Deimos integration: Pathfinder can run as Deimos' last step, reading an explicit,
versioned contract instead of guessing.

### Added
- **Upstream manifest.** Deimos workbooks carry a hidden `_manifest` sheet (schema
  `deimos-results` 1.x). Pathfinder takes the GO sheet → contrast mapping from it: Deimos
  shortens sheet names (common affixes stripped, 12 characters per side, `#N` on
  collision) and name parsing recognised only 1/6 to 0/12 of them on real naming
  schemes. A contrast declared by the manifest but absent from the statistics is
  reported, never re-guessed. The organism is read from the manifest before any
  free-text scan; an unknown major schema version is reported loudly. Without a
  manifest (older workbooks, Phobos) nothing changes.
- `--go none`: disables Meta-GO **and** the GO auto-detection, which could otherwise
  pick up the GO workbook of another project in the folder.
- `--version`, and `__version__` in the provenance (was "Pathfinder v1").
- The summary states the significance rule (`p.adj < pthr`, no fold-change criterion)
  and, when the manifest provides it, the upstream volcano rule — so a report listing 3
  proteins and a dashboard listing 128 no longer look contradictory.
- Contract test on a real Deimos workbook (`tests/fixtures/`).

### Fixed
- **The replicate index is no longer tested as a design factor.** `replicate` / `rep`
  (1, 2, 3 inside each condition, always present in Deimos designs) was tested against
  every PC, doubling the BH family: on the Deimos example, condition ~ PC1 went from
  p.adj 0.039 to 0.078. It is dropped and recorded in the provenance.
- `normalise_org()` accepts `NCBI taxid 1304` / `taxid:1304` (Deimos STRING mode), so the
  species is found in workbooks without a manifest. A taxid quoted inside free text (a
  strain in a Methods detail) is deliberately not taken.
- README: CI badge pointed to a repository that does not exist; test count.

## v1.0 — 2026-07

First release.

### Fixed
- **Text columns are detected whatever pandas calls them** ('object' in pandas 2,
  'str' in pandas 3): the WGCNA protein-to-module assignment was silently lost on
  recent pandas.
- **Unnamed columns can no longer pose as samples.** GO sheets exported without a
  header row carry 'Unnamed: 3' columns; they were counted as sample candidates, so a
  GO sheet could compete against the real intensity matrix for the PCA.
- **GO sheet mapping tolerates spelling differences between the design and the stats
  columns.** Conditions are read from the design ('CAM203_CTR') while contrasts come
  from the workbook columns ('CAM203CTR_vs_...'); the reconstructed contrast name was
  compared character for character, so an underscore, a space or a different case
  rejected every sheet. Comparison is now normalised.
- `--diagnose` explains *why* a sheet is unmatched: it prints the known condition
  spellings and the contrasts, instead of leaving the mismatch invisible.
- **Design labels are matched to columns by containment** before falling back on row
  order: a design listing 'Sple-15' now resolves a column named
  'LFQ.intensity.Sple-15' (and vice versa). Ambiguity is settled by the closest
  length, so 'Sple-1' cannot steal 'Sple-10'. The positional fallback - which silently
  depends on row order - is now a last resort rather than the common case.
- **GO sheets that drop the token shared by every condition are now recognised.**
  When all conditions end (or start) with the same string, upstream naming may strip
  it: 'EPHPS11' becomes 'EPHP', 'CS11' becomes 'C'. Every sheet was then rejected and
  the whole Meta-GO axis silently skipped. Conditions are matched through unambiguous
  aliases; an alias that could designate two conditions is discarded, not guessed.
- **The design verdict now follows the p-value, not just the effect size.** With many
  groups and few replicates a large epsilon squared can coexist with a non-significant
  test; the summary used to print "group separation is clear" at p=0.069. It now
  reports a large-but-not-significant effect as a trend, and states the power limit.
- **QC outliers require a material gap** (>= 0.05 below the median), not just a
  statistical one: with replicates tightly grouped, 3*MAD flagged 0.02 differences.
- **Duplicate gene symbols no longer crash the export.** Real proteomes repeat gene
  symbols (isoforms, protein groups) and leave some unannotated; sorting the PCA
  loadings by reindexing that duplicated index raised `ValueError: cannot reindex on
  an axis with duplicate labels`, and silently mis-sorted otherwise. Sorting is now
  positional, and the test fixtures carry repeated and empty symbols like real data.
- **The PCA no longer runs on a pre-selected subset.** The loader used to grab the
  z-scored heatmap sheet, which typically holds only the significant or most-variable
  proteins: assessing group separation on proteins picked *because* they separate the
  groups is circular and inflates the variance explained and the effect sizes. Every
  sheet is now profiled and the best per-sample matrix wins — design-matching columns
  first, then raw (non z-scored) log2 intensities, then protein count. When only a
  z-scored matrix exists, the report says so and warns that magnitudes are inflated.
- **Figures are much larger**: wrapping each figure for the export buttons made the
  SVGs collapse to their default intrinsic size (~300 px). They now carry explicit
  dimensions and scale to the panel; heatmap cells, scatters and the report width
  were enlarged too.
- **PNG export now works**: the serialised SVG kept `height:auto` (resolving to a
  zero-height document) and declared `xmlns` twice (unparseable). The image is also
  passed as a data: URL instead of a blob:, whose null origin taints the canvas when
  the report is opened from file://.
- **GO sheet mapping** now matches contrasts exactly instead of scoring substrings.
  With nested condition names (e.g. 'HPS7' inside 'EPHPS7'/'IPHPS7'), sheets were
  silently assigned to the wrong contrast and merged, losing whole contrasts.
  Collisions, unmatched sheets and contrasts without GO are now reported.
- **Publications axis** no longer runs enrichment on tiny gene sets: contrasts with
  fewer than 5 significant proteins are skipped and reported (a 1-gene "enrichment"
  is just a list of papers citing that gene). Gene-list capping and PMID truncation
  are now visible, and the significance criterion is written into the report.
- **Contrast similarity** now requires at least 10 shared terms before making a
  claim: a correlation computed over a handful of terms is unstable and would
  otherwise put tiny signatures at the top of the ranking. The number of shared
  terms is reported alongside every r.
- **Contrast similarity** no longer presents mechanical correlation as biology:
  pairs sharing a condition (common or reversed reference) are flagged and excluded
  from the highlighted statements, which now also require |r| >= 0.3.

### Added
- **Summary redesigned around the key facts**: colour-coded verdict cards (design,
  quality, core proteins, core module, top function) and a visual chain answering the
  question the report is for - do the core proteins fall in one module, and which
  comparison does that module follow? A missing link is drawn as a break, never
  implied.
- **Module vs contrast test**: the eigengene is compared between the two groups of
  every contrast (Mann-Whitney + rank-biserial). The exact test has a hard floor with
  few replicates (3 vs 3 cannot go below p=0.1), which is now stated rather than
  reported as 'not significant'.
- **WGCNA axis**: the co-expression modules computed upstream are read (never
  recomputed) and crossed with the rest of the report - module eigengene vs condition
  (Kruskal-Wallis + epsilon squared), Fisher exact test locating the core protein
  signature across modules, and module hubs. New WGCNA tab and three sheets. 'grey',
  WGCNA's unassigned bin, is flagged and excluded from module counts.
- **Protein axis**: protein x contrast matrix of significant fold-changes, core
  signature (recurrent proteins), contrast-specific proteins and opposite regulation,
  ranked by recurrence then Pi-score, with the robustness score as a confidence
  column. New Proteins tab, `Proteins_core` / `Proteins_matrix` sheets, a summary
  section, and contrast clustering computed on proteins - to be compared with the
  GO-level one. The Pi-score orders proteins but never selects them: inclusion stays
  on the adjusted p-value.
- **Species resolution against STRING's own list** (~12k organisms) when the built-in
  table does not know the organism: scientific names ('Camelina sativa'), genus alone
  ('cannabis') and gProfiler codes all work, and the list is cached on disk. An
  ambiguous query lists its candidates and asks - 'csativa' fits Cannabis, Camelina
  and Castanea sativa - rather than picking one. Resolving through STRING also proves
  upfront that the Publications axis can run for that organism.
- **Test suite (78 tests) and CI** on Python 3.9-3.12. Each fixed bug is pinned by a
  test: nested condition names in GO sheet names, matrix selection, mechanical
  contrast correlation, thin signatures, tiny publication gene sets, figure sizing,
  SVG/PNG export, French leftovers, cross-version syntax.
- **Meta-GO axis** — term × contrast matrix of signed z-scores, cross-directionality,
  Jaccard redundancy grouping, split by source (BP / MF / CC / KEGG / Reactome).
- **Exploration axis** — PCA (scree, loadings), PC ↔ condition association
  (Kruskal-Wallis + ε², BH-corrected), inter-sample correlation.
- **Multi-factor design testing** — every design column tested against every PC;
  categorical (Kruskal-Wallis + ε²) and numeric (Spearman) factors; technical
  factors driving a high-variance PC are flagged as batch effects.
- **Protein UMAP** — interactive, colourable by abundance, per-contrast LFC and
  GO term, with source / term-count selectors.
- **Publications axis** — STRING PMID enrichment per contrast with a dated cache
  snapshot. The species is auto-detected in the workbook (Deimos `go_organism`)
  and confirmed by an interactive prompt at launch; `--org` skips it, `--no-online`
  disables the axis. Accepts 'hsapiens', 'Homo sapiens', 'human' or an NCBI taxid.
- **Interpretive summary** — deterministic, exported to HTML, Excel and Markdown.
- **Contrast clustering** on the GO signature (Pearson on shared terms, hierarchical
  order) with a contrast x contrast heatmap; Meta-GO heatmaps ordered accordingly.
- Graduated z-score colour bars under every Meta-GO heatmap, and a
  **Representatives only / All terms** switch collapsing redundant GO terms.
- Adaptive inter-sample correlation scale with automatic detection of z-scored
  input (absolute r is low by construction) and a robust outlier rule.
- **SVG / PNG export buttons on every figure** (including the GO UMAP canvas),
  with filenames derived from the figure title.
- Self-contained HTML report and multi-sheet synthesis workbook.
- Content-based detection of input workbooks and columns; experimental design
  read by name with positional fallback; design template generator.
