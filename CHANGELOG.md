# Changelog

## v1.0 — 2026-07

First release.

### Fixed
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
