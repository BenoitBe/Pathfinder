# Changelog

## v1.0 — 2026-07

First release.

### Added
- **Meta-GO axis** — term × contrast matrix of signed z-scores, cross-directionality,
  Jaccard redundancy grouping, split by source (BP / MF / CC / KEGG / Reactome).
- **Exploration axis** — PCA (scree, loadings), PC ↔ condition association
  (Kruskal-Wallis + ε², BH-corrected), inter-sample correlation.
- **Multi-factor design testing** — every design column tested against every PC;
  categorical (Kruskal-Wallis + ε²) and numeric (Spearman) factors; technical
  factors driving a high-variance PC are flagged as batch effects.
- **Protein UMAP** — interactive, colourable by abundance, per-contrast LFC and
  GO term, with source / term-count selectors.
- **Publications axis** (opt-in) — STRING PMID enrichment per contrast with a
  dated cache snapshot.
- **Interpretive summary** — deterministic, exported to HTML, Excel and Markdown.
- Graduated z-score colour bars under every Meta-GO heatmap, and a
  **Representatives only / All terms** switch collapsing redundant GO terms.
- Adaptive inter-sample correlation scale with automatic detection of z-scored
  input (absolute r is low by construction) and a robust outlier rule.
- Self-contained HTML report and multi-sheet synthesis workbook.
- Content-based detection of input workbooks and columns; experimental design
  read by name with positional fallback; design template generator.
