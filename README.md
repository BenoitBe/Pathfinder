# Pathfinder — Pathway And Term Harmonisation For INtegrated Design Exploration & Reporting

![Pathfinder](pathfinder_icon.svg)

[![tests](https://github.com/BenoitBe/pathfinder-synthesis-pipeline/actions/workflows/tests.yml/badge.svg)](https://github.com/BenoitBe/pathfinder-synthesis-pipeline/actions/workflows/tests.yml)

**Downstream companion to Deimos and Phobos.**
Cross-comparison synthesis and unsupervised exploration — from the Excel report
of an upstream pipeline to a self-contained interactive HTML report.

Pathfinder **recomputes no statistics**: it reads the workbooks produced upstream,
reshapes them and visualises them to deliver a whole-project view across all
contrasts, plus an automatic interpretive summary.

---

## Pathfinder vs Deimos / Phobos

| Aspect      | Deimos (DIA/DIA-NN)         | Phobos (DDA/PEAKS)          | Pathfinder (downstream)                 |
| ----------- | --------------------------- | --------------------------- | --------------------------------------- |
| Input       | `report.pg_matrix.tsv`      | `protein-peptides.csv`      | Excel report from Deimos/Phobos         |
| Unit        | Protein (LFQ)               | Precursor peptide           | Contrast × GO term, sample × PC         |
| Statistics  | limma + DEqMS               | limma eBayes                | None — reads upstream results           |
| Focus       | Differential analysis       | Peptide-level analysis      | Cross-contrast synthesis, design QC     |
| Network     | Offline                     | Offline                     | Offline (Publications axis = optional)  |

---

## Features

- **Meta-GO** — term × contrast matrix of signed z-scores, cross-directionality
  (shared up/down vs opposite regulation), Jaccard redundancy grouping with a
  **Non-redundant / All terms** switch. Split by source (BP / MF / CC / KEGG / Reactome).
- **Contrast clustering** — contrasts correlated on their GO signature and ordered
  hierarchically: which comparisons move the same functions the same way, which are
  opposed. Pairs sharing a condition are flagged (a common reference inflates r
  mechanically) and excluded from the headline statements. The Meta-GO heatmaps
  follow that order.
- **Exploration** — PCA with scree, loadings, and **PC ↔ condition association**
  (Kruskal-Wallis + ε², BH-corrected); inter-sample correlation.
- **Multi-factor design testing** — every design column (condition, batch,
  injection order…) is tested against every PC. Categorical → Kruskal-Wallis + ε²;
  numeric → Spearman. Technical factors driving a high-variance PC are flagged.
- **Protein UMAP** — interactive, colourable by mean abundance, per-contrast LFC,
  and **GO term** (source and term-count selectors).
- **Publications** — STRING **PMID enrichment** per contrast: which papers are
  over-represented in the signature, with a dated cache snapshot. Contrasts with
  fewer than 5 significant proteins are skipped (enrichment describes a set, not a
  gene); the criterion, gene-list cap and PMID truncation are reported. The species is
  detected in the workbook and confirmed interactively at launch (`--org` to skip
  the prompt, `--no-online` to stay fully offline).
- **Automatic interpretive summary** — design capture, batch suspicion, outlier
  samples, functional synthesis per source. Deterministic and reproducible.

---

## Input files

| File                       | Required | Description                                                       |
| -------------------------- | -------- | ----------------------------------------------------------------- |
| Upstream stats workbook    | ✅        | Any `.xlsx` with `<contrast>_p.val` columns + per-sample intensities |
| GO enrichment workbook     | ⬜        | `GO_*` sheets with `term_name`, `source`, `p_value`, `intersection` |
| `ExperimentalDesign.csv`   | ⬜        | `sample` + `condition` (+ any extra factors)                        |

Files and columns are detected **by content**, not by filename — the GO sheets may
live in the same workbook as the statistics.

The design is optional but **required to colour the PCA and test PC ↔ condition**
when sample columns do not encode the condition (e.g. `Sple-1 … Sple-36`).
Labels are matched by name, falling back to positional order if the design lists
raw filenames. Generate a pre-filled template:

```
python pathfinder.py --write-design-template
# -> ExperimentalDesign_template.csv (fill in the condition column)
```

---

## Installation

Requires **Python 3.9+**.

```
pip install -r requirements.txt
```

Or manually:

```
pip install pandas numpy scipy scikit-learn openpyxl umap-learn requests
```

`umap-learn` (UMAP views) and `requests` (Publications axis) are optional —
without them the corresponding modules are skipped cleanly.

---

## Usage

```
# Offline — auto-detects the workbooks in the current folder
python pathfinder.py

# Explicit files
python pathfinder.py --stats ProteomicAnalysis_Results.xlsx --go GO_Results.xlsx

# With an experimental design (multi-factor supported)
python pathfinder.py --design ExperimentalDesign.csv

# Publications axis — Pathfinder detects the species in the workbook and asks to
# confirm at launch (Enter to accept, 'n' to skip). To answer up-front:
python pathfinder.py --org rnorvegicus        # gProfiler code
python pathfinder.py --org "Camelina sativa"  # scientific name (resolved via STRING)
python pathfinder.py --org 90675              # NCBI taxid, always unambiguous
python pathfinder.py --org 10116 --refresh # bypass the STRING cache
python pathfinder.py --no-online           # fully offline: no prompt, no query

# Inspect what was detected, then exit
python pathfinder.py --diagnose
```

### Options

| Option                    | Role                                                       |
| ------------------------- | ---------------------------------------------------------- |
| `--stats`                 | Statistics workbook (auto-detected by content if omitted)   |
| `--go`                    | GO enrichment workbook (may be the same file)               |
| `--design`                | Sample → condition mapping (+ extra factors)                |
| `--outdir`                | Output directory                                            |
| `--pthr`                  | p threshold masking the GO matrix (default 0.05)            |
| `--jaccard`               | Cut distance for GO redundancy grouping (default 0.5)       |
| `--org`                   | Species: gProfiler code, scientific name, or NCBI taxid     |
| `--no-online`             | Fully offline: no species prompt, no STRING query           |
| `--refresh`               | Bypass the STRING cache                                     |
| `--no-umap`               | Disable UMAP views                                          |
| `--diagnose`              | Report detected sheets/columns and group sizes              |
| `--write-design-template` | Write a pre-filled design template                          |

---

## Tests

```
pip install pytest
pytest
```

78 tests covering the loader (GO sheet mapping, per-sample matrix selection, design
resolution), the three axes, the interpretive summary and the HTML report, plus
repository invariants (Python 3.9-3.12 syntax, no untranslated strings). Most of them
encode a regression that actually happened on production data — nested condition names
silently merging contrasts, a PCA running on a pre-selected subset, correlation
reported as biology when it was mechanical, figures collapsing to 300 px. CI runs the
suite on Python 3.9 to 3.12 plus a smoke test on the synthetic dataset.

## Repository structure

```
pathfinder/
├── pathfinder.py              # Main orchestrator + CLI
├── io_layer.py                # Excel + design readers → shared objects
├── axis_meta_go.py            # Meta-GO axis (term × contrast)
├── axis_explore.py            # Exploration axis (PCA, multi-factor, UMAP)
├── axis_publications.py       # Publications axis (STRING PMID, online)
├── summary.py                 # Deterministic interpretive summary
├── species.py                 # gProfiler organism → NCBI taxid
├── cache.py                   # Disk cache for online calls
├── build_report_html.py       # Self-contained HTML report generator
├── build_report_xlsx.py       # Synthesis workbook generator
├── generate_icon.py           # Icon generator
├── pathfinder_icon.svg        # Pathfinder icon (README banner)
├── pathfinder_mark.svg        # Pathfinder mark (report navigation)
├── requirements.txt
├── LICENSE
├── example_data/              # Minimal synthetic dataset generator
│   └── make_synthetic.py
└── tests/                     # pytest suite (regression-driven)
    ├── conftest.py            # synthetic workbooks mirroring real-world shapes
    ├── test_io_layer.py
    ├── test_axes.py
    ├── test_summary.py
    └── test_report.py
```

---

## Output

```
pathfinder_output/
├── pathfinder_report.html      # Self-contained interactive report
├── pathfinder_synthesis.xlsx   # Multi-sheet synthesis workbook
└── pathfinder_summary.md       # Interpretive summary (standalone)
```

**HTML tabs**: Summary · Overview · Meta-GO · Exploration · Protein UMAP · Publications

Every figure carries **SVG** and **PNG** export buttons (hover over it): vector SVG for
Illustrator/Inkscape, 2× PNG for slides. Filenames are derived from the figure title.

**Excel sheets**: provenance · Summary · MetaGO_matrix · MetaGO_directionality ·
Contrast_similarity · PCA_scores · PCA_loadings · PC_factor_assoc · PC_factor_multi ·
Sample_correlation · PCA_scree · Publications

---

## Example

```
cd example_data && python make_synthetic.py
python ../pathfinder.py --design ExperimentalDesign.csv
```

---

## Methodological notes

**No recomputation.** Every value derives from the upstream workbooks. The GO
matrix stores the signed z-score masked by significance — `NaN` (not tested /
not significant), never `0`, which would distort clustering.

**GO redundancy without the ontology.** Terms are grouped by Jaccard distance on
their `intersection` accession sets: functionally identical terms collapse
without downloading the ontology, keeping the tool fully offline.

**PC ↔ factor association.** Kruskal-Wallis with ε² effect size (categorical) or
Spearman (numeric, e.g. injection order), BH-corrected across all tests. A
high-variance PC associated with a technical factor points to a batch effect —
consider `limma::removeBatchEffect`, ComBat, or adding it as a covariate.
Without replicates the test is reported as *not computable* rather than
misinterpreted.

**Matrix selection.** Every sheet is profiled and the best per-sample matrix is used
for the PCA: columns matching the design first, then raw (non z-scored) log2
intensities, then protein count. A z-scored heatmap is usually a *selected* subset
(significant or most-variable proteins), and running the PCA on it is circular — it
inflates both the variance explained and the effect sizes. When no raw matrix is
available Pathfinder still runs, but says so and flags the magnitudes as unreliable.

**QC correlation** follows the same matrix: per-protein z-scores are decorrelated by
construction (absolute r is meaningless, only relative structure is), so the colour
scale adapts and the outlier rule switches to a robust, relative criterion.

**UMAP guard-rail.** Sample-level UMAP requires `n ≥ 4 × n_conditions`; below
that it fabricates convincing but arbitrary clusters, so PCA remains the
reference for design assessment. Protein-level UMAP is always available.

**Publications reproducibility.** The literature corpus evolves, so every STRING
call is cached with its fetch date and logged in the provenance sheet. Known
bias: heavily studied proteins are over-represented — input genes are ranked by
|LFC| to attenuate it.

---

## Name

**Pathfinder** — the Mars lander that explored the surface those two moons orbit.
Acronym: **P**athway **A**nd **T**erm **H**armonisation **F**or **IN**tegrated
**D**esign **E**xploration & **R**eporting.
Downstream synthesis layer for the Deimos (DIA) and Phobos (DDA) suites.

---

## Related

- [deimos-proteomics-pipeline](https://github.com/BenoitBe/deimos-proteomics-pipeline) — DIA Expression Integrated Multi-Omics Suite
- [phobos-peptidomics-pipeline](https://github.com/BenoitBe/phobos-peptidomics-pipeline) — Peptide-level Heuristics for Bottom-up Omics Suite

---

## License

MIT — see [LICENSE](LICENSE).
