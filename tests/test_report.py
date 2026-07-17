"""HTML report and end-to-end runs, plus repository-wide invariants."""
import ast
import glob
import os
import re
import sys
from pathlib import Path

import pandas as pd
import pytest

import io_layer
import pathfinder

ROOT = Path(__file__).resolve().parents[1]


def _run(tmp_path, *args, umap=False):
    """UMAP is slow and orthogonal to what these tests check: off unless asked."""
    out = tmp_path / "out"
    extra = [] if umap else ["--no-umap"]
    pathfinder.main(["--outdir", str(out), "--no-online", *extra, *args])
    return out


class TestEndToEnd:
    def test_full_run_produces_every_output(self, tmp_path, workbook_full, go_workbook):
        out = _run(tmp_path, "--stats", workbook_full["path"], "--go", go_workbook["path"],
                   "--design", workbook_full["design"])
        for name in ("pathfinder_report.html", "pathfinder_synthesis.xlsx",
                     "pathfinder_summary.md"):
            assert (out / name).exists(), f"{name} missing"

    def test_workbook_sheets(self, tmp_path, workbook_full, go_workbook):
        out = _run(tmp_path, "--stats", workbook_full["path"], "--go", go_workbook["path"],
                   "--design", workbook_full["design"])
        sheets = pd.read_excel(out / "pathfinder_synthesis.xlsx", sheet_name=None)
        for s in ("provenance", "Summary", "MetaGO_matrix", "MetaGO_directionality",
                  "Contrast_similarity", "PCA_scores", "PC_factor_assoc",
                  "PC_factor_multi", "Sample_correlation", "PCA_scree"):
            assert s in sheets, f"sheet {s} missing"

    def test_provenance_records_the_matrix_and_mapping(self, tmp_path, workbook_full, go_workbook):
        out = _run(tmp_path, "--stats", workbook_full["path"], "--go", go_workbook["path"],
                   "--design", workbook_full["design"])
        prov = pd.read_excel(out / "pathfinder_synthesis.xlsx",
                             sheet_name="provenance").set_index("key")["value"]
        assert "Log2_Impute" in prov["pca_matrix"]
        assert prov["go_sheets_mapped"] == "28/28"
        assert prov["conditions_resolved"] == "yes"

    def test_runs_without_go(self, tmp_path, workbook_full):
        out = _run(tmp_path, "--stats", workbook_full["path"], "--design", workbook_full["design"])
        assert (out / "pathfinder_report.html").exists()

    def test_runs_without_design(self, tmp_path, workbook_zscore_only):
        out = _run(tmp_path, "--stats", workbook_zscore_only["path"])
        assert (out / "pathfinder_report.html").exists()

    def test_unreadable_go_does_not_crash(self, tmp_path, workbook_full):
        bad = tmp_path / "bad.xlsx"
        bad.write_text("not an excel file")
        out = _run(tmp_path, "--stats", workbook_full["path"], "--go", str(bad))
        assert (out / "pathfinder_report.html").exists()

    def test_duplicate_gene_symbols_do_not_break_loadings(self, tmp_path, workbook_full):
        """REGRESSION: real proteomes repeat gene symbols (isoforms, protein groups);
        sorting the loadings by reindexing a duplicated index raises ValueError."""
        out = _run(tmp_path, "--stats", workbook_full["path"],
                   "--design", workbook_full["design"])
        load = pd.read_excel(out / "pathfinder_synthesis.xlsx", sheet_name="PCA_loadings",
                             index_col=0)
        assert len(load) > 0
        assert load.index.duplicated().any(), "fixture must contain duplicated gene symbols"
        # and the sort must be correct despite the duplicates
        pc1 = load["PC1"].abs().to_numpy()
        assert (pc1[:-1] >= pc1[1:] - 1e-12).all(), "loadings are not sorted by |PC1|"

    def test_design_template(self, tmp_path, workbook_full, capsys):
        out = tmp_path / "tpl"
        pathfinder.main(["--stats", workbook_full["path"], "--outdir", str(out),
                         "--no-online", "--write-design-template"])
        tpl = out / "ExperimentalDesign_template.csv"
        assert tpl.exists()
        d = pd.read_csv(tpl)
        assert list(d.columns) == ["sample", "condition"]
        assert len(d) == 24


class TestHtmlReport:
    @pytest.fixture
    def html(self, tmp_path, workbook_full, go_workbook):
        out = _run(tmp_path, "--stats", workbook_full["path"], "--go", go_workbook["path"],
                   "--design", workbook_full["design"])
        return (out / "pathfinder_report.html").read_text(encoding="utf-8")

    def test_tabs_present(self, html):
        for tab in ("Summary", "Overview", "Meta-GO", "Exploration", "Protein UMAP"):
            assert f">{tab}<" in html

    def test_tab_switching_uses_safe_ids(self, html):
        """REGRESSION: onclick=show('Vue d'ensemble') broke on the apostrophe."""
        assert 'onclick="show(0)"' in html
        assert "id='sec-1'" in html

    def test_figures_carry_intrinsic_size(self, html):
        """REGRESSION: wrapping figures collapsed them to ~300 px."""
        svgs = re.findall(r"<svg viewBox='0 0 (\d+) (\d+)' width='(\d+)' height='(\d+)'", html)
        assert svgs, "no sized SVG found"
        for vw, vh, w, h in svgs:
            assert vw == w and vh == h

    def test_no_percentage_width_on_svg_root(self, html):
        assert "<svg viewBox='0 0" in html
        assert "width='100%'" not in html

    def test_export_script_present(self, html):
        assert "exportSVG" in html and "exportPNG" in html
        assert "svgDataUrl" in html          # data: URL, not blob: (canvas tainting)

    def test_export_removes_style_before_serialising(self, html):
        """REGRESSION: 'height:auto' serialised to a zero-height document."""
        assert "removeAttribute('style')" in html

    def test_export_does_not_duplicate_xmlns(self, html):
        """REGRESSION: a hand-set xmlns plus the serializer's own = invalid XML."""
        script = html.split("EXPORT")[-1]
        assert "setAttribute('xmlns'" not in html

    def test_go_widget_is_interactive(self, tmp_path, workbook_full, go_workbook):
        """The GO colouring widget lives in the protein UMAP tab, so this one needs UMAP."""
        out = _run(tmp_path, "--stats", workbook_full["path"], "--go", go_workbook["path"],
                   "--design", workbook_full["design"], umap=True)
        h = (out / "pathfinder_report.html").read_text(encoding="utf-8")
        assert "id='go_src'" in h and "id='go_k'" in h
        assert "id='go_cv'" in h                       # canvas
        assert '"coords"' in h and '"terms"' in h      # embedded data

    def test_meta_go_view_switch(self, html):
        assert "id='go_view'" in html
        assert "Non-redundant" in html and "All terms" in html

    def test_colourbars_have_unique_ids(self, html):
        ids = re.findall(r"linearGradient id='(\w+)'", html)
        assert len(ids) == len(set(ids)), "duplicate gradient ids collide"

    def test_report_is_self_contained(self, html):
        assert "http://" not in html.replace("http://www.w3.org", "")
        assert "cdn" not in html.lower()

    def test_lang_is_english(self, html):
        assert "lang='en'" in html


class TestRepositoryInvariants:
    def test_syntax_is_compatible_with_python_39_to_312(self):
        """REGRESSION: an f-string backslash only parses from 3.12 on."""
        for f in glob.glob(str(ROOT / "*.py")) + glob.glob(str(ROOT / "example_data/*.py")):
            src = Path(f).read_text(encoding="utf-8")
            for ver in (9, 10, 11, 12):
                try:
                    ast.parse(src, feature_version=(3, ver))
                except SyntaxError as e:
                    pytest.fail(f"{os.path.basename(f)} is not valid on Python 3.{ver}: {e.msg}")

    def test_no_french_strings_in_user_facing_code(self):
        accents = re.compile(r"[éèêàùçôîœ]")
        allowed = re.compile(r"Université|Proteogen|epsilon|ε|²")
        for f in glob.glob(str(ROOT / "*.py")):
            for i, line in enumerate(Path(f).read_text(encoding="utf-8").splitlines(), 1):
                if accents.search(line) and not allowed.search(line):
                    pytest.fail(f"{os.path.basename(f)}:{i} still in French: {line.strip()[:70]}")

    def test_species_normalisation(self):
        import species
        assert species.normalise_org("hsapiens")[0] == 9606
        assert species.normalise_org("Homo sapiens")[0] == 9606
        assert species.normalise_org("human")[0] == 9606
        assert species.normalise_org("10116")[0] == 10116
        assert species.normalise_org("bnapus")[0] == 3708
        assert species.normalise_org("not-a-species")[0] is None
