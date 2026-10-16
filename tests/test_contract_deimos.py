"""Contract test against a REAL Deimos workbook (tests/fixtures/deimos_contract.xlsx).

Synthetic fixtures check Pathfinder's logic; this one checks that Pathfinder still
reads what Deimos actually writes. When Deimos renames a sheet or a column it bumps
its manifest schema version, and the fixture is regenerated with
tests/fixtures/make_deimos_fixture.py.
"""
from pathlib import Path

import pandas as pd
import pytest

import axis_wgcna
import io_layer
import pathfinder
import species

FIX = Path(__file__).with_name("fixtures")
WB = str(FIX / "deimos_contract.xlsx")
DESIGN = str(FIX / "deimos_contract_design.csv")


@pytest.fixture(scope="module")
def manifest():
    return io_layer.read_manifest(WB)


@pytest.fixture(scope="module")
def stats():
    return io_layer.load_stats(WB, design_map=io_layer.load_design(DESIGN))


def test_manifest_schema_is_supported(manifest):
    assert manifest["schema"]["name"] == "deimos-results"
    assert io_layer.manifest_version_warning(manifest) is None


def test_declared_sheets_and_columns_exist(manifest):
    sheets = pd.read_excel(WB, sheet_name=None)
    comp = sheets[manifest["stats"]["sheet"]]
    for c in manifest["contrast"]:
        for col in (f"{c}_diff", f"{c}_p.val", f"{c}_p.adj",
                    f"Pi_Score_{c}", f"Robustness_Score_{c}"):
            assert col in comp.columns, col
    assert manifest["matrix"]["sheet"] in sheets
    for sheet in manifest["go_sheet"]:
        assert sheet in sheets
        assert {"term_id", "term_name", "source", "p_value", "intersection"} <= set(sheets[sheet].columns)


def test_stats_read_as_declared(stats, manifest):
    assert stats.contrasts == list(manifest["contrast"])
    assert stats.pca_source == manifest["matrix"]["sheet"]
    assert stats.conditions_resolved
    assert sorted(stats.conditions) == sorted({"Live_K12", "Heat_killed_K12",
                                               "Live_M18", "Heat_killed_M18"})


def test_every_go_sheet_mapped_through_the_manifest(stats):
    go = io_layer.load_go(WB, stats.contrasts, conditions=stats.conditions, verbose=False)
    assert go.unmatched == [] and go.duplicates == {} and go.missing == []
    assert set(go.mapping_how.values()) == {"manifest"}


def test_wgcna_sheets_are_read(stats):
    res = axis_wgcna.run(pd.read_excel(WB, sheet_name=None), stats)
    assert res and res["n_modules"] >= 2


def test_species_from_manifest():
    assert species.detect_org_in_workbook(WB) == "hsapiens"


def test_end_to_end(tmp_path):
    out = tmp_path / "out"
    pathfinder.main(["--stats", WB, "--go", WB, "--design", DESIGN, "--outdir", str(out),
                     "--no-online", "--no-umap"])
    syn = pd.read_excel(out / "pathfinder_synthesis.xlsx", sheet_name=None)
    prov = syn["provenance"].set_index("key")["value"]
    assert prov["go_sheets_mapped"] == "6/6 (6 via the upstream manifest)"
    assert prov["upstream_manifest"].startswith("deimos-results 1.")
    assert "volcano_rule_upstream" in prov.index
    contrast_cols = [c for c in syn["MetaGO_matrix"].columns if "_vs_" in str(c)]
    assert len(contrast_cols) == 6
    summary = (out / "pathfinder_summary.md").read_text(encoding="utf-8")
    assert "Significance rule: p.adj < 0.05" in summary
