"""Upstream manifest (hidden '_manifest' sheet written by Deimos >= schema 1.0).

The manifest is the contract between the Deimos workbook and Pathfinder: GO sheet ->
contrast is read from it instead of being reconstructed from shortened sheet names,
and the organism the analysis ran on comes from it before any free-text scan.
"""
import openpyxl
import pandas as pd
import pytest

import io_layer
import species

# Deimos shortens GO sheet names (common affixes stripped, 12 characters per side,
# '#N' on collision): names Pathfinder cannot parse back. Real output of
# go_enrichment._make_go_sheet_names() (Deimos 8e2faff) for three intra-strain contrasts.
STRAINS = ["K12", "M18", "JIM8777"]
CONTRASTS = [f"Live_{s}_vs_Heat_killed_{s}" for s in STRAINS]
SHEETS = {"Live_K12_vs_Heat_killed_K12": "GO_Live_K12vHeat_kil",
          "Live_M18_vs_Heat_killed_M18": "GO_Live_M18vHeat_k#2",
          "Live_JIM8777_vs_Heat_killed_JIM8777": "GO_Live_JIMvHeat_k#3"}
CONDITIONS = ["Live_K12", "Heat_killed_K12", "Live_M18", "Heat_killed_M18",
              "Live_JIM8777", "Heat_killed_JIM8777"]


def _write(path, sheets=SHEETS, manifest=True, version="1.0", extra_rows=()):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Differential_Expression"
    ws.append(["Protein.Group"] + [f"{c}_{s}" for c in CONTRASTS for s in ("diff", "p.val", "p.adj")])
    for k, (c, sheet) in enumerate(sheets.items()):
        g = wb.create_sheet(sheet)
        # Deimos writes the table from column U (21), plots on the left
        g.cell(row=1, column=21, value="source")
        for j, col in enumerate(["term_id", "term_name", "p_value", "intersection"], start=22):
            g.cell(row=1, column=j, value=col)
        for i in range(2 + k):                              # distinct count per contrast
            g.cell(row=2 + i, column=21, value="GO:BP")
            g.cell(row=2 + i, column=22, value=f"GO:{k}{i:06d}")
            g.cell(row=2 + i, column=23, value=f"term {k}.{i}")
            g.cell(row=2 + i, column=24, value=1e-4)
            g.cell(row=2 + i, column=25, value="P1,P2")
    if manifest:
        m = wb.create_sheet("_manifest")
        m.append(("section", "key", "value"))
        m.append(("schema", "name", "deimos-results"))
        m.append(("schema", "version", version))
        for c in CONTRASTS:
            m.append(("contrast", c, c.replace("_vs_", "|")))
        for c, sheet in sheets.items():
            m.append(("go_sheet", sheet, c))
        for r in extra_rows:
            m.append(r)
        m.sheet_state = "hidden"
    wb.save(path)
    return str(path)


class TestReadManifest:
    def test_absent_is_empty(self, tmp_path):
        assert io_layer.read_manifest(_write(tmp_path / "a.xlsx", manifest=False)) == {}

    def test_from_path_and_from_sheet_dict(self, tmp_path):
        p = _write(tmp_path / "a.xlsx")
        by_path = io_layer.read_manifest(p)
        by_dict = io_layer.read_manifest(pd.read_excel(p, sheet_name=None))
        assert by_path == by_dict
        assert by_path["schema"] == {"name": "deimos-results", "version": "1.0"}
        assert by_path["go_sheet"]["GO_Live_M18vHeat_k#2"] == "Live_M18_vs_Heat_killed_M18"

    def test_hidden_sheet_is_read(self, tmp_path):
        p = _write(tmp_path / "a.xlsx")
        assert openpyxl.load_workbook(p)["_manifest"].sheet_state == "hidden"
        assert "_manifest" in pd.read_excel(p, sheet_name=None)

    @pytest.mark.parametrize("version,warns", [("1.0", False), ("1.7", False),
                                               ("2.0", True), ("", True), ("beta", True)])
    def test_version_check(self, tmp_path, version, warns):
        m = io_layer.read_manifest(_write(tmp_path / "a.xlsx", version=version))
        assert (io_layer.manifest_version_warning(m) is not None) == warns

    def test_no_manifest_no_warning(self):
        assert io_layer.manifest_version_warning({}) is None


class TestGoMappingFromManifest:
    def test_name_parsing_alone_fails_on_these_sheets(self, tmp_path):
        """The problem the manifest solves: shortened names are not parseable."""
        p = _write(tmp_path / "a.xlsx", manifest=False)
        go = io_layer.load_go(p, CONTRASTS, conditions=CONDITIONS, verbose=False)
        assert len(go.unmatched) >= 2

    def test_every_sheet_mapped_via_manifest(self, tmp_path):
        p = _write(tmp_path / "a.xlsx")
        go = io_layer.load_go(p, CONTRASTS, conditions=CONDITIONS, verbose=False)
        assert go.unmatched == [] and go.duplicates == {} and go.missing == []
        assert set(go.mapping_how.values()) == {"manifest"}
        assert go.from_manifest == 3
        counts = go.long.groupby("contrast").size().to_dict()
        assert counts == {c: 2 + k for k, c in enumerate(SHEETS)}

    def test_declared_contrast_absent_from_stats_is_not_guessed(self, tmp_path):
        bad = dict(SHEETS)
        p = _write(tmp_path / "a.xlsx", sheets=bad)
        go = io_layer.load_go(p, CONTRASTS[:2], conditions=CONDITIONS, verbose=False)
        sheet = SHEETS[CONTRASTS[2]]
        assert go.mapping[sheet] is None
        assert "manifest declares" in go.mapping_how[sheet]

    def test_sheet_missing_from_manifest_falls_back_on_its_name(self, tmp_path):
        p = _write(tmp_path / "a.xlsx", sheets={**SHEETS})
        wb = openpyxl.load_workbook(p)
        g = wb.create_sheet("GO_Live_K12_vs_Heat_killed_K12")     # unlisted, parseable
        g.append(["source", "term_id", "term_name", "p_value", "intersection"])
        g.append(["GO:BP", "GO:9", "t", 1e-3, "P1"])
        del wb["GO_Live_K12vHeat_kil"]
        wb.save(p)
        go = io_layer.load_go(p, CONTRASTS, conditions=CONDITIONS, verbose=False)
        assert go.mapping["GO_Live_K12_vs_Heat_killed_K12"] == CONTRASTS[0]
        assert go.mapping_how["GO_Live_K12_vs_Heat_killed_K12"] == "exact"


class TestSpeciesFromManifest:
    def test_manifest_taxid_wins_over_methods(self, tmp_path):
        p = _write(tmp_path / "a.xlsx", extra_rows=[("organism", "taxid", "1304"),
                                                    ("organism", "string_mode", "ortholog")])
        wb = openpyxl.load_workbook(p)
        ws = wb.create_sheet("Methods", 0)
        ws.append(["Item", "Choice / Value", "Method / Detail"])
        ws.append(["Species", "hsapiens", "g:Profiler organism code"])   # stale/other value
        wb.save(p)
        assert species.detect_org_in_workbook(p) == "1304"

    def test_manifest_gprofiler_code(self, tmp_path):
        p = _write(tmp_path / "a.xlsx", extra_rows=[("organism", "gprofiler", "bnapus")])
        assert species.detect_org_in_workbook(p) == "bnapus"

    def test_no_organism_in_manifest_falls_back_on_methods(self, tmp_path):
        p = _write(tmp_path / "a.xlsx")
        wb = openpyxl.load_workbook(p)
        ws = wb.create_sheet("Methods", 0)
        ws.append(["Item", "Choice / Value", "Method / Detail"])
        ws.append(["Species", "NCBI taxid 1304", "Streptococcus salivarius (taxid 1200793)"])
        wb.save(p)
        assert species.detect_org_in_workbook(p) == "1304"
