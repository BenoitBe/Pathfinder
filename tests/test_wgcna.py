"""WGCNA axis: modules read from upstream, crossed with the design and the core signature."""
import itertools

import numpy as np
import pandas as pd
import pytest

import axis_proteins
import axis_wgcna
import io_layer


@pytest.fixture
def workbook_wgcna(tmp_path):
    """Known structure: module 'blue' both tracks the condition and holds the core
    signature; the other modules are noise; 'grey' is the unassigned bin."""
    rng = np.random.default_rng(41)
    conds = ["A", "B", "C", "D"]
    contrasts = [f"{a}_vs_{b}" for a, b in itertools.combinations(conds, 2)]
    samples = [f"Sple-{i + 1}" for i in range(12)]
    n = 500

    comp = pd.DataFrame({"Protein.Group": [f"P{i:04d}" for i in range(n)],
                         "Genes": [f"G{i:04d}" for i in range(n)],
                         "First.Protein.Description": [f"Synthetic protein {i} (test)"
                                                        for i in range(n)]})
    for c in contrasts:
        p, d, pi = np.ones(n), np.zeros(n), np.zeros(n)
        p[:30] = 1e-8                                   # core = P0000..P0029
        d[:30] = rng.normal(2, 0.3, 30)
        pi[:30] = rng.uniform(10, 20, 30)
        comp[f"{c}_p.adj"] = p
        comp[f"{c}_p.val"] = p
        comp[f"{c}_diff"] = d
        comp[f"Pi_Score_{c}"] = pi
    comp["significant"] = True

    log2 = pd.DataFrame(rng.normal(22, 3, (n, 12)), columns=samples)
    log2.insert(0, "Protein.Group", comp["Protein.Group"])
    log2.insert(1, "Genes", comp["Genes"])

    mods = np.array(["grey"] * n, dtype=object)
    mods[:60] = "blue"                                   # holds the 30 core proteins
    mods[60:180] = "brown"
    mods[180:330] = "turquoise"
    mods[330:420] = "yellow"
    wg = pd.DataFrame({"Protein.Group": comp["Protein.Group"], "Genes": comp["Genes"],
                       "Module": mods, "HubScore": rng.uniform(0, 1, n)})

    eg = pd.DataFrame({"Sample": samples})
    per_sample = [c for c in conds for _ in range(3)]
    eg["ME_blue"] = [{"A": -0.4, "B": -0.2, "C": 0.2, "D": 0.45}[c] + rng.normal(0, 0.03)
                     for c in per_sample]
    for m in ("brown", "turquoise", "yellow", "grey"):
        eg[f"ME_{m}"] = rng.normal(0, 0.3, 12)

    kme = pd.DataFrame({"Protein.Group": comp["Protein.Group"]})
    for m in ("blue", "brown", "turquoise", "yellow", "grey"):
        kme[f"kME_{m}"] = rng.uniform(-1, 1, n)

    path = tmp_path / "wgcna.xlsx"
    with pd.ExcelWriter(path) as xw:
        log2.to_excel(xw, sheet_name="Log2_Impute", index=False)
        comp.to_excel(xw, sheet_name="Differential_Expression", index=False)
        wg.to_excel(xw, sheet_name="WGCNA_Results", index=False)
        eg.to_excel(xw, sheet_name="WGCNA_Eigengenes", index=False)
        kme.to_excel(xw, sheet_name="WGCNA_kME", index=False)
    design = tmp_path / "design.csv"
    pd.DataFrame({"sample": samples, "condition": per_sample}).to_csv(design, index=False)
    return {"path": str(path), "design": str(design)}


def _run(wb):
    xls = pd.read_excel(wb["path"], sheet_name=None)
    stats = io_layer.load_stats(wb["path"], design_map=io_layer.load_design(wb["design"]))
    prots = axis_proteins.run(stats)
    core = list(prots["table"].head(30).index)
    return axis_wgcna.run(xls, stats, core_ids=core), stats


class TestLoading:
    def test_sheets_are_detected(self, workbook_wgcna):
        res, _ = _run(workbook_wgcna)
        assert res is not None
        assert res["assignment"] is not None and len(res["assignment"]) == 500
        assert res["eigengenes"] is not None and res["eigengenes"].shape == (12, 5)
        assert res["kme"] is not None

    def test_text_columns_are_found_whatever_pandas_calls_them(self, workbook_wgcna):
        """REGRESSION: pandas 2 says 'object', pandas 3 says 'str'; testing for object
        silently lost the protein->module assignment."""
        res, _ = _run(workbook_wgcna)
        assert res["assignment"] is not None
        assert set(res["assignment"].columns) >= {"protein", "module"}

    def test_grey_is_not_counted_as_a_module(self, workbook_wgcna):
        res, _ = _run(workbook_wgcna)
        assert res["n_modules"] == 4                     # blue, brown, turquoise, yellow
        assert res["core_enrichment"].loc[
            res["core_enrichment"]["module"] == "grey", "is_grey"].all()

    def test_missing_wgcna_returns_none(self, tmp_path, workbook_wgcna):
        src = pd.read_excel(workbook_wgcna["path"], sheet_name=None)
        p = tmp_path / "nowgcna.xlsx"
        with pd.ExcelWriter(p) as xw:
            for n, d in src.items():
                if not n.startswith("WGCNA"):
                    d.to_excel(xw, sheet_name=n, index=False)
        xls = pd.read_excel(p, sheet_name=None)
        stats = io_layer.load_stats(str(p))
        assert axis_wgcna.run(xls, stats, core_ids=[]) is None


class TestModuleCondition:
    def test_the_responding_module_is_found(self, workbook_wgcna):
        res, _ = _run(workbook_wgcna)
        a = res["association"].set_index("module")
        assert a.loc["ME_blue", "effect"] > 0.8
        others = a.drop(index="ME_blue")["effect"]
        assert (others < 0.5).all()

    def test_bh_column_present(self, workbook_wgcna):
        res, _ = _run(workbook_wgcna)
        assert "p_adj_BH" in res["association"].columns


class TestCoreEnrichment:
    def test_core_signature_is_located(self, workbook_wgcna):
        res, _ = _run(workbook_wgcna)
        e = res["core_enrichment"].set_index("module")
        assert e.loc["blue", "n_core"] == 30
        assert e.loc["blue", "p_adj_BH"] < 0.01
        assert (e.drop(index="blue")["n_core"] == 0).all()

    def test_no_core_yields_empty(self, workbook_wgcna):
        xls = pd.read_excel(workbook_wgcna["path"], sheet_name=None)
        stats = io_layer.load_stats(workbook_wgcna["path"])
        res = axis_wgcna.run(xls, stats, core_ids=[])
        assert res["core_enrichment"].empty


class TestModuleContrast:
    """Which comparison does a module follow - not just 'does it follow the design'."""

    def test_module_contrast_table(self, workbook_wgcna):
        res, _ = _run(workbook_wgcna)
        ca = res["contrast_association"]
        assert not ca.empty
        assert set(ca.columns) >= {"module", "contrast", "effect", "p_value", "p_adj_BH"}

    def test_the_responding_module_separates_contrasts(self, workbook_wgcna):
        res, _ = _run(workbook_wgcna)
        ca = res["contrast_association"]
        best = ca.iloc[0]
        assert best["module"] == "ME_blue"
        assert best["effect"] > 0.9

    def test_p_floor_is_reported_when_groups_are_small(self, workbook_wgcna):
        """REGRESSION: with 3 vs 3 the exact test cannot go below p=0.1, so a
        'not significant' verdict would be an artefact of group size."""
        res, _ = _run(workbook_wgcna)
        ca = res["contrast_association"]
        assert ca.attrs.get("n_per_group") == 3
        assert abs(ca.attrs.get("p_floor") - 0.1) < 1e-9
        assert (ca["p_value"].dropna() >= 0.1 - 1e-9).all()


class TestHubs:
    def test_hubs_per_module(self, workbook_wgcna):
        res, _ = _run(workbook_wgcna)
        h = res["hubs"]
        assert not h.empty
        assert h.groupby("module").size().max() <= 5
        for _, grp in h.groupby("module"):
            v = grp["hub_score"].values
            assert (v[:-1] >= v[1:]).all()

    def test_hub_description_is_joined_from_the_stats_workbook(self, workbook_wgcna):
        res, _ = _run(workbook_wgcna)
        h = res["hubs"]
        assert "description" in h.columns
        assert (h["description"] != "").all()
        for prot, desc in zip(h["protein"], h["description"]):
            assert desc == f"Synthetic protein {int(prot[1:])} (test)"

    def test_missing_description_column_does_not_crash(self, workbook_wgcna, tmp_path):
        """No First.Protein.Description / Protein.Names / Description in the workbook:
        the hub table must still be produced, just without that column."""
        xls = pd.read_excel(workbook_wgcna["path"], sheet_name=None)
        stats = io_layer.load_stats(workbook_wgcna["path"],
                                    design_map=io_layer.load_design(workbook_wgcna["design"]))
        stats.desc_col = None
        prots = axis_proteins.run(stats)
        core = list(prots["table"].head(30).index)
        res = axis_wgcna.run(xls, stats, core_ids=core)
        assert not res["hubs"].empty
        assert "description" not in res["hubs"].columns
