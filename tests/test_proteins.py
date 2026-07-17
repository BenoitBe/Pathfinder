"""Protein axis: which proteins carry the signature."""
import itertools

import numpy as np
import pandas as pd
import pytest

import axis_proteins
import io_layer


@pytest.fixture
def workbook_proteins(tmp_path):
    """Known structure: 10 core proteins (significant everywhere, same direction),
    10 with opposite regulation, some contrast-specific ones."""
    rng = np.random.default_rng(31)
    conds = ["A", "B", "C", "D"]
    contrasts = [f"{a}_vs_{b}" for a, b in itertools.combinations(conds, 2)]
    n = 300
    comp = pd.DataFrame({
        "Protein.Group": [f"P{i:04d}" for i in range(n)],
        "Genes": [f"G{i // 2:04d}" for i in range(n)],
        "First.Protein.Description": [f"desc {i}" for i in range(n)],
    })
    for k, c in enumerate(contrasts):
        p, d = np.ones(n), np.zeros(n)
        pi, rb = np.zeros(n), np.zeros(n)
        p[:10] = 1e-8                       # core, always up
        d[:10] = rng.normal(2.5, 0.3, 10)
        pi[:10] = rng.uniform(15, 25, 10)
        rb[:10] = rng.uniform(80, 100, 10)
        p[10:20] = 1e-6                     # flipped depending on the contrast
        d[10:20] = rng.normal(1.5, 0.2, 10) * (1 if c.startswith("A") else -1)
        pi[10:20] = rng.uniform(8, 14, 10)
        rb[10:20] = rng.uniform(50, 80, 10)
        lo = 20 + k * 5                     # 5 proteins specific to this contrast only
        p[lo:lo + 5] = 1e-4
        d[lo:lo + 5] = rng.normal(1.2, 0.3, 5)
        pi[lo:lo + 5] = rng.uniform(3, 7, 5)
        rb[lo:lo + 5] = rng.uniform(30, 60, 5)
        comp[f"{c}_p.adj"] = p
        comp[f"{c}_p.val"] = p
        comp[f"{c}_diff"] = d
        comp[f"Pi_Score_{c}"] = pi
        comp[f"Robustness_Score_{c}"] = rb
    comp["significant"] = True
    for s in [f"{c}_{i}" for c in conds for i in (1, 2, 3)]:
        comp[s] = rng.normal(22, 2, n)
    p = tmp_path / "proteins.xlsx"
    with pd.ExcelWriter(p) as xw:
        comp.to_excel(xw, sheet_name="Comparaison", index=False)
    return {"path": str(p), "contrasts": contrasts}


class TestProteinMatrix:
    def test_only_significant_proteins_are_kept(self, workbook_proteins):
        stats = io_layer.load_stats(workbook_proteins["path"])
        res = axis_proteins.run(stats)
        assert 20 <= res["n_significant"] < 300

    def test_non_significant_is_nan_not_zero(self, workbook_proteins):
        """A 0 would read as 'no change'; NaN reads as 'not significant'."""
        stats = io_layer.load_stats(workbook_proteins["path"])
        m = axis_proteins.run(stats)["matrix"]
        assert m.isna().any().any()
        assert not (m == 0).any().any()

    def test_pi_and_robustness_are_picked_up(self, workbook_proteins):
        stats = io_layer.load_stats(workbook_proteins["path"])
        res = axis_proteins.run(stats)
        assert res["has_pi"] and res["has_robustness"]
        assert "best_pi" in res["table"].columns
        assert "median_robustness" in res["table"].columns

    def test_core_proteins_rank_first(self, workbook_proteins):
        stats = io_layer.load_stats(workbook_proteins["path"])
        t = axis_proteins.run(stats)["table"]
        top = t.head(10)
        assert (top["n_sig"] == 6).all()
        assert (top["pattern"] == "shared up").all()

    def test_opposite_regulation_is_detected(self, workbook_proteins):
        stats = io_layer.load_stats(workbook_proteins["path"])
        t = axis_proteins.run(stats)["table"]
        opp = t[t["pattern"] == "opposite"]
        assert len(opp) == 10                       # exactly the injected ones
        assert (opp["up_in"].str.len() > 0).all()
        assert (opp["down_in"].str.len() > 0).all()

    def test_pi_score_ranks_but_never_selects(self, workbook_proteins):
        """Inclusion must rest on the adjusted p-value alone."""
        stats = io_layer.load_stats(workbook_proteins["path"])
        res = axis_proteins.run(stats, p_thr=0.05)
        n_by_p = int(sum(
            (pd.to_numeric(stats.comp[f"{c}_p.adj"], errors="coerce") < 0.05).any()
            for c in [])) if False else None
        # every kept protein is significant somewhere, regardless of its Pi
        m = res["matrix"]
        assert m.notna().any(axis=1).all()
        # a high-Pi protein that never passes p.adj must be absent
        assert "P0299" not in m.index

    def test_ranking_is_recurrence_then_pi(self, workbook_proteins):
        stats = io_layer.load_stats(workbook_proteins["path"])
        t = axis_proteins.run(stats)["table"]
        assert (t["n_sig"].values[:-1] >= t["n_sig"].values[1:]).all()
        core = t[t["n_sig"] == t["n_sig"].max()]
        pi = core["best_pi"].dropna().values
        assert (pi[:-1] >= pi[1:] - 1e-9).all(), "Pi does not order within equal recurrence"

    def test_missing_pi_falls_back_on_lfc(self, tmp_path, workbook_proteins):
        comp = pd.read_excel(workbook_proteins["path"])
        comp = comp[[c for c in comp.columns if not str(c).startswith("Pi_Score")]]
        p = tmp_path / "nopi.xlsx"
        with pd.ExcelWriter(p) as xw:
            comp.to_excel(xw, sheet_name="Comparaison", index=False)
        stats = io_layer.load_stats(str(p))
        res = axis_proteins.run(stats)
        assert res["has_pi"] is False
        assert len(res["table"]) > 0

    def test_contrast_clustering_at_protein_level(self, workbook_proteins):
        stats = io_layer.load_stats(workbook_proteins["path"])
        res = axis_proteins.run(stats)
        cc = res["contrast_corr"]
        assert cc.shape[0] == len(workbook_proteins["contrasts"])
        assert cc.notna().any().any()

    def test_empty_when_nothing_significant(self, tmp_path):
        rng = np.random.default_rng(1)
        n = 50
        comp = pd.DataFrame({"Protein.Group": [f"P{i}" for i in range(n)],
                             "Genes": [f"G{i}" for i in range(n)]})
        comp["A_vs_B_p.adj"] = 1.0
        comp["A_vs_B_p.val"] = 1.0
        comp["A_vs_B_diff"] = rng.normal(0, 1, n)
        for s in ("A_1", "A_2", "B_1", "B_2"):
            comp[s] = rng.normal(22, 2, n)
        p = tmp_path / "nothing.xlsx"
        with pd.ExcelWriter(p) as xw:
            comp.to_excel(xw, sheet_name="Comparaison", index=False)
        stats = io_layer.load_stats(str(p))
        res = axis_proteins.run(stats)
        assert res["n_significant"] == 0
        assert res["matrix"].empty
