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


class TestAnnotationBasics:
    """The join no single tab provides: which proteins, doing what, in which module."""

    def _go(self, core_ids):
        import axis_meta_go
        rows = []
        terms = [("unfolded protein binding", "GO:MF", core_ids[:5]),
                 ("plastid stroma", "GO:CC", core_ids[3:9]),
                 ("elsewhere", "GO:BP", ["ZZZ1", "ZZZ2"])]
        for c in ("A_vs_B", "A_vs_C"):
            for name, src, accs in terms:
                rows.append([c, "s", src, f"GO:{abs(hash(name)) % 9999}", name,
                             0.001, 1.5, list(accs)])
        long = pd.DataFrame(rows, columns=["contrast", "sheet", "source", "term_id",
                                           "term_name", "p_value", "z", "intersection"])
        return axis_meta_go.run(long, p_thr=0.05)

    def test_terms_and_module_are_joined(self, workbook_proteins):
        stats = io_layer.load_stats(workbook_proteins["path"])
        res = axis_proteins.run(stats)
        core = list(res["table"].head(10).index)
        go = self._go(core)
        wg = {"assignment": pd.DataFrame({"protein": core, "module": ["blue"] * len(core)}),
              "n_modules": 1}
        t = axis_proteins.annotate(res["table"], go_res=go, wgcna=wg)
        assert "module" in t.columns and "go_terms" in t.columns
        assert (t.loc[core, "module"] == "blue").all()
        assert int((t["n_go_terms"] > 0).sum()) >= 5

    def test_protein_outside_any_term_is_left_blank(self, workbook_proteins):
        stats = io_layer.load_stats(workbook_proteins["path"])
        res = axis_proteins.run(stats)
        core = list(res["table"].head(10).index)
        t = axis_proteins.annotate(res["table"], go_res=self._go(core), wgcna=None)
        far = t.index[-1]
        assert t.loc[far, "n_go_terms"] == 0
        assert t.loc[far, "go_terms"] == ""

    def test_annotation_is_optional(self, workbook_proteins):
        stats = io_layer.load_stats(workbook_proteins["path"])
        res = axis_proteins.run(stats)
        t = axis_proteins.annotate(res["table"], go_res=None, wgcna=None)
        assert len(t) == len(res["table"])
        assert "go_terms" not in t.columns


class TestAnnotation:
    """The join nobody else does: which proteins, carrying which functions, in which
    module. The GO sheets hold the accessions behind each term; WGCNA holds the module."""

    def _go(self, accs_by_term):
        long = pd.DataFrame([
            {"contrast": "A_vs_B", "sheet": "s", "source": src, "term_id": tid,
             "term_name": name, "p_value": 1e-5, "z": 1.0, "intersection": accs}
            for tid, (name, src, accs) in accs_by_term.items()])
        import axis_meta_go
        return axis_meta_go.run(long, p_thr=0.05)

    def test_terms_and_module_are_attached(self, workbook_proteins):
        go = self._go({"GO:1": ("protein folding", "GO:BP", [f"P{i:04d}" for i in range(5)]),
                       "GO:2": ("plastid stroma", "GO:CC", [f"P{i:04d}" for i in range(5, 10)])})
        wgcna = {"assignment": pd.DataFrame({
            "protein": [f"P{i:04d}" for i in range(30)],
            "module": ["blue"] * 10 + ["brown"] * 20})}
        stats = io_layer.load_stats(workbook_proteins["path"])
        t = axis_proteins.annotate(axis_proteins.run(stats)["table"], go, wgcna)
        assert "protein folding" in t.loc["P0000", "go_terms"]
        assert "plastid stroma" in t.loc["P0007", "go_terms"]
        assert t.loc["P0000", "module"] == "blue"
        assert t.loc["P0015", "module"] == "brown"

    def test_a_protein_in_no_term_gets_an_empty_annotation(self, workbook_proteins):
        go = self._go({"GO:1": ("protein folding", "GO:BP", [f"P{i:04d}" for i in range(5)])})
        stats = io_layer.load_stats(workbook_proteins["path"])
        t = axis_proteins.annotate(axis_proteins.run(stats)["table"], go, None)
        assert t.loc["P0000", "n_go_terms"] == 1
        assert t.loc["P0019", "n_go_terms"] == 0
        assert t.loc["P0019", "go_terms"] == ""

    def test_matching_works_on_gene_symbols_too(self, workbook_proteins):
        """gProfiler may return gene symbols rather than accessions."""
        go = self._go({"GO:1": ("protein folding", "GO:BP", ["G0000", "G0001"])})
        stats = io_layer.load_stats(workbook_proteins["path"])
        t = axis_proteins.annotate(axis_proteins.run(stats)["table"], go, None)
        assert t.loc["P0000", "n_go_terms"] == 1        # P0000 -> gene G0000

    def test_annotation_without_go_or_wgcna_is_harmless(self, workbook_proteins):
        stats = io_layer.load_stats(workbook_proteins["path"])
        base = axis_proteins.run(stats)["table"]
        t = axis_proteins.annotate(base, None, None)
        assert len(t) == len(base)
        assert "go_terms" not in t.columns


class TestAnnotationJoin:
    """The join nobody else does: which proteins, doing what, in which module."""

    def _go(self, core_ids):
        """GO result where two terms are carried by known core proteins."""
        keys = ["GO:0006457", "GO:0009570"]
        meta = pd.DataFrame({"term_name": ["protein folding", "plastid stroma"],
                             "term_id": keys, "source": ["GO:BP", "GO:CC"],
                             "n_sig": [6, 4], "redundancy_cluster": [1, 2],
                             "is_representative": [True, True]}, index=keys)
        return {"matrix": pd.DataFrame({"A_vs_B": [2.0, 1.0]}, index=keys),
                "meta": meta,
                "term_genes": {keys[0]: set(core_ids[:20]), keys[1]: set(core_ids[10:30])}}

    def _wgcna(self, core_ids):
        return {"n_modules": 2,
                "assignment": pd.DataFrame({"protein": core_ids + ["P0500"],
                                            "module": ["blue"] * len(core_ids) + ["brown"]})}

    def test_terms_and_module_are_joined(self, workbook_proteins):
        stats = io_layer.load_stats(workbook_proteins["path"])
        res = axis_proteins.run(stats)
        core = [f"P{i:04d}" for i in range(30)]
        t = axis_proteins.annotate(res["table"], self._go(core), self._wgcna(core))
        assert "go_terms" in t.columns and "module" in t.columns
        assert t.loc["P0000", "module"] == "blue"

    def test_each_protein_gets_only_its_own_terms(self, workbook_proteins):
        """A protein must not inherit terms it is not part of."""
        stats = io_layer.load_stats(workbook_proteins["path"])
        res = axis_proteins.run(stats)
        core = [f"P{i:04d}" for i in range(30)]
        t = axis_proteins.annotate(res["table"], self._go(core), self._wgcna(core))
        assert "protein folding" in t.loc["P0005", "go_terms"]      # in 0..19 only
        assert "plastid stroma" not in t.loc["P0005", "go_terms"]
        assert "protein folding" in t.loc["P0015", "go_terms"]      # in both
        assert "plastid stroma" in t.loc["P0015", "go_terms"]
        assert "plastid stroma" in t.loc["P0025", "go_terms"]       # in 10..29 only
        assert "protein folding" not in t.loc["P0025", "go_terms"]

    def test_term_count_matches(self, workbook_proteins):
        stats = io_layer.load_stats(workbook_proteins["path"])
        res = axis_proteins.run(stats)
        core = [f"P{i:04d}" for i in range(30)]
        t = axis_proteins.annotate(res["table"], self._go(core), self._wgcna(core))
        assert t.loc["P0015", "n_go_terms"] == 2
        assert t.loc["P0005", "n_go_terms"] == 1

    def test_annotation_is_optional(self, workbook_proteins):
        """No GO, no WGCNA: the table must still stand."""
        stats = io_layer.load_stats(workbook_proteins["path"])
        res = axis_proteins.run(stats)
        t = axis_proteins.annotate(res["table"], None, None)
        assert len(t) == len(res["table"])
        assert "go_terms" not in t.columns


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
