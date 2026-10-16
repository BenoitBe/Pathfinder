"""Reversed regulation only between comparable contrasts; heterogeneous conditions flagged.

Regressions from a real time-course x cultivar run: a protein up with time in both
cultivars and lower in cultivar A than B at 24 h was reported as 'opposite regulation'
(A24_vs_A0 up, A24_vs_B24 down share a condition: an ordering, not a reversal), and a
condition whose replicates sat at different stages of the response passed the
correlation-based quality check (r = 0.95 over all proteins).
"""
import numpy as np
import pandas as pd

import axis_meta_go
import axis_proteins
import build_report_html
import summary as summary_mod


class _Stats:
    def __init__(self, comp, contrasts):
        self.comp, self.contrasts = comp, contrasts
        self.id_col, self.gene_col, self.desc_col = "Protein.Group", "Genes", "Desc"


def _design_stats(seed=0, n=400):
    """Two strata (X, Y) x two levels (0, 1): X1_vs_X0 and Y1_vs_Y0 are parallel,
    X1_vs_Y1 is a stratum comparison sharing X1 with X1_vs_X0."""
    rng = np.random.default_rng(seed)
    contrasts = ["X1_vs_X0", "Y1_vs_Y0", "X1_vs_Y1"]
    resp = rng.normal(0, 1.5, n)                        # shared response to level 1
    lfc = {"X1_vs_X0": resp + rng.normal(0, .3, n),
           "Y1_vs_Y0": resp + rng.normal(0, .3, n),
           "X1_vs_Y1": rng.normal(0, 1, n)}
    lfc["Y1_vs_Y0"][0] = -2.0                           # P0000: genuine reversal
    lfc["X1_vs_X0"][0] = 2.0
    lfc["X1_vs_X0"][1], lfc["X1_vs_Y1"][1] = 2.0, -1.5  # P0001: up with level, lower than Y1
    lfc["Y1_vs_Y0"][1] = 2.5
    comp = pd.DataFrame({"Protein.Group": [f"P{i:04d}" for i in range(n)],
                         "Genes": [f"G{i}" for i in range(n)], "Desc": ""})
    for c in contrasts:
        comp[f"{c}_diff"] = lfc[c]
        comp[f"{c}_p.adj"] = np.where(np.abs(lfc[c]) > 1.4, 1e-4, 0.5)
    return _Stats(comp, contrasts)


def test_pairs_sharing_a_condition_are_never_comparable():
    pairs = axis_proteins.comparable_pairs(_design_stats())
    assert ("X1_vs_X0", "Y1_vs_Y0") in pairs
    assert not any("X1_vs_Y1" in p for p in pairs)


def test_reversal_vs_ordering():
    t = axis_proteins.run(_design_stats())["table"]
    assert t.loc["P0000", "pattern"] == "opposite"
    assert "Y1_vs_Y0" in t.loc["P0000", "opposite_in"]
    assert t.loc["P0001", "pattern"] == "mixed"          # B < A < C, not a contradiction
    assert t.loc["P0001", "opposite_in"] == ""


def test_mirrored_contrast_is_aligned():
    """Y0_vs_Y1 written the other way round: same signs are then the reversal."""
    st = _design_stats()
    st.comp["Y0_vs_Y1_diff"] = -st.comp.pop("Y1_vs_Y0_diff")
    st.comp["Y0_vs_Y1_p.adj"] = st.comp.pop("Y1_vs_Y0_p.adj")
    st.contrasts = ["X1_vs_X0", "Y0_vs_Y1", "X1_vs_Y1"]
    pairs = axis_proteins.comparable_pairs(st)
    assert pairs[("X1_vs_X0", "Y0_vs_Y1")] == -1
    t = axis_proteins.run(st, pairs=pairs)["table"]
    assert t.loc["P0000", "pattern"] == "opposite"


def test_without_fold_changes_legacy_rule_is_kept():
    st = _design_stats()
    matrix, meta, _, _ = axis_proteins.build_matrix(st)
    d = axis_proteins.directionality(matrix, meta, pairs=None)
    assert d.loc["P0001", "pattern"] == "opposite"


def test_go_terms_follow_the_same_rule():
    long = pd.DataFrame([
        ["X1_vs_X0", "s", "GO:BP", "GO:1", "t1", 1e-3, 2.0, ["P1"]],
        ["X1_vs_Y1", "s", "GO:BP", "GO:1", "t1", 1e-3, -2.0, ["P1"]],
        ["X1_vs_X0", "s", "GO:BP", "GO:2", "t2", 1e-3, 2.0, ["P2"]],
        ["Y1_vs_Y0", "s", "GO:BP", "GO:2", "t2", 1e-3, -2.0, ["P2"]],
    ], columns=["contrast", "sheet", "source", "term_id", "term_name", "p_value", "z",
                "intersection"])
    pairs = {("X1_vs_X0", "Y1_vs_Y0"): 1}
    d = axis_meta_go.run(long, pairs=pairs)["directionality"]
    assert d.loc["GO:1", "pattern"] == "mixed"
    assert d.loc["GO:2", "pattern"] == "opposite"


def _explore(spread_b=1.0, seed=1):
    rng = np.random.default_rng(seed)
    rows, centres = [], {"A": (0, 0), "B": (10, 0), "C": (0, 10), "D": (10, 10)}
    for cond, (x, y) in centres.items():
        for r in range(3):
            s = spread_b if cond == "B" else 1.0
            rows.append({"sample": f"{cond}{r}", "condition": cond,
                         "PC1": x + rng.normal(0, s), "PC2": y + rng.normal(0, s)})
    sc = pd.DataFrame(rows).set_index("sample")
    if spread_b > 1:                                   # one B replicate still at A's stage
        sc.loc["B0", ["PC1", "PC2"]] = [1.0, 0.5]
    assoc = pd.DataFrame({"PC": ["PC1", "PC2"], "H": [10, 10], "p_value": [0.01, 0.01],
                          "p_adj_BH": [0.01, 0.01], "epsilon2": [0.8, 0.8]})
    return {"pca": {"scores": sc, "var_exp": np.array([40.0, 30.0])}, "pc_factor": assoc}


def test_heterogeneous_condition_is_flagged():
    flags = summary_mod._dispersion_flags(_explore(spread_b=3.0))
    assert len(flags) == 1 and "Condition B" in flags[0]
    assert "B0 lies closer to A" in flags[0]


def test_homogeneous_conditions_are_not_flagged():
    assert summary_mod._dispersion_flags(_explore(spread_b=1.0)) == []


def test_long_contrast_name_is_not_cut_in_the_chain():
    svg = build_report_html._svg_chain({
        "core_n": 200, "core_title": "top-ranked signature",
        "core_detail": "recurrence, then Pi-score", "module": "blue",
        "module_detail": "92/676 proteins", "contrast": "embAviso24_vs_embAviso0",
        "contrast_detail": "separation=1.00", "verdict": "v"})
    assert "embAviso24_vs_embAviso0<" in svg
    assert "top-ranked signature".upper() in svg
