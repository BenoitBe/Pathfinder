"""Protein selection = the upstream (Deimos) volcano rule, read from the manifest."""
import numpy as np
import pandas as pd

import io_layer
from io_layer import SelectionRule, selection_rule


def test_default_without_manifest_is_historical():
    r = selection_rule({})
    assert (r.p_column, r.p_thr, r.ratio_min, r.lfc_min) == ("p.adj", 0.05, 1.0, 0.0)
    assert r.text() == "p.adj < 0.05, no fold-change criterion"
    assert not r.source.startswith("upstream")


def test_structured_keys_win_over_text():
    m = {"stats": {"volcano_rule": "p.adj < 0.05 & ratio >= 1.5",
                   "volcano_p_column": "p.val", "volcano_p_thresh": "0.01",
                   "volcano_ratio_min": "2.0"}}
    r = selection_rule(m)
    assert (r.p_column, r.p_thr, r.ratio_min) == ("p.val", 0.01, 2.0)
    assert r.lfc_min == 1.0 and r.source == "upstream manifest (volcano rule)"


def test_text_rule_of_older_manifests():
    for text, exp in [("p.adj < 0.05 & ratio >= 2.0", ("p.adj", 0.05, 2.0)),
                      ("p.value < 0.01 & ratio >= 1.5", ("p.val", 0.01, 1.5)),
                      ("p.adj < 0.1", ("p.adj", 0.1, 1.0))]:
        r = selection_rule({"stats": {"volcano_rule": text}})
        assert (r.p_column, r.p_thr, r.ratio_min) == exp, text


def test_cli_overrides_are_flagged():
    m = {"stats": {"volcano_rule": "p.adj < 0.05 & ratio >= 2.0"}}
    r = selection_rule(m, p_thr=0.01)
    assert r.p_thr == 0.01 and r.ratio_min == 2.0 and "overridden by --pthr 0.01" in r.source
    r = selection_rule(m, min_ratio=1)
    assert r.lfc_min == 0 and "--min-ratio 1" in r.source
    r = selection_rule(m, p_thr=0.05, min_ratio=2.0)       # same values: not an override
    assert "overridden" not in r.source


def test_mask_matches_deimos_volcano_status():
    """Same comparisons as Deimos: p < thr strictly, |log2FC| >= log2(ratio), NaN never."""
    r = SelectionRule("p.adj", 0.05, 2.0)
    p = pd.Series([0.01, 0.05, 0.049, 0.01, np.nan, 0.01])
    d = pd.Series([1.0, 3.0, -1.0, 0.999, 2.0, np.nan])
    assert r.mask(p, d).tolist() == [True, False, True, False, False, False]
    r1 = SelectionRule("p.adj", 0.05, 1.0)
    assert r1.mask(p, pd.Series([0.0] * 6)).tolist() == [True, False, True, True, False, True]


def test_protein_axis_uses_the_rule():
    import axis_proteins

    class S:
        pass
    s = S()
    s.comp = pd.DataFrame({"name": ["a", "b", "c", "d"], "Genes": ["A", "B", "C", "D"],
                           "X_vs_Y_diff": [2.0, 0.5, -1.5, 1.2],
                           "X_vs_Y_p.adj": [0.01, 0.001, 0.2, 0.04],
                           "X_vs_Y_p.val": [0.001, 0.0001, 0.01, 0.004]})
    s.id_col, s.gene_col, s.desc_col, s.contrasts = "name", "Genes", "desc", ["X_vs_Y"]
    m, *_ = axis_proteins.build_matrix(s, rule=SelectionRule("p.adj", 0.05, 2.0))
    assert sorted(m.index) == ["a", "d"]                    # b: |LFC| < 1 ; c: p.adj 0.2
    m, *_ = axis_proteins.build_matrix(s, rule=SelectionRule("p.val", 0.05, 2.0))
    assert sorted(m.index) == ["a", "c", "d"]               # raw p : c passes
    m, *_ = axis_proteins.build_matrix(s, p_thr=0.05)       # historical call unchanged
    assert sorted(m.index) == ["a", "b", "d"]
