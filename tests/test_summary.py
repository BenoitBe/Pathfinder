"""The interpretive summary must never claim more than the data supports."""
import numpy as np
import pandas as pd

import axis_explore
import axis_meta_go
import io_layer
import summary as summary_mod


def _md(stats, meta_go, explore):
    return summary_mod.generate(stats, meta_go, explore, {"date": "t"})["markdown"]


def _go_long(rows):
    return pd.DataFrame(rows, columns=["contrast", "sheet", "source", "term_id",
                                       "term_name", "p_value", "z", "intersection"])


def _empty_go():
    return {"matrix": pd.DataFrame(),
            "meta": pd.DataFrame(columns=["term_name", "term_id", "source", "n_sig",
                                          "redundancy_cluster", "is_representative"]),
            "directionality": pd.DataFrame(), "term_genes": {},
            "contrast_corr": pd.DataFrame(), "contrast_order": [],
            "contrast_shared": pd.DataFrame(), "contrast_relation": None}


class TestDesignSection:
    def test_reports_capture_when_design_is_clean(self, workbook_full):
        stats = io_layer.load_stats(workbook_full["path"],
                                    design_map=io_layer.load_design(workbook_full["design"]))
        ex = axis_explore.run(stats, do_umap=False)
        md = _md(stats, _empty_go(), ex)
        assert "Condition structures the variance" in md
        assert "batch effect" not in md.split("## Quality")[0]

    def test_warns_when_matrix_is_a_zscored_subset(self, workbook_zscore_only):
        """REGRESSION: assessing separation on selected proteins is circular."""
        stats = io_layer.load_stats(workbook_zscore_only["path"],
                                    design_map=io_layer.load_design(workbook_zscore_only["design"]))
        ex = axis_explore.run(stats, do_umap=False)
        md = _md(stats, _empty_go(), ex)
        assert "z-scored matrix" in md
        assert "inflated" in md

    def test_no_caveat_on_raw_log2(self, workbook_full):
        stats = io_layer.load_stats(workbook_full["path"],
                                    design_map=io_layer.load_design(workbook_full["design"]))
        ex = axis_explore.run(stats, do_umap=False)
        assert "z-scored matrix" not in _md(stats, _empty_go(), ex)


class TestQualitySection:
    def test_only_material_gaps_are_flagged(self):
        """REGRESSION: with replicates tightly grouped (median r=0.92), the robust
        rule alone flagged a 0.02 drop as an outlier."""
        samples = [f"S{i}" for i in range(9)]
        rs = {"S0": 0.92, "S1": 0.92, "S2": 0.93,
              "S3": 0.90, "S4": 0.92, "S5": 0.86,      # only 0.86 is material
              "S6": 0.92, "S7": 0.90, "S8": 0.93}
        corr = pd.DataFrame(np.eye(9), index=samples, columns=samples)
        for i, a in enumerate(samples):
            for j, b in enumerate(samples):
                if i != j:
                    corr.loc[a, b] = rs[a] if (i // 3 == j // 3) else 0.5
        corr.attrs["zscored"] = False
        cond = {s: f"C{i // 3}" for i, s in enumerate(samples)}
        lines = summary_mod._qc_section({"correlation": corr}, cond)
        flagged = [l for l in lines if l.startswith("⚠")]
        assert len(flagged) == 1
        assert "S5" in flagged[0]

    def test_verdict_follows_the_p_value_not_only_the_effect_size(self):
        """REGRESSION: 'group separation is clear' was printed at p=0.069 because the
        verdict only looked at epsilon squared."""
        def fake(p, eps2):
            scores = pd.DataFrame({"PC1": range(24), "PC2": range(24),
                                   "condition": [f"C{i // 3}" for i in range(24)]})
            assoc = pd.DataFrame({"PC": ["PC1", "PC2"], "H": [1, 1],
                                  "p_value": [p, p], "p_adj_BH": [p, p],
                                  "epsilon2": [eps2, eps2]})
            return {"pca": {"scores": scores, "var_exp": [26, 15]}, "pc_factor": assoc}

        marginal = " ".join(summary_mod._design_section(fake(0.069, 0.59)))
        assert "does not reach significance" in marginal
        assert "consistent trend, not a proof" in marginal
        assert "separation is clear" not in marginal

        clear = " ".join(summary_mod._design_section(fake(0.004, 0.87)))
        assert "separation is clear" in clear

    def test_global_pattern_is_summarised_not_listed(self, workbook_zscore_only):
        """REGRESSION: 24 identical warnings are unreadable; say it once."""
        stats = io_layer.load_stats(workbook_zscore_only["path"],
                                    design_map=io_layer.load_design(workbook_zscore_only["design"]))
        ex = axis_explore.run(stats, do_umap=False)
        md = _md(stats, _empty_go(), ex)
        qc = md.split("## Quality & homogeneity")[1]
        assert qc.count("correlate poorly with their replicates") <= 1

    def test_zscore_note_is_present(self, workbook_zscore_only):
        stats = io_layer.load_stats(workbook_zscore_only["path"],
                                    design_map=io_layer.load_design(workbook_zscore_only["design"]))
        ex = axis_explore.run(stats, do_umap=False)
        assert "per-protein z-scores" in _md(stats, _empty_go(), ex)


class TestContrastSection:
    """REGRESSIONS: mechanical correlation reported as biology; and a claim of
    'same direction' backed by r = -0.07 over 3 shared terms."""

    def _build(self, rows, workbook_full):
        stats = io_layer.load_stats(workbook_full["path"])
        go = axis_meta_go.run(_go_long(rows), p_thr=0.05)
        ex = axis_explore.run(stats, do_umap=False)
        return _md(stats, go, ex)

    def test_pairs_sharing_a_condition_are_excluded(self, workbook_full):
        rng = np.random.default_rng(2)
        prof = rng.normal(0, 1.5, 40)
        rows = []
        # A_vs_C and B_vs_C share reference C and are identical -> would top the list
        for c, z in [("A_vs_C", prof), ("B_vs_C", prof),
                     ("D_vs_E", prof), ("F_vs_G", prof * 0.9)]:
            for i, zz in enumerate(z):
                rows.append([c, "s", "GO:BP", f"GO:{i}", f"t{i}", 0.001, zz, ["P1"]])
        md = self._build(rows, workbook_full)
        assert "pairs share a condition" in md
        # the reported pair must be one of the disjoint ones
        line = [l for l in md.splitlines() if "Closest independent" in l][0]
        assert "A_vs_C and B_vs_C" not in line

    def test_thin_pairs_are_excluded_from_claims(self, workbook_full):
        rng = np.random.default_rng(4)
        prof = rng.normal(0, 1.5, 40)
        rows = []
        for i, zz in enumerate(prof):
            rows.append(["A_vs_B", "s", "GO:BP", f"GO:{i}", f"t{i}", 0.001, zz, ["P1"]])
            rows.append(["C_vs_D", "s", "GO:BP", f"GO:{i}", f"t{i}", 0.001,
                         zz + rng.normal(0, .2), ["P1"]])
        # tiny signature, perfectly anti-correlated over 4 terms -> must NOT win
        for i in range(4):
            rows.append(["E_vs_F", "s", "GO:BP", f"GO:{i}", f"t{i}", 0.001, -prof[i], ["P1"]])
        md = self._build(rows, workbook_full)
        assert "fewer than 10 shared terms" in md
        assert "E_vs_F" not in md.split("## Per-contrast detail")[1].split("- A_vs_B")[0]

    def test_shared_term_count_is_always_reported(self, workbook_full):
        rng = np.random.default_rng(6)
        prof = rng.normal(0, 1.5, 40)
        rows = []
        for c, z in [("A_vs_B", prof), ("C_vs_D", prof)]:
            for i, zz in enumerate(z):
                rows.append([c, "s", "GO:BP", f"GO:{i}", f"t{i}", 0.001, zz, ["P1"]])
        md = self._build(rows, workbook_full)
        assert "shared terms" in md

    def test_unrelated_signatures_produce_no_claim(self, workbook_full):
        rng = np.random.default_rng(8)
        rows = []
        for c in ["A_vs_B", "C_vs_D"]:
            for i in range(40):
                rows.append([c, "s", "GO:BP", f"GO:{i}", f"t{i}", 0.001,
                             rng.normal(0, 1.5), ["P1"]])
        md = self._build(rows, workbook_full)
        assert "largely unrelated" in md or "Closest independent" not in md

    def test_few_contrasts_stay_one_line_each(self, workbook_full):
        """Below the grouping threshold, per-contrast detail must stay literal."""
        rng = np.random.default_rng(9)
        rows = []
        for c in ["A_vs_B", "C_vs_D", "E_vs_F"]:
            for i in range(20):
                rows.append([c, "s", "GO:BP", f"GO:{i}", f"t{i}", 0.001,
                             rng.normal(0, 1.5), ["P1"]])
        md = self._build(rows, workbook_full)
        detail = md.split("## Per-contrast detail")[1]
        assert "- A_vs_B: 20 significantly enriched terms." in detail
        assert "Family 1" not in detail

    def test_many_contrasts_are_grouped_into_families(self, workbook_full):
        """REGRESSION: 28 unread per-contrast lines -> a handful of families instead.

        Three blocks of contrasts share a near-identical GO profile (r ~ 1 within
        a block); the fourth is pure noise. With 12 contrasts (> MIN_CONTRASTS_FOR_FAMILIES),
        the detail section must read as families, not 12 individual bullets.
        """
        rng = np.random.default_rng(10)
        prof = {0: rng.normal(0, 1.5, 60), 1: rng.normal(0, 1.5, 60), 2: rng.normal(0, 1.5, 60)}
        rows = []
        contrasts = [f"X{i}_vs_Y{i}" for i in range(12)]
        for i, c in enumerate(contrasts):
            block = i % 4
            if block < 3:
                z = prof[block] + rng.normal(0, 0.05, 60)   # near-identical within block
            else:
                z = rng.normal(0, 1.5, 60)                  # unrelated singleton
            for k, zz in enumerate(z):
                rows.append([c, "s", "GO:BP", f"GO:{k}", f"t{k}", 0.001, zz, ["P1"]])
        md = self._build(rows, workbook_full)
        detail = md.split("## Per-contrast detail")[1]
        assert "Family 1" in detail
        assert "families summarising 12 contrasts" in detail
        # a tightly correlated block must not be spelled out as 3 separate bullets
        assert detail.count(": 60 significantly enriched terms.") < len(contrasts)


class TestFunctionalSection:
    def test_absent_go_is_stated(self, workbook_full):
        stats = io_layer.load_stats(workbook_full["path"])
        ex = axis_explore.run(stats, do_umap=False)
        assert "GO enrichment not provided" in _md(stats, _empty_go(), ex)

    def test_counts_are_non_redundant(self, workbook_full):
        genes = [f"P{i}" for i in range(10)]
        rows = [["A_vs_B", "s", "GO:BP", f"GO:{i}", f"term {i}", 0.001, 1.0 + i * 0.1,
                 genes[:9] + [f"X{i}"]] for i in range(4)]
        stats = io_layer.load_stats(workbook_full["path"])
        go = axis_meta_go.run(_go_long(rows), p_thr=0.05, jaccard_cut=0.5)
        ex = axis_explore.run(stats, do_umap=False)
        md = _md(stats, go, ex)
        assert "non-redundant" in md
