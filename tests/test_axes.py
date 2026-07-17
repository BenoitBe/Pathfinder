"""Analysis axes: Meta-GO, Exploration (PCA / factors / correlation), Publications."""
import numpy as np
import pandas as pd
import pytest

import axis_explore
import axis_meta_go
import axis_publications
import cache
import io_layer


def _go_long(rows):
    return pd.DataFrame(rows, columns=["contrast", "sheet", "source", "term_id",
                                       "term_name", "p_value", "z", "intersection"])


# --------------------------------------------------------------------------- #
# Meta-GO                                                                      #
# --------------------------------------------------------------------------- #
class TestMetaGoMatrix:
    def test_non_significant_is_nan_not_zero(self):
        """A 0 would read as 'no change' and distort clustering; NaN reads as 'no data'."""
        long = _go_long([
            ["A_vs_B", "s", "GO:BP", "GO:1", "t1", 0.001, 2.0, ["P1"]],
            ["A_vs_C", "s", "GO:BP", "GO:1", "t1", 0.900, 2.0, ["P1"]],
        ])
        matrix, meta, _ = axis_meta_go.build_matrix(long, p_thr=0.05)
        assert matrix.loc["GO:1", "A_vs_B"] == 2.0
        assert pd.isna(matrix.loc["GO:1", "A_vs_C"])

    def test_directionality_patterns(self):
        long = _go_long([
            ["A_vs_B", "s", "GO:BP", "GO:1", "shared", 0.001, 2.0, ["P1"]],
            ["A_vs_C", "s", "GO:BP", "GO:1", "shared", 0.001, 1.5, ["P1"]],
            ["A_vs_B", "s", "GO:BP", "GO:2", "mixed", 0.001, 2.0, ["P2"]],
            ["A_vs_C", "s", "GO:BP", "GO:2", "mixed", 0.001, -2.0, ["P2"]],
            ["A_vs_B", "s", "GO:BP", "GO:3", "only", 0.001, -1.0, ["P3"]],
        ])
        res = axis_meta_go.run(long, p_thr=0.05)
        d = res["directionality"]
        assert d.loc["GO:1", "pattern"] == "shared up"
        assert d.loc["GO:2", "pattern"] == "mixed"
        assert d.loc["GO:3", "pattern"] == "specific down"

    def test_redundant_terms_collapse_to_one_representative(self):
        """Nested GO terms hitting the same proteins must not fill the heatmap."""
        genes = [f"P{i}" for i in range(10)]
        rows = []
        for i in range(4):                       # 4 terms, ~90% shared genes
            rows.append(["A_vs_B", "s", "GO:BP", f"GO:{i}", f"purine term {i}",
                         0.001, 1.0 + i * 0.1, genes[:9] + [f"X{i}"]])
        rows.append(["A_vs_B", "s", "GO:BP", "GO:9", "unrelated", 0.001, 1.0,
                     [f"Z{i}" for i in range(10)]])
        res = axis_meta_go.run(_go_long(rows), p_thr=0.05, jaccard_cut=0.5)
        meta = res["meta"]
        assert int(meta["is_representative"].sum()) == 2      # one per family
        assert meta["redundancy_cluster"].nunique() == 2


class TestContrastSimilarity:
    def test_pair_relation_detects_shared_reference(self):
        rel, why = axis_meta_go._pair_relation("A_vs_C", "B_vs_C")
        assert rel == "shared reference" and "C" in why

    def test_pair_relation_detects_flipped_reference(self):
        rel, _ = axis_meta_go._pair_relation("C_vs_A", "B_vs_C")
        assert rel == "flipped reference"

    def test_pair_relation_disjoint(self):
        rel, _ = axis_meta_go._pair_relation("A_vs_B", "C_vs_D")
        assert rel == "disjoint"

    def test_similarity_recovers_injected_structure(self):
        rng = np.random.default_rng(3)
        profile = rng.normal(0, 1.5, 40)
        rows = []
        for c, z in [("A_vs_B", profile), ("C_vs_D", profile),      # identical
                     ("E_vs_F", -profile)]:                          # opposite
            for i, zz in enumerate(z):
                rows.append([c, "s", "GO:BP", f"GO:{i}", f"t{i}", 0.001, zz, ["P1"]])
        res = axis_meta_go.run(_go_long(rows), p_thr=0.05)
        cc = res["contrast_corr"]
        assert cc.loc["A_vs_B", "C_vs_D"] > 0.95
        assert cc.loc["A_vs_B", "E_vs_F"] < -0.95

    def test_too_few_shared_terms_gives_nan(self):
        rows = [["A_vs_B", "s", "GO:BP", "GO:1", "t1", 0.001, 1.0, ["P1"]],
                ["C_vs_D", "s", "GO:BP", "GO:2", "t2", 0.001, 1.0, ["P1"]]]
        res = axis_meta_go.run(_go_long(rows), p_thr=0.05)
        assert pd.isna(res["contrast_corr"].loc["A_vs_B", "C_vs_D"])


# --------------------------------------------------------------------------- #
# Exploration                                                                  #
# --------------------------------------------------------------------------- #
class TestExploration:
    def test_pca_and_association_on_real_design(self, workbook_full):
        dmap = io_layer.load_design(workbook_full["design"])
        stats = io_layer.load_stats(workbook_full["path"], design_map=dmap)
        res = axis_explore.run(stats, do_umap=False)
        assoc = res["pc_factor"].set_index("PC")
        assert assoc.loc["PC1", "epsilon2"] > 0.5
        assert assoc.loc["PC1", "p_adj_BH"] < 0.1

    def test_epsilon2_never_negative(self, workbook_full):
        stats = io_layer.load_stats(workbook_full["path"],
                                    design_map=io_layer.load_design(workbook_full["design"]))
        res = axis_explore.run(stats, do_umap=False)
        assert (res["pc_factor"]["epsilon2"] >= 0).all()
        assert (res["pc_factor"]["epsilon2"] <= 1).all()

    def test_bh_adjusted_p_is_present_and_monotone(self, workbook_full):
        stats = io_layer.load_stats(workbook_full["path"],
                                    design_map=io_layer.load_design(workbook_full["design"]))
        res = axis_explore.run(stats, do_umap=False)
        t = res["pc_factor"].dropna(subset=["p_value"])
        assert (t["p_adj_BH"] >= t["p_value"] - 1e-12).all()

    def test_no_replicates_yields_uncomputable_test(self, tmp_path):
        rng = np.random.default_rng(1)
        conds = ["A", "B", "C", "D"]
        n = 100
        comp = pd.DataFrame({"Protein.Group": [f"P{i}" for i in range(n)]})
        for c in ["A_vs_B"]:
            comp[f"{c}_p.val"] = rng.uniform(0, 1, n)
            comp[f"{c}_p.adj"] = rng.uniform(0, 1, n)
            comp[f"{c}_diff"] = rng.normal(0, 1, n)
        for s in conds:                       # one sample per condition
            comp[s] = rng.normal(22, 2, n)
        p = tmp_path / "norep.xlsx"
        with pd.ExcelWriter(p) as xw:
            comp.to_excel(xw, sheet_name="Comparaison", index=False)
        stats = io_layer.load_stats(str(p))
        res = axis_explore.run(stats, do_umap=False)
        assert res["pc_factor"]["p_value"].isna().all()

    def test_correlation_flags_zscored_input(self, workbook_zscore_only):
        stats = io_layer.load_stats(workbook_zscore_only["path"])
        res = axis_explore.run(stats, do_umap=False)
        assert res["correlation"].attrs.get("zscored") is True

    def test_correlation_on_log2_is_not_flagged(self, workbook_full):
        stats = io_layer.load_stats(workbook_full["path"])
        res = axis_explore.run(stats, do_umap=False)
        assert res["correlation"].attrs.get("zscored") is False

    def test_pca_refuses_a_matrix_without_samples(self):
        df = pd.DataFrame({"a": [1.0, 2.0], "b": [3.0, 4.0]})
        with pytest.raises(ValueError, match="sample column"):
            axis_explore.pca(df, ["a"], {"a": "X"})

    def test_multifactor_separates_condition_and_technical_drift(self, workbook_full, tmp_path):
        d = pd.read_csv(workbook_full["design"])
        d["injection_order"] = range(1, len(d) + 1)
        alt = tmp_path / "design_multi.csv"
        d.to_csv(alt, index=False)
        stats = io_layer.load_stats(workbook_full["path"],
                                    design_map=io_layer.load_design(str(alt)))
        frame, label_col = io_layer.load_design_frame(str(alt))
        factors = io_layer.map_design_factors(stats.sample_int_cols, frame, label_col)
        res = axis_explore.run(stats, do_umap=False, factors_df=factors)
        multi = res["pc_factor_multi"]
        assert set(multi["factor"]) == {"condition", "replicate", "injection_order"}
        num = multi[multi["factor"] == "injection_order"]
        assert (num["test"] == "Spearman").all()
        cat = multi[multi["factor"] == "condition"]
        assert (cat["test"] == "Kruskal-Wallis").all()


# --------------------------------------------------------------------------- #
# Publications                                                                 #
# --------------------------------------------------------------------------- #
class TestPublications:
    TSV = ("category\tterm\tnumber_of_genes\tinputGenes\tpreferredNames\tp_value\tfdr\tdescription\n"
           "PMID\tPMID:111\t8\tG1\tGene1\t1e-9\t3e-7\tA relevant paper\n"
           "KEGG\thsa04110\t5\tG2\tGene2\t1e-4\t2e-3\tCell cycle\n"
           "PMID\tPMID:222\t6\tG3\tGene3\t5e-6\t8e-5\tAnother paper\n")

    @pytest.fixture
    def stats_edge(self, tmp_path):
        """Contrasts with 600 / 1 / 0 / 50 significant proteins."""
        rng = np.random.default_rng(7)
        n = 1000
        comp = pd.DataFrame({"Protein.Group": [f"P{i:04d}" for i in range(n)],
                             "Genes": [f"G{i:04d}" for i in range(n)]})
        for c, k in {"A_vs_B": 600, "A_vs_C": 1, "A_vs_D": 0, "B_vs_C": 50}.items():
            p = np.ones(n)
            p[:k] = 1e-6
            comp[f"{c}_p.adj"] = p
            comp[f"{c}_p.val"] = p
            comp[f"{c}_diff"] = rng.normal(0, 2, n)
        comp["significant"] = True
        for s in [f"{c}_{i}" for c in "ABCD" for i in (1, 2, 3)]:
            comp[s] = rng.normal(22, 2, n)
        p = tmp_path / "pubs.xlsx"
        with pd.ExcelWriter(p) as xw:
            comp.to_excel(xw, sheet_name="Comparaison", index=False)
        return io_layer.load_stats(str(p))

    def test_only_pmid_rows_are_kept(self):
        rows = axis_publications._parse_pmid(self.TSV, top=15)
        assert [r["pmid"] for r in rows] == ["111", "222"]      # KEGG row dropped
        assert rows[0]["url"].endswith("/111/")

    def test_rows_are_sorted_by_fdr(self):
        rows = axis_publications._parse_pmid(self.TSV, top=15)
        assert rows[0]["fdr"] <= rows[1]["fdr"]

    def test_tiny_gene_sets_are_skipped(self, stats_edge, monkeypatch, tmp_path):
        """REGRESSION: a 1-gene 'enrichment' is just a list of papers citing it."""
        monkeypatch.setattr(cache, "CACHE_DIR", tmp_path / "cache")
        calls = []
        pub = axis_publications.run(stats_edge, taxid=3708,
                                    fetcher=lambda ids, tax: calls.append(len(ids)) or self.TSV)
        assert "A_vs_C" not in pub["results"]
        assert "skipped" in pub["provenance"]["A_vs_C"]
        assert all(n >= axis_publications.MIN_GENES for n in calls)

    def test_no_significant_protein_is_reported(self, stats_edge, monkeypatch, tmp_path):
        monkeypatch.setattr(cache, "CACHE_DIR", tmp_path / "cache")
        pub = axis_publications.run(stats_edge, taxid=3708, fetcher=lambda i, t: self.TSV)
        assert "no protein" in pub["provenance"]["A_vs_D"]

    def test_capping_is_visible(self, stats_edge, monkeypatch, tmp_path):
        monkeypatch.setattr(cache, "CACHE_DIR", tmp_path / "cache")
        pub = axis_publications.run(stats_edge, taxid=3708, max_ids=400,
                                    fetcher=lambda i, t: self.TSV)
        assert "400/600" in pub["provenance"]["A_vs_B"]
        assert "capped" in pub["provenance"]["A_vs_B"]

    def test_network_failure_is_graceful(self, stats_edge, monkeypatch, tmp_path):
        monkeypatch.setattr(cache, "CACHE_DIR", tmp_path / "cache")
        def boom(ids, tax):
            raise ConnectionError("no network")
        pub = axis_publications.run(stats_edge, taxid=3708, fetcher=boom)
        assert pub["results"] == {}
        assert "network failure" in pub["provenance"]["A_vs_B"]

    def test_cache_prevents_a_second_call(self, stats_edge, monkeypatch, tmp_path):
        monkeypatch.setattr(cache, "CACHE_DIR", tmp_path / "cache")
        calls = []
        f = lambda ids, tax: calls.append(1) or self.TSV
        axis_publications.run(stats_edge, taxid=3708, fetcher=f)
        first = len(calls)
        axis_publications.run(stats_edge, taxid=3708, fetcher=f)
        assert len(calls) == first, "cache did not prevent refetching"

    def test_refresh_bypasses_the_cache(self, stats_edge, monkeypatch, tmp_path):
        monkeypatch.setattr(cache, "CACHE_DIR", tmp_path / "cache")
        calls = []
        f = lambda ids, tax: calls.append(1) or self.TSV
        axis_publications.run(stats_edge, taxid=3708, fetcher=f)
        first = len(calls)
        axis_publications.run(stats_edge, taxid=3708, fetcher=f, refresh=True)
        assert len(calls) > first
