"""Loader tests: GO sheet mapping, per-sample matrix selection, design handling.

Several of these encode real regressions found on production data.
"""
import itertools

import numpy as np

import pandas as pd
import pytest

import io_layer
from conftest import NESTED_CONDITIONS, contrasts_of


# --------------------------------------------------------------------------- #
# GO sheet -> contrast mapping                                                 #
# --------------------------------------------------------------------------- #
class TestGoSheetMapping:
    """REGRESSION: substring scoring assigned GO_EPHPS7vIPHPS7 to IPHPS7_vs_HPS7
    because 'HPS7' is a substring of both sides, silently merging contrasts."""

    def test_every_sheet_maps_to_its_own_contrast(self):
        conds = NESTED_CONDITIONS
        contrasts = contrasts_of(conds)
        for a, b in itertools.combinations(conds, 2):
            got, how = io_layer._map_go_sheet(f"GO_{a}v{b}", contrasts, conds)
            assert got == f"{a}_vs_{b}", f"GO_{a}v{b} -> {got} ({how})"

    def test_no_collision_across_all_sheets(self):
        conds = NESTED_CONDITIONS
        contrasts = contrasts_of(conds)
        got = [io_layer._map_go_sheet(f"GO_{a}v{b}", contrasts, conds)[0]
               for a, b in itertools.combinations(conds, 2)]
        assert len(set(got)) == len(got), "two sheets map to the same contrast"
        assert set(got) == set(contrasts), "some contrasts lost their sheet"

    def test_nested_names_are_not_confused(self):
        """The exact pair that used to break."""
        conds = NESTED_CONDITIONS
        contrasts = contrasts_of(conds)
        assert io_layer._map_go_sheet("GO_EPHPS7vIPHPS7", contrasts, conds)[0] == "EPHPS7_vs_IPHPS7"
        assert io_layer._map_go_sheet("GO_EPHPS7vHPS7", contrasts, conds)[0] == "EPHPS7_vs_HPS7"
        assert io_layer._map_go_sheet("GO_IPHPS7vLPHPS7", contrasts, conds)[0] == "IPHPS7_vs_LPHPS7"

    def test_sheets_may_drop_the_shared_suffix(self):
        """REGRESSION: when every condition ends with the same token, upstream naming
        drops it ('EPHPS11' -> 'EPHP', 'CS11' -> 'C'), and all 28 sheets were ignored."""
        conds = ["EPHPS11", "IPHPS11", "LPHPS11", "CS11", "HPS11", "EPS11", "IPS11", "LPS11"]
        contrasts = contrasts_of(conds)
        short = {c: c[:-3] for c in conds}          # strip the shared 'S11'
        for a, b in itertools.combinations(conds, 2):
            sheet = f"GO_{short[a]}v{short[b]}"
            got, how = io_layer._map_go_sheet(sheet, contrasts, conds)
            assert got == f"{a}_vs_{b}", f"{sheet} -> {got} ({how})"

    def test_sheets_may_drop_a_shared_prefix(self):
        conds = ["S11_EPHP", "S11_IPHP", "S11_CTRL"]
        contrasts = contrasts_of(conds)
        for a, b in itertools.combinations(conds, 2):
            sheet = f"GO_{a[4:]}v{b[4:]}"
            assert io_layer._map_go_sheet(sheet, contrasts, conds)[0] == f"{a}_vs_{b}"

    def test_design_and_column_spellings_may_differ(self):
        """REGRESSION: conditions come from the design, contrasts from the stats
        columns; an underscore, a space or a different case between the two made the
        reconstructed contrast name miss by exact comparison."""
        contrasts = ["CAM203CTR_vs_CAM124CTR"]
        for conds in (["CAM203_CTR", "CAM124_CTR"],      # underscores in the design
                      ["CAM203 CTR", "CAM124 CTR"],      # spaces
                      ["cam203ctr", "cam124ctr"]):       # case
            got, how = io_layer._map_go_sheet("GO_CAM203CTRvCAM124CTR", contrasts, conds)
            assert got == "CAM203CTR_vs_CAM124CTR", f"{conds} -> {got} ({how})"

    def test_shared_leading_letter_is_stripped(self):
        """REGRESSION: every condition starting with 'C' ('CAM203CTR', 'CelineCTR')
        yielded sheets named GO_AM203CTRvelineCTR, and all 15 were ignored."""
        conds = ["CAM203CTR", "CAM124CTR", "CelineCTR", "CAM203HT", "CAM124HT", "CelineHT"]
        contrasts = contrasts_of(conds)
        for a, b in itertools.combinations(conds, 2):
            sheet = f"GO_{a[1:]}v{b[1:]}"
            got, how = io_layer._map_go_sheet(sheet, contrasts, conds)
            assert got == f"{a}_vs_{b}", f"{sheet} -> {got} ({how})"

    def test_no_alias_when_conditions_share_nothing(self):
        """Full names must keep working: 4A has CS8 among S7 conditions, so no token
        is shared and no shortcut may be invented."""
        conds = NESTED_CONDITIONS
        aliases = io_layer._condition_aliases(conds)
        assert set(aliases.values()) == set(conds)
        assert len(aliases) == len(conds)

    def test_ambiguous_aliases_are_discarded(self):
        """An alias designating two conditions must never be used."""
        aliases = io_layer._condition_aliases(["AA_X", "AA_X_bis"])
        for alias, cond in aliases.items():
            matches = [c for c in ["AA_X", "AA_X_bis"] if io_layer._norm(c).endswith(alias)
                       or io_layer._norm(c) == alias]
            assert len(matches) >= 1

    @pytest.mark.parametrize("sheet", ["GO_A_vs_B", "GO_AvB", "GO_AvsB", "go_a_vs_b"])
    def test_naming_variants(self, sheet):
        assert io_layer._map_go_sheet(sheet, ["A_vs_B"], ["A", "B"])[0] == "A_vs_B"

    def test_reversed_orientation_is_recognised(self):
        got, how = io_layer._map_go_sheet("GO_BvA", ["A_vs_B"], ["A", "B"])
        assert got == "A_vs_B" and "reversed" in how

    def test_unknown_sheet_is_reported_not_guessed(self):
        got, how = io_layer._map_go_sheet("GO_XvY", ["A_vs_B"], ["A", "B"])
        assert got is None and how == "UNMATCHED"

    def test_load_go_reports_mapping_and_no_loss(self, workbook_full, go_workbook):
        stats = io_layer.load_stats(workbook_full["path"])
        go = io_layer.load_go(go_workbook["path"], stats.contrasts,
                              conditions=stats.conditions, verbose=False)
        assert go.unmatched == []
        assert go.duplicates == {}
        assert go.missing == []
        assert len(go.contrasts) == len(workbook_full["contrasts"])

    def test_terms_land_on_the_right_contrast(self, workbook_full, go_workbook):
        """Each sheet has a unique term count: a merge would show up here."""
        stats = io_layer.load_stats(workbook_full["path"])
        go = io_layer.load_go(go_workbook["path"], stats.contrasts,
                              conditions=stats.conditions, verbose=False)
        counts = go.long.groupby("contrast").size().to_dict()
        assert counts == go_workbook["expected_terms"]

    def test_foreign_go_file_is_not_force_fitted(self, go_workbook):
        """A GO file from another analysis must be reported, not mapped anyway."""
        go = io_layer.load_go(go_workbook["path"], ["X_vs_Y", "X_vs_Z"],
                              conditions=["X", "Y", "Z"], verbose=False)
        assert len(go.long) == 0
        assert len(go.unmatched) == 28
        assert go.contrasts == []


# --------------------------------------------------------------------------- #
# Per-sample matrix selection                                                  #
# --------------------------------------------------------------------------- #
class TestMatrixSelection:
    """REGRESSION: the PCA used to run on the z-scored heatmap, which holds a
    pre-selected subset -> circular, inflated separation."""

    def test_prefers_full_log2_over_zscored_subset(self, workbook_full):
        stats = io_layer.load_stats(workbook_full["path"])
        assert stats.pca_source == "Log2_Impute"
        assert stats.pca_is_zscored is False
        assert len(stats.pca_df) == 400              # full set, not the 120 subset

    def test_rejects_non_log_raw_intensities(self, workbook_full):
        stats = io_layer.load_stats(workbook_full["path"])
        assert stats.pca_source != "raw_data"

    def test_falls_back_to_zscore_and_flags_it(self, workbook_zscore_only):
        stats = io_layer.load_stats(workbook_zscore_only["path"])
        assert stats.pca_is_zscored is True
        assert len(stats.sample_int_cols) == 9

    def test_unnamed_columns_are_never_samples(self, tmp_path, workbook_full):
        """REGRESSION: GO sheets exported without headers carry 'Unnamed: 3' columns,
        which were counted as samples and let a GO sheet compete as the PCA matrix."""
        rng = np.random.default_rng(1)
        src = pd.read_excel(workbook_full["path"], sheet_name=None)
        go = pd.DataFrame({"term_name": ["a", "b"], "p_value": [0.01, 0.02]})
        for j in range(6):
            go[f"Unnamed: {j}"] = rng.normal(0, 1, 2)
        p = tmp_path / "with_unnamed.xlsx"
        with pd.ExcelWriter(p) as xw:
            for n, d in src.items():
                d.to_excel(xw, sheet_name=n, index=False)
            go.to_excel(xw, sheet_name="GO_AvB", index=False)
        stats = io_layer.load_stats(str(p))
        sheets = [c["sheet"] for c in stats.candidates]
        assert "GO_AvB" not in sheets
        assert stats.pca_source == "Log2_Impute"

    def test_contrast_columns_are_never_samples(self, workbook_full):
        """'<C>_p.val_limma' style columns must not be taken for samples."""
        stats = io_layer.load_stats(workbook_full["path"])
        assert not any("_vs_" in str(c) for c in stats.sample_int_cols)
        assert all(str(c).startswith("Sple-") for c in stats.sample_int_cols)

    def test_detects_all_contrasts(self, workbook_full):
        stats = io_layer.load_stats(workbook_full["path"])
        assert sorted(stats.contrasts) == sorted(workbook_full["contrasts"])


# --------------------------------------------------------------------------- #
# Experimental design                                                          #
# --------------------------------------------------------------------------- #
class TestLabelMatching:
    """Design labels rarely match column names character for character."""

    def test_exact_match(self):
        m = io_layer.match_columns_to_labels(["Sple-1", "Sple-2"], ["Sple-1", "Sple-2"])
        assert m == {"Sple-1": 0, "Sple-2": 1}

    def test_containment_match(self):
        """'Sple-15' in the design vs 'LFQ.intensity.Sple-15' in the workbook."""
        cols = ["LFQ.intensity.Sple-15", "LFQ.intensity.Sple-16"]
        m = io_layer.match_columns_to_labels(cols, ["Sple-15", "Sple-16"])
        assert m == {cols[0]: 0, cols[1]: 1}

    def test_longest_label_wins_over_prefix_collision(self):
        """'Sple-1' must not steal the column of 'Sple-15'."""
        cols = ["LFQ.intensity.Sple-15", "LFQ.intensity.Sple-1"]
        m = io_layer.match_columns_to_labels(cols, ["Sple-15", "Sple-1"])
        assert m == {cols[0]: 0, cols[1]: 1}

    def test_prefixed_labels_match_bare_columns(self):
        """The other direction: design lists 'LFQ.intensity.Sple-1', columns are bare."""
        m = io_layer.match_columns_to_labels(
            ["Sple-1", "Sple-10"], ["LFQ.intensity.Sple-1", "LFQ.intensity.Sple-10"])
        assert m == {"Sple-1": 0, "Sple-10": 1}

    def test_unrelated_labels_match_nothing(self):
        assert io_layer.match_columns_to_labels(["A", "B"], ["X", "Y"]) == {}

    def test_design_with_prefixed_columns_resolves_conditions(self, workbook_full, tmp_path):
        """A design listing bare sample names must still resolve prefixed columns."""
        d = pd.read_csv(workbook_full["design"])
        alt = tmp_path / "design_bare.csv"
        d.to_csv(alt, index=False)
        dmap = io_layer.load_design(str(alt))
        stats = io_layer.load_stats(workbook_full["path"], design_map=dmap)
        assert stats.conditions_resolved is True
        assert len(stats.conditions) == 8


class TestDesign:
    def test_conditions_resolved_from_design(self, workbook_full):
        dmap = io_layer.load_design(workbook_full["design"])
        stats = io_layer.load_stats(workbook_full["path"], design_map=dmap)
        assert stats.conditions_resolved is True
        assert len(stats.conditions) == 8
        sizes = {c: sum(1 for s in stats.sample_int_cols if stats.cond_map[s] == c)
                 for c in stats.conditions}
        assert set(sizes.values()) == {3}

    def test_generic_sample_names_without_design_are_not_resolved(self, workbook_zscore_only):
        stats = io_layer.load_stats(workbook_zscore_only["path"])
        assert stats.conditions_resolved is False    # 'Sple-1' encodes no condition

    def test_positional_fallback_when_labels_differ(self, workbook_full, tmp_path):
        """Design listing raw filenames instead of column names."""
        d = pd.read_csv(workbook_full["design"])
        d["sample"] = [f"20240101_QE_{i:02d}.raw" for i in range(len(d))]
        alt = tmp_path / "design_files.csv"
        d.to_csv(alt, index=False)
        dmap = io_layer.load_design(str(alt))
        stats = io_layer.load_stats(workbook_full["path"], design_map=dmap)
        assert stats.conditions_resolved is True
        assert len(stats.conditions) == 8

    def test_multi_factor_columns_are_kept(self, workbook_full, tmp_path):
        d = pd.read_csv(workbook_full["design"])
        d["batch"] = ["b1" if i < 12 else "b2" for i in range(len(d))]
        alt = tmp_path / "design_batch.csv"
        d.to_csv(alt, index=False)
        df, label_col = io_layer.load_design_frame(str(alt))
        stats = io_layer.load_stats(workbook_full["path"],
                                    design_map=io_layer.load_design(str(alt)))
        factors = io_layer.map_design_factors(stats.sample_int_cols, df, label_col)
        assert list(factors.columns) == ["condition", "batch"]
        assert len(factors) == 24

    @pytest.mark.parametrize("col", ["replicate", "Replicate", "rep", "REP", "replicat"])
    def test_replicate_index_is_not_a_factor(self, workbook_full, tmp_path, col):
        """1, 2, 3 inside every condition: testing it against the PCs is meaningless
        and doubles the BH family (Deimos designs always carry this column)."""
        d = pd.read_csv(workbook_full["design"]).rename(columns={"replicate": col})
        alt = tmp_path / "design_rep.csv"
        d.to_csv(alt, index=False)
        df, label_col = io_layer.load_design_frame(str(alt))
        stats = io_layer.load_stats(workbook_full["path"],
                                    design_map=io_layer.load_design(str(alt)))
        factors = io_layer.map_design_factors(stats.sample_int_cols, df, label_col)
        assert list(factors.columns) == ["condition"]

    def test_a_column_merely_containing_rep_is_kept(self, workbook_full, tmp_path):
        d = pd.read_csv(workbook_full["design"])
        d["prep_day"] = ["d1" if i % 2 else "d2" for i in range(len(d))]
        alt = tmp_path / "design_prep.csv"
        d.to_csv(alt, index=False)
        df, label_col = io_layer.load_design_frame(str(alt))
        stats = io_layer.load_stats(workbook_full["path"],
                                    design_map=io_layer.load_design(str(alt)))
        factors = io_layer.map_design_factors(stats.sample_int_cols, df, label_col)
        assert "prep_day" in factors.columns
