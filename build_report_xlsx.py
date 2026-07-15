"""Pathfinder synthesis workbook export (matrices + tables)."""
from __future__ import annotations
import pandas as pd


def build(path, meta_go, explore, provenance: dict, resume=None, pubs=None):
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        # provenance first
        pd.DataFrame(list(provenance.items()), columns=["key", "value"]).to_excel(
            xw, sheet_name="provenance", index=False)

        # publications (online axis)
        if pubs is not None:
            recs = []
            for c, rows in pubs.get("results", {}).items():
                for r in rows:
                    recs.append({"contrast": c, **r})
            if recs:
                pd.DataFrame(recs).to_excel(xw, sheet_name="Publications", index=False)

        # interpretive summary
        if resume is not None and not resume["lines"].empty:
            resume["lines"].to_excel(xw, sheet_name="Summary", index=False)

        # --- Meta-GO ---
        mat = meta_go["matrix"].copy()
        mmeta = meta_go["meta"]
        mat_out = mat.copy()
        mat_out.insert(0, "term_name", mmeta["term_name"])
        mat_out.insert(1, "source", mmeta["source"])
        mat_out.insert(2, "n_sig", mmeta["n_sig"])
        mat_out.insert(3, "redundancy_cluster", mmeta["redundancy_cluster"])
        mat_out.insert(4, "is_representative", mmeta["is_representative"])
        mat_out = mat_out.sort_values(["source", "n_sig"], ascending=[True, False])
        mat_out.to_excel(xw, sheet_name="MetaGO_matrix")
        direction = meta_go["directionality"]
        if "source" in direction.columns:
            direction = direction.sort_values(["source", "n_sig"], ascending=[True, False])
        direction.to_excel(xw, sheet_name="MetaGO_directionality")

        # --- Exploration ---
        explore["pca"]["scores"].to_excel(xw, sheet_name="PCA_scores")
        # loadings: top |PC1| contributors
        load = explore["pca"]["loadings"]
        if "PC1" in load.columns:
            load = load.reindex(load["PC1"].abs().sort_values(ascending=False).index)
        load.head(200).to_excel(xw, sheet_name="PCA_loadings")
        explore["pc_factor"].to_excel(xw, sheet_name="PC_factor_assoc", index=False)
        multi = explore.get("pc_factor_multi")
        if multi is not None and not multi.empty:
            multi.to_excel(xw, sheet_name="PC_factor_multi", index=False)
        explore["correlation"].to_excel(xw, sheet_name="Sample_correlation")

        # explained variance
        ve = explore["pca"]["var_exp"]
        pd.DataFrame({"PC": [f"PC{i+1}" for i in range(len(ve))],
                      "explained_variance_%": [round(float(v), 2) for v in ve]}).to_excel(
            xw, sheet_name="PCA_scree", index=False)
