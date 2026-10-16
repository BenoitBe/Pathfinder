"""Pathfinder synthesis workbook export (matrices + tables)."""
from __future__ import annotations
import numpy as np
import pandas as pd


def build(path, meta_go, explore, provenance: dict, resume=None, pubs=None, prots=None, wgcna=None):
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
        cc = meta_go.get("contrast_corr")
        if cc is not None and not cc.empty:
            cc.to_excel(xw, sheet_name="Contrast_similarity")

        # --- Proteins (which proteins carry the signature) ---
        if prots is not None and not prots["matrix"].empty:
            t = prots["table"]
            cols = [c for c in ("gene", "description", "n_sig", "pattern", "recurrence",
                                "best_pi", "median_pi", "median_robustness", "max_abs_lfc",
                                "module", "n_go_terms", "go_terms",
                                "up_in", "down_in", "opposite_in") if c in t.columns]
            t[cols].to_excel(xw, sheet_name="Proteins_core")
            m = prots["matrix"].copy()
            for c in ("gene", "n_sig", "pattern"):
                if c in t.columns:
                    m.insert(0, c, t[c].reindex(m.index))
            m.to_excel(xw, sheet_name="Proteins_matrix")
            cc = prots.get("contrast_corr")
            if cc is not None and not cc.empty:
                cc.to_excel(xw, sheet_name="Contrast_similarity_prot")

        # --- WGCNA (read from upstream) ---
        if wgcna:
            if wgcna["association"] is not None and not wgcna["association"].empty:
                wgcna["association"].to_excel(xw, sheet_name="WGCNA_module_condition", index=False)
            if wgcna["core_enrichment"] is not None and not wgcna["core_enrichment"].empty:
                wgcna["core_enrichment"].to_excel(xw, sheet_name="WGCNA_core_enrichment", index=False)
            if wgcna.get("contrast_association") is not None and not wgcna["contrast_association"].empty:
                wgcna["contrast_association"].to_excel(xw, sheet_name="WGCNA_module_contrast", index=False)
            if wgcna["hubs"] is not None and not wgcna["hubs"].empty:
                wgcna["hubs"].to_excel(xw, sheet_name="WGCNA_hubs", index=False)

        # --- Exploration ---
        explore["pca"]["scores"].to_excel(xw, sheet_name="PCA_scores")
        # loadings: top |PC1| contributors
        load = explore["pca"]["loadings"]
        # sort positionally: gene symbols repeat (isoforms, protein groups), and
        # reindexing on a duplicated index raises.
        if "PC1" in load.columns:
            order = np.argsort(-load["PC1"].abs().to_numpy())
            load = load.iloc[order]
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
