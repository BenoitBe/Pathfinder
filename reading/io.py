"""Data of the biological reading, taken from what Pathfinder has already loaded:
differential results (comparison sheet), contrasts, Pathfinder's selection rule (= the
Deimos volcano rule), and, next to the workbook, what Deimos writes in its output
folder: matrices/ (DPC or imputed log2 abundances, samples.tsv), PCA_Scores,
last_config.yaml.

Everything is optional except the differential results: each figure checks what it
needs and switches itself off cleanly otherwise.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np
import pandas as pd


@dataclass
class Run:
    out_dir: str                           # folder of the workbook (= Deimos output)
    xlsx: str
    de: pd.DataFrame                       # one row per protein, index = protein id
    contrasts: list                        # manifest order
    pairs: dict                            # contrast -> (numerator, denominator)
    rule: object = None                    # Pathfinder's io_layer.SelectionRule
    y_axis: str = "pvalue"                 # Deimos volcano_y_axis
    samples: pd.DataFrame | None = None    # label, condition, (factors)
    conditions: list = field(default_factory=list)
    means: pd.DataFrame | None = None      # proteins x conditions (log2)
    abundance: pd.DataFrame | None = None  # proteins x samples (log2)
    observed: pd.DataFrame | None = None   # proteins x samples (1/0)
    abundance_kind: str = ""
    scores: pd.DataFrame | None = None     # PCA_Scores
    labels: pd.Series | None = None        # protein id -> readable name
    label_source: str = ""
    config: dict = field(default_factory=dict)

    # --- selection rule: Pathfinder's, hence the one of the Deimos volcanos ---
    @property
    def use_padj(self) -> bool:
        return self.rule.p_column == "p.adj"

    @property
    def p_thr(self) -> float:
        return float(self.rule.p_thr)

    @property
    def ratio_min(self) -> float:
        return float(self.rule.ratio_min)

    @property
    def lfc_min(self) -> float:
        return float(self.rule.lfc_min)

    @property
    def rule_source(self) -> str:
        return str(self.rule.source)

    def sig(self, c: str) -> pd.Series:
        """Proteins significant in a contrast: SelectionRule.mask, the same selection as
        the Proteins tab and as Deimos' volcano_status."""
        col = c + ("_p.adj" if self.use_padj else "_p.val")
        if col not in self.de.columns:
            return pd.Series(False, index=self.de.index)
        return pd.Series(self.rule.mask(self.de[col].values, self.de[c + "_diff"].values),
                         index=self.de.index)

    def rule_text(self) -> str:
        return self.rule.text()

    def short(self, c: str) -> str:
        num, den = self.pairs.get(c, (None, None))
        return f"{num} vs {den}" if num else str(c).replace("_vs_", " vs ")


def _read_tsv(path: str) -> pd.DataFrame | None:
    return pd.read_csv(path, sep="\t") if os.path.exists(path) else None


def _config(out_dir: str) -> dict:
    p = os.path.join(out_dir, "last_config.yaml")
    if not os.path.exists(p):
        return {}
    try:
        import yaml
        with open(p, encoding="utf-8") as fh:
            return yaml.safe_load(fh) or {}
    except Exception:                                                 # noqa: BLE001
        return {}


def _de_ids(comp: pd.DataFrame) -> pd.Series:
    """Protein identifier: Deimos' `name` (the one of the GO intersections and of
    matrices/), else the head of the Protein.Group."""
    if "name" in comp.columns:
        ids = comp["name"].astype(str)
    else:
        col = "Protein.Group" if "Protein.Group" in comp.columns else comp.columns[0]
        ids = comp[col].astype(str).str.split(";").str[0].str.strip()
    return ids


def id_lookup(de: pd.DataFrame) -> dict:
    """Token -> identifier of `de`: the identifier itself and every accession of the
    Protein.Group. Links the GO intersections and the matrices to the proteins."""
    out = {str(i): str(i) for i in de.index}
    if "Protein.Group" in de.columns:
        for pid, pg in zip(de.index, de["Protein.Group"].astype(str)):
            for tok in pg.replace(",", ";").split(";"):
                tok = tok.strip()
                if tok and tok not in out:
                    out[tok] = str(pid)
    return out


def _index_on_de(mat: pd.DataFrame, id_col: str, lookup: dict) -> pd.DataFrame:
    ids = mat[id_col].astype(str)
    mapped = ids.map(lambda s: lookup.get(s) or lookup.get(s.split(";")[0].strip()))
    out = mat.loc[mapped.notna()].copy()
    out.index = mapped[mapped.notna()].values
    return out[~out.index.duplicated()]


def load_run(stats, stats_path: str, rule, manifest: dict | None = None,
             reference_fasta: str | None = None, verbose: bool = True) -> Run:
    """stats: io_layer.StatsData; rule: io_layer.SelectionRule; manifest:
    io_layer.read_manifest(...)."""
    out_dir = os.path.dirname(os.path.abspath(stats_path))
    manifest = manifest or {}
    de = stats.comp.copy()
    de["protein"] = _de_ids(de).values
    de = de.drop_duplicates("protein").set_index("protein", drop=False)

    contrasts = [c for c in stats.contrasts if c + "_diff" in de.columns]
    declared = [c for c in manifest.get("contrast", {}) if c in contrasts]
    contrasts = declared + [c for c in contrasts if c not in declared]
    pairs = {}
    for c, v in manifest.get("contrast", {}).items():
        if "|" in str(v):
            pairs[c] = tuple(str(v).split("|", 1))
    for c in contrasts:
        if c not in pairs and "_vs_" in c:
            pairs[c] = tuple(c.split("_vs_", 1))

    cfg = _config(out_dir)
    run = Run(out_dir=out_dir, xlsx=stats_path, de=de, contrasts=contrasts, pairs=pairs,
              rule=rule, config=cfg,
              y_axis=str(cfg.get("volcano_y_axis", "pvalue")).lower() if cfg else "pvalue")
    xl = pd.ExcelFile(stats_path)
    lookup = id_lookup(de)

    # --- samples / conditions: samples.tsv, else PCA_Scores, else Pathfinder's mapping ---
    mdir = os.path.join(out_dir, "matrices")
    samples = _read_tsv(os.path.join(mdir, "samples.tsv"))
    if samples is None and "PCA_Scores" in xl.sheet_names:
        s = pd.read_excel(xl, sheet_name="PCA_Scores")
        if {"label", "condition"} <= set(s.columns):
            samples = s[["label", "condition"]].copy()
    if samples is None and stats.conditions_resolved and stats.sample_int_cols:
        samples = pd.DataFrame({"label": [str(c) for c in stats.sample_int_cols],
                                "condition": [str(stats.cond_map[c]) for c in stats.sample_int_cols]})
    if samples is not None:
        samples["label"] = samples["label"].astype(str)
        samples["condition"] = samples["condition"].astype(str)
        run.samples = samples
        run.conditions = list(dict.fromkeys(samples["condition"]))

    # --- abundances: DPC, else imputed log2 (file, then sheet), else the matrix of
    #     Pathfinder's PCA when it is not z-scored ---
    ab, kind, pid = _read_tsv(os.path.join(mdir, "dpc_abundance.tsv")), "DPC", "protein"
    if ab is None:
        ab, kind = _read_tsv(os.path.join(mdir, "log2_imputed.tsv")), "log2 imputed"
    if ab is None and "Log2_Impute" in xl.sheet_names:
        ab, kind, pid = pd.read_excel(xl, sheet_name="Log2_Impute"), "log2 imputed", "Protein.Group"
    if ab is None and not stats.pca_is_zscored and stats.pca_df is not None:
        ab, kind, pid = stats.pca_df.copy(), f"log2 ({stats.pca_source})", stats.id_col
    if ab is not None and samples is not None:
        if pid not in ab.columns:
            pid = ab.columns[0]
        cols = [c for c in samples["label"] if c in ab.columns]
        if cols:
            A = _index_on_de(ab, pid, lookup)[cols].apply(pd.to_numeric, errors="coerce")
            run.abundance, run.abundance_kind = A, kind
            ob = _read_tsv(os.path.join(mdir, "dpc_observed.tsv")) if kind == "DPC" else None
            if ob is not None and "protein" in ob.columns:
                run.observed = _index_on_de(ob, "protein", lookup)[cols].fillna(0).astype(int)
            cond = samples.drop_duplicates("label").set_index("label").loc[cols, "condition"]
            run.means = A.T.groupby(cond.values).mean().T.reindex(columns=run.conditions)
            run.means = run.means.loc[:, run.means.notna().any()]

    if "PCA_Scores" in xl.sheet_names:
        run.scores = pd.read_excel(xl, sheet_name="PCA_Scores")

    from .labels import best_labels
    run.labels, run.label_source = best_labels(run, xl, reference_fasta)
    if verbose:
        print(f"[Pathfinder] reading: {len(de)} proteins, {len(run.contrasts)} contrasts, "
              f"{len(run.conditions)} conditions, abundance: {run.abundance_kind or 'none'}, "
              f"names: {run.label_source}")
    return run
