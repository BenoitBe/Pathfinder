"""
Pathfinder - shared input layer.

Reads the Excel workbooks produced by Deimos (DEP module) and the GO module, and
rebuilds the objects feeding every axis. No statistic is recomputed: read, extract,
reshape (the pipeline's founding rule).

Detection is CONTENT-based, not filename-based:
  - Stats : a sheet carrying '<C>_p.val' columns + a per-sample intensity matrix
  - GO    : 'GO_*' sheets carrying a 'term_name' column
"""
from __future__ import annotations
import os
import re
import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

warnings.simplefilter("ignore", category=UserWarning)

# Known metadata columns from the DEP module (never sample intensities)
_META = {
    "name", "Protein.Group", "Protein.Names", "Genes", "First.Protein.Description",
    "N.Sequences", "N.Proteotypic.Sequences", "ID", "imputed", "num_NAs", "significant",
}


# --------------------------------------------------------------------------- #
# Data structures                                                                   #
# --------------------------------------------------------------------------- #
@dataclass
class StatsData:
    comp: pd.DataFrame                    # comparison sheet (stats + accessions)
    zscore: Optional[pd.DataFrame]        # Z-score_Heatmap sheet (may be None)
    pca_df: pd.DataFrame                  # sheet selected for the PCA (per-sample matrix)
    pca_source: str                       # name of the PCA sheet
    pca_is_zscored: bool                  # True if that matrix is per-protein z-scored
    contrasts: list                       # ['A_vs_B', ...]
    sample_int_cols: list                 # per-sample intensity columns (within pca_df)
    comp_sample_cols: list                 # raw intensity columns in comp (for abundance)
    conditions: list                      # order of appearance
    cond_map: dict                        # {sample_column: condition}
    conditions_resolved: bool             # True when the condition mapping is reliable
    id_col: str                           # protein identifier column
    gene_col: Optional[str]
    desc_col: Optional[str]


@dataclass
class GOData:
    long: pd.DataFrame                    # contrast, source, term_id, term_name, p_value, z, intersection(list)
    p_col: str
    z_col: str
    contrasts: list = field(default_factory=list)


# --------------------------------------------------------------------------- #
# Helpers                                                                       #
# --------------------------------------------------------------------------- #
def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def _find_col(cols, *cands, contains_all=None):
    """Tolerant column lookup: normalised equality, then substring match."""
    ncols = {_norm(c): c for c in cols}
    for cand in cands:
        if _norm(cand) in ncols:
            return ncols[_norm(cand)]
    if contains_all:
        toks = [t.lower() for t in contains_all]
        for c in cols:
            cl = str(c).lower()
            if all(t in cl for t in toks):
                return c
    return None


def map_sample_to_condition(sample_cols, conditions):
    """Map each sample column to its condition.

    Prefers the longest matching prefix among known conditions (robust to names
    containing underscores), falling back to rsplit('_', 1).
    """
    conds_sorted = sorted(conditions, key=len, reverse=True)
    cond_map = {}
    for col in sample_cols:
        hit = None
        for cond in conds_sorted:
            if col == cond or col.startswith(cond + "_") or _norm(col).startswith(_norm(cond)):
                hit = cond
                break
        cond_map[col] = hit if hit is not None else col.rsplit("_", 1)[0]
    return cond_map


# --------------------------------------------------------------------------- #
# Statistics                                                                         #
# --------------------------------------------------------------------------- #
def load_design_frame(path_or_df):
    """Read a design and return (df, label_col) - every column is kept for
    multi-factor use (condition, treatment, batch, injection order...)."""
    if isinstance(path_or_df, pd.DataFrame):
        df = path_or_df.copy()
    elif str(path_or_df).lower().endswith((".csv", ".tsv", ".txt")):
        df = None
        for sep in (None, ",", ";", "\t"):
            try:
                d = pd.read_csv(path_or_df, sep=sep, engine="python")
                if d.shape[1] >= 2:
                    df = d
                    break
            except Exception:
                continue
        if df is None:
            df = pd.read_csv(path_or_df)
    else:
        df = pd.read_excel(path_or_df, sheet_name=0)
    label_col = _find_col(df.columns, "sample", "label", "run", "name", "file",
                          "file.name", "filename", "sample_name", "id", "raw") or df.columns[0]
    return df, label_col


def map_design_factors(sample_cols, df, label_col, cond_map=None):
    """Align sample columns with design rows (by name, else positionally) and return
    a DataFrame indexed by sample column, one column per factor."""
    labels = [str(x).strip() for x in df[label_col].tolist()]
    matched = match_columns_to_labels(sample_cols, labels)
    if len(matched) == len(sample_cols):
        rows = [matched[c] for c in sample_cols]
    elif len(labels) == len(sample_cols):
        rows = list(range(len(sample_cols)))          # positional fallback
    else:
        return None
    out = df.iloc[rows].copy()
    out.index = sample_cols
    keep = [c for c in out.columns if c != label_col]
    return out[keep]


def load_design(path_or_df):
    """Read an ExperimentalDesign (CSV/xlsx/DataFrame) -> {sample_label: condition}.

    Tolerant detection of the label and condition columns. Keys are indexed raw AND
    normalised to tolerate separator variants (Sple-1 / Sple_1).
    """
    if isinstance(path_or_df, pd.DataFrame):
        df = path_or_df
    elif str(path_or_df).lower().endswith((".csv", ".tsv", ".txt")):
        for sep in (None, ",", ";", "\t"):
            try:
                df = pd.read_csv(path_or_df, sep=sep, engine="python")
                if df.shape[1] >= 2:
                    break
            except Exception:
                continue
    else:
        df = pd.read_excel(path_or_df, sheet_name=0)

    label_col = _find_col(df.columns, "sample", "label", "run", "name", "file",
                          "file.name", "filename", "sample_name", "id", "raw")
    cond_col = _find_col(df.columns, "condition", "group", "groupe", "class",
                         "treatment", "type", "category", "cond")
    if label_col is None or cond_col is None:
        # heuristic: first column = label, least-unique column = condition
        label_col = label_col or df.columns[0]
        if cond_col is None:
            nunq = {c: df[c].nunique() for c in df.columns if c != label_col}
            cond_col = min(nunq, key=nunq.get) if nunq else df.columns[-1]

    mp = {}
    for _, r in df.iterrows():
        lab, cond = str(r[label_col]).strip(), str(r[cond_col]).strip()
        mp[lab] = cond
        mp[_norm(lab)] = cond
    mp["__order__"] = [str(x).strip() for x in df[cond_col].tolist()]   # positional fallback
    mp["__labels__"] = [str(x).strip() for x in df[label_col].tolist()]
    return mp


def load_stats(path: str, design_map=None) -> StatsData:
    xls = pd.read_excel(path, sheet_name=None)  # dict {sheet: df}

    # 1) sheet de comparaison = celle qui porte le plus de colonnes '_p.val'
    comp_sheet, comp_df, best = None, None, -1
    for name, df in xls.items():
        n = sum(1 for c in df.columns if str(c).endswith("_p.val"))
        if n > best:
            best, comp_sheet, comp_df = n, name, df
    if best <= 0:
        raise ValueError("Aucune sheet avec colonnes '<C>_p.val' : fichier Stats non reconnu.")

    contrasts = [str(c)[:-6] for c in comp_df.columns if str(c).endswith("_p.val")]

    # 2) sheet Z-score (facultative)
    z_df = None
    for name, df in xls.items():
        if "zscore" in _norm(name) or "heatmap" in _norm(name):
            z_df = df
            break

    # 3) POSITIVE identification of sample columns (per sheet)
    #    Excludes: metadata, contrast columns (containing '_vs_'), annotations.
    stat_cols = set()
    for c in contrasts:
        for suf in ("_diff", "_p.val", "_p.adj"):
            stat_cols.add(c + suf)
        stat_cols.add("Robustness_Score_" + c)
        stat_cols.add("Pi_Score_" + c)
    ANNOT = ("count", "coverage", "sequences", "proteotypic", "razor", "unique",
             "qvalue", "q.value", "imputed", "num_", "nas", "peptide", "score",
             "robust", "_pi", "pi_", "fdr", "logfc", "adj.p", "p.val", "p.adj",
             "intensity_mean", "mean_intensity", "cv", "rank")

    def _samples(df):
        out = []
        for c in df.columns:
            cs = str(c)
            if cs in _META or cs in stat_cols:
                continue
            # 'Unnamed: 3' is an artefact of a sheet exported without headers (GO sheets
            # do this): never a sample, and it must not compete as a PCA matrix.
            if re.match(r"^unnamed:?\s*\d*$", cs.strip().lower()):
                continue
            if not pd.api.types.is_numeric_dtype(df[c]):
                continue
            cl = cs.lower()
            if "_vs_" in cl or any(a in cl for a in ANNOT):
                continue
            out.append(c)
        return out

    comp_samples = _samples(comp_df)

    # 4) PCA matrix: evaluate EVERY sheet and pick the best per-sample matrix.
    #    Priority: columns matching the design > raw (non z-scored) > log-scaled >
    #    more proteins. A z-scored heatmap is usually a SELECTED subset (significant
    #    or most-variable proteins): running the PCA on it is circular - it inflates
    #    group separation - and its correlations are decorrelated by construction.
    design_labels = set()
    if design_map:
        design_labels = {k for k in design_map.keys() if k != "__order__"}

    def _profile(df, cols):
        X = df[cols].apply(pd.to_numeric, errors="coerce")
        rm, rs = X.mean(axis=1), X.std(axis=1)
        zscored = bool((rm.abs().median() < 0.05) and (abs(rs.median() - 1.0) < 0.15))
        med = float(np.nanmedian(X.values)) if X.size else 0.0
        log_scale = bool(3.0 < med < 45.0)          # log2 intensities land here; raw LFQ does not
        if design_labels:
            hit = sum(1 for c in cols if str(c) in design_labels or _norm(str(c)) in design_labels)
            match = hit / max(len(cols), 1)
        else:
            match = 0.0
        return {"zscored": zscored, "log_scale": log_scale, "match": match,
                "n_rows": len(df), "n_cols": len(cols)}

    candidates = []
    for name, df in xls.items():
        cols = _samples(df)
        if len(cols) < 3:
            continue
        prof = _profile(df, cols)
        prof.update({"sheet": name, "df": df, "cols": cols})
        candidates.append(prof)

    if candidates:
        best = max(candidates, key=lambda c: (round(c["match"], 2),
                                              0 if c["zscored"] else 1,
                                              1 if c["log_scale"] else 0,
                                              c["n_rows"]))
        pca_df, pca_source, sample_int_cols = best["df"], best["sheet"], best["cols"]
        pca_is_zscored = best["zscored"]
    else:
        pca_df, pca_source, sample_int_cols = comp_df, comp_sheet, comp_samples
        pca_is_zscored = False

    # 5) conditions: from the design when provided, else parsed from column names
    def _cond_of(col):
        return str(col).rsplit("_", 1)[0]

    if design_map:
        labels = design_map.get("__labels__") or []
        order = design_map.get("__order__") or []
        cond_map, how = {}, "exact"
        matched = match_columns_to_labels(sample_int_cols, labels) if labels else {}
        hits = len(matched)
        if hits == len(sample_int_cols) and labels:
            cond_map = {c: order[i] for c, i in matched.items()}
            if any(_norm(str(c)) != _norm(str(labels[i])) for c, i in matched.items()):
                how = "by containment"
        else:
            # fall back on the old key lookup, then on row order
            for c in sample_int_cols:
                v = design_map.get(str(c), design_map.get(_norm(str(c))))
                cond_map[c] = v if v is not None else _cond_of(c)
            hits = sum(1 for c in sample_int_cols
                       if design_map.get(str(c), design_map.get(_norm(str(c)))) is not None)
            if hits < 0.5 * len(sample_int_cols) and order and len(order) == len(sample_int_cols):
                cond_map = {c: order[i] for i, c in enumerate(sample_int_cols)}
                hits, how = len(sample_int_cols), "positional"
                print("[Pathfinder] \u26a0 design labels not recognised among the columns: "
                      "POSITIONAL mapping applied (same count). Check the row order of the design.")
        if how == "by containment":
            print("[Pathfinder] design labels matched to columns by containment "
                  "(e.g. 'Sple-15' inside 'LFQ.intensity.Sple-15').")
        conditions_resolved = hits >= 0.5 * len(sample_int_cols) and \
            len(set(cond_map.values())) < len(sample_int_cols)
    else:
        conditions = list(dict.fromkeys(_cond_of(c) for c in sample_int_cols))
        cond_map = map_sample_to_condition(sample_int_cols, conditions)
        # resolved only if parsing actually grouped (fewer conditions than samples)
        conditions_resolved = len(set(cond_map.values())) < len(sample_int_cols)

    seen = []
    for c in sample_int_cols:
        cd = cond_map[c]
        if cd not in seen:
            seen.append(cd)
    conditions = seen

    id_col = _find_col(comp_df.columns, "Protein.Group", "name", "ID") or comp_df.columns[0]
    gene_col = _find_col(comp_df.columns, "Genes", "Gene_Name", "Gene")
    desc_col = _find_col(comp_df.columns, "First.Protein.Description", "Protein.Names", "Description")

    out = StatsData(
        comp=comp_df, zscore=z_df, pca_df=pca_df, pca_source=pca_source,
        pca_is_zscored=pca_is_zscored,
        contrasts=contrasts, sample_int_cols=sample_int_cols, comp_sample_cols=comp_samples,
        conditions=conditions, cond_map=cond_map, conditions_resolved=conditions_resolved,
        id_col=id_col, gene_col=gene_col, desc_col=desc_col,
    )
    out.candidates = [{k: v for k, v in c.items() if k != "df"} for c in candidates]
    return out


# --------------------------------------------------------------------------- #
# GO                                                                            #
# --------------------------------------------------------------------------- #
def _sheet_variants(a, b):
    """Naming variants an upstream module may use for '<a> vs <b>'."""
    return {_norm(a + "v" + b), _norm(a + "vs" + b), _norm(a + "_vs_" + b),
            _norm(a + "_v_" + b), _norm(a + "-vs-" + b), _norm(a + b)}


def _condition_aliases(conditions):
    """Map every unambiguous spelling of a condition to that condition.

    Upstream sheet names sometimes drop the token shared by every condition
    ('EPHPS11' -> 'EPHP', 'CS11' -> 'C' when all conditions end in 'S11'), so exact
    matching on the full names alone would reject every sheet. Aliases that would
    designate two different conditions are discarded rather than guessed.
    """
    conds = [str(c) for c in conditions]
    norms = [_norm(c) for c in conds]
    aliases = {}

    def add(alias, cond):
        a = _norm(alias)
        if a:
            aliases.setdefault(a, set()).add(cond)

    for c in conds:
        add(c, c)
    if len(norms) > 1:
        pre = os.path.commonprefix(norms)
        suf = os.path.commonprefix([n[::-1] for n in norms])[::-1]
        for c, n in zip(conds, norms):
            if suf and len(n) > len(suf):
                add(n[: len(n) - len(suf)], c)
            if pre and len(n) > len(pre):
                add(n[len(pre):], c)
            if pre and suf and len(n) > len(pre) + len(suf):
                add(n[len(pre): len(n) - len(suf)], c)
    return {a: next(iter(s)) for a, s in aliases.items() if len(s) == 1}


def match_columns_to_labels(cols, labels):
    """Match sample columns to design labels. Returns {column: label_index}.

    Three levels, from safest to loosest:
      1. exact (normalised) equality;
      2. containment, when unambiguous - a design listing 'Sple-15' still matches a
         column named 'LFQ.intensity.Sple-15', and vice versa;
      3. nothing - the caller decides whether to fall back on row order.
    """
    nlabels = [_norm(str(l)) for l in labels]
    exact = {}
    for i, nl in enumerate(nlabels):
        exact.setdefault(nl, i)
    out = {}
    for c in cols:
        nc = _norm(str(c))
        if nc in exact:
            out[c] = exact[nc]
            continue
        hits = [i for i, nl in enumerate(nlabels) if nl and (nl in nc or nc in nl)]
        # Tie-break on the closest length, which works whichever side contains the
        # other: 'Sple-1' must not grab 'LFQ.intensity.Sple-10', and 'S1' must not
        # grab a column named after 'S11'.
        if hits:
            hits.sort(key=lambda i: abs(len(nlabels[i]) - len(nc)))
            best_gap = abs(len(nlabels[hits[0]]) - len(nc))
            best = [i for i in hits if abs(len(nlabels[i]) - len(nc)) == best_gap]
            if len(best) == 1:
                out[c] = best[0]
    return out


def _map_go_sheet(sheet, contrasts, conditions=None):
    """Map a 'GO_<abbrev>' sheet to its contrast. Returns (contrast, how).

    EXACT matching only. Substring scoring is unsafe when condition names nest:
    'HPS7' is a substring of 'EPHPS7' and 'IPHPS7', so a fuzzy score would happily
    assign GO_EPHPS7vIPHPS7 to IPHPS7_vs_HPS7 and silently merge two contrasts.
    """
    s = re.sub(r"^GO[_-]?", "", str(sheet), flags=re.I)
    ns = _norm(s)

    # 1) exact match against the contrast name and its usual abbreviations
    for c in contrasts:
        if "_vs_" in str(c):
            a, b = str(c).split("_vs_", 1)
            if ns == _norm(c) or ns in _sheet_variants(a, b):
                return c, "exact"

    # 2) split the sheet name into two KNOWN conditions (exact match on both sides,
    #    tolerating the shared prefix/suffix that upstream naming may drop)
    if conditions:
        conds = _condition_aliases(conditions)
        # contrasts come from the stats columns, conditions often from the design:
        # the two spellings may differ (underscores, spaces, case), so compare the
        # reconstructed contrast in normalised form rather than character for character.
        by_norm = {_norm(c): c for c in contrasts}
        for i in range(1, len(ns)):
            left = ns[:i]
            if left not in conds:
                continue
            for sep in ("vs", "v", ""):
                j = i + len(sep)
                if ns[i:j] != sep:
                    continue
                right = ns[j:]
                if right in conds and right != left:
                    a, b = conds[left], conds[right]
                    hit = by_norm.get(_norm(f"{a}_vs_{b}"))
                    if hit is not None:
                        return hit, "parsed"
                    hit = by_norm.get(_norm(f"{b}_vs_{a}"))
                    if hit is not None:
                        return hit, "parsed (reversed)"
    return None, "UNMATCHED"


def load_go(path: str, contrasts, conditions=None, verbose=True) -> GOData:
    xls = pd.read_excel(path, sheet_name=None)
    go_sheets = {n: d for n, d in xls.items() if str(n).lower().startswith("go")
                 and _find_col(d.columns, "term_name")}
    if not go_sheets:
        raise ValueError("Aucune sheet 'GO_*' avec colonne 'term_name' : fichier GO non reconnu.")

    # columns detected on the first sheet, assumed homogeneous
    sample_df = next(iter(go_sheets.values()))
    c_source = _find_col(sample_df.columns, "source")
    c_termid = _find_col(sample_df.columns, "term_id", "native")
    c_termnm = _find_col(sample_df.columns, "term_name")
    c_pval = _find_col(sample_df.columns, "p_value", "p.value", "adjusted_p_value", "p_val", "p.adj")
    c_inter = _find_col(sample_df.columns, "intersection", "intersections")
    c_z = _find_col(sample_df.columns, "z_score", "zscore", "z.score", "z-score", contains_all=["z", "score"])

    # map every GO sheet to a contrast, exactly, and report any problem loudly
    mapping, how = {}, {}
    for sheet in go_sheets:
        if contrasts:
            c, h = _map_go_sheet(sheet, contrasts, conditions)
        else:
            c, h = sheet, "no contrast list"
        mapping[sheet] = c
        how[sheet] = h

    collisions = {}
    for sheet, c in mapping.items():
        if c is not None:
            collisions.setdefault(c, []).append(sheet)
    dupes = {c: sh for c, sh in collisions.items() if len(sh) > 1}
    unmatched = [sh for sh, c in mapping.items() if c is None]
    missing = [c for c in (contrasts or []) if c not in collisions]

    if verbose:
        if dupes:
            for c, sh in dupes.items():
                print(f"[Pathfinder] \u26a0 GO mapping collision: {', '.join(sh)} all map to '{c}' "
                      "- their terms would be merged. Check the sheet names.")
        if unmatched:
            print(f"[Pathfinder] \u26a0 {len(unmatched)} GO sheet(s) could not be matched to a "
                  f"contrast and are IGNORED: {', '.join(unmatched)}")
        if missing:
            print(f"[Pathfinder] note: {len(missing)} contrast(s) have no GO sheet: "
                  f"{', '.join(map(str, missing[:6]))}{' ...' if len(missing) > 6 else ''}")

    rows = []
    for sheet, df in go_sheets.items():
        contrast = mapping[sheet]
        if contrast is None:
            continue
        for _, r in df.iterrows():
            inter = r.get(c_inter, "") if c_inter else ""
            acc = [a.strip() for a in re.split(r"[,;/ ]+", str(inter)) if a.strip()] if pd.notna(inter) else []
            rows.append({
                "contrast": contrast,
                "sheet": sheet,
                "source": r.get(c_source, "NA") if c_source else "NA",
                "term_id": r.get(c_termid, "") if c_termid else "",
                "term_name": r.get(c_termnm, ""),
                "p_value": pd.to_numeric(r.get(c_pval, np.nan), errors="coerce") if c_pval else np.nan,
                "z": pd.to_numeric(r.get(c_z, np.nan), errors="coerce") if c_z else np.nan,
                "intersection": acc,
            })
    long = pd.DataFrame(rows, columns=["contrast", "sheet", "source", "term_id", "term_name",
                                       "p_value", "z", "intersection"])
    go = GOData(long=long, p_col=str(c_pval), z_col=str(c_z),
                contrasts=sorted(long["contrast"].dropna().unique().tolist()) if len(long) else [])
    go.mapping = mapping
    go.mapping_how = how
    go.unmatched = unmatched
    go.duplicates = dupes
    go.missing = missing
    return go
