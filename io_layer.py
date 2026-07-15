"""
Pathfinder — couche d'entrée commune.

Lit les Excel produits par Deimos (module DEP) et le module GO, et reconstruit
les trois objets qui alimentent tous les axes. Aucune statistique n'est
recalculée : on lit, on extrait, on reshape (règle fondamentale du pipeline).

Reconnaissance par CONTENU, pas par nom de fichier :
  - Stats  : une sheet contenant des colonnes '<C>_p.val' + une sheet 'Z-score_Heatmap'
  - GO     : des sheets 'GO_*' avec une colonne 'term_name'
"""
from __future__ import annotations
import re
import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

warnings.simplefilter("ignore", category=UserWarning)

# Colonnes méta connues du module DEP (jamais des intensités échantillon)
_META = {
    "name", "Protein.Group", "Protein.Names", "Genes", "First.Protein.Description",
    "N.Sequences", "N.Proteotypic.Sequences", "ID", "imputed", "num_NAs", "significant",
}


# --------------------------------------------------------------------------- #
# Structures                                                                   #
# --------------------------------------------------------------------------- #
@dataclass
class StatsData:
    comp: pd.DataFrame                    # sheet de comparaison (stats + accessions)
    zscore: Optional[pd.DataFrame]        # sheet Z-score_Heatmap (peut être None)
    pca_df: pd.DataFrame                  # sheet retenue pour la PCA (matrice par échantillon)
    pca_source: str                       # nom de la sheet PCA
    contrasts: list                       # ['A_vs_B', ...]
    sample_int_cols: list                 # colonnes intensité par échantillon (dans pca_df)
    comp_sample_cols: list                 # colonnes intensité brutes dans comp (pour abondance)
    conditions: list                      # ordre d'apparition
    cond_map: dict                        # {colonne_echantillon: condition}
    conditions_resolved: bool             # True si le mapping condition est fiable
    id_col: str                           # colonne identifiant protéine
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
    """Recherche tolérante d'une colonne : égalité normalisée puis sous-chaînes."""
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
    """Associe chaque colonne échantillon à sa condition.

    Priorité au préfixe le plus long parmi les conditions connues (robuste aux
    noms contenant des underscores), fallback sur rsplit('_', 1).
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
# Stats                                                                         #
# --------------------------------------------------------------------------- #
def load_design_frame(path_or_df):
    """Lit un design et renvoie (df, label_col) — toutes les colonnes conservées
    pour un usage multi-facteurs (condition, traitement, batch, ordre d'injection…)."""
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
    """Aligne les colonnes échantillon avec les lignes du design (par nom, sinon
    positionnellement) et renvoie un DataFrame indexé par colonne échantillon,
    une colonne par facteur."""
    labels = [str(x).strip() for x in df[label_col].tolist()]
    byname = {l: i for i, l in enumerate(labels)}
    byname.update({_norm(l): i for i, l in enumerate(labels)})
    rows, hits = [], 0
    for c in sample_cols:
        i = byname.get(str(c), byname.get(_norm(str(c))))
        if i is not None:
            hits += 1
        rows.append(i)
    if hits < 0.5 * len(sample_cols) and len(labels) == len(sample_cols):
        rows = list(range(len(sample_cols)))          # repli positionnel
    if any(i is None for i in rows):
        return None
    out = df.iloc[rows].copy()
    out.index = sample_cols
    keep = [c for c in out.columns if c != label_col]
    return out[keep]


def load_design(path_or_df):
    """Lit un ExperimentalDesign (CSV/xlsx/DataFrame) → {label_échantillon: condition}.

    Détection tolérante des colonnes label et condition. Clé indexée en brut ET
    normalisé pour tolérer les variantes de séparateur (Sple-1 / Sple_1).
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
        # heuristique : 1re colonne = label, colonne la moins unique = condition
        label_col = label_col or df.columns[0]
        if cond_col is None:
            nunq = {c: df[c].nunique() for c in df.columns if c != label_col}
            cond_col = min(nunq, key=nunq.get) if nunq else df.columns[-1]

    mp = {}
    for _, r in df.iterrows():
        lab, cond = str(r[label_col]).strip(), str(r[cond_col]).strip()
        mp[lab] = cond
        mp[_norm(lab)] = cond
    mp["__order__"] = [str(x).strip() for x in df[cond_col].tolist()]   # repli positionnel
    return mp


def find_design_in_workbook(path, sample_cols, contrasts):
    """Cherche un mapping échantillon→condition dans N'IMPORTE quelle sheet, par CONTENU.

    Critère : une colonne dont les valeurs recouvrent les noms d'échantillons, plus
    une autre colonne peu cardinale (candidate condition). Bonus si ses valeurs
    matchent les tokens de conditions extraits des contrastes.
    """
    want = {_norm(c) for c in sample_cols}
    cond_tokens = set()
    for c in contrasts or []:
        for t in str(c).split("_vs_"):
            cond_tokens.add(_norm(t))

    best = None
    try:
        sheets = pd.read_excel(path, sheet_name=None)
    except Exception:
        return None, None
    for name, df in sheets.items():
        if df.empty or df.shape[1] < 2:
            continue
        for lab_c in df.columns:
            vals = {_norm(v) for v in df[lab_c].dropna().astype(str)}
            cover = len(want & vals) / max(len(want), 1)
            if cover < 0.8:                       # doit couvrir la plupart des échantillons
                continue
            for cond_c in df.columns:
                if cond_c == lab_c:
                    continue
                nun = df[cond_c].nunique(dropna=True)
                if not (1 < nun < len(df)):       # ni constant ni unique par ligne
                    continue
                cvals = {_norm(v) for v in df[cond_c].dropna().astype(str)}
                match = len(cvals & cond_tokens) / max(len(cvals), 1) if cond_tokens else 0
                score = cover + match
                if best is None or score > best[0]:
                    best = (score, name, lab_c, cond_c, df)
    if not best:
        return None, None
    _, name, lab_c, cond_c, df = best
    mp = {}
    for _, r in df.iterrows():
        lab, cond = str(r[lab_c]).strip(), str(r[cond_c]).strip()
        mp[lab] = cond
        mp[_norm(lab)] = cond
    return mp, f"{name}#{lab_c}→{cond_c}"


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

    # 3) identification POSITIVE des colonnes échantillon (par sheet)
    #    Exclut : méta, colonnes de contraste (contiennent '_vs_'), annotations.
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
            if not pd.api.types.is_numeric_dtype(df[c]):
                continue
            cl = cs.lower()
            if "_vs_" in cl or any(a in cl for a in ANNOT):
                continue
            out.append(c)
        return out

    comp_samples = _samples(comp_df)
    z_samples = _samples(z_df) if z_df is not None else []

    # 4) matrice PCA = sheet avec le plus de colonnes échantillon (préférence Z-score
    #    si elle en contient autant/plus : c'est la matrice par échantillon la plus propre)
    if z_df is not None and len(z_samples) >= max(len(comp_samples), 1):
        pca_df, pca_source, sample_int_cols = z_df, "Z-score_Heatmap", z_samples
    else:
        pca_df, pca_source, sample_int_cols = comp_df, comp_sheet, comp_samples

    # 5) conditions : via design si fourni, sinon déduites des noms de colonnes
    def _cond_of(col):
        return str(col).rsplit("_", 1)[0]

    if design_map:
        cond_map, hits = {}, 0
        for c in sample_int_cols:
            v = design_map.get(str(c), design_map.get(_norm(str(c))))
            if v is not None:
                hits += 1
            cond_map[c] = v if v is not None else _cond_of(c)
        # repli positionnel : libellés du design != noms de colonnes (ex. fichiers .raw
        # côté design vs Sple-N côté matrice) mais même effectif → mapping par ordre.
        order = design_map.get("__order__")
        if hits < 0.5 * len(sample_int_cols) and order and len(order) == len(sample_int_cols):
            cond_map = {c: order[i] for i, c in enumerate(sample_int_cols)}
            hits = len(sample_int_cols)
            print("[Pathfinder] ⚠ libellés du design non reconnus dans les colonnes : "
                  "mapping POSITIONNEL appliqué (même effectif). Vérifie l'ordre des lignes du design.")
        conditions_resolved = hits >= 0.5 * len(sample_int_cols) and \
            len(set(cond_map.values())) < len(sample_int_cols)
    else:
        conditions = list(dict.fromkeys(_cond_of(c) for c in sample_int_cols))
        cond_map = map_sample_to_condition(sample_int_cols, conditions)
        # résolu seulement si le parsing a réellement groupé (moins de conditions que d'échantillons)
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

    return StatsData(
        comp=comp_df, zscore=z_df, pca_df=pca_df, pca_source=pca_source,
        contrasts=contrasts, sample_int_cols=sample_int_cols, comp_sample_cols=comp_samples,
        conditions=conditions, cond_map=cond_map, conditions_resolved=conditions_resolved,
        id_col=id_col, gene_col=gene_col, desc_col=desc_col,
    )


# --------------------------------------------------------------------------- #
# GO                                                                            #
# --------------------------------------------------------------------------- #
def _map_go_sheet(sheet: str, contrasts):
    """Associe une sheet 'GO_<abbrev>' au contraste '<A>_vs_<B>' correspondant."""
    s = re.sub(r"^GO[_-]?", "", str(sheet), flags=re.I)
    # tokens de la sheet : coupe sur '_vs_' ou un 'v' isolé entre deux mots
    s2 = re.sub(r"([A-Za-z0-9])v([A-Z])", r"\1_vs_\2", s)
    sheet_toks = {_norm(t) for t in re.split(r"_vs_|_", s2) if t}

    best, best_score = None, -1.0
    for c in contrasts:
        c_toks = {_norm(t) for t in c.split("_vs_")}
        if sheet_toks and c_toks:
            inter = sum(1 for st in sheet_toks for ct in c_toks if st and (st in ct or ct in st))
            score = inter / max(len(sheet_toks), len(c_toks))
        else:
            score = 0.0
        # départage par similarité de caractères
        score += 0.001 * sum(1 for a, b in zip(sorted(_norm(s)), sorted(_norm(c))) if a == b)
        if score > best_score:
            best, best_score = c, score
    return best


def load_go(path: str, contrasts) -> GOData:
    xls = pd.read_excel(path, sheet_name=None)
    go_sheets = {n: d for n, d in xls.items() if str(n).lower().startswith("go")
                 and _find_col(d.columns, "term_name")}
    if not go_sheets:
        raise ValueError("Aucune sheet 'GO_*' avec colonne 'term_name' : fichier GO non reconnu.")

    # colonnes détectées sur la première sheet, supposées homogènes
    sample_df = next(iter(go_sheets.values()))
    c_source = _find_col(sample_df.columns, "source")
    c_termid = _find_col(sample_df.columns, "term_id", "native")
    c_termnm = _find_col(sample_df.columns, "term_name")
    c_pval = _find_col(sample_df.columns, "p_value", "p.value", "adjusted_p_value", "p_val", "p.adj")
    c_inter = _find_col(sample_df.columns, "intersection", "intersections")
    c_z = _find_col(sample_df.columns, "z_score", "zscore", "z.score", "z-score", contains_all=["z", "score"])

    rows = []
    for sheet, df in go_sheets.items():
        contrast = _map_go_sheet(sheet, contrasts) if contrasts else sheet
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
    long = pd.DataFrame(rows)
    return GOData(long=long, p_col=str(c_pval), z_col=str(c_z),
                  contrasts=sorted(long["contrast"].dropna().unique().tolist()))
