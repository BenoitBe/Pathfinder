"""Response programmes: the proteins significant in at least one contrast are
grouped by the shape of their profile across conditions (per-condition means,
z-scored per protein). Independent of any notion of time: a programme is a group of
proteins that move together.

Each programme is then described by:
  - the conditions where it is high / low;
  - the functional groups over-represented among its proteins (hypergeometric
    test, universe = selected proteins carrying at least one term, BH q < 0.05 like
    the upstream GO);
  - its most marked proteins with a readable name.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

THEME_Q_THR = 0.05      # same threshold as the upstream GO (gProfiler / STRING, q < 5 %)


@dataclass
class Programmes:
    z: pd.DataFrame            # proteins x conditions (z-scores)
    labels: pd.Series          # protein -> programme (1..k, ordered)
    k: int
    themes: dict               # programme -> [(functional group name, q, n)]
    names: dict                # programme -> [readable names]
    high: dict                 # programme -> condition where it is highest
    low: dict
    stability: dict = None     # programme -> mean Jaccard over the replicate resamplings
    protein_stability: pd.Series = None   # protein -> fraction of resamplings kept
    n_boot: int = 0
    seed: int = 42


def _choose_k(Z, kmin=3, kmax=8, seed=42):
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score
    best, best_s = kmin, -1
    n = len(Z)
    for k in range(kmin, min(kmax, max(kmin, n // 15)) + 1):
        lab = KMeans(n_clusters=k, n_init=10, random_state=seed).fit_predict(Z)
        s = silhouette_score(Z, lab, sample_size=min(n, 2000), random_state=seed)
        if s > best_s + 0.01:          # a larger k must really separate better
            best, best_s = k, s
    return best


def build_programmes(run, design, themes=None, k=None, seed=42) -> Programmes | None:
    from sklearn.cluster import KMeans
    if run.means is None or len(run.conditions) < 3:
        return None
    sig = np.zeros(len(run.de), bool)
    for c in run.contrasts:
        sig |= run.sig(c).fillna(False).values
    ids = [p for p in run.de.index[sig] if p in run.means.index]
    if len(ids) < 20:
        return None
    conds = [c for c in design.conditions if c in run.means.columns]
    if len(conds) < 3:
        return None
    P = run.means.loc[ids, conds]
    sd = P.std(axis=1).replace(0, np.nan)
    Z = P.sub(P.mean(axis=1), axis=0).div(sd, axis=0).dropna()
    k = k or _choose_k(Z.values, seed=seed)
    km = KMeans(n_clusters=k, n_init=30, random_state=seed).fit(Z.values)
    raw = pd.Series(km.labels_, index=Z.index)
    # order: most populated programmes first
    order = raw.value_counts().index.tolist()
    lab = raw.map({c: i + 1 for i, c in enumerate(order)})

    lfc_cols = [c + "_diff" for c in run.contrasts]
    peak = run.de.loc[Z.index, lfc_cols].abs().max(axis=1)
    names, high, low, th = {}, {}, {}, {}
    from .labels import informative
    for p in range(1, k + 1):
        members = lab.index[lab == p]
        prof = Z.loc[members].median()
        high[p] = prof.idxmax(); low[p] = prof.idxmin()
        cand = peak.loc[members].sort_values(ascending=False).index
        nm = [run.labels.get(i, i) for i in cand if informative(run.labels.get(i, ""))]
        names[p] = list(dict.fromkeys(nm))[:6]
    if themes is not None:
        th = _theme_enrichment(lab, themes)
    return Programmes(z=Z, labels=lab, k=k, themes=th, names=names, high=high, low=low,
                      seed=seed)


# Jaccard thresholds of Hennig's clusterboot (2007): >= 0.75 stable, 0.6-0.75 a pattern
# that holds only partly, < 0.6 a programme that dissolves under resampling.
STABLE_J, PARTIAL_J = 0.75, 0.6


def stability_class(j) -> str:
    if j is None or not np.isfinite(j):
        return "n/a"
    return "stable" if j >= STABLE_J else ("partial" if j >= PARTIAL_J else "unstable")


def programme_stability(run, prog, n_boot=50, seed=None):
    """Stability of the programmes under replicate resampling (clusterboot-like).

    Each resampling draws the replicates of every condition with replacement, recomputes
    the condition means and their per-protein z-scores, and re-runs the k-means (same k).
    Its clusters are matched one-to-one to the original programmes (Hungarian assignment
    on the Jaccard index of their members). Per programme: mean Jaccard with its match
    (Hennig 2007: >= 0.75 stable, < 0.6 dissolved). Per protein: fraction of the
    resamplings where it stays in the match of its programme. Sets prog.stability /
    prog.protein_stability in place and returns prog; needs per-sample abundances."""
    from scipy.optimize import linear_sum_assignment
    from sklearn.cluster import KMeans
    if prog is None or not n_boot or run.abundance is None or run.samples is None:
        return prog
    seed = prog.seed if seed is None else seed
    conds = list(prog.z.columns)
    smp = run.samples.drop_duplicates("label").set_index("label")["condition"]
    cols = {c: [l for l in run.abundance.columns if smp.get(l) == c] for c in conds}
    if any(not v for v in cols.values()) or max(len(v) for v in cols.values()) < 2:
        return prog                      # no replicate to resample
    ids = list(prog.z.index)
    A = run.abundance.reindex(ids)
    lab0 = prog.labels.reindex(ids).to_numpy()
    k = prog.k
    members0 = [set(np.flatnonzero(lab0 == p)) for p in range(1, k + 1)]
    rng = np.random.default_rng(seed)
    jac = []
    kept = np.zeros(len(ids))
    seen = np.zeros(len(ids))
    for b in range(n_boot):
        M = np.column_stack([A[rng.choice(cols[c], size=len(cols[c]), replace=True)]
                             .mean(axis=1).to_numpy() for c in conds])
        sd = M.std(axis=1, ddof=1)
        ok = np.isfinite(sd) & (sd > 0)
        Z = (M[ok] - M[ok].mean(axis=1, keepdims=True)) / sd[ok][:, None]
        if len(Z) <= k:
            continue
        lab = np.full(len(ids), -1)
        lab[ok] = KMeans(n_clusters=k, n_init=10, random_state=seed + 1 + b).fit_predict(Z)
        members = [set(np.flatnonzero(lab == q)) for q in range(k)]
        J = np.array([[len(a & m) / len(a | m) if (a | m) else 0.0 for m in members]
                      for a in members0])
        r, c_ = linear_sum_assignment(-J)
        match = dict(zip(r, c_))
        jac.append([J[p, match[p]] for p in range(k)])
        for p in range(k):
            idx = np.fromiter(members0[p], int)
            kept[idx] += (lab[idx] == match[p])
        seen[ok] += 1
    if not jac:
        return prog
    jac = np.array(jac)
    prog.stability = {p + 1: float(jac[:, p].mean()) for p in range(k)}
    prog.protein_stability = pd.Series(np.where(seen > 0, kept / np.maximum(seen, 1), np.nan),
                                       index=ids)
    prog.n_boot = len(jac)
    return prog


def _theme_enrichment(lab, themes, max_themes=3):
    """Functional groups over-represented in each programme (hypergeometric + BH)."""
    from scipy.stats import hypergeom
    universe = set(lab.index)
    tg = {t: themes.genes(t) & universe for t in themes.names}
    tg = {t: g for t, g in tg.items() if len(g) >= 3}
    annotated = set().union(*tg.values()) if tg else set()
    if not annotated:
        return {}
    N = len(annotated)
    out = {}
    for p in sorted(lab.unique()):
        members = set(lab.index[lab == p]) & annotated
        n = len(members)
        rows = []
        for t, g in tg.items():
            k_ = len(g & members)
            if k_ >= 2:
                rows.append((t, hypergeom.sf(k_ - 1, N, len(g), n), k_))
        if rows:
            ps = np.array([r[1] for r in rows])
            order = np.argsort(ps)
            q = np.minimum.accumulate((ps[order] * len(ps) / (np.arange(len(ps)) + 1))[::-1])[::-1]
            sel = [(themes.names[rows[i][0]], float(qq), rows[i][2])
                   for i, qq in zip(order, q) if qq < THEME_Q_THR]
            out[p] = sel[:max_themes]
        else:
            out[p] = []
    return out
