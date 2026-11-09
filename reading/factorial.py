r"""Factorial designs (Hordeum-like: genotype x light x time): how much of each protein's
variance each factor and each two-factor interaction carries.

Factors, in order of preference:
  1. --factor-regex: named groups read in the condition names, e.g.
     '(?P<genotype>[^_]+)_T(?P<time>\d+)(?:_(?P<light>Light|Dark))?' (absent group = '—');
  2. extra columns of Deimos' matrices/samples.tsv (e.g. genotype, light, time);
  3. the design file given to Pathfinder (--design), extra columns;
  4. the condition names, split on '_' / '-' / '.' / ' ' (Tron_T3_Light -> Tron | T3 |
     Light), then, if needed, letters glued to a number (Grp_Male4 -> Grp | Male | 4).
     Time-like tokens (T3, 24h, D7...) form the 'time' factor wherever they sit;
     conditions with fewer tokens (Tron_T0 next to Tron_T3_Light) are aligned on the
     vocabulary of each position, the missing factor reading '—'.
A factor must be constant within a condition, have 2 levels or more, and not simply
re-label the conditions. At least two factors make the design factorial.

Roles for the displays: x = the ordered factor (time, else the factor with most levels),
lines = the factor with most levels among the others (genotype), facets = the rest
(light). A level of x present in a single facet (T0 in the light only) is the shared
starting point of every facet.

Variance partition, per protein, on the per-sample log2 abundances: linear model with the
main effects and the two-factor interactions (sum-to-zero coding), type-II sums of
squares (each main effect adjusted for the other main effects and the interactions that
do not contain it; each interaction adjusted for everything else), then
  - "higher-order": what the condition means carry beyond that model (3-way and above,
    or combinations the model cannot separate);
  - "replicates": the within-condition residual.
Shares are divided by the total sum of squares; in an unbalanced design (missing
combinations, e.g. T0 in the light only) correlated terms make them sum to slightly
more or less than 1, which the tables say.
"""
from __future__ import annotations

import itertools
import re

import numpy as np
import pandas as pd

_NOT_FACTORS = {"label", "condition", "replicate", "rep", "sample", "file", "files", "run",
                "raw", "path", "order", "injection", "injection_order", "run_order",
                "date", "batch_order", "name", "id", "sample_id"}
_TIME = re.compile(r"^(?:[TDJ]|day|time)?\d+(?:\.\d+)?\s*(?:h|hr|hrs|min|d|dpi|hai|hpi|dap|w|wk)?$",
                   re.I)


def _factor_ok(series_by_cond, n_cond):
    lv = series_by_cond.astype(str).nunique()
    return 2 <= lv < n_cond


def _from_columns(df, conds):
    """Columns of a sample table that are constant within each condition."""
    out = {}
    if df is None or "condition" not in df.columns:
        return out
    d = df[df["condition"].astype(str).isin(conds)].copy()
    d["condition"] = d["condition"].astype(str)
    for col in d.columns:
        if str(col).strip().lower() in _NOT_FACTORS or col == "condition":
            continue
        per = d.groupby("condition")[col].nunique(dropna=False)
        if (per != 1).any():
            continue
        s = d.drop_duplicates("condition").set_index("condition")[col].reindex(conds)
        if s.isna().any():
            continue
        if pd.api.types.is_numeric_dtype(s) and s.nunique() > 12:
            continue                      # continuous covariate, not a design factor
        if _factor_ok(s, len(conds)):
            out[str(col)] = s.astype(str)
    return out


MISSING = "—"           # level of a factor absent from a condition name (Tron_T0: no light)


def _name_of(lv, j):
    lv = [v for v in lv if v != MISSING]
    # 'time' only with a unit or a prefix somewhere (T3, 24h); bare numbers (4 / 14:
    # an age? a dose?) keep their levels as name
    if lv and all(_TIME.match(v) for v in lv) and not all(v.replace(".", "").isdigit() for v in lv):
        return "time"
    if len(lv) <= 4:
        return "/".join(lv)
    return f"factor {j + 1}"


_GLUED = re.compile(r"^([A-Za-z]+?)(\d+(?:\.\d+)?)$")


def _from_names(conds):
    """Factors read in the condition names: separators first; if that gives fewer than two
    factors, letters glued to a number are split too (Grp_Male4 -> Grp | Male | 4,
    embAviso24 -> embAviso | 24), time-like tokens (T3, D7, 24h) staying whole."""
    toks = [[t for t in re.split(r"[_\-. ]+", c) if t] for c in conds]
    out = _parse_tokens(conds, toks)
    if len(out) >= 2:
        return out
    split = []
    for t in toks:
        row = []
        for x in t:
            m = _GLUED.match(x)
            if m and not (_TIME.match(x) and not x.isdigit() and len(m.group(1)) <= 1):
                row += [m.group(1), m.group(2)]
            else:
                row.append(x)
        split.append(row)
    if split != toks:
        out2 = _parse_tokens(conds, split)
        if len(out2) >= 2:
            return out2
    return out


def _parse_tokens(conds, toks):
    if min(len(t) for t in toks) < 1 or max(len(t) for t in toks) < 2:
        return {}
    # time-like tokens -> one 'time' factor, wherever they sit in the name
    # (bare numbers are not taken for time: 4 / 14 may be an age, a dose, a line)
    tt = [[x for x in t if _TIME.match(x) and not x.replace(".", "").isdigit()] for t in toks]
    use_time = all(len(x) == 1 for x in tt) and len({x[0] for x in tt}) >= 2
    rest = [[x for x in t if not (use_time and x == tt[i][0])] for i, t in enumerate(toks)]
    L = max(len(r) for r in rest)
    slots = [[None] * len(conds) for _ in range(L)]
    full = [r for r in rest if len(r) == L]
    vocab = [{r[j] for r in full} for j in range(L)]
    for i, r in enumerate(rest):
        if len(r) == L:
            for j, x in enumerate(r):
                slots[j][i] = x
            continue
        free = list(range(L))
        for x in r:                   # shorter name: each token goes to the slot that knows it
            cand = [j for j in free if x in vocab[j]]
            if len(cand) != 1:
                return {}             # unknown or ambiguous token: no reliable reading
            slots[cand[0]][i] = x
            free.remove(cand[0])
    out = {}
    if use_time:
        out["time"] = pd.Series([x[0] for x in tt], index=conds)
    for j in range(L):
        s = pd.Series([v if v is not None else MISSING for v in slots[j]], index=conds)
        if not _factor_ok(s, len(conds)):
            continue
        name = _name_of(list(dict.fromkeys(s)), j)
        while name in out:
            name += "'"
        out[name] = s
    return out


def _from_regex(conds, rx):
    """Named groups of a regex read in every condition name ('—' for an absent group)."""
    pat = re.compile(rx)
    if not pat.groupindex:
        raise ValueError("--factor-regex needs named groups, e.g. (?P<genotype>[^_]+)")
    rows = {}
    for c in conds:
        m = pat.search(c)
        if m is None:
            raise ValueError(f"--factor-regex does not match the condition '{c}'")
        rows[c] = {k: (v if v is not None else MISSING) for k, v in m.groupdict().items()}
    F = pd.DataFrame.from_dict(rows, orient="index")
    return {k: F[k] for k in F.columns if F[k].nunique() >= 2}


def design_factors(run, conds, factors_df=None, regex=None):
    """(DataFrame condition x factor, source) with >= 2 factors, else (None, reason)."""
    conds = [str(c) for c in conds]
    src = ""
    F = {}
    if regex:
        F = _from_regex(conds, regex)
        src = "--factor-regex" if len(F) >= 2 else ""
    if len(F) < 2:
        F = _from_columns(run.samples, conds)
        src = "samples.tsv columns" if len(F) >= 2 else ""
    if len(F) < 2 and factors_df is not None and run.samples is not None:
        fd = factors_df.copy()
        fd = fd.drop(columns=[c for c in fd.columns if str(c).strip().lower() == "condition"])
        fd.index = fd.index.astype(str)
        smp = run.samples.drop_duplicates("label").set_index("label")
        fd = fd.join(smp[["condition"]], how="inner")
        F = _from_columns(fd.reset_index(), conds)
        src = "design file columns" if len(F) >= 2 else ""
    if len(F) < 2:
        F = _from_names(conds)
        src = "condition names" if len(F) >= 2 else ""
    if len(F) < 2:
        return None, "fewer than two design factors"
    F = pd.DataFrame(F).reindex(conds)
    # two factors that are the same partition (aliases) carry no separate information
    keep = []
    for col in F.columns:
        if not any(F.groupby(col)[k].nunique().max() == 1 and F.groupby(k)[col].nunique().max() == 1
                   for k in keep):
            keep.append(col)
    F = F[keep]
    if F.shape[1] < 2:
        return None, "the factors found are aliases of each other"
    return F, src


def _coding(levels_per_sample):
    """Sum-to-zero (deviation) coding of one factor -> n x (L-1) matrix."""
    s = pd.Series(levels_per_sample).astype(str)
    lv = list(dict.fromkeys(s))
    X = np.zeros((len(s), len(lv) - 1))
    for j, l in enumerate(lv[1:]):
        X[:, j] = (s == l).astype(float) - (s == lv[0]).astype(float)
    return X


def _rss(X, Y):
    """Residual sum of squares of every column of Y after projection on X."""
    if X.shape[1] == 0:
        R = Y
    else:
        B, *_ = np.linalg.lstsq(X, Y, rcond=None)
        R = Y - X @ B
    return (R ** 2).sum(axis=0)


def variance_partition(run, ids, F):
    """DataFrame protein x [factor terms..., 'higher-order', 'replicates', 'sum'] of
    variance shares (0-1), plus the list of model terms."""
    if run.abundance is None or run.samples is None:
        return None, []
    smp = run.samples.drop_duplicates("label").set_index("label")["condition"].astype(str)
    labs = [l for l in run.abundance.columns if smp.get(l) in F.index]
    if len(labs) < len(F) + 2:
        return None, []
    A = run.abundance.reindex(ids)[labs].astype(float)
    A = A[A.notna().all(axis=1)]
    if A.empty:
        return None, []
    Y = (A.to_numpy() - A.to_numpy().mean(axis=1, keepdims=True)).T      # n x p
    cond = smp.loc[labs].to_numpy()
    n = len(labs)
    blocks = {f: _coding(F.loc[cond, f].to_numpy()) for f in F.columns}
    terms = list(F.columns)
    for a, b in itertools.combinations(F.columns, 2):
        P = np.column_stack([blocks[a][:, i] * blocks[b][:, j]
                             for i in range(blocks[a].shape[1]) for j in range(blocks[b].shape[1])])
        if np.linalg.matrix_rank(P) > 0:
            blocks[f"{a} × {b}"] = P
            terms.append(f"{a} × {b}")
    one = np.ones((n, 1))

    def X_of(ts):
        return np.column_stack([one] + [blocks[t] for t in ts]) if ts else one

    sst = (Y ** 2).sum(axis=0)
    rss_full = _rss(X_of(terms), Y)
    out = {}
    for t in terms:
        if " × " in t:
            others = [u for u in terms if u != t]
            out[t] = _rss(X_of(others), Y) - rss_full
        else:
            base = [u for u in terms if u != t and t not in u.split(" × ")]
            out[t] = _rss(X_of(base), Y) - _rss(X_of(base + [t]), Y)
    C = pd.get_dummies(pd.Series(cond)).to_numpy(float)
    rss_cond = _rss(C, Y)
    out["higher-order"] = np.clip(rss_full - rss_cond, 0, None)
    out["replicates"] = rss_cond
    with np.errstate(invalid="ignore", divide="ignore"):
        S = pd.DataFrame({k: np.clip(v, 0, None) / sst for k, v in out.items()}, index=A.index)
    S["sum"] = S.sum(axis=1)
    return S, terms


def _num(v):
    m = re.search(r"\d+(?:\.\d+)?", str(v))
    return float(m.group()) if m else float("inf")


def is_ordered(s) -> bool:
    lv = [v for v in dict.fromkeys(s.astype(str)) if v != MISSING]
    if len(lv) < 3:
        return False
    return (all(_TIME.match(v) for v in lv)
            or all(re.fullmatch(r"-?\d+(?:\.\d+)?", v) for v in lv))


def levels(F, f):
    """Levels of a factor in display order: numeric order for an ordered factor, else
    order of first appearance; '—' last."""
    lv = list(dict.fromkeys(F[f].astype(str)))
    if is_ordered(F[f]) or all(re.fullmatch(r"-?\d+(?:\.\d+)?", v) for v in lv if v != MISSING):
        lv = sorted(lv, key=lambda v: (v == MISSING, _num(v)))
    return [v for v in lv if v != MISSING] + ([MISSING] if MISSING in lv else [])


def roles(F):
    """{'x', 'lines', 'facets'}: which factor goes on the x axis, which one draws the
    lines, which ones split the panels."""
    nl = {f: len([v for v in set(F[f]) if v != MISSING]) for f in F.columns}
    ordered = [f for f in F.columns if is_ordered(F[f])]
    if not ordered:                     # numeric levels (age 4 / 14, dose) go on x
        ordered = [f for f in F.columns
                   if all(re.fullmatch(r"-?\d+(?:\.\d+)?", str(v))
                          for v in set(F[f]) if v != MISSING)]
    x = ordered[0] if ordered else max(F.columns, key=lambda f: nl[f])
    rest = [f for f in F.columns if f != x]
    lines = max(rest, key=lambda f: nl[f]) if rest else None
    facets = [f for f in rest if f != lines]
    return {"x": x, "lines": lines, "facets": facets}


def condition_order(F, r):
    """Conditions ordered lines > facets > x (genotype blocks, then light, then time)."""
    keys = [r["lines"]] + r["facets"] + [r["x"]]
    # '—' (factor absent: the shared baseline, e.g. T0 without light) comes first
    rank = {f: {v: (-1 if v == MISSING else i) for i, v in enumerate(levels(F, f))}
            for f in keys if f}
    return sorted(F.index, key=lambda c: tuple(rank[f].get(F.at[c, f], 99) for f in keys if f))


def contrast_factors(pairs, F, contrasts):
    """contrast -> (factors that differ, fixed context {factor: level}, numerator,
    denominator). Contrasts whose conditions are not both in F are left out."""
    out = {}
    # contrasts name the conditions as in the design (Grp_Male4); samples.tsv / PCA_Scores
    # may write them with '-' (Grp-Male4): match on the name without separators
    canon = {re.sub(r"[_\-. ]", "", str(c)).lower(): c for c in F.index}

    def _find(x):
        return x if x in F.index else canon.get(re.sub(r"[_\-. ]", "", str(x)).lower())
    for c in contrasts:
        num, den = pairs.get(c, (None, None))
        num, den = _find(num), _find(den)
        if num is None or den is None:
            continue
        diff = [f for f in F.columns if F.at[num, f] != F.at[den, f]]
        ctx = {f: F.at[num, f] for f in F.columns if f not in diff}
        out[c] = (diff, ctx, num, den)
    return out


def factor_blocks(cf, F):
    """Group the one-factor contrasts by the factor they vary: [(factor, [contrast...])],
    largest blocks first; contrasts varying 2+ factors under the key 'several'."""
    blocks = {}
    for c, (diff, ctx, num, den) in cf.items():
        key = diff[0] if len(diff) == 1 else ("several" if diff else None)
        if key:
            blocks.setdefault(key, []).append(c)
    return sorted(blocks.items(), key=lambda kv: (kv[0] == "several", -len(kv[1])))


def summarise(S, prog=None):
    """Median share per term, for all proteins and per programme."""
    cols = [c for c in S.columns if c != "sum"]
    rows = [("all selected", len(S), S[cols].median())]
    if prog is not None:
        lab = prog.labels.reindex(S.index)
        for p in range(1, prog.k + 1):
            sub = S[lab == p]
            if len(sub):
                rows.append((f"P{p}", len(sub), sub[cols].median()))
    out = pd.DataFrame([dict(set=r[0], n_proteins=r[1], **r[2].round(3).to_dict()) for r in rows])
    return out
