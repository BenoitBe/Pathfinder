"""Design structure, read from the condition names or from the samples.tsv columns.

Three cases, from the richest to the simplest:
  - group x ordered level (time, stage, dose): trajectories, line profiles;
  - ordered level alone (one series);
  - categorical conditions (no order): everything stays readable, no trajectory.

Nothing is assumed. A number in a condition name is read as an order only with
evidence:
  - a unit glued to the number (24h, 7d, 15DAP, 48hai, 30min...);
  - an explicit prefix (T3, D7, J14, Day7, E14.5, Dose10...);
  - or a crossed design: at least two groups sharing at least two levels
    (embAviso0/2/8/24 and embMajor0/2/8/24).
Otherwise (Line1/Line2/Line3, Grp_Male4, Clone2...) the conditions stay categorical.
At least 3 distinct levels and 2 levels per group are also required.
Overrides (order_regex, order_column, categorical) always win.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import pandas as pd

_NUM = r"(?P<n>\d+(?:\.\d+)?)"
_UNIT = (r"(?P<u>hours?|hrs?|hpi|hai|hpf|hps|hpa|h|min|days?|dpi|dpa|dap|daa|das|dai|dpf|"
         r"d|weeks?|wks?|w|months?|mo)")
# explicit prefixes; no ambiguous single letter (S1, P2, E3 = sites, patients...)
_KW = r"(?:day|jour|dose|conc|stage|time|wk|week|t|d|j)"
_PREFIX = re.compile(rf"(?:^|[^A-Za-z])(?P<k>{_KW})[_\-. ]*$", re.I)

_PATTERNS = [
    # number (+ unit) at the end: embAviso24, WT_24h, Ctrl_D7, Stage3
    re.compile(rf"^(?P<a>.*?){_NUM}\s*{_UNIT}?$", re.I),
    # number + unit in the middle: WT_24h_Light
    re.compile(rf"^(?P<a>(?:.*?[^0-9.])?){_NUM}\s*{_UNIT}(?P<b>[^A-Za-z0-9].*)$", re.I),
    # T/D/J prefix + number in the middle: Tron_T0_Light
    re.compile(rf"^(?P<a>(?:.*?[^A-Za-z])?)(?P<k>T|D|J|Day){_NUM}(?P<b>[^0-9.].*)$", re.I),
    # number + unit at the start: 15DAP_WT, 24h_Ctrl
    re.compile(rf"^{_NUM}\s*{_UNIT}(?P<b>(?:[^A-Za-z].*)?)$", re.I),
]

_LEVEL_NAME = [  # (units / prefixes, axis name)
    ({"h", "hr", "hrs", "hour", "hours", "hpi", "hai", "hpf", "hps", "hpa"}, "time (h)"),
    ({"d", "day", "days", "dpi", "dpa", "dap", "daa", "das", "dai", "dpf", "j", "jour"},
     "time (d)"),
    ({"min"}, "time (min)"),
    ({"w", "wk", "wks", "week", "weeks"}, "time (wk)"),
    ({"mo", "month", "months"}, "time (months)"),
    ({"t", "time"}, "time"),
    ({"dose", "conc"}, "dose"),
    ({"stage"}, "stage"),
]
_UNIT_LABEL = {"time (h)": " h", "time (d)": " d", "time (min)": " min", "time (wk)": " wk",
               "time (months)": " mo"}


@dataclass
class Design:
    conditions: list            # display order
    group: dict | None = None   # condition -> group (None: no groups)
    level: dict | None = None   # condition -> ordered level (float) or None
    level_name: str = "level"
    source: str = "categorical"

    @property
    def ordered(self) -> bool:
        return self.level is not None

    @property
    def groups(self) -> list:
        return list(dict.fromkeys(self.group.values())) if self.group else []

    @property
    def levels(self) -> list:
        return sorted(set(self.level.values())) if self.level else []

    def display(self, c) -> str:
        """Short display name: common prefix / suffix removed
        (Tissue_Liver, Tissue_Brain -> Liver, Brain)."""
        return _display_map(tuple(self.conditions)).get(c, str(c))

    def label(self, lv) -> str:
        if self.level_name == "time":
            return f"T{lv:g}"
        return f"{lv:g}{_UNIT_LABEL.get(self.level_name, '')}"


_SEP = "_-. "


def _display_map(conds: tuple) -> dict:
    import os
    cs = [str(c) for c in conds]
    if len(cs) < 2:
        return {c: c for c in cs}
    pre = os.path.commonprefix(cs)
    cut = max((i + 1 for i, ch in enumerate(pre) if ch in _SEP), default=0)
    suf = os.path.commonprefix([c[::-1] for c in cs])
    cut_s = max((i + 1 for i, ch in enumerate(suf) if ch in _SEP), default=0)
    out = {}
    for c in cs:
        d = c[cut:len(c) - cut_s if cut_s else None].strip(_SEP)
        out[c] = d or c
    if len(set(out.values())) < len(cs):          # never ambiguous
        return {c: c for c in cs}
    return out


def _clean(*parts) -> str:
    parts = [str(p or "").strip("_-. ·") for p in parts]
    return " · ".join(p for p in parts if p)


def _parse(c, pat):
    """-> (group, level, unit/prefix hint or None) or None."""
    m = pat.match(c)
    if not m:
        return None
    d = m.groupdict()
    a, b, u, k = d.get("a") or "", d.get("b") or "", d.get("u"), d.get("k")
    hint = (u or k or "").lower() or None
    if not hint:
        mk = _PREFIX.search(a)
        if mk:
            hint = mk.group("k").lower()
            a = a[:mk.start("k")]
    return _clean(a, b) or "all", float(d["n"]), hint


def _try(conds, pat):
    g, lv, hints = {}, {}, []
    for c in conds:
        r = _parse(c, pat)
        if r is None:
            return None
        g[c], lv[c], h = r
        hints.append(h)
    levels = set(lv.values())
    per_group = pd.Series(lv).groupby(pd.Series(g)).agg(lambda s: set(s))
    if len(levels) < 3 or (per_group.map(len) < 2).any():
        return None
    evidence = all(hints)
    if not evidence and len(per_group) >= 2:           # crossed design
        shared = set.intersection(*per_group.tolist())
        evidence = len(shared) >= 2
    if not evidence:
        return None
    name = "level"
    if all(hints):
        for keys, nm in _LEVEL_NAME:
            if all(h in keys for h in hints):
                name = nm
                break
    return g, lv, name


def detect_design(conditions: list, samples: pd.DataFrame | None = None,
                  overrides: dict | None = None) -> Design:
    """overrides: {"categorical": True}; {"order_column": ..., "group_column": ...}
    (samples.tsv columns); or {"order_regex": r"T(\\d+)"} to force the reading."""
    conds = list(conditions)
    ov = overrides or {}
    if ov.get("categorical"):
        return Design(conditions=conds, source="forced categorical")
    if samples is not None and ov.get("order_column") in samples.columns:
        s = samples.drop_duplicates("condition").set_index("condition")
        lv = {c: float(s.at[c, ov["order_column"]]) for c in conds}
        gcol = ov.get("group_column")
        g = {c: str(s.at[c, gcol]) for c in conds} if gcol in s.columns else None
        return _ordered(conds, g, lv, str(ov["order_column"]), "samples.tsv")
    if ov.get("order_regex"):
        pat = re.compile(ov["order_regex"])
        g, lv = {}, {}
        for c in conds:
            m = pat.search(c)
            if not m:
                return Design(conditions=conds)
            lv[c] = float(m.group(1))
            g[c] = _clean(c[:m.start()], c[m.end():]) or "all"
        return _ordered(conds, g, lv, "level", "order_regex")
    for pat in _PATTERNS:
        res = _try(conds, pat)
        if res:
            g, lv, name = res
            return _ordered(conds, g, lv, name, "condition names")
    return Design(conditions=conds)


def _ordered(conds, g, lv, name, source):
    groups = list(dict.fromkeys(g.values())) if g else ["all"]
    g = g if g and len(groups) > 1 else None
    order = sorted(conds, key=lambda c: ((groups.index(g[c]) if g else 0), lv[c]))
    return Design(conditions=order, group=g, level=lv, level_name=name, source=source)
