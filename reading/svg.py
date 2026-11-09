"""Figures of the biological reading, as interactive SVG in Pathfinder's dark style.

Each function returns (list of SVG strings, facts for the text) or None when the data
it needs are missing. Hover: every mark carries a <title> (name, log2FC, p, programme,
group...); legend entries tagged `data-t` hide / show their series on click (READING_JS
in build_report_html). The SVGs export as they are (SVG / PNG buttons of the report).

Palette: fixed-order categorical (identity of the design groups), red / blue diverging
(direction), grey = not significant. Same colours as the rest of the report.
"""
from __future__ import annotations

import html
import re

import numpy as np
import pandas as pd

from build_report_html import (BG, COND_PALETTE, DOWN, FG, GRID, MUT, PANEL, UP,  # noqa: F401
                               _diverging)

FONT = "-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif"
WHITE = "#ffffff"
FAINT = "#3d444d"
GROUP_COLS = ["#3fb950", "#bc8cff", "#d29922"]     # functional groups of a volcano
CAT = COND_PALETTE


# --------------------------------------------------------------------------- helpers

def _esc(s) -> str:
    return html.escape(str(s), quote=True)


def _n(v) -> str:
    return f"{float(v):.1f}"


def _short(s, n=34) -> str:
    s = str(s)
    return s if len(s) <= n else s[:n - 1] + "…"


def _tw(s, fs) -> float:
    """Approximate text width (px)."""
    return len(str(s)) * fs * 0.56


def _p(p) -> str:
    p = float(p)
    if p == 0:
        return "0"
    return f"{p:.2e}" if p < 1e-3 else f"{p:.3g}"


def _seq(t) -> str:
    """Blue sequential scale (ordered level), t in [0, 1]."""
    a, b = np.array([13, 65, 157]), np.array([165, 214, 255])
    r, g, bl = (a + (b - a) * max(0.0, min(1.0, t))).astype(int)
    return f"rgb({r},{g},{bl})"


def nice_ticks(lo, hi, n=5):
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return [lo]
    raw = (hi - lo) / max(1, n)
    mag = 10 ** np.floor(np.log10(raw))
    step = min((m * mag for m in (1, 2, 2.5, 5, 10) if m * mag >= raw), default=raw)
    start = np.ceil(lo / step) * step
    out = list(np.arange(start, hi + step * 0.001, step))
    return [0.0 if abs(t) < step * 1e-9 else float(t) for t in out]


def _fmt_tick(t) -> str:
    return f"{t:g}" if abs(t) >= 1e-3 or t == 0 else f"{t:.0e}"


def svg_open(W, H, label="") -> str:
    return (f"<svg class='rd' viewBox='0 0 {W:.0f} {H:.0f}' width='{W:.0f}' height='{H:.0f}' "
            f"style='max-width:100%;height:auto' font-family='{FONT}' role='img' "
            f"aria-label='{_esc(label)}'>")


def title_block(x, y, title, subtitle=None) -> str:
    s = (f"<text x='{_n(x)}' y='{_n(y)}' fill='{WHITE}' font-size='13' "
         f"font-weight='600'>{_esc(title)}</text>")
    if subtitle:
        s += (f"<text x='{_n(x)}' y='{_n(y + 16)}' fill='{MUT}' font-size='10'>"
              f"{_esc(subtitle)}</text>")
    return s


class Ax:
    """Axes: pixel box (x0, y0, w, h) and data limits."""

    def __init__(self, x0, y0, w, h, xlim, ylim):
        self.x0, self.y0, self.w, self.h = x0, y0, w, h
        self.xl = (float(xlim[0]), float(xlim[1]) if xlim[1] != xlim[0] else float(xlim[0]) + 1)
        self.yl = (float(ylim[0]), float(ylim[1]) if ylim[1] != ylim[0] else float(ylim[0]) + 1)

    def X(self, v):
        return self.x0 + (np.asarray(v, float) - self.xl[0]) / (self.xl[1] - self.xl[0]) * self.w

    def Y(self, v):
        return self.y0 + self.h - (np.asarray(v, float) - self.yl[0]) / (self.yl[1] - self.yl[0]) * self.h

    def frame(self, xlabel="", ylabel="", xticks=None, yticks=None, grid=True,
              xfmt=_fmt_tick, yfmt=_fmt_tick, xtick_labels=None) -> str:
        s = []
        xt = nice_ticks(*self.xl) if xticks is None else xticks
        yt = nice_ticks(*self.yl) if yticks is None else yticks
        for i, t in enumerate(xt):
            if not (self.xl[0] - 1e-9 <= t <= self.xl[1] + 1e-9):
                continue
            x = float(self.X(t))
            if grid:
                s.append(f"<line x1='{_n(x)}' y1='{_n(self.y0)}' x2='{_n(x)}' "
                         f"y2='{_n(self.y0 + self.h)}' stroke='{GRID}' stroke-width='0.6'/>")
            lab = xtick_labels[i] if xtick_labels else xfmt(t)
            s.append(f"<text x='{_n(x)}' y='{_n(self.y0 + self.h + 14)}' fill='{MUT}' "
                     f"font-size='10' text-anchor='middle'>{_esc(lab)}</text>")
        for t in yt:
            if not (self.yl[0] - 1e-9 <= t <= self.yl[1] + 1e-9):
                continue
            y = float(self.Y(t))
            if grid:
                s.append(f"<line x1='{_n(self.x0)}' y1='{_n(y)}' x2='{_n(self.x0 + self.w)}' "
                         f"y2='{_n(y)}' stroke='{GRID}' stroke-width='0.6'/>")
            s.append(f"<text x='{_n(self.x0 - 6)}' y='{_n(y + 3.5)}' fill='{MUT}' "
                     f"font-size='10' text-anchor='end'>{_esc(yfmt(t))}</text>")
        s.append(f"<line x1='{_n(self.x0)}' y1='{_n(self.y0 + self.h)}' x2='{_n(self.x0 + self.w)}' "
                 f"y2='{_n(self.y0 + self.h)}' stroke='{MUT}' stroke-width='0.8'/>")
        s.append(f"<line x1='{_n(self.x0)}' y1='{_n(self.y0)}' x2='{_n(self.x0)}' "
                 f"y2='{_n(self.y0 + self.h)}' stroke='{MUT}' stroke-width='0.8'/>")
        if xlabel:
            s.append(f"<text x='{_n(self.x0 + self.w / 2)}' y='{_n(self.y0 + self.h + 32)}' "
                     f"fill='{FG}' font-size='11' text-anchor='middle'>{_esc(xlabel)}</text>")
        if ylabel:
            yc = self.y0 + self.h / 2
            s.append(f"<text x='{_n(self.x0 - 40)}' y='{_n(yc)}' fill='{FG}' font-size='11' "
                     f"text-anchor='middle' transform='rotate(-90 {_n(self.x0 - 40)} {_n(yc)})'>"
                     f"{_esc(ylabel)}</text>")
        return "".join(s)


class Placer:
    """Greedy non-overlapping label placement (replaces adjustText): candidate positions
    around the point, further and further away; a leader line beyond the first ring."""

    def __init__(self, bounds, fs=10):
        self.b = bounds              # x0, y0, x1, y1
        self.boxes = []
        self.fs = fs

    def _free(self, bx):
        x0, y0, x1, y1 = bx
        if x0 < self.b[0] or y0 < self.b[1] or x1 > self.b[2] or y1 > self.b[3]:
            return False
        return all(x1 < a or x0 > c or y1 < b or y0 > d for a, b, c, d in self.boxes)

    def block(self, bx):
        self.boxes.append(bx)

    def place(self, x, y, text, r=4):
        w, h = _tw(text, self.fs), self.fs + 2
        for k, ring in enumerate((1.0, 2.2, 3.6, 5.2)):
            d = (r + 3) * ring
            cands = [(x + d, y - d / 2 - h / 2), (x - d - w, y - d / 2 - h / 2),
                     (x + d, y + d / 2), (x - d - w, y + d / 2),
                     (x - w / 2, y - d - h), (x - w / 2, y + d)]
            for cx, cy in cands:
                bx = (cx - 1, cy - 1, cx + w + 1, cy + h + 1)
                if self._free(bx):
                    self.boxes.append(bx)
                    return cx, cy + h - 3, k > 0
        return None


def _label(x, y, pos, text, color=FG, italic=True, fs=10) -> str:
    tx, ty, lead = pos
    s = ""
    if lead:
        ex = tx if tx > x else tx + _tw(text, fs)
        s += (f"<line x1='{_n(x)}' y1='{_n(y)}' x2='{_n(ex)}' y2='{_n(ty - fs / 2 + 1)}' "
              f"stroke='{MUT}' stroke-width='0.6'/>")
    st = " font-style='italic'" if italic else ""
    return s + (f"<text x='{_n(tx)}' y='{_n(ty)}' fill='{color}' font-size='{fs}'{st} "
                f"paint-order='stroke' stroke='{BG}' stroke-width='3' "
                f"stroke-linejoin='round'>{_esc(text)}</text>")


def _legend_item(x, y, color, text, toggle=None, shape="circle") -> str:
    attr = f" class='lg' data-t='{_esc(toggle)}' style='cursor:pointer'" if toggle else ""
    mark = (f"<circle cx='{_n(x + 5)}' cy='{_n(y - 4)}' r='4.5' fill='{color}'/>" if shape == "circle"
            else f"<rect x='{_n(x)}' y='{_n(y - 9)}' width='10' height='10' rx='2' fill='{color}'/>")
    return (f"<g{attr}>{mark}<text x='{_n(x + 15)}' y='{_n(y)}' fill='{FG}' font-size='10.5'>"
            f"{_esc(text)}</text></g>")


def zbar(x, y, vmax, uid, label="z", width=170, height=9) -> str:
    """Diverging colour bar (blue -> grey -> red), without a NaN swatch."""
    stops = "".join(f"<stop offset='{i / 24 * 100:.0f}%' stop-color='"
                    f"{_diverging(-vmax + i / 24 * 2 * vmax, vmax)}'/>" for i in range(25))
    s = [f"<defs><linearGradient id='{uid}' x1='0' y1='0' x2='1' y2='0'>{stops}"
         f"</linearGradient></defs>",
         f"<rect x='{x}' y='{y}' width='{width}' height='{height}' fill='url(#{uid})' rx='2' "
         f"stroke='{GRID}'/>"]
    for i, v in enumerate((-vmax, -vmax / 2, 0, vmax / 2, vmax)):
        lx = x + width * i / 4
        s.append(f"<line x1='{lx:.1f}' y1='{y + height}' x2='{lx:.1f}' y2='{y + height + 3}' "
                 f"stroke='{MUT}'/><text x='{lx:.1f}' y='{y + height + 12}' fill='{MUT}' "
                 f"font-size='8.5' text-anchor='middle'>{v:+.1f}</text>")
    s.append(f"<text x='{x + width + 8}' y='{y + height - 0.5}' fill='{MUT}' font-size='9'>"
             f"{_esc(label)}</text>")
    return "".join(s)


def group_colors(design):
    """Identity colour per design group, fixed order; beyond the palette, a single hue
    (the groups are then named in clear)."""
    if not design.group:
        return {}
    gs = design.groups
    if len(gs) > len(CAT):
        return {g: CAT[0] for g in gs}
    return {g: CAT[i] for i, g in enumerate(gs)}


def _prog_of(prog):
    return prog.labels.to_dict() if prog is not None else {}


# --------------------------------------------------------------------------- 1. landscape

def fig_landscape(run):
    rows = []
    for c in run.contrasts:
        s = run.sig(c).fillna(False)
        lf = run.de[c + "_diff"]
        rows.append((c, int((s & (lf > 0)).sum()), int((s & (lf < 0)).sum())))
    if not rows:
        return None
    df = pd.DataFrame(rows, columns=["contrast", "up", "down"])
    labs = [_short(run.short(c), 44) for c in df.contrast]
    left = min(300, max(110, max(_tw(l, 10.5) for l in labs) + 18))
    rh, top = 22, 62
    W = left + 560
    H = top + rh * len(df) + 48
    m = max(int(df.up.max()), int(df.down.max()), 1)
    ax = Ax(left, top, W - left - 30, rh * len(df), (-m * 1.2, m * 1.2), (0, len(df)))
    rule = run.rule_text()
    s = [svg_open(W, H, "How much changes in each comparison"),
         title_block(14, 20, "How much changes in each comparison",
                     f"Differential proteins ({rule})")]
    xt = [t for t in nice_ticks(-m, m, 6)]
    s.append(ax.frame("significant proteins (down | up)", "", xticks=xt, yticks=[],
                      xfmt=lambda t: f"{abs(t):g}"))
    x0 = float(ax.X(0))
    for i, (r, lab) in enumerate(zip(df.itertuples(), labs)):
        yc = top + rh * i + rh / 2
        tip = f"<title>{_esc(run.short(r.contrast))}: {r.up} up, {r.down} down</title>"
        xu, xd = float(ax.X(r.up)), float(ax.X(-r.down))
        s.append(f"<g>{tip}<rect x='{_n(x0)}' y='{_n(yc - 7)}' width='{_n(max(xu - x0, 0))}' "
                 f"height='14' fill='{UP}' rx='1.5'/><rect x='{_n(xd)}' y='{_n(yc - 7)}' "
                 f"width='{_n(max(x0 - xd, 0))}' height='14' fill='{DOWN}' rx='1.5'/>"
                 f"<rect x='{left}' y='{_n(yc - rh / 2)}' width='{_n(ax.w)}' height='{rh}' "
                 f"fill='transparent'/></g>")
        s.append(f"<text x='{_n(xu + 4)}' y='{_n(yc + 3.5)}' fill='{FG}' font-size='9.5'>{r.up}</text>")
        s.append(f"<text x='{_n(xd - 4)}' y='{_n(yc + 3.5)}' fill='{FG}' font-size='9.5' "
                 f"text-anchor='end'>{r.down}</text>")
        s.append(f"<text x='{_n(left - 8)}' y='{_n(yc + 3.5)}' fill='{FG}' font-size='10.5' "
                 f"text-anchor='end'>{_esc(lab)}</text>")
    s.append(f"<line x1='{_n(x0)}' y1='{top}' x2='{_n(x0)}' y2='{_n(top + ax.h)}' stroke='{MUT}'/>")
    s.append(_legend_item(W - 150, 44, DOWN, "down", shape="rect")
             + _legend_item(W - 90, 44, UP, "up", shape="rect"))
    s.append("</svg>")
    top3 = df.assign(tot=df.up + df.down).sort_values("tot", ascending=False)
    facts = {"n_contrasts": len(df), "empty": int(((df.up + df.down) == 0).sum()),
             "top": [(run.short(r.contrast), r.up, r.down) for r in top3.head(3).itertuples()]}
    return ["".join(s)], facts


# --------------------------------------------------------------------------- 2. condition map

def _ellipse(ax, pts, color):
    if len(pts) < 3:
        return ""
    C = np.cov(pts.T)
    w, V = np.linalg.eigh(C)
    cx, cy = pts.mean(0)
    sx = ax.w / (ax.xl[1] - ax.xl[0])
    sy = ax.h / (ax.yl[1] - ax.yl[0])
    # semi-axes in pixels (1 SD); angle in screen coordinates (y pointing down)
    v = V[:, 1]
    ang = np.degrees(np.arctan2(-v[1] * sy, v[0] * sx))
    a = np.sqrt(max(w[1], 0)) * np.hypot(v[0] * sx, v[1] * sy)
    u = V[:, 0]
    b = np.sqrt(max(w[0], 0)) * np.hypot(u[0] * sx, u[1] * sy)
    X, Y = float(ax.X(cx)), float(ax.Y(cy))
    return (f"<ellipse cx='{_n(X)}' cy='{_n(Y)}' rx='{_n(a)}' ry='{_n(b)}' "
            f"transform='rotate({ang:.1f} {_n(X)} {_n(Y)})' fill='{color}' fill-opacity='0.09' "
            f"stroke='{color}' stroke-opacity='0.25'/>")


SHAPES = ("circle", "triangle", "square", "diamond", "tri-down", "cross")


def _marker(shape, X, Y, r, fill, title="", stroke=BG):
    """One mark of a given shape (factor-coded maps and legends)."""
    t = f"<title>{_esc(title)}</title>" if title else ""
    st = f"stroke='{stroke}' stroke-width='1.2'"
    if shape == "circle":
        return f"<circle class='pt' cx='{_n(X)}' cy='{_n(Y)}' r='{r:.1f}' fill='{fill}' {st}>{t}</circle>"
    if shape == "square":
        return (f"<rect class='pt' x='{_n(X - r * .9)}' y='{_n(Y - r * .9)}' width='{_n(1.8 * r)}' "
                f"height='{_n(1.8 * r)}' fill='{fill}' {st}>{t}</rect>")
    pts = {"triangle": [(0, -1.25), (1.1, .75), (-1.1, .75)],
           "tri-down": [(0, 1.25), (1.1, -.75), (-1.1, -.75)],
           "diamond": [(0, -1.3), (1.1, 0), (0, 1.3), (-1.1, 0)],
           "cross": [(-.4, -1.1), (.4, -1.1), (.4, -.4), (1.1, -.4), (1.1, .4), (.4, .4), (.4, 1.1),
                     (-.4, 1.1), (-.4, .4), (-1.1, .4), (-1.1, -.4), (-.4, -.4)]}[shape]
    pp = " ".join(f"{_n(X + dx * r)},{_n(Y + dy * r)}" for dx, dy in pts)
    return f"<polygon class='pt' points='{pp}' fill='{fill}' {st}>{t}</polygon>"


def fig_map(run, design, F=None, roles=None):
    sc = run.scores
    if sc is None or "PC1" not in sc.columns or "PC2" not in sc.columns:
        return None
    sc = sc.copy()
    sc["condition"] = sc["condition"].astype(str)
    lab_col = "label" if "label" in sc.columns else None
    v1 = sc.filter(like="PC1_var").iloc[0, 0] if len(sc.filter(like="PC1_var").columns) else None
    v2 = sc.filter(like="PC2_var").iloc[0, 0] if len(sc.filter(like="PC2_var").columns) else None
    cen = sc.groupby("condition")[["PC1", "PC2"]].mean()
    nrep = sc.groupby("condition").size()
    gcol = group_colors(design)
    W, H = 780, 600
    xs, ys = sc.PC1.values.astype(float), sc.PC2.values.astype(float)
    px, py = np.ptp(xs) or 1, np.ptp(ys) or 1
    ax = Ax(70, 60, 560, 470, (xs.min() - 0.14 * px, xs.max() + 0.14 * px),
            (ys.min() - 0.14 * py, ys.max() + 0.14 * py))
    traj = design.ordered and len(design.groups or ["all"]) <= 6
    facts = {"var": (v1, v2)}
    if F is not None and roles and roles.get("lines") and not traj:
        return _factor_map(run, sc, cen, nrep, ax, F, roles, v1, v2, facts)
    if traj:
        sub = (f"Condition centroids linked in order of {design.level_name} · replicates "
               "faint · ellipses = 1 SD")
    else:
        sub = "Condition centroids · replicates faint · ellipses = 1 SD"
        if design.ordered:
            sub += (f" · colour = {design.level_name} ({len(design.groups)} groups, too many "
                    "for trajectories)")
    s = [svg_open(W, H, "Where the conditions sit in protein space"),
         title_block(14, 20, "Where the conditions sit in protein space", sub),
         ax.frame(f"PC1 ({v1:.1f} %)" if v1 is not None else "PC1",
                  f"PC2 ({v2:.1f} %)" if v2 is not None else "PC2")]
    if ax.xl[0] < 0 < ax.xl[1]:
        s.append(f"<line x1='{_n(ax.X(0))}' y1='{ax.y0}' x2='{_n(ax.X(0))}' y2='{ax.y0 + ax.h}' "
                 f"stroke='{MUT}' stroke-width='0.5' stroke-dasharray='3 3'/>")
    if ax.yl[0] < 0 < ax.yl[1]:
        s.append(f"<line x1='{ax.x0}' y1='{_n(ax.Y(0))}' x2='{ax.x0 + ax.w}' y2='{_n(ax.Y(0))}' "
                 f"stroke='{MUT}' stroke-width='0.5' stroke-dasharray='3 3'/>")
    lv_cols = {}
    if design.ordered and not traj:
        lv = design.levels
        lv_cols = {l: _seq(0.25 + 0.75 * i / max(1, len(lv) - 1)) for i, l in enumerate(lv)}

    def col_of(c):
        if lv_cols:
            return lv_cols[design.level[c]]
        if design.group:
            return gcol.get(design.group[c], CAT[0])
        return CAT[0]

    for c in design.conditions:
        if c not in cen.index:
            continue
        col = col_of(c)
        r = sc[sc.condition == c]
        s.append(_ellipse(ax, r[["PC1", "PC2"]].values.astype(float), col))
        for _, row in r.iterrows():
            lab = row[lab_col] if lab_col else ""
            s.append(f"<circle class='pt' cx='{_n(ax.X(row.PC1))}' cy='{_n(ax.Y(row.PC2))}' r='3.2' "
                     f"fill='{col}' fill-opacity='0.45'><title>{_esc(lab)} · {_esc(c)}\n"
                     f"PC1 {row.PC1:.2f} · PC2 {row.PC2:.2f}</title></circle>")
    placer = Placer((ax.x0, ax.y0, ax.x0 + ax.w, ax.y0 + ax.h), fs=10)
    labels = []
    if traj:
        groups = design.groups or ["all"]
        steps = []
        defs = []
        for gi, g in enumerate(groups):
            cs = [c for c in design.conditions if (design.group or {}).get(c, "all") == g
                  and c in cen.index]
            cs = sorted(cs, key=lambda c: design.level[c])
            col = gcol.get(g, CAT[0])
            mid = f"rdah{gi}"
            defs.append(f"<marker id='{mid}' viewBox='0 0 10 10' refX='9' refY='5' markerWidth='7' "
                        f"markerHeight='7' orient='auto-start-reverse'><path d='M0,0 L10,5 L0,10 z' "
                        f"fill='{col}'/></marker>")
            for a, b in zip(cs[:-1], cs[1:]):
                xa, ya = float(ax.X(cen.loc[a, "PC1"])), float(ax.Y(cen.loc[a, "PC2"]))
                xb, yb = float(ax.X(cen.loc[b, "PC1"])), float(ax.Y(cen.loc[b, "PC2"]))
                d = np.hypot(xb - xa, yb - ya) or 1
                sh = min(8, d / 3)
                xa2, ya2 = xa + (xb - xa) * sh / d, ya + (yb - ya) * sh / d
                xb2, yb2 = xb - (xb - xa) * sh / d, yb - (yb - ya) * sh / d
                dist = float(np.linalg.norm(cen.loc[b] - cen.loc[a]))
                steps.append((g, a, b, dist))
                s.append(f"<line x1='{_n(xa2)}' y1='{_n(ya2)}' x2='{_n(xb2)}' y2='{_n(yb2)}' "
                         f"stroke='{col}' stroke-width='2' marker-end='url(#{mid})'>"
                         f"<title>{_esc(a)} → {_esc(b)} · distance {dist:.2f}</title></line>")
            for i, c in enumerate(cs):
                rr = 4 + 4 * i / max(1, len(cs) - 1)
                X, Y = float(ax.X(cen.loc[c, "PC1"])), float(ax.Y(cen.loc[c, "PC2"]))
                placer.block((X - rr, Y - rr, X + rr, Y + rr))
                s.append(f"<circle class='pt' cx='{_n(X)}' cy='{_n(Y)}' r='{rr:.1f}' fill='{col}' "
                         f"stroke='{BG}' stroke-width='1.4'><title>{_esc(c)} · centroid of "
                         f"{int(nrep.get(c, 0))} replicate(s)</title></circle>")
                labels.append((X, Y, design.label(design.level[c]), MUT, False, rr))
            if cs and design.group:
                X, Y = float(ax.X(cen.loc[cs[-1], "PC1"])), float(ax.Y(cen.loc[cs[-1], "PC2"]))
                labels.insert(0, (X, Y, _short(g, 22), col, False, 9))
        s.insert(2, "<defs>" + "".join(defs) + "</defs>")
        facts["steps"] = sorted(steps, key=lambda t: -t[3])[:3]
    else:
        present = [c for c in design.conditions if c in cen.index]
        if len(present) <= 16:
            named = set(present)
        else:
            dist = np.linalg.norm(cen.loc[present] - cen.loc[present].mean(), axis=1)
            named = set(np.array(present)[np.argsort(-dist)[:12]])
        for c in present:
            X, Y = float(ax.X(cen.loc[c, "PC1"])), float(ax.Y(cen.loc[c, "PC2"]))
            placer.block((X - 6, Y - 6, X + 6, Y + 6))
            s.append(f"<circle class='pt' cx='{_n(X)}' cy='{_n(Y)}' r='6' fill='{col_of(c)}' "
                     f"stroke='{BG}' stroke-width='1.4'><title>{_esc(c)} · centroid of "
                     f"{int(nrep.get(c, 0))} replicate(s)</title></circle>")
            if c in named:
                labels.append((X, Y, _short(design.display(c), 24), FG, False, 6))
        if lv_cols:
            y = 80
            s.append(f"<text x='650' y='{y - 14}' fill='{MUT}' font-size='10'>"
                     f"{_esc(design.level_name)}</text>")
            for l, colr in lv_cols.items():
                s.append(_legend_item(650, y, colr, design.label(l)))
                y += 17
        elif design.group and len(design.groups) <= len(CAT):
            y = 80
            for g in design.groups:
                s.append(_legend_item(650, y, gcol[g], _short(g, 18)))
                y += 17
    for X, Y, text, colr, it, rr in labels:
        pos = placer.place(X, Y, text, r=rr)
        if pos:
            bold = colr not in (FG, MUT)
            lab = _label(X, Y, pos, text, colr, italic=it, fs=11 if bold else 10)
            if bold:
                lab = lab.replace("<text ", "<text font-weight='600' ", 1)
            s.append(lab)
    if traj and design.group:
        y = 80
        for g in design.groups:
            s.append(_legend_item(650, y, gcol.get(g, CAT[0]), _short(g, 18)))
            y += 17
    d = cen.values
    if len(d) >= 2:
        D = np.linalg.norm(d[:, None] - d[None], axis=2)
        i, j = np.unravel_index(D.argmax(), D.shape)
        facts["farthest"] = (cen.index[i], cen.index[j])
    s.append("</svg>")
    return ["".join(s)], facts


def _factor_map(run, sc, cen, nrep, ax, F, r, v1, v2, facts):
    """Condition map of a factorial design: colour = lines factor (genotype), shape =
    first facet factor (light), size = x level (time)."""
    from .factorial import MISSING, levels
    W, H = 780, 600
    lv_line = levels(F, r["lines"])
    colf = {l: CAT[i % len(CAT)] for i, l in enumerate(lv_line)}
    fac = r["facets"][0] if r["facets"] else None
    lv_fac = levels(F, fac) if fac else []
    shp = {l: SHAPES[i % len(SHAPES)] for i, l in enumerate(lv_fac)}
    lv_x = levels(F, r["x"])
    size = {l: 4.5 + 4.5 * i / max(1, len([v for v in lv_x if v != MISSING]) - 1)
            for i, l in enumerate(v for v in lv_x if v != MISSING)}
    size[MISSING] = 5
    sub = (f"Condition centroids · colour = {r['lines']}" + (f" · shape = {fac}" if fac else "")
           + f" · size = {r['x']} · replicates faint · hover a mark")
    s = [svg_open(W, H, "Where the conditions sit in protein space"),
         title_block(14, 20, "Where the conditions sit in protein space", sub),
         ax.frame(f"PC1 ({v1:.1f} %)" if v1 is not None else "PC1",
                  f"PC2 ({v2:.1f} %)" if v2 is not None else "PC2")]
    present = [c for c in F.index if c in cen.index]
    for c in present:
        col = colf.get(F.at[c, r["lines"]], CAT[0])
        for _, row in sc[sc.condition == c].iterrows():
            s.append(f"<circle cx='{_n(ax.X(row.PC1))}' cy='{_n(ax.Y(row.PC2))}' r='2.2' "
                     f"fill='{col}' fill-opacity='0.35'/>")
    gl = {}
    for c in present:
        line = F.at[c, r["lines"]]
        X, Y = float(ax.X(cen.loc[c, "PC1"])), float(ax.Y(cen.loc[c, "PC2"]))
        desc = " · ".join(f"{f} {F.at[c, f]}" for f in F.columns if F.at[c, f] != MISSING)
        mk = _marker(shp.get(F.at[c, fac], "circle") if fac else "circle", X, Y,
                     size.get(F.at[c, r["x"]], 6), colf.get(line, CAT[0]),
                     f"{c}\n{desc}\ncentroid of {int(nrep.get(c, 0))} replicate(s)")
        gl.setdefault(line, []).append(mk)
    for i, l in enumerate(lv_line):
        s.append(f"<g data-g='l{i}'>" + "".join(gl.get(l, [])) + "</g>")
    y = 80
    s.append(f"<text x='650' y='{y - 14}' fill='{MUT}' font-size='10'>{_esc(_short(r['lines'], 20))}</text>")
    for i, l in enumerate(lv_line):
        s.append(_legend_item(650, y, colf[l], _short(l, 16), toggle=f"l{i}"))
        y += 16
    if fac:
        y += 10
        s.append(f"<text x='650' y='{y}' fill='{MUT}' font-size='10'>{_esc(_short(fac, 20))}</text>")
        y += 16
        for l in lv_fac:
            s.append(_marker(shp[l], 655, y - 4, 5, MUT, stroke="none")
                     + f"<text x='665' y='{y}' fill='{FG}' font-size='10.5'>{_esc(_short(l, 14))}</text>")
            y += 16
    y += 10
    s.append(f"<text x='650' y='{y}' fill='{MUT}' font-size='10'>size = {_esc(_short(r['x'], 14))}</text>")
    y += 14
    s.append(f"<text x='650' y='{y}' fill='{FG}' font-size='10'>"
             f"{_esc(' < '.join(v for v in lv_x if v != MISSING))}</text>")
    s.append("</svg>")
    d = cen.loc[present].values
    if len(d) >= 2:
        D = np.linalg.norm(d[:, None] - d[None], axis=2)
        i, j = np.unravel_index(D.argmax(), D.shape)
        facts["farthest"] = (present[i], present[j])
    # how far apart the levels of each factor sit on the map (centroid spread)
    spread = {}
    for f in F.columns:
        g = cen.loc[present].groupby(F.loc[present, f].values).mean()
        g = g[g.index != MISSING]
        if len(g) >= 2:
            spread[f] = float(np.linalg.norm(g - g.mean(), axis=1).mean())
    facts["factor_spread"] = spread
    return ["".join(s)], facts


# --------------------------------------------------------------------------- 3. programmes

def _top_names(run, prog, p, n=10):
    members = prog.labels.index[prog.labels == p]
    lfc = run.de.loc[members, [c + "_diff" for c in run.contrasts]].abs().max(axis=1)
    return [str(run.labels.get(i, i)) for i in lfc.sort_values(ascending=False).index[:n]]


def fig_programmes(run, design, prog, order=None, blocks=None):
    """order: condition order (factorial designs: genotype > light > time); blocks:
    condition -> block label for the heatmap separators."""
    if prog is None:
        return None
    k = prog.k
    conds = [c for c in design.conditions if c in prog.z.columns]
    gcol = group_colors(design)
    lines = design.ordered and len(design.levels) >= 3 and len(design.groups or ["all"]) <= 6
    if not lines and len(conds) > 14:
        return _programme_heatmap(run, design, prog, order, blocks)
    ncol = 3 if k > 4 else 2
    nrow = int(np.ceil(k / ncol))
    pw, ph = 300, 190
    gx, gy = 70, 92
    top = 70
    W = 20 + ncol * (pw + gx)
    H = top + nrow * (ph + gy) + 10
    q1 = min(float(prog.z[prog.labels == p][conds].quantile(0.25).min()) for p in range(1, k + 1))
    q3 = max(float(prog.z[prog.labels == p][conds].quantile(0.75).max()) for p in range(1, k + 1))
    lo, hi = min(q1, -0.5) - 0.15, max(q3, 0.5) + 0.15
    s = [svg_open(W, H, "Response programmes"),
         title_block(14, 20, "Response programmes: proteins that move together",
                     f"{len(prog.z)} differential proteins grouped by profile shape (k-means, "
                     f"k={k}) · " + ("median ± IQR" if lines else "bar = median, whisker = IQR")
                     + " · subtitle = enriched functional groups, or proteins · hover for values")]
    if lines and design.group:
        x = 14
        for g in design.groups:
            s.append(_legend_item(x, 56, gcol.get(g, CAT[0]), _short(g, 18)))
            x += 30 + _tw(_short(g, 18), 10.5)
    for p in range(1, k + 1):
        r, cidx = divmod(p - 1, ncol)
        x0 = 60 + cidx * (pw + gx)
        y0 = top + 26 + r * (ph + gy)
        Z = prog.z[prog.labels == p]
        names = _top_names(run, prog, p)
        tip = (f"P{p} · {len(Z)} proteins · highest in {prog.high[p]}, lowest in {prog.low[p]}\n"
               f"most marked: {', '.join(names)}")
        s.append(f"<g><title>{_esc(tip)}</title><rect x='{x0}' y='{y0}' width='{pw}' height='{ph}' "
                 f"fill='{PANEL}' opacity='0.35'/></g>")
        if lines:
            lv = design.levels
            ax = Ax(x0, y0, pw, ph, (-0.4, len(lv) - 0.6), (lo, hi))
            s.append(ax.frame("", "z-score" if cidx == 0 else "", xticks=list(range(len(lv))),
                              xtick_labels=[design.label(l) for l in lv]))
            for g in (design.groups or ["all"]):
                cs = {design.level[c]: c for c in conds if (design.group or {}).get(c, "all") == g}
                cols = [cs[l] for l in lv if l in cs]
                xs = [lv.index(l) for l in lv if l in cs]
                m = Z[cols]
                col = gcol.get(g, CAT[0])
                med, qa, qb = m.median(), m.quantile(0.25), m.quantile(0.75)
                band = ([f"{_n(ax.X(x))},{_n(ax.Y(v))}" for x, v in zip(xs, qb)]
                        + [f"{_n(ax.X(x))},{_n(ax.Y(v))}" for x, v in zip(xs[::-1], qa[::-1])])
                gl = f" data-g='{_esc(g)}'" if design.group else ""
                s.append(f"<g{gl}><polygon points='{' '.join(band)}' fill='{col}' fill-opacity='0.16'/>"
                         f"<polyline points='{' '.join(f'{_n(ax.X(x))},{_n(ax.Y(v))}' for x, v in zip(xs, med))}' "
                         f"fill='none' stroke='{col}' stroke-width='2'/>")
                for x, c, v, a, b in zip(xs, cols, med, qa, qb):
                    s.append(f"<circle class='pt' cx='{_n(ax.X(x))}' cy='{_n(ax.Y(v))}' r='3.5' "
                             f"fill='{col}' stroke='{BG}'><title>P{p} · {_esc(c)}\nmedian z "
                             f"{v:+.2f} (IQR {a:+.2f} to {b:+.2f})</title></circle>")
                s.append("</g>")
        else:
            ax = Ax(x0, y0, pw, ph, (-0.6, len(conds) - 0.4), (lo, hi))
            s.append(ax.frame("", "z-score" if cidx == 0 else "", xticks=[]))
            med, qa, qb = Z[conds].median(), Z[conds].quantile(0.25), Z[conds].quantile(0.75)
            bw = min(26, pw / max(1, len(conds)) * 0.62)
            for i, c in enumerate(conds):
                X = float(ax.X(i))
                y0b, y1b = float(ax.Y(0)), float(ax.Y(med[c]))
                col = UP if med[c] > 0 else DOWN
                s.append(f"<g><title>P{p} · {_esc(c)}\nmedian z {med[c]:+.2f} (IQR {qa[c]:+.2f} "
                         f"to {qb[c]:+.2f})</title><rect class='pt' x='{_n(X - bw / 2)}' "
                         f"y='{_n(min(y0b, y1b))}' width='{_n(bw)}' height='{_n(abs(y1b - y0b))}' "
                         f"fill='{col}' fill-opacity='0.85'/><line x1='{_n(X)}' y1='{_n(ax.Y(qa[c]))}' "
                         f"x2='{_n(X)}' y2='{_n(ax.Y(qb[c]))}' stroke='{FG}' stroke-width='1'/></g>")
                lab = _short(design.display(c), 14)
                s.append(f"<text x='{_n(X)}' y='{_n(y0 + ph + 12)}' fill='{MUT}' font-size='9' "
                         f"text-anchor='end' transform='rotate(-45 {_n(X)} {_n(y0 + ph + 12)})'>"
                         f"{_esc(lab)}</text>")
        s.append(f"<line x1='{x0}' y1='{_n(ax.Y(0))}' x2='{x0 + pw}' y2='{_n(ax.Y(0))}' "
                 f"stroke='{MUT}' stroke-width='0.7' stroke-dasharray='2 3'/>")
        th = prog.themes.get(p, [])
        line2 = ("; ".join(_short(t[0], 30) for t in th[:2]) if th
                 else ", ".join(prog.names.get(p, [])[:4]))
        stab = ""
        if prog.stability:
            stab = f" · stability {prog.stability[p]:.2f}"
        s.append(f"<text x='{x0}' y='{y0 - 20}' fill='{WHITE}' font-size='12' font-weight='600'>"
                 f"<title>Stability = mean Jaccard between this programme and its match over "
                 f"{prog.n_boot} replicate resamplings (≥ 0.75 stable, &lt; 0.6 dissolves)</title>"
                 f"P{p} · {len(Z)} proteins{stab}</text>")
        s.append(f"<text x='{x0}' y='{y0 - 6}' fill='{MUT}' font-size='9.5'"
                 f"{'' if th else ' font-style=' + chr(39) + 'italic' + chr(39)}>"
                 f"{_esc(_short(line2, 58))}</text>")
    s.append("</svg>")
    return ["".join(s)], {"layout": "lines" if lines else "bars"}


def _programme_heatmap(run, design, prog, order=None, blocks=None):
    conds = [c for c in (order or design.conditions) if c in prog.z.columns]
    grp = blocks or design.group
    M = pd.DataFrame({p: prog.z[prog.labels == p][conds].median()
                      for p in range(1, prog.k + 1)}).T
    ylab = []
    for p in M.index:
        th = prog.themes.get(p, [])
        extra = _short(th[0][0], 34) if th else ", ".join(prog.names.get(p, [])[:3])
        st_ = f", J={prog.stability[p]:.2f}" if prog.stability else ""
        ylab.append(f"P{p} ({int((prog.labels == p).sum())}{st_}) · {extra}")
    left = min(330, max(_tw(l, 10) for l in ylab) + 20)
    cw = max(12, min(30, 900 / max(1, len(conds))))
    ch = 26
    top = 92 if grp else 70
    W = left + cw * len(conds) + 120
    H = top + ch * len(M) + 120
    lim = max(1.0, float(np.nanmax(np.abs(M.values))))
    s = [svg_open(W, H, "Response programmes"),
         title_block(14, 20, "Response programmes: proteins that move together",
                     f"{len(prog.z)} differential proteins grouped by profile shape (k={prog.k}) · "
                     "cell = median z-score of the programme · hover for values")]
    for i, p in enumerate(M.index):
        y = top + i * ch
        names = _top_names(run, prog, p, 8)
        s.append(f"<text x='{left - 8}' y='{y + ch * 0.65:.1f}' fill='{FG}' font-size='10' "
                 f"text-anchor='end'><title>{_esc(', '.join(names))}</title>{_esc(ylab[i])}</text>")
        for j, c in enumerate(conds):
            v = M.loc[p, c]
            s.append(f"<rect x='{left + j * cw:.1f}' y='{y}' width='{cw - 1:.1f}' height='{ch - 1}' "
                     f"fill='{_diverging(v, lim)}'><title>P{p} · {_esc(c)}: median z {v:+.2f}"
                     f"</title></rect>")
    for j, c in enumerate(conds):
        x = left + j * cw + cw / 2
        y = top + ch * len(M) + 8
        s.append(f"<text x='{x:.1f}' y='{y}' fill='{MUT}' font-size='9' text-anchor='end' "
                 f"transform='rotate(-60 {x:.1f} {y})'>{_esc(_short(design.display(c), 16))}</text>")
    if grp:
        gs = [grp[c] for c in conds]
        starts = [i for i in range(len(gs)) if i == 0 or gs[i] != gs[i - 1]]
        for i in starts[1:]:
            x = left + i * cw - 0.5
            s.append(f"<line x1='{x:.1f}' y1='{top}' x2='{x:.1f}' y2='{top + ch * len(M)}' "
                     f"stroke='{BG}' stroke-width='3'/>")
        for a, b in zip(starts, starts[1:] + [len(gs)]):
            s.append(f"<text x='{left + (a + b) / 2 * cw:.1f}' y='{top - 6}' fill='{FG}' "
                     f"font-size='9.5' text-anchor='middle'>"
                     f"{_esc(_short(gs[a], max(6, int((b - a) * cw / 6))))}</text>")
    s.append(zbar(left, H - 30, round(lim, 1), "rdcbprog", "median z-score"))
    s.append("</svg>")
    return ["".join(s)], {"layout": "heatmap"}


# --------------------------------------------------------------------------- 4. functional groups

def _spring(n, Wm, iters=250, seed=7):
    """Force-directed layout (Fruchterman-Reingold), weights = Jaccard."""
    rng = np.random.default_rng(seed)
    P = rng.uniform(-1, 1, (n, 2))
    if n == 1:
        return np.zeros((1, 2))
    k = 1.2 / np.sqrt(n)
    t = 0.12
    for _ in range(iters):
        d = P[:, None, :] - P[None, :, :]
        dist = np.linalg.norm(d, axis=2) + 1e-6
        rep = (k * k / dist ** 2)[..., None] * d
        att = (Wm * dist / k)[..., None] * d
        disp = (rep - att).sum(1)
        ln = np.linalg.norm(disp, axis=1, keepdims=True) + 1e-9
        P += disp / ln * np.minimum(ln, t)
        t *= 0.985
    P -= P.mean(0)
    return P / (np.abs(P).max() or 1)


def _packed(groups, keep, edge_j):
    """One block per group (force-directed on its own links), blocks packed on shelves,
    most significant first."""
    from .groups import jaccard
    terms = groups.terms
    blocks = []
    for t in keep:
        sub = terms[terms.theme == t].sort_values("p")
        ids = list(sub.term)
        gs = list(sub.genes)
        n = len(ids)
        Wm = np.zeros((n, n))
        edges = []
        for i in range(n):
            for j in range(i + 1, n):
                jj = jaccard(gs[i], gs[j])
                if jj >= edge_j:
                    Wm[i, j] = Wm[j, i] = jj
                    edges.append((ids[i], ids[j], jj))
        P = _spring(n, Wm)
        r = 0.35 + 0.30 * np.sqrt(n)
        blocks.append((t, dict(zip(ids, P * r)), r, edges))
    pad = 1.1
    total = sum((2 * r + pad) ** 2 for _, _, r, _ in blocks)
    # shelf width aiming at a width / height ratio of ~1.6 (wide panel)
    width = max(max(2 * r + pad for _, _, r, _ in blocks), np.sqrt(total * 1.6) * 1.1) \
        if blocks else 1
    pos, x, y, row_h, out = {}, 0.0, 0.0, 0.0, []
    for t, p, r, edges in blocks:
        d = 2 * r + pad
        if x > 0 and x + d > width:
            x, y, row_h = 0.0, y - row_h, 0.0
        c = np.array([x + d / 2, y - d / 2])
        for k_, v in p.items():
            pos[k_] = v + c
        out.append((t, c, r, edges))
        x += d
        row_h = max(row_h, d)
    return pos, out


EDGE_J = 0.35          # line between two terms of a group: Jaccard >= 0.35


def fig_groups(run, groups, n_groups=16):
    if groups is None:
        return None
    keep = groups.top(n_groups)
    if not keep:
        return None
    order = [c for c in run.contrasts if c in set(groups.table.contrast)]
    tab = groups.table[groups.table.theme.isin(keep)]
    vmax = 2.5
    size_t = groups.terms.groupby("theme").term.nunique()
    gsz = dict(zip(groups.terms.term, groups.terms.genes.map(len)))
    zt = groups.hits.groupby("term").z.median()
    info = groups.terms.set_index("term")

    # A. enrichment map
    pos, blocks = _packed(groups, keep, EDGE_J)
    x0 = min(c[0] - r for _, c, r, _ in blocks) - 0.3
    x1 = max(c[0] + r for _, c, r, _ in blocks) + 0.3
    y0 = min(c[1] - r for _, c, r, _ in blocks) - 0.6
    y1 = max(c[1] + r for _, c, r, _ in blocks) + 0.3
    WA = 760
    # one x / y scale (no distortion), capped for small maps
    scale = min((WA - 40) / (x1 - x0), 900 / (y1 - y0), 120)
    HA = 70 + (y1 - y0) * scale + 60
    xoff = 20 + ((WA - 40) - (x1 - x0) * scale) / 2

    def px(v):
        return xoff + (v[0] - x0) * scale, 64 + (y1 - v[1]) * scale

    a = [svg_open(WA, HA, "Enrichment map"),
         title_block(14, 20, "A  Enrichment map",
                     (f"{len(keep)} most significant groups of GSEA pathways · pathways linked "
                      "by shared leading-edge proteins · size = proteins · colour = NES"
                      if groups.origin == "GSEA" else
                      f"{len(keep)} most significant functional groups (Meta-GO redundancy groups) "
                      "· terms linked by shared proteins · size = proteins · colour = direction"))]
    for t, c, r, edges in blocks:
        for u, v, jj in edges:
            (xa, ya), (xb, yb) = px(pos[u]), px(pos[v])
            a.append(f"<line x1='{_n(xa)}' y1='{_n(ya)}' x2='{_n(xb)}' y2='{_n(yb)}' "
                     f"stroke='#57606a' stroke-width='{0.4 + 1.6 * jj:.2f}'/>")
    for t, c, r, edges in blocks:
        sub = groups.terms[groups.terms.theme == t]
        for term in sub.term:
            X, Y = px(pos[term])
            ng = gsz.get(term, 0)
            rad = np.sqrt((14 + 5 * ng) / np.pi) * 1.25
            z = float(zt.get(term, 0))
            row = info.loc[term]
            tip = (f"{row['name']} ({term})\n{row['source']} · {ng} proteins · best FDR "
                   f"{_p(row['p'])} · median z {z:+.2f}\ngroup G{t + 1}: {groups.names[t]}")
            a.append(f"<circle class='pt' cx='{_n(X)}' cy='{_n(Y)}' r='{rad:.1f}' "
                     f"fill='{_diverging(z, vmax)}' stroke='{BG}' stroke-width='0.8'>"
                     f"<title>{_esc(tip)}</title></circle>")
        pts = np.array([px(pos[term]) for term in sub.term])
        rads = np.array([np.sqrt((14 + 5 * gsz.get(term, 0)) / np.pi) * 1.25 for term in sub.term])
        ylow = float((pts[:, 1] + rads).max()) + 12
        med = (px(c)[0], 0)
        room = int((2 * r + 1.0) * scale / (10 * 0.56)) - 5       # characters fitting the block
        name = _short(groups.names[t], max(10, min(34, room)))
        a.append(f"<text x='{_n(med[0])}' y='{_n(ylow)}' fill='{FG}' font-size='10' "
                 f"text-anchor='middle' paint-order='stroke' stroke='{BG}' stroke-width='3'>"
                 f"<title>G{t + 1}: {_esc(groups.names[t])} · {int(size_t.get(t, 0))} terms</title>"
                 f"G{t + 1} · {_esc(name)}</text>")
    dname = "NES" if groups.origin == "GSEA" else "median z"
    a.append(zbar(20, HA - 34, vmax, "rdcbgrpA", f"direction ({dname}): down ← → up"))
    a.append("</svg>")

    # B. group x contrast matrix
    mat = tab.pivot(index="theme", columns="contrast", values="p").reindex(columns=order)
    mat = -np.log10(mat.astype(float)).fillna(0)
    if len(mat) > 2 and mat.shape[1] >= 1:
        from scipy.cluster.hierarchy import leaves_list, linkage
        ordr = list(mat.index[leaves_list(linkage(mat.values, "average"))])
    else:
        ordr = list(mat.index)
    rlabs = [f"G{t + 1} · {_short(groups.names[t], 42)} ({int(size_t.get(t, 0))})" for t in ordr]
    cw, ch = 34, 24
    top = 150
    left = 20
    lab_w = max(_tw(l, 10) for l in rlabs) + 14
    WB = max(left + cw * len(order) + lab_w + 20, left + 230 + 140 + 130)
    HB = top + ch * len(ordr) + 70
    b = [svg_open(WB, HB, "Functional groups across comparisons"),
         title_block(14, 20, "B  Functional groups across comparisons",
                     "size = −log10(FDR) of the best term · colour = direction (median "
                     + ("NES" if groups.origin == "GSEA" else "z") + ") · (n terms)")]
    for j, c in enumerate(order):
        x = left + j * cw + cw / 2
        b.append(f"<text x='{_n(x)}' y='{top - 10}' fill='{MUT}' font-size='9.5' "
                 f"transform='rotate(-60 {_n(x)} {top - 10})'>{_esc(_short(run.short(c), 26))}</text>")
        b.append(f"<line x1='{_n(x)}' y1='{top - 4}' x2='{_n(x)}' y2='{top + ch * len(ordr)}' "
                 f"stroke='{GRID}' stroke-width='0.5'/>")
    for i, t in enumerate(ordr):
        y = top + i * ch + ch / 2
        b.append(f"<line x1='{left}' y1='{_n(y)}' x2='{left + cw * len(order)}' y2='{_n(y)}' "
                 f"stroke='{GRID}' stroke-width='0.5'/>")
        b.append(f"<text x='{left + cw * len(order) + 10}' y='{_n(y + 3.5)}' fill='{FG}' "
                 f"font-size='10'><title>{_esc(groups.names[t])}</title>{_esc(rlabs[i])}</text>")
        sub = tab[tab.theme == t].set_index("contrast")
        for j, c in enumerate(order):
            if c not in sub.index:
                continue
            p, z, nt = float(sub.at[c, "p"]), float(sub.at[c, "z"]), int(sub.at[c, "n"])
            lp = -np.log10(max(p, 1e-300))
            rad = np.sqrt(30 + 70 * min(lp, 6)) * 0.62
            b.append(f"<circle class='pt' cx='{_n(left + j * cw + cw / 2)}' cy='{_n(y)}' "
                     f"r='{rad:.1f}' fill='{_diverging(z, vmax)}' stroke='{FG}' stroke-width='0.4'>"
                     f"<title>G{t + 1} {_esc(groups.names[t])} · {_esc(run.short(c))}\nbest FDR "
                     f"{_p(p)} · median z {z:+.2f} · {nt} term(s)</title></circle>")
    ly = HB - 30
    b.append(f"<text x='{left}' y='{ly - 14}' fill='{MUT}' font-size='9.5'>best-term FDR</text>")
    for k_, (v, lab) in enumerate(((2, "1e-2"), (4, "1e-4"), (6, "1e-6"))):
        x = left + 14 + k_ * 62
        b.append(f"<circle cx='{x}' cy='{ly}' r='{np.sqrt(30 + 70 * v) * 0.62:.1f}' fill='none' "
                 f"stroke='{FG}' stroke-width='0.6'/><text x='{x + 14}' y='{ly + 3.5}' fill='{MUT}' "
                 f"font-size='9.5'>{lab}</text>")
    b.append(zbar(left + 230, ly - 6, vmax, "rdcbgrpB",
                  "direction (median NES)" if groups.origin == "GSEA" else "direction (median z)",
                  width=140))
    b.append("</svg>")
    facts = {"n_groups": len(groups.names), "n_terms": len(groups.terms), "shown": len(keep),
             "top": [(groups.names[t], [run.short(c) for c in
                                         tab[tab.theme == t].sort_values("p").contrast.head(3)])
                     for t in keep[:5]]}
    return ["".join(a), "".join(b)], facts


# --------------------------------------------------------------------------- 5. volcanos

def fig_volcano(run, c, groups, prog=None, idx=1):
    from .labels import informative
    de = run.de
    lfc = de[c + "_diff"].astype(float)
    y_padj = run.use_padj and run.y_axis == "padj"
    p = de[c + ("_p.adj" if y_padj else "_p.val")].astype(float).clip(lower=1e-300)
    lp = -np.log10(p)
    ok = lfc.notna() & lp.notna()
    s = run.sig(c).fillna(False) & ok
    up, dn = s & (lfc > 0), s & (lfc < 0)
    grp = []
    if groups is not None:
        t = groups.table[groups.table.contrast == c].sort_values("p")
        used = set()
        for th in t.theme:
            genes = (groups.genes(th) & set(de.index[s])) - used
            if len(genes) >= 3:
                grp.append((int(th), _short(groups.names[th], 34), GROUP_COLS[len(grp)],
                            de.index.isin(list(genes))))
                used |= genes
            if len(grp) == 3:
                break
    W, H = 780, 640
    xm = float(np.nanmax(np.abs(lfc[ok].values))) if ok.any() else 1.0
    ymax = float(np.nanmax(lp[ok].values)) if ok.any() else 1.0
    ax = Ax(70, 60, 680, 440, (-xm * 1.15, xm * 1.15), (-0.02 * ymax, ymax * 1.08))
    sub = (f"{int(up.sum())} up · {int(dn.sum())} down · "
           + ("colour = functional group" if grp else "colour = direction")
           + " · names = most marked proteins · hover a coloured point · click the legend")
    sv = [svg_open(W, H, run.short(c)), title_block(14, 20, run.short(c), sub),
          ax.frame("log2 fold change", "−log10(p.adj)" if y_padj else "−log10(p)")]
    stat = de[c + ("_p.adj" if run.use_padj else "_p.val")].astype(float)
    cut = None
    if y_padj or not run.use_padj:
        cut = run.p_thr
    elif (stat < run.p_thr).any():
        hi = float(p[stat < run.p_thr].max())
        above = p[(stat >= run.p_thr) & (p > hi)]
        cut = float(np.sqrt(hi * above.min())) if len(above) else hi
    if cut is not None and ax.yl[0] < -np.log10(cut) < ax.yl[1]:
        y = float(ax.Y(-np.log10(cut)))
        sv.append(f"<line x1='{ax.x0}' y1='{_n(y)}' x2='{ax.x0 + ax.w}' y2='{_n(y)}' stroke='{MUT}' "
                  f"stroke-dasharray='4 3' stroke-width='0.8'/><text x='{ax.x0 + ax.w - 4}' "
                  f"y='{_n(y - 4)}' fill='{MUT}' font-size='9.5' text-anchor='end'>"
                  f"{'p.adj' if run.use_padj else 'p'} &lt; {run.p_thr:g}</text>")
    for v in (run.lfc_min, -run.lfc_min):
        if v:
            x = float(ax.X(v))
            sv.append(f"<line x1='{_n(x)}' y1='{ax.y0}' x2='{_n(x)}' y2='{ax.y0 + ax.h}' "
                      f"stroke='{MUT}' stroke-dasharray='4 3' stroke-width='0.8'/>")
    # not significant: a single light path, no tooltip
    ns = ok & ~s
    X, Y = ax.X(lfc[ns].values), ax.Y(lp[ns].values)
    path = "".join(f"M{x:.1f} {y:.1f}h0" for x, y in zip(X, Y))
    sv.append(f"<g data-g='ns'><path d='{path}' stroke='{FAINT}' stroke-width='3.4' "
              f"stroke-linecap='round' fill='none'/></g>")
    prog_of = _prog_of(prog)
    grp_of = {}
    for th, name, col, m in grp:
        for pid in de.index[m]:
            grp_of[pid] = (th, name)
    padj = de[c + "_p.adj"] if c + "_p.adj" in de.columns else None
    pval = de[c + "_p.val"] if c + "_p.val" in de.columns else None

    def pt(pid, col, r, key):
        nm = str(run.labels.get(pid, pid))
        tip = [f"{nm}" + (f" ({pid})" if nm != pid else ""),
               f"log2FC {lfc[pid]:+.2f}"
               + (f" · p {_p(pval[pid])}" if pval is not None else "")
               + (f" · p.adj {_p(padj[pid])}" if padj is not None else "")]
        if pid in prog_of:
            tip.append(f"programme P{int(prog_of[pid])}")
        if pid in grp_of:
            tip.append(f"group G{grp_of[pid][0] + 1}: {grp_of[pid][1]}")
        return (f"<circle class='pt' cx='{_n(ax.X(lfc[pid]))}' cy='{_n(ax.Y(lp[pid]))}' r='{r}' "
                f"fill='{col}' stroke='{BG}' stroke-width='0.6'><title>{_esc(chr(10).join(tip))}"
                f"</title></circle>")

    legend = [(FAINT, f"not significant ({int(ns.sum())})", "ns")]
    if grp:
        other = s.copy()
        for _, _, _, m in grp:
            other &= ~m
        sv.append("<g data-g='other'>" + "".join(pt(i, MUT, 2.6, "other")
                                                 for i in de.index[other]) + "</g>")
        legend.append((MUT, f"other significant ({int(other.sum())})", "other"))
        for th, name, col, m in grp:
            sv.append(f"<g data-g='g{th}'>" + "".join(pt(i, col, 4, f"g{th}")
                                                     for i in de.index[m]) + "</g>")
            legend.append((col, f"{name} ({int(m.sum())})", f"g{th}"))
    else:
        sv.append("<g data-g='up'>" + "".join(pt(i, UP, 3.2, "up") for i in de.index[up]) + "</g>")
        sv.append("<g data-g='down'>" + "".join(pt(i, DOWN, 3.2, "down") for i in de.index[dn]) + "</g>")
        legend += [(UP, f"up ({int(up.sum())})", "up"), (DOWN, f"down ({int(dn.sum())})", "down")]
    named = pd.Series([informative(run.labels.get(i, "")) for i in de.index], index=de.index)
    score = (lp * lfc.abs()).where(ok, 0)
    cand = list(score[s & named].sort_values(ascending=False).index[:10])
    for _, _, _, m in grp:
        cand += [i for i in score[m & named].sort_values(ascending=False).index[:3] if i not in cand]
    placer = Placer((ax.x0 + 2, ax.y0 + 2, ax.x0 + ax.w - 2, ax.y0 + ax.h - 2), fs=10)
    for i in cand:
        placer.block((float(ax.X(lfc[i])) - 4, float(ax.Y(lp[i])) - 4,
                      float(ax.X(lfc[i])) + 4, float(ax.Y(lp[i])) + 4))
    for i in cand:
        x, y = float(ax.X(lfc[i])), float(ax.Y(lp[i]))
        text = _short(run.labels.get(i, i), 16)
        pos = placer.place(x, y, text)
        if pos:
            sv.append(_label(x, y, pos, text))
    lx, ly = 70, ax.y0 + ax.h + 56
    for col, text, key in legend:
        w = 30 + _tw(text, 10.5)
        if lx + w > W - 10:
            lx, ly = 70, ly + 18
        sv.append(_legend_item(lx, ly, col, text, toggle=key))
        lx += w + 6
    sv.append("</svg>")
    up_names = [run.labels.get(i, i) for i in score[up & named].sort_values(ascending=False).index[:5]]
    dn_names = [run.labels.get(i, i) for i in score[dn & named].sort_values(ascending=False).index[:5]]
    facts = {"contrast": run.short(c), "up": int(up.sum()), "down": int(dn.sum()),
             "up_names": up_names, "down_names": dn_names, "groups": [g[1] for g in grp]}
    return ["".join(sv)], facts


# --------------------------------------------------------------------------- 6. spotlight

def spotlight_picks(run, prog, n=8):
    from .labels import informative
    if run.abundance is None or run.samples is None:
        return []
    lfc = run.de[[c + "_diff" for c in run.contrasts]].abs().max(axis=1)
    picks = []
    if prog is not None:
        for rnd in range(3):
            for p in range(1, prog.k + 1):
                members = prog.labels.index[prog.labels == p]
                cand = [i for i in lfc.loc[members].sort_values(ascending=False).index
                        if informative(run.labels.get(i, "")) and i in run.abundance.index
                        and i not in [x[0] for x in picks]]
                if len(cand) > rnd:
                    picks.append((cand[rnd], p))
                if len(picks) >= n:
                    break
            if len(picks) >= n:
                break
    if len(picks) < n:
        sig = np.zeros(len(run.de), bool)
        for c in run.contrasts:
            sig |= run.sig(c).fillna(False).values
        for i in lfc[sig].sort_values(ascending=False).index:
            if len(picks) >= n:
                break
            if i in run.abundance.index and informative(run.labels.get(i, "")) \
                    and i not in [x[0] for x in picks]:
                picks.append((i, None))
    return picks


def fig_spotlight(run, design, prog, n=8):
    picks = spotlight_picks(run, prog, n)
    if not picks:
        return None
    samples = run.samples.drop_duplicates("label").set_index("label")
    labs = [l for l in run.abundance.columns if l in samples.index]
    cond = samples.loc[labs, "condition"].values
    gcol = group_colors(design)
    lines = design.ordered and len(design.groups or ["all"]) <= 6
    many = (not lines) and len(design.conditions) > 12
    ncol = 2 if many else (4 if len(picks) > 4 else len(picks))
    nrow = int(np.ceil(len(picks) / ncol))
    pw = min(560, 120 + 22 * len(design.conditions)) if many else 210
    ph = 150
    gx, gy = 70, 100 if not lines else 76
    top = 64
    W = 30 + ncol * (pw + gx)
    H = top + nrow * (ph + gy) + 20
    hollow = " (hollow = not detected, DPC estimate)" if run.observed is not None else ""
    s = [svg_open(W, H, "Spotlight"),
         title_block(14, 20, "Spotlight: the proteins behind the programmes",
                     f"dots = replicates{hollow} · {'line' if lines else 'dash'} = mean · most "
                     "marked named protein of each programme · hover a dot")]
    prog_of = _prog_of(prog)
    for k, (pid, p) in enumerate(picks):
        r, ci = divmod(k, ncol)
        x0 = 60 + ci * (pw + gx)
        y0 = top + 36 + r * (ph + gy)
        v = run.abundance.loc[pid, labs].values.astype(float)
        seen = (run.observed.loc[pid, labs].values == 1) if (
            run.observed is not None and pid in run.observed.index) else np.isfinite(v)
        fin = np.isfinite(v)
        lo, hi = (np.nanmin(v), np.nanmax(v)) if fin.any() else (0, 1)
        padv = (hi - lo) * 0.12 or 0.5
        if lines:
            lv = design.levels
            ax = Ax(x0, y0, pw, ph, (-0.4, len(lv) - 0.6), (lo - padv, hi + padv))
            s.append(ax.frame("", "log2" if ci == 0 else "", xticks=list(range(len(lv))),
                              xtick_labels=[design.label(l) for l in lv]))
            gl = design.groups or ["all"]
            for gi, g in enumerate(gl):
                col = gcol.get(g, CAT[0])
                m = np.array([(design.group or {}).get(c, "all") == g for c in cond])
                xs = np.array([lv.index(design.level[c]) for c in cond[m]], float)
                off = (gi - (len(gl) - 1) / 2) * 0.10 if design.group else 0
                mean = pd.Series(v[m]).groupby(xs).mean()
                s.append(f"<polyline points='{' '.join(f'{_n(ax.X(a + off))},{_n(ax.Y(b))}' for a, b in mean.items() if np.isfinite(b))}' "
                         f"fill='none' stroke='{col}' stroke-width='1.8'/>")
                for xx, vv, ll, cc, ob in zip(xs, v[m], np.array(labs)[m], cond[m], seen[m]):
                    if not np.isfinite(vv):
                        continue
                    fill = col if ob else BG
                    s.append(f"<circle class='pt' cx='{_n(ax.X(xx + off))}' cy='{_n(ax.Y(vv))}' "
                             f"r='2.8' fill='{fill}' fill-opacity='{0.65 if ob else 1}' stroke='{col}' "
                             f"stroke-width='{0 if ob else 1}'><title>{_esc(ll)} · {_esc(cc)}\n"
                             f"log2 {vv:.2f}{'' if ob else ' · not detected (estimate)'}"
                             f"</title></circle>")
        else:
            conds = design.conditions
            ax = Ax(x0, y0, pw, ph, (-0.6, len(conds) - 0.4), (lo - padv, hi + padv))
            s.append(ax.frame("", "log2" if ci == 0 else "", xticks=[]))
            xs = np.array([conds.index(c) if c in conds else -1 for c in cond], float)
            for xx, vv, ll, cc, ob in zip(xs, v, labs, cond, seen):
                if xx < 0 or not np.isfinite(vv):
                    continue
                col = (gcol.get(design.group[cc], CAT[0]) if design.group and len(design.groups) <= 3
                       else CAT[0])
                fill = col if ob else BG
                s.append(f"<circle class='pt' cx='{_n(ax.X(xx))}' cy='{_n(ax.Y(vv))}' r='2.8' "
                         f"fill='{fill}' fill-opacity='{0.7 if ob else 1}' stroke='{col}' "
                         f"stroke-width='{0 if ob else 1}'><title>{_esc(ll)} · {_esc(cc)}\n"
                         f"log2 {vv:.2f}{'' if ob else ' · not detected (estimate)'}</title></circle>")
            mean = pd.Series(v).groupby(xs).mean()
            for xx, mm in mean.items():
                if xx < 0 or not np.isfinite(mm):
                    continue
                s.append(f"<line x1='{_n(ax.X(xx) - 7)}' y1='{_n(ax.Y(mm))}' x2='{_n(ax.X(xx) + 7)}' "
                         f"y2='{_n(ax.Y(mm))}' stroke='{WHITE}' stroke-width='1.6'>"
                         f"<title>{_esc(conds[int(xx)])}: mean {mm:.2f}</title></line>")
            for i, c in enumerate(conds):
                X = float(ax.X(i))
                yy = y0 + ph + 10
                s.append(f"<text x='{_n(X)}' y='{yy}' fill='{MUT}' font-size='8.5' text-anchor='end' "
                         f"transform='rotate(-60 {_n(X)} {yy})'>"
                         f"{_esc(_short(design.display(c), 16 if many else 12))}</text>")
        name = str(run.labels.get(pid, pid))
        desc = ""
        if "First.Protein.Description" in run.de.columns:
            desc = re.sub(r"\s*(OS=|\[).*$", "", str(run.de.at[pid, "First.Protein.Description"]))
        if not desc or desc.lower().startswith(name.rstrip("…").lower()[:12]) or desc == "nan":
            desc = str(pid)
        pp = prog_of.get(pid, p)
        s.append(f"<text x='{x0}' y='{y0 - 20}' fill='{WHITE}' font-size='12' font-style='italic' "
                 f"font-weight='600'><title>{_esc(name)} ({_esc(pid)})</title>{_esc(_short(name, 22))}</text>")
        s.append(f"<text x='{x0}' y='{y0 - 7}' fill='{MUT}' font-size='9'>"
                 f"{_esc(_short((f'P{int(pp)} · ' if pp else '') + desc, 44))}</text>")
    if lines and design.group:
        x = 14
        for g in design.groups:
            s.append(_legend_item(x, 56, gcol.get(g, CAT[0]), _short(g, 18)))
            x += 30 + _tw(_short(g, 18), 10.5)
    s.append("</svg>")
    return ["".join(s)], {"proteins": [run.labels.get(i, i) for i, _ in picks]}


# --------------------------------------------------------------------------- 7. factors

def fig_factors(summary, terms, F, src):
    """Stacked bars: median variance share of each factor / interaction, for all selected
    proteins and for each programme."""
    if summary is None or summary.empty:
        return None
    parts = list(terms) + ["higher-order", "replicates"]
    cols = {t: CAT[i % len(CAT)] for i, t in enumerate(terms)}
    cols["higher-order"] = "#6e7681"
    cols["replicates"] = FAINT
    labs = [f"{r['set']} ({int(r['n_proteins'])})" for _, r in summary.iterrows()]
    left = max(130, max(_tw(l, 10.5) for l in labs) + 20)
    W = left + 660
    leg, x, y = [], 14, 62
    for t in parts:                            # legend, wrapped
        w = 30 + _tw(t, 10.5)
        if x + w > W - 10:
            x, y = 14, y + 17
        leg.append(_legend_item(x, y, cols[t], t, shape="rect"))
        x += w
    rh, top = 26, y + 26
    H = top + rh * len(summary) + 56
    ax = Ax(left, top, 600, rh * len(summary), (0, 1.0), (0, len(summary)))
    from .factorial import levels
    lv = "; ".join(f"{c} ({', '.join(map(str, levels(F, c)))})" for c in F.columns)
    s = [svg_open(W, H, "What drives the response"),
         title_block(14, 20, "What drives the response: variance share per design factor",
                     f"Factors from {src}: {_short(lv, 110)} · median type-II share per "
                     "protein")]
    s.append(ax.frame("median share of the protein's variance", "", xticks=[0, .25, .5, .75, 1],
                      yticks=[], xfmt=lambda t: f"{t * 100:.0f} %"))
    s += leg
    for i, (_, rd) in enumerate(summary.iterrows()):
        yc = top + rh * i + rh / 2
        x0 = 0.0
        for t in parts:
            v = float(rd.get(t, 0) or 0) if t in summary.columns else 0.0
            if v <= 0:
                continue
            xa, xb = float(ax.X(x0)), float(ax.X(min(x0 + v, 1.05)))
            s.append(f"<rect class='pt' x='{_n(xa)}' y='{_n(yc - 8)}' width='{_n(max(xb - xa, 0.5))}' "
                     f"height='16' fill='{cols[t]}'><title>{_esc(rd['set'])} · {_esc(t)}: "
                     f"{v * 100:.0f} % (median share)</title></rect>")
            x0 += v
        s.append(f"<text x='{_n(left - 8)}' y='{_n(yc + 3.5)}' fill='{FG}' font-size='10.5' "
                 f"text-anchor='end'>{_esc(labs[i])}</text>")
    s.append("</svg>")
    return ["".join(s)], {}


# --------------------------------------------------------------------------- 8. factor effects

def _lv(f, v, always=False):
    """A level as text: a bare number gets its factor name ('time 3'), unless the factor
    is named after its levels ('4/14')."""
    v = str(v)
    if "/" in str(f):
        return v
    return f"{f} {v}" if always or re.fullmatch(r"-?\d+(?:\.\d+)?", v) else v


def fig_factor_effects(run, F, cf, blocks):
    """Comparisons arranged by the design factor they vary: one block per factor, rows =
    the context level with most values (genotype), columns = the other context levels
    (time), cell = up / down counts under the selection rule."""
    from .factorial import MISSING, levels
    if not blocks:
        return None
    counts = {}
    for c in cf:
        sg = run.sig(c).fillna(False)
        lf = run.de[c + "_diff"]
        counts[c] = (int((sg & (lf > 0)).sum()), int((sg & (lf < 0)).sum()))
    m = max([max(u, d) for u, d in counts.values()] + [1])
    cw, rh, left = 132, 22, 170
    parts, facts, y = [], [], 58
    width = left
    for fct, cs in blocks:
        if fct == "several":
            rows = [(c, None) for c in cs]
            colkeys = [None]
            title = f"Comparisons varying several factors ({len(cs)})"
            head = ""
        else:
            ctxf = [f for f in F.columns if f != fct]
            nval = {f: len({cf[c][1].get(f) for c in cs}) for f in ctxf}
            var = [f for f in ctxf if nval[f] > 1]
            rowf = max(var, key=lambda f: nval[f]) if var else None
            colf = [f for f in var if f != rowf]
            nums = {F.at[cf[c][2], fct] for c in cs}
            dens = {F.at[cf[c][3], fct] for c in cs}
            if rowf is None:              # nothing varies in the context: one row per level
                rows = [(c, f"{F.at[cf[c][2], fct]} vs {F.at[cf[c][3], fct]}") for c in cs]
                colkeys = [None]
            else:
                rk = levels(F, rowf)
                rows = sorted({cf[c][1][rowf] for c in cs}, key=lambda v: rk.index(v) if v in rk else 99)
                rows = [(None, v) for v in rows]
                colkeys = sorted({tuple(cf[c][1][f] for f in colf) for c in cs},
                                 key=lambda t: tuple(levels(F, f).index(v) if v in levels(F, f) else 99
                                                     for f, v in zip(colf, t))) if colf else [()]
            fixed = {f: next(iter({cf[c][1].get(f) for c in cs})) for f in ctxf if nval[f] == 1}
            fx = " · ".join(_lv(f, v, True) for f, v in fixed.items() if v != MISSING)
            title = f"{fct} varied — {len(cs)} comparison(s)" + (f" · at {fx}" if fx else "")
            head = (f"up = higher in {next(iter(nums))} than in {next(iter(dens))}"
                    if len(nums) == 1 and len(dens) == 1 else
                    "up = higher in the first condition of each comparison")
        parts.append(f"<text x='14' y='{y}' fill='{WHITE}' font-size='12' font-weight='600'>"
                     f"{_esc(title)}</text>")
        if head:
            parts.append(f"<text x='14' y='{y + 14}' fill='{MUT}' font-size='10'>{_esc(head)}</text>")
        y += 34
        if fct != "several" and rowf is not None:
            parts.append(f"<text x='{left - 10}' y='{y + (0 if colf else -4)}' fill='{MUT}' "
                         f"font-size='10' text-anchor='end'>{_esc(_short(rowf, 24))}</text>")
        if fct != "several" and rowf is not None and colf:
            for j, ck in enumerate(colkeys):
                lab = " · ".join(_lv(f, v) for f, v in zip(colf, ck) if v != MISSING) or "—"
                parts.append(f"<text x='{left + j * cw + cw / 2:.0f}' y='{y}' fill='{MUT}' font-size='10' "
                             f"text-anchor='middle'>{_esc(lab)}</text>")
            y += 8
        width = max(width, left + cw * len(colkeys) + 20)
        best = []
        for i, (c0, rlab) in enumerate(rows):
            yy = y + i * rh
            label = rlab if rlab is not None else run.short(c0)
            parts.append(f"<text x='{left - 10}' y='{yy + 15}' fill='{FG}' font-size='10.5' "
                         f"text-anchor='end'>{_esc(_short(label, 26))}</text>")
            for j, ck in enumerate(colkeys):
                if c0 is not None:
                    c = c0
                else:
                    c = next((x for x in cs if cf[x][1][rowf] == rlab
                              and tuple(cf[x][1][f] for f in colf) == tuple(ck)), None)
                x0 = left + j * cw
                parts.append(f"<rect x='{x0 + 2}' y='{yy + 2}' width='{cw - 4}' height='{rh - 4}' "
                             f"fill='{PANEL}' rx='3'/>")
                if c is None:
                    parts.append(f"<text x='{x0 + cw / 2:.0f}' y='{yy + 15}' fill='{FAINT}' "
                                 f"font-size='10' text-anchor='middle'>not compared</text>")
                    continue
                u, d = counts[c]
                best.append((u + d, c, u, d))
                xc, half = x0 + cw / 2, (cw - 72) / 2
                wu, wd = half * u / m, half * d / m
                parts.append(f"<g><title>{_esc(run.short(c))}: {u} up, {d} down</title>"
                             f"<rect x='{x0 + 2}' y='{yy + 2}' width='{cw - 4}' height='{rh - 4}' "
                             f"fill='transparent'/>"
                             f"<rect x='{xc:.1f}' y='{yy + 6}' width='{wu:.1f}' height='{rh - 12}' fill='{UP}'/>"
                             f"<rect x='{xc - wd:.1f}' y='{yy + 6}' width='{wd:.1f}' height='{rh - 12}' fill='{DOWN}'/>"
                             f"<text x='{x0 + 8}' y='{yy + 15}' fill='{FG}' font-size='9.5'>{d}</text>"
                             f"<text x='{x0 + cw - 8}' y='{yy + 15}' fill='{FG}' font-size='9.5' "
                             f"text-anchor='end'>{u}</text></g>")
        y += rh * len(rows) + 22
        if best:
            best.sort()
            tot = [b[0] for b in best]
            facts.append({"factor": fct, "n": len(best), "min": best[0], "max": best[-1],
                          "median": float(np.median(tot))})
    W = max(width, 780)
    H = y + 10
    out = [svg_open(W, H, "Effects of each design factor"),
           title_block(14, 20, "Effects of each design factor",
                       "Comparisons grouped by the factor they vary · bars = down | up differential "
                       f"proteins ({run.rule_text()}) · shared scale · hover a cell")]
    out += parts
    out.append("</svg>")
    return ["".join(out)], facts


# --------------------------------------------------------------------------- 9. programme profiles

def fig_factor_profiles(run, prog, F, r, cf=None):
    """Interaction plots: one row per programme, one column per facet level (light),
    x = ordered factor (time), one line per level of the lines factor (genotype);
    y = median z of the programme. A level of x seen in a single facet (T0 in the light
    only) is drawn in every facet as the shared start (hollow mark)."""
    from .factorial import MISSING, levels
    if prog is None or not r.get("lines"):
        return None
    x, lines = r["x"], r["lines"]
    facets = r["facets"][:1]
    lv_x = [v for v in levels(F, x) if v != MISSING]
    lv_l = [v for v in levels(F, lines) if v != MISSING]
    lv_f = [v for v in levels(F, facets[0]) if v != MISSING] if facets else [None]
    if len(lv_x) < 2:
        return None
    conds = [c for c in F.index if c in prog.z.columns]
    Fc = F.loc[conds]
    # shared x levels: present in only one facet (or with the facet absent)
    shared = set()
    if facets:
        for v in lv_x:
            fv = set(Fc.loc[Fc[x] == v, facets[0]])
            if len(fv - {MISSING}) <= 1:
                shared.add(v)
    # reference level of the lines factor: the usual denominator of its comparisons
    ref = None
    if cf:
        dens = [F.at[den, lines] for c, (diff, ctx, num, den) in cf.items() if diff == [lines]]
        if dens:
            ref = max(set(dens), key=dens.count)
    colf = {l: CAT[i % len(CAT)] for i, l in enumerate(lv_l)}
    med = {}
    for p in range(1, prog.k + 1):
        Z = prog.z[prog.labels == p]
        med[p] = (Z[conds].median(), len(Z))
    allv = np.concatenate([m[0].values for m in med.values()])
    lo, hi = float(np.nanmin(allv)), float(np.nanmax(allv))
    pad = (hi - lo) * 0.1 or 0.5
    pw, ph, gx, gy = 220, 130, 26, 44
    left, top = 190, 92
    ncol = len(lv_f)
    W = max(left + ncol * (pw + gx) + 10, 780)
    H = top + prog.k * (ph + gy) + 10

    def cond_of(xv, lv, fv):
        sel = (Fc[x] == xv) & (Fc[lines] == lv)
        if facets:
            sel &= (Fc[facets[0]] == fv) if xv not in shared else True
        hit = Fc.index[sel]
        return hit[0] if len(hit) else None

    s = [svg_open(W, H, "Programmes by design factor"),
         title_block(14, 20, "Programmes by design factor",
                     f"median z of each programme · x = {x} · one line per {lines}"
                     + (f" · panels = {facets[0]}" if facets else "")
                     + (f" · {ref} in white" if ref else "")
                     + (" · hollow = shared start" if shared else "") + " · click the legend")]
    xl, yl = 14, 60
    for i, l in enumerate(lv_l):
        colr = WHITE if l == ref else colf[l]
        w = 30 + _tw(_short(l, 14), 10.5)
        if xl + w > W - 10:
            xl, yl = 14, yl + 16
        s.append(_legend_item(xl, yl, colr, _short(l, 14), toggle=f"l{i}"))
        xl += w
    for j, fv in enumerate(lv_f):
        if fv is not None:
            s.append(f"<text x='{left + j * (pw + gx) + pw / 2:.0f}' y='{top - 8}' fill='{FG}' "
                     f"font-size='11' font-weight='600' text-anchor='middle'>{_esc(fv)}</text>")
    for p in range(1, prog.k + 1):
        m, n = med[p]
        y0 = top + (p - 1) * (ph + gy)
        st = f" · J={prog.stability[p]:.2f}" if prog.stability else ""
        s.append(f"<text x='{left - 44}' y='{y0 + ph / 2:.0f}' fill='{WHITE}' font-size='12' "
                 f"font-weight='600' text-anchor='end'>P{p}</text>")
        s.append(f"<text x='{left - 44}' y='{y0 + ph / 2 + 14:.0f}' fill='{MUT}' font-size='9.5' "
                 f"text-anchor='end'>{n} proteins{st}</text>")
        for j, fv in enumerate(lv_f):
            x0 = left + j * (pw + gx)
            ax = Ax(x0, y0, pw, ph, (-0.3, len(lv_x) - 0.7), (lo - pad, hi + pad))
            s.append(ax.frame("", "", xticks=list(range(len(lv_x))),
                              xtick_labels=[_lv(x, v) for v in lv_x],
                              yticks=nice_ticks(lo - pad, hi + pad, 3)))
            s.append(f"<line x1='{x0}' y1='{_n(ax.Y(0))}' x2='{x0 + pw}' y2='{_n(ax.Y(0))}' "
                     f"stroke='{MUT}' stroke-width='0.6' stroke-dasharray='2 3'/>")
            order = [l for l in lv_l if l != ref] + ([ref] if ref in lv_l else [])
            for l in order:
                i = lv_l.index(l)
                colr = WHITE if l == ref else colf[l]
                pts = []
                for k_, xv in enumerate(lv_x):
                    c = cond_of(xv, l, fv)
                    if c is not None and np.isfinite(m.get(c, np.nan)):
                        pts.append((k_, float(m[c]), c, xv in shared))
                if not pts:
                    continue
                g = [f"<g data-g='l{i}'>"]
                if len(pts) > 1:
                    g.append(f"<polyline points='{' '.join(f'{_n(ax.X(a))},{_n(ax.Y(b))}' for a, b, _, _ in pts)}' "
                             f"fill='none' stroke='{colr}' stroke-width='{2.4 if l == ref else 1.5}' "
                             f"stroke-opacity='0.9'/>")
                for a, b, c, sh in pts:
                    g.append(f"<circle class='pt' cx='{_n(ax.X(a))}' cy='{_n(ax.Y(b))}' r='3' "
                             f"fill='{BG if sh else colr}' stroke='{colr}' stroke-width='1.2'>"
                             f"<title>P{p} · {_esc(c)}\nmedian z {b:+.2f} ({n} proteins)"
                             f"{' · shared start' if sh else ''}</title></circle>")
                g.append("</g>")
                s.append("".join(g))
    s.append("</svg>")
    # facts: per programme, the line / facet / x combination at its highest and lowest
    facts = {}
    for p, (m, n) in med.items():
        mm = m.dropna()
        if len(mm):
            facts[p] = (mm.idxmax(), float(mm.max()), mm.idxmin(), float(mm.min()))
    return ["".join(s)], {"extremes": facts, "ref": ref, "shared": sorted(shared)}
