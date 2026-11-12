"""
Self-contained Pathfinder HTML report.

Fully offline rendering: figures are generated as SVG on the Python side (no CDN,
no charting JS). Palette consistent with the pipeline dashboard.
"""
from __future__ import annotations
import html
import json
import numpy as np
import pandas as pd

BG = "#0d1117"; PANEL = "#161b22"; GRID = "#30363d"; FG = "#c9d1d9"; MUT = "#8b949e"
BLUE = "#388bfd"; UP = "#f85149"; DOWN = "#388bfd"; GREEN = "#3fb950"; PURP = "#bc8cff"; AMBER = "#d29922"
COND_PALETTE = ["#388bfd", "#3fb950", "#bc8cff", "#d29922", "#f85149", "#39c5cf", "#db61a2", "#a5d6ff",
                "#7ee787", "#ffa657", "#d2a8ff", "#f2cc60", "#79c0ff", "#ff7b72", "#56d364", "#e3b341"]
# Day / night theme (button in the report header). The night palette above stays the
# source of every figure; the day theme is applied by the browser:
#   - THEME_VARS: page colours as CSS variables (night = the constants above);
#   - LIGHT_SWAP: neutral colours swapped on every figure element (fill and stroke);
#   - LIGHT_TEXT: light hues darkened on text only (contrast on white); marks keep them.
# The same maps recolour an exported SVG / PNG when the day theme is on.
THEME_VARS = {
    "dark": {"bg": BG, "panel": PANEL, "grid": GRID, "fg": FG, "mut": MUT, "hi": "#ffffff",
             "stripe": "#11161d", "blue": BLUE, "amber": AMBER, "green": GREEN, "purp": PURP,
             "red": UP},
    "light": {"bg": "#ffffff", "panel": "#f6f8fa", "grid": "#d0d7de", "fg": "#1f2328", "mut": "#59636e",
              "hi": "#1f2328", "stripe": "#f6f8fa", "blue": "#0969da", "amber": "#9a6700",
              "green": "#1a7f37", "purp": "#8250df", "red": "#cf222e"},
}
LIGHT_SWAP = {
    "#0d1117": "#ffffff",   # BG
    "#161b22": "#f6f8fa",   # PANEL
    "#1c2128": "#eaeef2",   # not significant / missing cells
    "#30363d": "#d0d7de",   # GRID
    "#3d444d": "#c8d1da",   # faint guides
    "#c9d1d9": "#1f2328",   # FG
    "#8b949e": "#59636e",   # MUT
    "#6e7681": "#6e7781",
    "#57606a": "#8c959f",   # edges of the functional-group network
    "#e6edf3": "#8c959f",   # diagonal of the correlation heatmaps
    "#ffffff": "#1f2328",   # white text and lines (titles, means, reference level)
    "#fff": "#1f2328",
}
LIGHT_TEXT = {
    "#388bfd": "#0969da", "#f85149": "#cf222e", "#3fb950": "#1a7f37", "#bc8cff": "#8250df",
    "#d29922": "#9a6700", "#39c5cf": "#1b7c83", "#db61a2": "#bf3989", "#a5d6ff": "#0550ae",
    "#7ee787": "#1a7f37", "#ffa657": "#bc4c00", "#d2a8ff": "#8250df", "#f2cc60": "#9a6700",
    "#79c0ff": "#0969da", "#ff7b72": "#cf222e", "#56d364": "#1a7f37", "#e3b341": "#9a6700",
}


def theme_css() -> str:
    """CSS of the day / night theme: variables and figure recolouring."""
    def _vars(d):
        return ";".join(f"--{k}:{v}" for k, v in d.items())
    out = [f":root,html[data-theme=dark]{{{_vars(THEME_VARS['dark'])}}}",
           f"html[data-theme=light]{{{_vars(THEME_VARS['light'])}}}"]
    L = "html[data-theme=light] main svg"
    for a, b in LIGHT_SWAP.items():
        out.append(f"{L} [fill='{a}' i]{{fill:{b}}} {L} [stroke='{a}' i]{{stroke:{b}}}")
    for a, b in LIGHT_TEXT.items():
        out.append(f"{L} text[fill='{a}' i],{L} tspan[fill='{a}' i]{{fill:{b}}}")
    out.append("html[data-theme=light] #go_cv{background:var(--panel)!important}")
    return "\n".join(out)


THEME_HEAD_JS = """(function(){try{var t=localStorage.getItem('pathfinder-theme');
if(t==='light'||t==='dark')document.documentElement.setAttribute('data-theme',t);}catch(e){}})();"""

THEME_JS = """
(function(){
  function cur(){ return document.documentElement.getAttribute('data-theme') === 'light' ? 'light' : 'dark'; }
  function label(){ var b = document.getElementById('themeBtn'); if(!b) return;
    var day = cur() === 'light';
    b.textContent = day ? '\\u263E Night' : '\\u2600 Day';
    b.title = day ? 'Switch to the night theme' : 'Switch to the day theme'; }
  window.toggleTheme = function(){
    var t = cur() === 'light' ? 'dark' : 'light';
    document.documentElement.setAttribute('data-theme', t);
    try{ localStorage.setItem('pathfinder-theme', t); }catch(e){}
    label();
  };
  if(document.readyState === 'loading') document.addEventListener('DOMContentLoaded', label); else label();
})();
"""
SRC_ORDER = ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"]
SRC_LABEL = {"GO:BP": "GO:BP — Biological process", "GO:MF": "GO:MF — Molecular function",
             "GO:CC": "GO:CC — Cellular component", "KEGG": "KEGG", "REAC": "Reactome"}


def _sources_present(meta):
    present = [s for s in SRC_ORDER if (meta["source"] == s).any()]
    others = [s for s in meta["source"].dropna().unique() if s not in SRC_ORDER]
    return present + sorted(others)


def _esc(s):
    return html.escape(str(s))


def _diverging(v, vmax):
    """signed z -> colour (red=up, blue=down), NaN -> dark grey."""
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "#1c2128"
    t = max(-1.0, min(1.0, v / vmax if vmax else 0.0))
    if t >= 0:
        r, g, b = 248, int(81 + (1 - t) * 120), int(73 + (1 - t) * 120)
    else:
        r, g, b = int(56 + (1 + t) * 120), int(139 + (1 + t) * 60), 253
    return f"rgb({r},{g},{b})"


# --------------------------------------------------------------------------- #
# Figures SVG                                                                   #
# --------------------------------------------------------------------------- #
def _colorbar(x, y, vmax, width=170, height=9, uid="cb"):
    """Graduated z-score colour bar (blue -> grey -> red) + NaN swatch."""
    s = []
    stops = []
    n = 24
    for i in range(n + 1):
        t = i / n
        v = -vmax + t * 2 * vmax
        stops.append(f"<stop offset='{t*100:.0f}%' stop-color='{_diverging(v, vmax)}'/>")
    s.append(f"<defs><linearGradient id='{uid}' x1='0' y1='0' x2='1' y2='0'>{''.join(stops)}</linearGradient></defs>")
    s.append(f"<rect x='{x}' y='{y}' width='{width}' height='{height}' fill='url(#{uid})' rx='2'/>")
    s.append(f"<rect x='{x}' y='{y}' width='{width}' height='{height}' fill='none' stroke='{GRID}' rx='2'/>")
    for i, lab in enumerate([-vmax, -vmax / 2, 0, vmax / 2, vmax]):
        lx = x + width * i / 4
        s.append(f"<line x1='{lx}' y1='{y+height}' x2='{lx}' y2='{y+height+3}' stroke='{MUT}'/>")
        s.append(f"<text x='{lx}' y='{y+height+11}' fill='{MUT}' font-size='8' text-anchor='middle'>{lab:+.1f}</text>")
    s.append(f"<text x='{x-6}' y='{y+height-1}' fill='{MUT}' font-size='8' text-anchor='end'>z</text>")
    # NaN swatch
    nx = x + width + 16
    s.append(f"<rect x='{nx}' y='{y}' width='{height}' height='{height}' fill='#1c2128' stroke='{GRID}' rx='2'/>")
    s.append(f"<text x='{nx+height+4}' y='{y+height-1}' fill='{MUT}' font-size='8'>not significant</text>")
    return "".join(s)


def svg_meta_heatmap(matrix: pd.DataFrame, meta: pd.DataFrame, max_terms=40, uid="cb"):
    if matrix.empty:
        return "<p class='mut'>No significant GO term.</p>"
    m = matrix.reindex(meta.sort_values("n_sig", ascending=False).index).head(max_terms)
    names = [meta.loc[k, "term_name"] for k in m.index]
    contrasts = list(m.columns)
    vmax = np.nanmax(np.abs(m.values)) or 1.0
    cw, ch = 54, 24
    left, top = 360, 150
    W = left + len(contrasts) * cw + 30
    H = top + len(m) * ch + 42
    s = [f"<svg viewBox='0 0 {W} {H}' width='{W}' height='{H}' style='max-width:100%;height:auto'>"]
    for j, c in enumerate(contrasts):
        x = left + j * cw + cw / 2
        s.append(f"<text x='{x}' y='{top-8}' fill='{MUT}' font-size='9' "
                 f"text-anchor='start' transform='rotate(-90 {x} {top-8})'>{_esc(c)}</text>")
    for i, k in enumerate(m.index):
        y = top + i * ch
        label = names[i][:52]
        s.append(f"<text x='{left-8}' y='{y+ch*0.7}' fill='{FG}' font-size='10' text-anchor='end'>{_esc(label)}</text>")
        for j, c in enumerate(contrasts):
            v = m.iloc[i, j]
            col = _diverging(v if pd.notna(v) else None, vmax)
            x = left + j * cw
            tip = "" if pd.isna(v) else f"<title>{_esc(names[i])} — {_esc(c)} : z={v:.2f}</title>"
            s.append(f"<rect x='{x}' y='{y}' width='{cw-2}' height='{ch-2}' fill='{col}' rx='2'>{tip}</rect>")
    # graduated colour bar
    s.append(_colorbar(left, H - 28, vmax, uid=uid))
    s.append("</svg>")
    return "".join(s)


def svg_scatter(df, xcol, ycol, cat_col, title_x, title_y, size=(680, 540)):
    W, Hh = size
    xs = df[xcol].astype(float).values
    ys = df[ycol].astype(float).values
    cats = list(dict.fromkeys(df[cat_col].tolist())) if cat_col else [None]
    cmap = {c: COND_PALETTE[i % len(COND_PALETTE)] for i, c in enumerate(cats)}

    # legend gutter: the plot shrinks so the legend never overlaps it
    if cat_col:
        maxlab = max((len(str(c)) for c in cats[:14]), default=0)
        gutter = min(190, 26 + int(maxlab * 5.6))
        W = W + gutter
    else:
        gutter = 0

    pad = 46
    right = pad + gutter
    xmin, xmax = xs.min(), xs.max(); ymin, ymax = ys.min(), ys.max()
    xr = (xmax - xmin) or 1; yr = (ymax - ymin) or 1
    xmin -= 0.08 * xr; xmax += 0.08 * xr; ymin -= 0.08 * yr; ymax += 0.08 * yr

    def px(x): return pad + (x - xmin) / (xmax - xmin) * (W - pad - right)
    def py(y): return Hh - pad - (y - ymin) / (ymax - ymin) * (Hh - 2 * pad)

    s = [f"<svg viewBox='0 0 {W} {Hh}' width='{W}' height='{Hh}' style='max-width:100%;height:auto'>"]
    s.append(f"<rect x='0' y='0' width='{W}' height='{Hh}' fill='{PANEL}' rx='6'/>")
    # axes 0
    if xmin < 0 < xmax:
        s.append(f"<line x1='{px(0)}' y1='{pad}' x2='{px(0)}' y2='{Hh-pad}' stroke='{GRID}' stroke-dasharray='3 3'/>")
    if ymin < 0 < ymax:
        s.append(f"<line x1='{pad}' y1='{py(0)}' x2='{W-right}' y2='{py(0)}' stroke='{GRID}' stroke-dasharray='3 3'/>")
    for _, r in df.iterrows():
        col = cmap[r[cat_col]] if cat_col else BLUE
        s.append(f"<circle cx='{px(float(r[xcol])):.1f}' cy='{py(float(r[ycol])):.1f}' r='5' "
                 f"fill='{col}' fill-opacity='0.85' stroke='{BG}' stroke-width='0.5'>"
                 f"<title>{_esc(r.name if cat_col else '')} {_esc(r[cat_col]) if cat_col else ''}</title></circle>")
    s.append(f"<text x='{(pad + (W-right))/2:.0f}' y='{Hh-12}' fill='{MUT}' font-size='11' text-anchor='middle'>{_esc(title_x)}</text>")
    s.append(f"<text x='14' y='{Hh/2}' fill='{MUT}' font-size='11' text-anchor='middle' transform='rotate(-90 14 {Hh/2})'>{_esc(title_y)}</text>")
    # legend inside the gutter
    if cat_col:
        lx = W - gutter + 8
        for i, c in enumerate(cats[:14]):
            yy = pad + i * 16
            lab = (str(c)[:22] + "…") if len(str(c)) > 23 else str(c)
            s.append(f"<rect x='{lx}' y='{yy}' width='10' height='10' fill='{cmap[c]}' rx='2'/>")
            s.append(f"<text x='{lx+15}' y='{yy+9}' fill='{FG}' font-size='10'>{_esc(lab)}</text>")
        if len(cats) > 14:
            s.append(f"<text x='{lx}' y='{pad+14*16+9}' fill='{MUT}' font-size='9'>+{len(cats)-14} autres</text>")
    s.append("</svg>")
    return "".join(s)


def svg_scatter_value(df, xcol, ycol, valcol, diverging, size=(720, 560), r=3):
    W, Hh = size
    xs = df[xcol].astype(float).values; ys = df[ycol].astype(float).values
    vals = pd.to_numeric(df[valcol], errors="coerce").values
    pad = 40
    xmin, xmax = np.nanmin(xs), np.nanmax(xs); ymin, ymax = np.nanmin(ys), np.nanmax(ys)
    xr = (xmax - xmin) or 1; yr = (ymax - ymin) or 1
    def px(x): return pad + (x - xmin) / xr * (W - 2 * pad)
    def py(y): return Hh - pad - (y - ymin) / yr * (Hh - 2 * pad)
    finite = vals[np.isfinite(vals)]
    if diverging:
        vmax = np.nanmax(np.abs(finite)) if finite.size else 1.0
        def color(v):
            return _diverging(v if np.isfinite(v) else None, vmax or 1.0)
    else:
        lo, hi = (np.nanpercentile(finite, 2), np.nanpercentile(finite, 98)) if finite.size else (0, 1)
        rng = (hi - lo) or 1
        def color(v):
            if not np.isfinite(v):
                return "#30363d"
            t = max(0.0, min(1.0, (v - lo) / rng))
            r_, g_, b_ = int(30 + t * 30), int(60 + t * 150), int(120 + t * 80)
            return f"rgb({r_},{g_},{b_})"
    s = [f"<svg viewBox='0 0 {W} {Hh}' width='{W}' height='{Hh}' style='max-width:100%;height:auto'>",
         f"<rect x='0' y='0' width='{W}' height='{Hh}' fill='{PANEL}' rx='6'/>"]
    for x, y, v in zip(xs, ys, vals):
        s.append(f"<circle cx='{px(x):.1f}' cy='{py(y):.1f}' r='{r}' fill='{color(v)}' fill-opacity='0.8'/>")
    # legend bar
    s.append(f"<text x='{W/2}' y='{Hh-10}' fill='{MUT}' font-size='11' text-anchor='middle'>{_esc(xcol)} / {_esc(ycol)}</text>")
    s.append("</svg>")
    return "".join(s)


def svg_scatter_colored(df, xcol, ycol, colors, legend, size=(560, 440), r=3):
    W, Hh = size
    xs = df[xcol].astype(float).values; ys = df[ycol].astype(float).values
    pad = 40
    xmin, xmax = np.nanmin(xs), np.nanmax(xs); ymin, ymax = np.nanmin(ys), np.nanmax(ys)
    xr = (xmax - xmin) or 1; yr = (ymax - ymin) or 1
    def px(x): return pad + (x - xmin) / xr * (W - 2 * pad)
    def py(y): return Hh - pad - (y - ymin) / yr * (Hh - 2 * pad)
    s = [f"<svg viewBox='0 0 {W} {Hh}' width='{W}' height='{Hh}' style='max-width:100%;height:auto'>",
         f"<rect x='0' y='0' width='{W}' height='{Hh}' fill='{PANEL}' rx='6'/>"]
    # grey points first, coloured on top
    order = np.argsort([0 if c == "#30363d" else 1 for c in colors])
    for i in order:
        s.append(f"<circle cx='{px(xs[i]):.1f}' cy='{py(ys[i]):.1f}' r='{r if colors[i]!='#30363d' else 2}' "
                 f"fill='{colors[i]}' fill-opacity='{0.9 if colors[i]!='#30363d' else 0.35}'/>")
    for i, (lab, col) in enumerate(legend[:10]):
        yy = pad + i * 15
        s.append(f"<rect x='{W-190}' y='{yy}' width='10' height='10' fill='{col}' rx='2'/>")
        lab2 = (str(lab)[:26] + "…") if len(str(lab)) > 27 else str(lab)
        s.append(f"<text x='{W-176}' y='{yy+9}' fill='{FG}' font-size='9'>{_esc(lab2)}</text>")
    s.append("</svg>")
    return "".join(s)


def _go_umap_widget(udf, meta, term_genes, max_terms=24):
    """Interactive widget: protein UMAP recoloured by GO term, with source and term-count
    selectors. Data embedded, rendered client-side on a canvas."""
    import json as _json

    def nrm(x):
        return re.sub(r"[^a-z0-9]", "", str(x).lower())

    m = meta.copy()
    if "is_representative" in m.columns:
        m = m[m["is_representative"]]
    m = m.sort_values("n_sig", ascending=False)

    terms, term_sets = [], []
    for key, row in m.iterrows():
        accs = {nrm(a) for a in term_genes.get(key, set())}
        if not accs:
            continue
        terms.append({"n": str(row["term_name"]), "s": str(row["source"]), "k": int(row["n_sig"])})
        term_sets.append(accs)
        if len(terms) >= max_terms:
            break
    if not terms:
        return ""

    # tokens per protein + membership (term indices)
    coords, memb = [], []
    for _, r in udf.iterrows():
        toks = {nrm(p) for p in re.split(r"[;,|]", str(r.get("pid", ""))) if p.strip()}
        toks.add(nrm(r.get("label", "")))
        coords.append([round(float(r["UMAP1"]), 2), round(float(r["UMAP2"]), 2)])
        memb.append([i for i, accs in enumerate(term_sets) if toks & accs])

    sources = [s for s in SRC_ORDER if any(t["s"] == s for t in terms)]
    data = {"coords": coords, "memb": memb, "terms": terms,
            "sources": sources, "palette": COND_PALETTE}

    src_opts = "".join(f"<option value='{_esc(s)}'>{_esc(SRC_LABEL.get(s, s))}</option>" for s in sources)
    tmpl = """
<div class='goctl'>
  <label>Source <select id='go_src'><option value='ALL'>All</option>__SRC__</select></label>
  <label>Terms <select id='go_k'><option>4</option><option selected>8</option><option>12</option><option>16</option><option>24</option></select></label>
  <span id='go_info' class='mut'></span>
</div>
<canvas id='go_cv' width='1520' height='1160' style='width:100%;max-width:560px;height:auto;background:#161b22;border-radius:6px'></canvas>
<div id='go_leg' class='golegend'></div>
<script>
(function(){
 var D=__DATA__, C=D.coords, M=D.memb, T=D.terms, PAL=D.palette;
 var cv=document.getElementById('go_cv'), ctx=cv.getContext('2d');
 var xs=C.map(function(p){return p[0];}), ys=C.map(function(p){return p[1];});
 var xmin=Math.min.apply(null,xs), xmax=Math.max.apply(null,xs);
 var ymin=Math.min.apply(null,ys), ymax=Math.max.apply(null,ys);
 var W=760,H=580,pad=34; ctx.scale(2,2);
 function px(x){return pad+(x-xmin)/((xmax-xmin)||1)*(W-2*pad);}
 function py(y){return H-pad-(y-ymin)/((ymax-ymin)||1)*(H-2*pad);}
 function selectedTerms(){
   var src=document.getElementById('go_src').value;
   var k=parseInt(document.getElementById('go_k').value,10);
   var idx=[]; for(var i=0;i<T.length;i++){ if(src==='ALL'||T[i].s===src) idx.push(i); if(idx.length>=k) break; }
   return idx;
 }
 function draw(){
   var sel=selectedTerms();
   var pos={}; sel.forEach(function(ti,j){pos[ti]=j;});
   ctx.clearRect(0,0,W,H);
   // gris d'abord
   for(var i=0;i<C.length;i++){ var hit=null; for(var m=0;m<M[i].length;m++){ if(pos.hasOwnProperty(M[i][m])){ if(hit===null||pos[M[i][m]]<pos[hit]) hit=M[i][m]; } }
     if(hit===null){ ctx.fillStyle='rgba(120,130,140,0.28)'; ctx.beginPath(); ctx.arc(px(C[i][0]),py(C[i][1]),1.6,0,6.283); ctx.fill(); } }
   var ncol=0;
   for(var i=0;i<C.length;i++){ var hit=null; for(var m=0;m<M[i].length;m++){ if(pos.hasOwnProperty(M[i][m])){ if(hit===null||pos[M[i][m]]<pos[hit]) hit=M[i][m]; } }
     if(hit!==null){ ctx.fillStyle=PAL[pos[hit]%PAL.length]; ctx.beginPath(); ctx.arc(px(C[i][0]),py(C[i][1]),2.6,0,6.283); ctx.fill(); ncol++; } }
   var leg=document.getElementById('go_leg'); leg.innerHTML='';
   sel.forEach(function(ti,j){ var d=document.createElement('div'); d.className='legit';
     d.innerHTML="<span class='sw' style='background:"+PAL[j%PAL.length]+"'></span>"+T[ti].n+" ["+T[ti].s+"] · "+T[ti].k+" contrasts"; leg.appendChild(d); });
   document.getElementById('go_info').textContent=ncol+' / '+C.length+' proteins coloured';
 }
 document.getElementById('go_src').addEventListener('change',draw);
 document.getElementById('go_k').addEventListener('change',draw);
 draw();
})();
</script>"""
    return tmpl.replace("__SRC__", src_opts).replace("__DATA__", _json.dumps(data))


import re  # noqa: E402 (used by the GO widget)


def svg_scree(var_exp):
    ve = [float(v) for v in var_exp]
    W, Hh = 660, 300; pad = 44
    n = len(ve); bw = (W - 2 * pad) / max(n, 1)
    vmax = max(ve) or 1
    s = [f"<svg viewBox='0 0 {W} {Hh}' width='{W}' height='{Hh}' style='max-width:100%;height:auto'>",
         f"<rect x='0' y='0' width='{W}' height='{Hh}' fill='{PANEL}' rx='6'/>"]
    cum = 0
    for i, v in enumerate(ve):
        x = pad + i * bw
        h = (v / vmax) * (Hh - 2 * pad)
        s.append(f"<rect x='{x+4}' y='{Hh-pad-h:.1f}' width='{bw-8:.1f}' height='{h:.1f}' fill='{BLUE}' rx='2'>"
                 f"<title>PC{i+1} : {v:.1f}%</title></rect>")
        s.append(f"<text x='{x+bw/2:.1f}' y='{Hh-pad+14}' fill='{MUT}' font-size='9' text-anchor='middle'>PC{i+1}</text>")
        s.append(f"<text x='{x+bw/2:.1f}' y='{Hh-pad-h-4:.1f}' fill='{FG}' font-size='9' text-anchor='middle'>{v:.0f}</text>")
    s.append(f"<text x='{W/2}' y='16' fill='{FG}' font-size='11' text-anchor='middle'>Explained variance (%)</text>")
    s.append("</svg>")
    return "".join(s)


def svg_protein_heatmap(matrix: pd.DataFrame, meta: pd.DataFrame, uid="cbprot"):
    """protein x contrast heatmap of significant log2 fold-changes."""
    if matrix.empty:
        return "<p class='mut'>No significant protein.</p>"
    contrasts = list(matrix.columns)
    vmax = float(np.nanmax(np.abs(matrix.values))) or 1.0
    cw, ch = 54, 22
    left, top = 260, 150
    W = left + len(contrasts) * cw + 30
    H = top + len(matrix) * ch + 42
    s = [f"<svg viewBox='0 0 {W} {H}' width='{W}' height='{H}' style='max-width:100%;height:auto'>"]
    for j, c in enumerate(contrasts):
        x = left + j * cw + cw / 2
        s.append(f"<text x='{x}' y='{top-8}' fill='{MUT}' font-size='9' text-anchor='start' "
                 f"transform='rotate(-90 {x} {top-8})'>{_esc(c)}</text>")
    for i, key in enumerate(matrix.index):
        y = top + i * ch
        g = meta.loc[key, "gene"] if "gene" in meta.columns else None
        label = str(g) if isinstance(g, str) and g not in ("", "nan") else str(key)
        n = int(meta.loc[key, "n_sig"]) if "n_sig" in meta.columns else 0
        s.append(f"<text x='{left-8}' y='{y+ch*0.7}' fill='{FG}' font-size='10' "
                 f"text-anchor='end'>{_esc(label[:30])} <tspan fill='{MUT}'>({n})</tspan></text>")
        for j, c in enumerate(contrasts):
            v = matrix.iloc[i, j]
            col = _diverging(v if pd.notna(v) else None, vmax)
            tip = "" if pd.isna(v) else f"<title>{_esc(label)} — {_esc(c)} : LFC={v:+.2f}</title>"
            s.append(f"<rect x='{left+j*cw}' y='{y}' width='{cw-2}' height='{ch-2}' "
                     f"fill='{col}' rx='2'>{tip}</rect>")
    s.append(_colorbar(left, H - 28, vmax, uid=uid))
    s.append(f"<text x='{left-8}' y='{H-19}' fill='{MUT}' font-size='8' text-anchor='end'>log2 FC</text>")
    return "".join(s) + "</svg>"


def svg_contrast_corr(corr: pd.DataFrame, shared: pd.DataFrame, relation=None, uid="cbctr"):
    """contrast x contrast similarity heatmap (diverging, hierarchical order)."""
    labels = list(corr.columns)
    n = len(labels)
    if n < 2:
        return ""
    cell = min(30, max(12, int(900 / n)))
    lab_w = min(190, 8 + max(len(str(l)) for l in labels) * 5)
    left, top = lab_w, lab_w
    W = left + n * cell + 24
    H = top + n * cell + 44
    s = [f"<svg viewBox='0 0 {W} {H}' width='{W}' height='{H}' style='max-width:100%;height:auto'>"]
    for i, li in enumerate(labels):
        s.append(f"<text x='{left-6}' y='{top+i*cell+cell*0.72}' fill='{FG}' font-size='9' "
                 f"text-anchor='end'>{_esc(li)}</text>")
        x = left + i * cell + cell / 2
        s.append(f"<text x='{x}' y='{top-6}' fill='{MUT}' font-size='9' text-anchor='start' "
                 f"transform='rotate(-90 {x} {top-6})'>{_esc(li)}</text>")
    for i, li in enumerate(labels):
        for j, lj in enumerate(labels):
            v = corr.iloc[i, j]
            sh = int(shared.iloc[i, j]) if shared is not None else 0
            if i == j:
                fill = "#e6edf3"
            elif pd.isna(v):
                fill = "#1c2128"
            else:
                fill = _diverging(float(v), 1.0)
            rl = ""
            if relation is not None and i != j:
                try:
                    rr = str(relation.loc[li, lj])
                    if rr and rr != "disjoint":
                        rl = f" · {rr.upper()} (mechanical)"
                except Exception:
                    rl = ""
            tip = (f"{_esc(li)} / {_esc(lj)} : "
                   + ("r=n/a" if pd.isna(v) else f"r={float(v):+.2f}")
                   + f" · {sh} shared terms{rl}")
            s.append(f"<rect x='{left+j*cell}' y='{top+i*cell}' width='{cell-1}' height='{cell-1}' "
                     f"fill='{fill}'><title>{tip}</title></rect>")
    s.append(_colorbar(left, H - 28, 1.0, width=150, uid=uid))
    return "".join(s) + "</svg>"


def svg_corr(corr: pd.DataFrame):
    labels = list(corr.columns)
    n = len(labels); cell = min(32, max(12, int(760 / max(n, 1))))
    left, top = 130, 130
    W = left + n * cell + 20; H = top + n * cell + 46
    V = corr.values.astype(float)
    off = V[~np.eye(n, dtype=bool)] if n > 1 else np.array([1.0])
    off = off[np.isfinite(off)]
    lo = float(np.nanmin(off)) if off.size else 0.0
    hi = float(np.nanmax(off)) if off.size else 1.0
    if hi - lo < 1e-6:
        hi = lo + 1e-6
    diverging = lo < 0 < hi
    amp = max(abs(lo), abs(hi))

    def col(v):
        if not np.isfinite(v):
            return "#1c2128"
        if diverging:                       # blue (neg) -> grey -> red (pos)
            return _diverging(v, amp)
        t = max(0.0, min(1.0, (v - lo) / (hi - lo)))     # sequential green
        return f"rgb({int(28 + t*20)},{int(45 + t*175)},{int(58 + t*40)})"

    s = [f"<svg viewBox='0 0 {W} {H}' width='{W}' height='{H}' style='max-width:100%;height:auto'>"]
    for i, ri in enumerate(labels):
        s.append(f"<text x='{left-6}' y='{top+i*cell+cell*0.7}' fill='{FG}' font-size='9' text-anchor='end'>{_esc(ri)}</text>")
        x = left + i * cell + cell / 2
        s.append(f"<text x='{x}' y='{top-6}' fill='{MUT}' font-size='9' text-anchor='end' transform='rotate(-50 {x} {top-6})'>{_esc(ri)}</text>")
    for i, ri in enumerate(labels):
        for j, cj in enumerate(labels):
            v = float(V[i, j])
            fill = "#e6edf3" if i == j else col(v)
            s.append(f"<rect x='{left+j*cell}' y='{top+i*cell}' width='{cell-1}' height='{cell-1}' fill='{fill}'>"
                     f"<title>{_esc(ri)} / {_esc(cj)} : r={v:.3f}</title></rect>")
    # graduated colour bar over the observed off-diagonal range
    barw = 150
    if diverging:
        s.append(_colorbar(left, H - 30, amp, width=barw, uid="cbcorr"))
    else:
        stops = "".join(f"<stop offset='{k*100/20:.0f}%' stop-color='{col(lo + (hi-lo)*k/20)}'/>" for k in range(21))
        s.append(f"<defs><linearGradient id='cbcorr' x1='0' y1='0' x2='1' y2='0'>{stops}</linearGradient></defs>")
        s.append(f"<rect x='{left}' y='{H-30}' width='{barw}' height='9' fill='url(#cbcorr)' rx='2' stroke='{GRID}'/>")
        for k, lab in enumerate([lo, (lo+hi)/2, hi]):
            lx = left + barw * k / 2
            s.append(f"<line x1='{lx}' y1='{H-21}' x2='{lx}' y2='{H-18}' stroke='{MUT}'/>")
            s.append(f"<text x='{lx}' y='{H-10}' fill='{MUT}' font-size='8' text-anchor='middle'>{lab:.2f}</text>")
        s.append(f"<text x='{left-6}' y='{H-22}' fill='{MUT}' font-size='8' text-anchor='end'>r</text>")
    s.append("</svg>")
    return "".join(s)


# --------------------------------------------------------------------------- #
# Tables                                                                        #
# --------------------------------------------------------------------------- #
def _table(df: pd.DataFrame, max_rows=60, floatfmt="{:.3g}"):
    df = df.head(max_rows)
    cols = list(df.columns)
    h = "".join(f"<th>{_esc(c)}</th>" for c in cols)
    rows = []
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, float):
                v = "" if pd.isna(v) else floatfmt.format(v)
            cells.append(f"<td>{_esc(v)}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return f"<table><thead><tr>{h}</tr></thead><tbody>{''.join(rows)}</tbody></table>"


# --------------------------------------------------------------------------- #
# Assemblage                                                                    #
# --------------------------------------------------------------------------- #
def _svg_chain(conv):
    """Core proteins -> module -> comparison, as one glance."""
    if not conv:
        return ""
    W, H = 980, 132
    box_w, box_h, gap = 260, 76, 60
    y = 26
    s = [f"<svg viewBox='0 0 {W} {H}' width='{W}' height='{H}' style='max-width:100%;height:auto'>"]

    def box(x, title, value, detail, colour, dim=False):
        op = "0.35" if dim else "1"
        s.append(f"<g opacity='{op}'>")
        s.append(f"<rect x='{x}' y='{y}' width='{box_w}' height='{box_h}' rx='8' "
                 f"fill='{PANEL}' stroke='{colour}' stroke-width='1.5'/>")
        s.append(f"<text x='{x+14}' y='{y+19}' fill='{MUT}' font-size='10' "
                 f"letter-spacing='0.06em'>{_esc(title.upper())}</text>")
        val = str(value)
        if len(val) > 40:
            val = val[:39] + "…"
        fs = 19 if len(val) <= 20 else max(11, int(19 * 20 / len(val)))   # fit, don't cut
        s.append(f"<text x='{x+14}' y='{y+43}' fill='#fff' font-size='{fs}' "
                 f"font-weight='700'>{_esc(val)}</text>")
        s.append(f"<text x='{x+14}' y='{y+63}' fill='{MUT}' font-size='10'>"
                 f"{_esc(str(detail)[:44])}</text>")
        s.append("</g>")

    def arrow(x, ok=True):
        col = GREEN if ok else GRID
        s.append(f"<line x1='{x}' y1='{y+box_h/2}' x2='{x+gap-14}' y2='{y+box_h/2}' "
                 f"stroke='{col}' stroke-width='2'/>")
        s.append(f"<path d='M{x+gap-14},{y+box_h/2} l-8,-5 l0,10 z' fill='{col}'/>")

    linked = bool(conv.get("module"))
    box(0, conv.get("core_title", "core signature"), f"{conv['core_n']} proteins",
        conv.get("core_detail", "recurrent across contrasts"), BLUE)
    arrow(box_w, linked)
    if linked:
        box(box_w + gap, "co-expression module", conv["module"],
            conv.get("module_detail", ""), GREEN)
        arrow(2 * box_w + gap, bool(conv.get("contrast")))
        if conv.get("contrast"):
            box(2 * (box_w + gap), "separates best", conv["contrast"],
                conv.get("contrast_detail", "")[:44], PURP)
        else:
            box(2 * (box_w + gap), "comparison", "not linked",
                "no eigengene available", GRID, dim=True)
    else:
        box(box_w + gap, "co-expression module", "scattered",
            "no module concentrates it", AMBER)
        box(2 * (box_w + gap), "comparison", "—", "", GRID, dim=True)
    s.append("</svg>")
    return "".join(s)


def _resume_html(resume):
    if resume is None or not resume.get("sections"):
        return "<p class='mut'>Summary unavailable.</p>"
    out = []

    # verdict cards: the five facts, colour-coded
    cards = resume.get("headline") or []
    if cards:
        cols = {"ok": GREEN, "warn": AMBER, "info": BLUE}
        chips = "".join(
            f"<div class='vcard' style='border-left:3px solid {cols.get(c.get('level'), BLUE)}'>"
            f"<div class='vl'>{_esc(c['label'])}</div>"
            f"<div class='vv'>{_esc(str(c['value']))}</div>"
            f"<div class='vd'>{_esc(str(c.get('detail', '')))}</div></div>"
            for c in cards)
        out.append(f"<div class='vcards'>{chips}</div>")

    # the chain everyone asks about: core proteins -> module -> comparison
    conv = resume.get("convergence")
    if conv:
        out.append("<h3>Where the signature lives</h3>" + _svg_chain(conv))
        if conv.get("verdict"):
            out.append(f"<p class='verdict'>{_esc(conv['verdict'])}</p>")

    for title, lines in resume["sections"]:
        items = []
        for l in lines:
            cls = "warn" if l.strip().startswith("⚠") else ""
            items.append(f"<li class='{cls}'>{_esc(l)}</li>")
        out.append(f"<h3>{_esc(title)}</h3><ul class='resume'>{''.join(items)}</ul>")
    return "".join(out)


EXPORT_JS = """
(function(){
  var BG_NIGHT = '#0d1117';
  function isDay(){ return document.documentElement.getAttribute('data-theme') === 'light'; }
  function bg(){ return isDay() && typeof PF_LIGHT_BG !== 'undefined' ? PF_LIGHT_BG : BG_NIGHT; }
  // day theme: the page recolours figures with CSS; an exported file has no such CSS,
  // so the same swaps are written into the copy's attributes
  function dayColours(root){
    if(!isDay() || typeof PF_LIGHT_SWAP === 'undefined') return;
    var els = [root].concat(Array.prototype.slice.call(root.querySelectorAll('*')));
    els.forEach(function(el){
      ['fill', 'stroke'].forEach(function(a){
        var v = el.getAttribute(a); if(!v) return; var k = v.trim().toLowerCase();
        if(PF_LIGHT_SWAP[k]) el.setAttribute(a, PF_LIGHT_SWAP[k]);
        else if(a === 'fill' && /^(text|tspan)$/i.test(el.tagName) && PF_LIGHT_TEXT[k]) el.setAttribute(a, PF_LIGHT_TEXT[k]);
      });
    });
  }
  function slug(s){ return (s||'figure').toString().trim().toLowerCase()
      .replace(/[^a-z0-9]+/g,'_').replace(/^_+|_+$/g,'').slice(0,60) || 'figure'; }
  function headText(h){
    var c = h.cloneNode(true);
    c.querySelectorAll('.mut').forEach(function(e){ e.remove(); });   // drop the counters
    return c.textContent;
  }
  function figName(el){
    var node = el, base = null;
    var lab = el.querySelector ? el.querySelector('svg.rd[aria-label]') : null;
    if(lab) return slug(lab.getAttribute('aria-label'));
    for(var up = 0; up < 5 && node; up++){
      var n = node.previousElementSibling, guard = 0;
      while(n && guard++ < 10){
        if(/^H[23]$/.test(n.tagName)){ base = slug(headText(n)); break; }
        n = n.previousElementSibling;
      }
      if(base) break;
      if(node.tagName === 'SECTION') break;
      node = node.parentElement;
    }
    if(!base){
      var sec = el.closest('section');
      base = slug(sec && sec.dataset.name ? sec.dataset.name : 'figure');
    }
    var view = el.closest('.hm-rep') ? '_non_redundant' : (el.closest('.hm-all') ? '_all_terms' : '');
    return base + view;
  }
  function dl(href, name){
    var a = document.createElement('a'); a.href = href; a.download = name;
    document.body.appendChild(a); a.click(); document.body.removeChild(a);
  }
  function svgString(svg){
    var c = svg.cloneNode(true);
    var vb = (svg.getAttribute('viewBox')||'').split(/\\s+/);
    var w = vb.length===4 ? +vb[2] : (svg.clientWidth||800);
    var h = vb.length===4 ? +vb[3] : (svg.clientHeight||600);
    // the inline style carries 'height:auto', which resolves to 0 once the SVG is
    // serialised as a standalone document: drop it before fixing explicit sizes.
    c.removeAttribute('style');
    c.setAttribute('width', w); c.setAttribute('height', h);
    dayColours(c);
    var back = document.createElementNS('http://www.w3.org/2000/svg','rect');
    back.setAttribute('x',0); back.setAttribute('y',0);
    back.setAttribute('width',w); back.setAttribute('height',h); back.setAttribute('fill',bg());
    c.insertBefore(back, c.firstChild);
    // XMLSerializer emits the SVG namespace itself; declaring it by hand would
    // produce a duplicate xmlns attribute and an unparseable document.
    var txt = new XMLSerializer().serializeToString(c);
    if (txt.indexOf('xmlns=') === -1){
      txt = txt.replace('<svg', '<svg xmlns="http://www.w3.org/2000/svg"');
    }
    return {txt:'<?xml version="1.0" encoding="UTF-8"?>\\n'+txt, w:w, h:h};
  }
  function exportSVG(svg, name){
    var o = svgString(svg);
    var blob = new Blob([o.txt], {type:'image/svg+xml;charset=utf-8'});
    var url = URL.createObjectURL(blob);
    dl(url, name+'.svg'); setTimeout(function(){URL.revokeObjectURL(url);}, 4000);
  }
  function svgDataUrl(txt){
    // data: URL, not blob: — a blob URL has a null origin when the report is opened
    // from file://, which taints the canvas and makes toDataURL throw.
    return 'data:image/svg+xml;base64,' + btoa(unescape(encodeURIComponent(txt)));
  }
  function exportPNG(svg, name, scale){
    var o = svgString(svg); scale = scale || 2;
    var img = new Image();
    img.onload = function(){
      try {
        var cv = document.createElement('canvas');
        cv.width = Math.round(o.w * scale); cv.height = Math.round(o.h * scale);
        var ctx = cv.getContext('2d');
        ctx.fillStyle = bg(); ctx.fillRect(0, 0, cv.width, cv.height);
        ctx.drawImage(img, 0, 0, cv.width, cv.height);
        dl(cv.toDataURL('image/png'), name + '.png');
      } catch(err){
        alert('PNG export blocked by the browser (' + err.name + '). Use the SVG button — '
              + 'it opens in Illustrator, Inkscape or any browser.');
      }
    };
    img.onerror = function(){
      alert('PNG export failed while rendering the figure. Use the SVG button instead.');
    };
    img.src = svgDataUrl(o.txt);
  }
  function canvasPNG(cv, name){
    var out = document.createElement('canvas');
    out.width = cv.width; out.height = cv.height;
    var ctx = out.getContext('2d');
    ctx.fillStyle = bg(); ctx.fillRect(0,0,out.width,out.height);
    ctx.drawImage(cv, 0, 0);
    dl(out.toDataURL('image/png'), name+'.png');
  }
  function addBtns(el, isCanvas){
    if(el.parentElement && el.parentElement.classList.contains('figwrap')) return;
    var wrap = document.createElement('div'); wrap.className = 'figwrap';
    el.parentNode.insertBefore(wrap, el); wrap.appendChild(el);
    var bar = document.createElement('div'); bar.className = 'figbtns';
    var name = 'pathfinder_' + figName(wrap);
    if(!isCanvas){
      var b1 = document.createElement('button'); b1.className='figbtn'; b1.textContent='SVG';
      b1.title='Download as vector SVG';
      b1.onclick = function(e){ e.stopPropagation(); exportSVG(el, name); };
      bar.appendChild(b1);
    }
    var b2 = document.createElement('button'); b2.className='figbtn'; b2.textContent='PNG';
    b2.title='Download as PNG (2x)';
    b2.onclick = function(e){ e.stopPropagation();
      isCanvas ? canvasPNG(el, name) : exportPNG(el, name, 2); };
    bar.appendChild(b2);
    wrap.appendChild(bar);
  }
  function init(){
    document.querySelectorAll('main svg').forEach(function(s){ addBtns(s, false); });
    document.querySelectorAll('main canvas').forEach(function(c){ addBtns(c, true); });
  }
  if(document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
  else init();
})();
"""


READING_JS = """
(function(){
  // clickable legends of the reading figures: hide / show the data-g series
  document.addEventListener('click', function(e){
    var it = e.target.closest ? e.target.closest('svg.rd .lg[data-t]') : null;
    if(!it) return;
    var svg = it.closest('svg'), k = it.getAttribute('data-t');
    var off = !it.classList.contains('off');
    it.classList.toggle('off', off);
    svg.querySelectorAll("[data-g='" + k + "']").forEach(function(g){
      g.style.display = off ? 'none' : ''; });
  });
})();
"""


def _reading_html(reading):
    """"Biological reading" tab (ex-InSight): numeric facts + SVG figures."""
    if reading and reading.get("error"):
        # The reading failed: keep the tab and say why. Dropping the tab leaves the
        # reason in pathfinder_log.txt alone, where nobody looks — a report missing
        # a whole section has to say so itself.
        e = reading["error"]
        out = ["<div class='failbox'><h4>\u26a0 Biological reading failed</h4>",
               "<p>This tab is empty: the biological reading raised an error, so its "
               "programmes, functional groups and figures could not be computed. "
               "Everything else in this report is unaffected.</p>",
               f"<p><code>{_esc(str(e))}</code></p>"]
        if "SyntaxError" in str(e) or "ImportError" in str(e):
            out.append("<p>A <code>SyntaxError</code> or <code>ImportError</code> here "
                       "means the <code>reading/</code> package could not be loaded at "
                       "all — usually a partial update, or a Python older than the one "
                       "the code was written for. Check that the whole "
                       "<code>reading/</code> folder matches <code>pathfinder.py</code>.</p>")
        out.append("<p class='mut'>Full traceback: set <code>PATHFINDER_DEBUG=1</code> "
                   "and run Pathfinder again.</p></div>")
        return "".join(out)
    if not reading or not reading.get("sections"):
        return ("<p class='mut'>Biological reading unavailable: it needs the differential "
                "results and, for the programmes and the spotlight, per-sample abundances "
                "(Deimos <code>matrices/</code> folder or a <code>Log2_Impute</code> sheet).</p>")
    out = []
    k = "".join(f"<div class='kpi'><div class='kv'>{_esc(v)}</div><div class='kl'>{_esc(n)}</div></div>"
                for n, v in reading.get("kpis", []))
    out.append(f"<div class='cards'>{k}</div>")
    m = reading.get("meta", {})
    lines = [f"Design: {m.get('reading_design', '')}",
             f"Proteins: {reading.get('rule_text', '')} — the selection of the Proteins tab",
             f"Abundances: {m.get('reading_abundance', 'n/a')} · protein names: "
             f"{m.get('reading_protein_names', '')}"]
    out.append("<p class='mut'>" + "<br>".join(_esc(l) for l in lines) + "</p>")
    chips = "".join(f"<a href='#rd-{_esc(s['id'])}'>{_esc(s['title'])}</a>"
                    for s in reading["sections"])
    out.append(f"<div class='rdnav'>{chips}</div>")
    for s in reading["sections"]:
        bl = "".join(f"<li>{b}</li>" for b in s.get("bullets", []))      # HTML already escaped
        figs = "".join(f"<div class='rdfig'>{svg}</div>" for svg in s.get("svgs", []))
        out.append(f"<h3 id='rd-{_esc(s['id'])}'>{_esc(s['title'])}</h3>"
                   f"<ul class='resume'>{bl}</ul>{figs}")
    out.append("<p class='mut'>Facts are computed from the results; interpretation is yours. "
               "Hover a mark for its values; click a legend entry of a volcano to hide or show "
               "that series. Tables: sheets <code>Programmes</code>, <code>Groups</code>, "
               "<code>Reading_facts</code> of <code>pathfinder_synthesis.xlsx</code>.</p>")
    return "".join(out)


def build(path, meta_go, explore, provenance: dict, kpis: dict, resume=None, branding=None, pubs=None, prots=None, wgcna=None, reading=None):
    b = {"title": "Pathfinder",
         "subtitle": "Pathway And Term Harmonisation For INtegrated Design Exploration & Reporting",
         "affiliation": "Proteogen · Université de Caen",
         "icon_b64": None, "icon_svg": None}
    if branding:
        b.update({k: v for k, v in branding.items() if v is not None})
    tabs = ["Summary", "Overview", "Proteins", "Meta-GO", "Exploration", "Protein UMAP"]
    # KPI cards
    kcards = "".join(
        f"<div class='kpi'><div class='kv'>{_esc(v)}</div><div class='kl'>{_esc(k)}</div></div>"
        for k, v in kpis.items())

    # Overview
    prov_tbl = _table(pd.DataFrame(list(provenance.items()), columns=["key", "value"]), max_rows=30)
    overview = f"<div class='cards'>{kcards}</div><h3>Provenance</h3>{prov_tbl}"

    # Meta-GO — split by source, with a representatives / all-terms switch
    matrix, meta, direction = meta_go["matrix"], meta_go["meta"], meta_go["directionality"]
    if matrix.empty:
        metago = "<p class='mut'>GO enrichment not provided — functional synthesis unavailable.</p>"
    else:
        has_rep = "is_representative" in meta.columns and bool(meta["is_representative"].any())
        # reorder contrasts by GO-signature similarity (falls back to given order)
        corder = [c for c in (meta_go.get("contrast_order") or []) if c in matrix.columns]
        if len(corder) == matrix.shape[1]:
            matrix = matrix[corder]
        n_clusters = meta["redundancy_cluster"].nunique() if (meta["redundancy_cluster"] >= 0).any() else 0
        n_rep = int(meta["is_representative"].sum()) if has_rep else matrix.shape[0]
        ctl = ("<div class='goctl'>"
               "<label>Terms shown <select id='go_view'>"
               "<option value='rep' selected>Non-redundant</option>"
               "<option value='all'>All terms</option></select></label>"
               f"<span class='mut'>{n_rep} non-redundant · {matrix.shape[0]} terms total · "
               f"{n_clusters} redundancy groups (Jaccard)</span></div>") if has_rep else ""
        parts = [ctl]
        ccorr = meta_go.get("contrast_corr")
        if ccorr is not None and not ccorr.empty and ccorr.shape[0] >= 2:
            cc = svg_contrast_corr(ccorr, meta_go.get("contrast_shared"), meta_go.get("contrast_relation"))
            if cc:
                off = ccorr.where(~np.eye(ccorr.shape[0], dtype=bool)).stack().dropna()
                rel = meta_go.get("contrast_relation")
                extra = ""
                if len(off):
                    if rel is not None:
                        disj = off[[rel.loc[a, b] == "disjoint" for a, b in off.index]]
                    else:
                        disj = off
                    R_MIN, N_MIN = 0.3, 10
                    shm = meta_go.get("contrast_shared")
                    n_thin = 0
                    if shm is not None and len(disj):
                        enough = [int(shm.loc[a, b]) >= N_MIN for a, b in disj.index]
                        n_thin = len(disj) - sum(enough)
                        disj = disj[enough]

                    def _nsh(a, b):
                        return int(shm.loc[a, b]) if shm is not None else 0

                    if len(disj):
                        tp = disj.idxmax(); tv = disj.max()
                        np_ = disj.idxmin(); nv = disj.min()
                        if tv >= R_MIN:
                            extra = (f" Most similar independent pair: <b>{_esc(tp[0])}</b> ↔ "
                                     f"<b>{_esc(tp[1])}</b> (r={tv:+.2f}, {_nsh(*tp)} shared terms).")
                        if nv <= -R_MIN:
                            extra += (f" Most opposed: <b>{_esc(np_[0])}</b> ↔ "
                                      f"<b>{_esc(np_[1])}</b> (r={nv:+.2f}, {_nsh(*np_)} shared terms).")
                        if tv < R_MIN and nv > -R_MIN:
                            extra = (f" No independent pair reaches |r| ≥ {R_MIN}: signatures are "
                                     "largely unrelated.")
                    else:
                        extra = (f" No pair is both independent and based on ≥ {N_MIN} shared terms.")
                    n_sh = int(sum(1 for a_, b_ in off.index if rel is not None
                                   and rel.loc[a_, b_] != "disjoint"))
                    caveats = []
                    if n_sh:
                        caveats.append(f"{n_sh}/{len(off)} pairs share a condition (a common "
                                       "reference inflates r, a reversed one inflates −r, "
                                       "regardless of biology — hover a cell to see the relation)")
                    if n_thin:
                        caveats.append(f"{n_thin} pair(s) rest on fewer than {N_MIN} shared terms, "
                                       "where correlation is unstable")
                    if caveats:
                        extra += (" <b>Excluded from the highlights:</b> " + "; ".join(caveats) + ".")
                parts.append(
                    "<h3>Contrast clustering <span class='mut'>· by GO signature</span></h3>" + cc +
                    "<p class='mut'>Pearson correlation between the signed z-score profiles of two "
                    "contrasts, over the terms significant in both. Red = they move the same functions "
                    "the same way; blue = opposite directions; grey = too few shared terms to compare. "
                    "Rows/columns — and the heatmaps below — are ordered by this clustering." + extra +
                    "</p>")
        parts.append(f"<p class='mut'>{matrix.shape[0]} significant terms across {matrix.shape[1]} "
                     f"contrast(s). Red = up, blue = down, grey = not significant. "
                     f"Redundant terms (hitting the same proteins) are collapsed into a single "
                     f"non-redundant entry; switch above to list every term.</p>")
        for src in _sources_present(meta):
            keys = meta.index[meta["source"] == src]
            sid = re.sub(r"[^A-Za-z0-9]", "", str(src))
            sub_meta_all = meta.loc[keys]
            rep_keys = keys[sub_meta_all["is_representative"]] if has_rep else keys
            blocks = []
            for view, kk in (("rep", rep_keys), ("all", keys)):
                if has_rep is False and view == "rep":
                    continue
                hide = "" if view == "rep" else " style='display:none'"
                if len(kk) == 0:
                    blocks.append(f"<div class='hm-{view}'{hide}>"
                                  f"<p class='mut'>No term in this view.</p></div>")
                    continue
                hm = svg_meta_heatmap(matrix.loc[kk], meta.loc[kk], uid=f"cb{sid}{view}")
                dt = direction.loc[kk].drop(columns=["source"], errors="ignore")
                tbl = _table(dt, max_rows=25)
                blocks.append(f"<div class='hm-{view}'{hide}>{hm}{tbl}</div>")
            label_rep = f" <span class='mut'>· {len(rep_keys)} non-redundant / {len(keys)} terms</span>" \
                if has_rep else f" <span class='mut'>· {len(keys)} terms</span>"
            parts.append(f"<h3>{_esc(SRC_LABEL.get(src, src))}{label_rep}</h3>" + "".join(blocks))
        if has_rep:
            parts.append("""<script>
(function(){var sel=document.getElementById('go_view');if(!sel)return;
 sel.addEventListener('change',function(){var rep=sel.value==='rep';
  document.querySelectorAll('.hm-rep').forEach(function(e){e.style.display=rep?'block':'none';});
  document.querySelectorAll('.hm-all').forEach(function(e){e.style.display=rep?'none':'block';});});
})();
</script>""")
        metago = "".join(parts)

    # Exploration
    pca = explore["pca"]; assoc = explore["pc_factor"]
    scatter = svg_scatter(pca["scores"], "PC1", "PC2", "condition",
                          f"PC1 ({pca['var_exp'][0]:.1f}%)",
                          f"PC2 ({pca['var_exp'][1]:.1f}%)" if len(pca['var_exp']) > 1 else "PC2")
    scree = svg_scree(pca["var_exp"])
    corr = svg_corr(explore["correlation"])
    assoc_tbl = _table(assoc, max_rows=10)
    # top PC1 loadings
    load = pca["loadings"]
    load_tbl = ""
    if "PC1" in load.columns:
        # positional sort: gene symbols repeat, and reindexing on duplicates raises
        top = load.iloc[np.argsort(-load["PC1"].abs().to_numpy())].head(15)
        load_tbl = _table(top[["PC1"] + [c for c in ["PC2", "PC3"] if c in top.columns]].reset_index()
                          .rename(columns={"index": "protein"}), max_rows=15)
    # empty association (design unresolved) -> note instead of an empty table
    if assoc["p_value"].isna().all():
        assoc_block = ("<p class='warn'>⚠ PC ↔ condition test not computable: conditions unresolved "
                       "(each sample isolated) or no replicates. Provide a design "
                       "(<code>--design</code>) to enable this test.</p>")
    else:
        assoc_block = _table(assoc, max_rows=10)

    # sample UMAP (design)
    umap_s_html = ""
    us = explore.get("umap_samples")
    if us:
        if us.get("df") is not None:
            umap_s_html = ("<h3>Sample UMAP</h3>"
                           + svg_scatter(us["df"], "UMAP1", "UMAP2", "condition", "UMAP1", "UMAP2"))
        elif us.get("warning"):
            umap_s_html = f"<p class='mut'>Sample UMAP: {_esc(us['warning'])}</p>"

    # multi-factor: which factor drives which axis
    multi = explore.get("pc_factor_multi")
    multi_block = ""
    if multi is not None and not multi.empty:
        sig = multi[(multi["p_adj_BH"] < 0.10) | (multi["effect"] >= 0.14)]
        show = (sig if len(sig) else multi).sort_values(["effect"], ascending=False)
        multi_block = ("<h3>Which factor drives which axis? <span class='mut'>· PC × design factors</span></h3>"
                       + _table(show, max_rows=24)
                       + "<p class='mut'>Categorical: Kruskal-Wallis, effect = ε². Numeric (e.g. injection "
                         "order): Spearman, effect = |rho|. A technical factor tied to a high-variance "
                         "PC ⇒ batch effect / drift to correct.</p>")

    corr_note = ("<p class='mut'>Pearson r on <b>per-protein z-scores</b> (no raw intensity matrix in "
                 "the workbook): rows are centred, so absolute r values are low by construction "
                 "(mean ≈ −1/(n−1)) — read the <b>relative</b> structure, not the magnitude. "
                 "Replicates should be relatively redder/greener than unrelated samples.</p>"
                 if explore["correlation"].attrs.get("zscored") else
                 "<p class='mut'>Pearson r on log₂ intensities, ordered by hierarchical clustering.</p>")

    explore_html = (
        f"<div class='grid2'><div><h3>PCA — PC1 vs PC2</h3>{scatter}</div>"
        f"<div><h3>Scree</h3>{scree}</div></div>"
        f"<h3>PC ↔ condition association (Kruskal-Wallis, ε² = effect size)</h3>{assoc_block}"
        f"<p class='mut'>High ε² + low p on a PC ⇒ that axis captures the biological factor; "
        f"on an unexpected PC, suspect a batch effect.</p>"
        f"{multi_block}"
        f"<div class='grid2'><div><h3>Inter-sample correlation</h3>{corr}{corr_note}</div>"
        f"<div><h3>Top PC1 contributors</h3>{load_tbl}</div></div>"
        f"{umap_s_html}")

    # Protein UMAP (dedicated tab)
    up = explore.get("umap_proteins")
    if up and up.get("df") is not None:
        udf = up["df"]
        blocks = ["<p class='mut'>Each point is a protein, positioned by the similarity of its "
                  "intensity profile across samples. Colour reveals structure: high/low abundance "
                  "zones, and regions regulated in each contrast.</p>"]
        blocks.append("<h3>Coloured by mean abundance (log2)</h3>"
                      + svg_scatter_value(udf, "UMAP1", "UMAP2", "abondance", diverging=False))
        for c in up.get("contrasts_lfc", [])[:3]:
            blocks.append(f"<h3>Coloured by LFC — {_esc(c)} <span class='mut'>(red=up, blue=down)</span></h3>"
                          + svg_scatter_value(udf, "UMAP1", "UMAP2", f"LFC:{c}", diverging=True))
        # colour by GO term — interactive widget (source + term count)
        if not meta_go["matrix"].empty and meta_go.get("term_genes"):
            widget = _go_umap_widget(udf, meta_go["meta"], meta_go["term_genes"])
            if widget:
                blocks.append("<h3>Coloured by GO term <span class='mut'>· interactive</span></h3>"
                              + widget
                              + "<p class='mut'>Pick the source (e.g. GO:CC for compartments) and the number of "
                                "terms. A function concentrated in one region ⇒ spatial co-regulation.</p>")
        umap_p_html = "".join(blocks)
    elif up and up.get("warning"):
        umap_p_html = f"<p class='mut'>{_esc(up['warning'])}</p>"
    else:
        umap_p_html = "<p class='mut'>Protein UMAP unavailable.</p>"

    # Proteins - which proteins carry the signature
    if prots is not None and not prots["matrix"].empty:
        t, matrix = prots["table"], prots["matrix"]
        corder = [c for c in (prots.get("contrast_order") or []) if c in matrix.columns]
        m = matrix[corder] if len(corder) == matrix.shape[1] else matrix
        top = t.head(40)
        hm = svg_protein_heatmap(m.reindex(top.index), top)
        cols = [c for c in ("gene", "n_sig", "pattern", "best_pi", "median_robustness",
                            "max_abs_lfc", "module", "go_terms", "up_in", "down_in")
                if c in t.columns]
        tbl = _table(t[cols].head(30).reset_index().rename(columns={"index": "protein"}),
                     max_rows=30)
        cc = prots.get("contrast_corr")
        cc_html = ""
        if cc is not None and not cc.empty and cc.shape[0] >= 2:
            fig = svg_contrast_corr(cc, prots.get("contrast_shared"), prots.get("contrast_relation"),
                                    uid="cbctrprot")
            if fig:
                cc_html = ("<h3>Contrast clustering <span class='mut'>· by protein signature</span></h3>"
                           + fig + "<p class='mut'>Same measure as the Meta-GO tab, computed on "
                           "protein fold-changes instead of term z-scores. Agreement between the "
                           "two views is a strong sign; disagreement means the shared functions "
                           "rest on different proteins.</p>")
        pi_note = ("Ranked by recurrence, then by Pi-score (|LFC| × −log10 p). The Pi-score orders "
                   "proteins; inclusion rests on the adjusted p-value alone."
                   if prots.get("has_pi") else
                   "No Pi-score found upstream: ranked by recurrence then |LFC|.")
        proteins_html = (
            f"<p class='mut'>{len(matrix)} proteins significant in at least one of the "
            f"{matrix.shape[1]} contrasts. {pi_note}</p>"
            f"<h3>Protein × contrast <span class='mut'>· top 40 by recurrence</span></h3>{hm}"
            f"<h3>Core signature</h3>{tbl}"
            f"<p class='mut'>Full table in the <code>Proteins_core</code> sheet, complete matrix "
            f"in <code>Proteins_matrix</code>.</p>{cc_html}")
    else:
        proteins_html = "<p class='mut'>No protein reaches significance in any contrast.</p>"

    panels = {"Summary": _resume_html(resume), "Overview": overview, "Proteins": proteins_html,
              "Meta-GO": metago, "Exploration": explore_html, "Protein UMAP": umap_p_html}
    if reading is not None:
        panels["Biological reading"] = _reading_html(reading)
        tabs = tabs[:1] + ["Biological reading"] + tabs[1:]

    # WGCNA - co-expression structure read from upstream
    if wgcna and wgcna.get("n_modules"):
        parts = [f"<p class='mut'>{wgcna['n_modules']} modules read from the upstream WGCNA "
                 "(Pathfinder does not recompute them). 'grey' is WGCNA's bin for unassigned "
                 "proteins, flagged as such and never treated as a module.</p>"]
        a = wgcna.get("association")
        if a is not None and not a.empty:
            parts.append("<h3>Module ↔ condition <span class='mut'>· eigengene per sample, "
                         "Kruskal-Wallis + ε²</span></h3>" + _table(a, max_rows=20)
                         + "<p class='mut'>The eigengene is the module's first principal "
                           "component: one value per sample, so it is testable exactly like a "
                           "PC score. A high ε² means the module responds to the design.</p>")
        e = wgcna.get("core_enrichment")
        if e is not None and not e.empty:
            parts.append("<h3>Where does the core signature sit? <span class='mut'>· Fisher exact"
                         "</span></h3>" + _table(e, max_rows=20)
                         + "<p class='mut'>Tests whether the recurrent proteins of the Proteins "
                           "tab concentrate in one module. Concentration means the signature is "
                           "co-expressed; dispersion means it cuts across the network.</p>")
        ca = wgcna.get("contrast_association")
        if ca is not None and not ca.empty:
            floor = ca.attrs.get("p_floor")
            note = ""
            if floor and floor > 0.05:
                note = (f" <b>With {ca.attrs.get('n_per_group')} replicates per group the exact "
                        f"test cannot go below p={floor:.2f}</b>: read the separation, not the "
                        "p-value.")
            parts.append("<h3>Which comparison does a module follow? <span class='mut'>· "
                         "eigengene, Mann-Whitney</span></h3>"
                         + _table(ca.head(20), max_rows=20)
                         + "<p class='mut'>Separation is the rank-biserial correlation: 1 means "
                           "the two groups of that contrast are perfectly split by the module "
                           "eigengene." + note + "</p>")
        h = wgcna.get("hubs")
        if h is not None and not h.empty:
            parts.append("<h3>Module hubs</h3>" + _table(h, max_rows=30))
        panels["WGCNA"] = "".join(parts)
        i_ = tabs.index("Proteins") + 1
        tabs = tabs[:i_] + ["WGCNA"] + tabs[i_:]

    # Publications (online axis) - tab added only when present
    if pubs is not None and any(pubs.get("results", {}).values()):
        crit = pubs.get("provenance", {}).get("_criterion", "")
        skip = pubs.get("provenance", {}).get("_skipped", "")
        pub_parts = ["<p class='mut'>Papers over-represented in each contrast signature (STRING "
                     "PMID enrichment, ranked by FDR). Known bias: heavily studied proteins surface "
                     "more readily — cross-check against the LFC.</p>"]
        if crit:
            pub_parts.append(f"<p class='mut'><b>Criterion:</b> {_esc(crit)}."
                             + (f" {_esc(skip)}." if skip else "")
                             + " Contrasts with fewer than 5 significant proteins are skipped: "
                               "enrichment is a statement about a set, not about one gene.</p>")
        for c, rows in pubs["results"].items():
            if not rows:
                continue
            trs = []
            for r in rows:
                link = f"<a href='{_esc(r.get('url',''))}' target='_blank' style='color:{BLUE}'>PMID:{_esc(r.get('pmid',''))}</a>"
                fdr = r.get("fdr")
                fdr_s = f"{fdr:.1e}" if isinstance(fdr, float) and fdr == fdr else ""
                trs.append(f"<tr><td>{link}</td><td>{_esc(r.get('title',''))}</td>"
                           f"<td>{_esc(r.get('n_genes',''))}</td><td>{fdr_s}</td></tr>")
            pub_parts.append(f"<h3>{_esc(c)}</h3><table><thead><tr><th>PMID</th><th>Title</th>"
                             f"<th>#genes</th><th>FDR</th></tr></thead><tbody>{''.join(trs)}</tbody></table>")
        panels["Publications"] = "".join(pub_parts)
        tabs = tabs + ["Publications"]
    nav = "".join(f"<button class='tab{' active' if i==0 else ''}' onclick=\"show({i})\">{_esc(t)}</button>"
                  for i, t in enumerate(tabs))
    sections = "".join(
        f"<section id='sec-{i}' data-name='{_esc(t)}' class='panel' style='display:{'block' if i==0 else 'none'}'>{panels[t]}</section>"
        for i, t in enumerate(tabs))

    doc = f"""<!doctype html><html lang='en' data-theme='dark'><head><meta charset='utf-8'>
<script>{THEME_HEAD_JS}</script>
<meta name='viewport' content='width=device-width, initial-scale=1'>
<title>Pathfinder — Exploration & Synthesis Report</title>
<style>
{theme_css()}
*{{box-sizing:border-box}} body{{margin:0;background:var(--bg);color:var(--fg);
 font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;font-size:14px}}
header{{padding:14px 22px;border-bottom:1px solid var(--grid);display:flex;align-items:center;gap:14px}}
header .logo{{width:46px;height:46px;flex:none;display:block}}
header .logo svg{{width:46px;height:46px;display:block}}
header img.logo{{width:46px;height:auto;border-radius:8px;flex:none}}
header .brand{{display:flex;flex-direction:column;gap:1px;line-height:1.25}}
header .brand .t{{font-size:20px;font-weight:700;color:var(--hi)}}
header .brand .s{{font-size:11px;color:var(--mut);letter-spacing:.02em}}
header .brand .a{{font-size:11px;color:var(--blue)}}
nav{{display:flex;gap:6px;padding:10px 24px;border-bottom:1px solid var(--grid);flex-wrap:wrap}}
.tab{{background:var(--panel);color:var(--fg);border:1px solid var(--grid);border-radius:6px;padding:7px 14px;cursor:pointer;font-size:13px}}
.tab.active{{background:var(--blue);border-color:var(--blue);color:#fff}}
main{{padding:20px 24px;max-width:1500px}}
h3{{color:var(--hi);font-size:14px;margin:22px 0 8px;font-weight:600}}
.cards{{display:flex;gap:12px;flex-wrap:wrap;margin-bottom:8px}}
.kpi{{background:var(--panel);border:1px solid var(--grid);border-radius:8px;padding:14px 18px;min-width:120px}}
.kv{{font-size:24px;font-weight:700;color:var(--hi)}} .kl{{color:var(--mut);font-size:11px;margin-top:2px}}
.grid2{{display:grid;grid-template-columns:1fr 1fr;gap:18px}}
@media(max-width:820px){{.grid2{{grid-template-columns:1fr}}}}
table{{border-collapse:collapse;width:100%;margin:6px 0 12px;font-size:12px}}
th,td{{border:1px solid var(--grid);padding:5px 8px;text-align:left}} th{{background:var(--panel);color:var(--fg)}}
td{{color:var(--fg)}} tbody tr:nth-child(even){{background:var(--stripe)}}
.mut{{color:var(--mut);font-size:12px}} .warn{{color:var(--amber);font-size:12px;border-left:3px solid var(--amber);padding-left:8px}}
.failbox{{border:1px solid var(--red);border-left:5px solid var(--red);border-radius:6px;padding:12px 14px;margin:4px 0 14px}}
.failbox h4{{margin:0 0 6px;color:var(--red);font-size:14px}}
.failbox code{{background:var(--panel);padding:1px 5px;border-radius:3px;font-size:12px;word-break:break-word}}
.failbox p{{margin:6px 0;font-size:13px}}
ul.resume{{margin:6px 0 4px;padding-left:20px}} ul.resume li{{margin:5px 0;line-height:1.5}}
ul.resume li.warn{{list-style:none;margin-left:-12px}}
.vcards{{display:flex;gap:10px;flex-wrap:wrap;margin:4px 0 18px}}
.vcard{{background:var(--panel);border:1px solid var(--grid);border-radius:8px;padding:10px 14px;min-width:170px;flex:1 1 170px}}
.vcard .vl{{color:var(--mut);font-size:10px;letter-spacing:.06em;text-transform:uppercase}}
.vcard .vv{{color:var(--hi);font-size:18px;font-weight:700;margin:3px 0 2px;line-height:1.2}}
.vcard .vd{{color:var(--mut);font-size:11px;line-height:1.35}}
.verdict{{background:var(--panel);border-left:3px solid var(--green);border-radius:6px;
 padding:10px 14px;margin:6px 0 16px;font-size:13px;color:var(--fg)}}
code{{background:var(--panel);padding:1px 5px;border-radius:4px;color:var(--purp)}}
.goctl{{display:flex;gap:16px;align-items:center;flex-wrap:wrap;margin:6px 0 10px;font-size:12px;color:var(--fg)}}
.goctl select{{background:var(--panel);color:var(--fg);border:1px solid var(--grid);border-radius:5px;padding:3px 6px;font-size:12px}}
.golegend{{display:flex;flex-wrap:wrap;gap:6px 16px;margin:8px 0 4px;font-size:11px;color:var(--fg)}}
.golegend .legit{{display:flex;align-items:center;gap:5px}}
.golegend .sw{{width:10px;height:10px;border-radius:2px;display:inline-block}}
.figwrap{{position:relative;display:block;width:fit-content;max-width:100%;margin:2px 0}}
.figbtns{{position:absolute;top:4px;right:4px;display:flex;gap:4px;opacity:0;transition:opacity .15s}}
.figwrap:hover .figbtns{{opacity:1}}
.figbtn{{background:var(--panel);color:var(--fg);border:1px solid var(--grid);border-radius:4px;
 padding:2px 7px;font-size:10px;cursor:pointer;line-height:1.5}}
.figbtn:hover{{background:var(--blue);border-color:var(--blue);color:#fff}}
.rdnav{{display:flex;flex-wrap:wrap;gap:6px;margin:4px 0 6px}}
.rdnav a{{color:var(--blue);text-decoration:none;font-size:12px;border:1px solid var(--grid);border-radius:14px;padding:3px 11px}}
.rdnav a:hover{{border-color:var(--blue)}}
.rdfig{{margin:6px 0 18px}}
svg.rd .pt:hover{{stroke:var(--hi);stroke-width:1.6}}
svg.rd .lg.off{{opacity:.35}}
.themebtn{{margin-left:auto;background:var(--panel);color:var(--fg);border:1px solid var(--grid);
 border-radius:16px;padding:5px 13px;font-size:12px;cursor:pointer;white-space:nowrap}}
.themebtn:hover{{border-color:var(--blue)}}
</style></head><body>
<header>{(f"<span class='logo'>{b['icon_svg']}</span>" if b.get('icon_svg') else (f"<img class='logo' src='data:image/png;base64,{b['icon_b64']}' alt=''/>" if b['icon_b64'] else ""))}<div class='brand'><span class='t'>{_esc(b['title'])}</span><span class='s'>{_esc(b['subtitle'])}</span><span class='a'>{_esc(b['affiliation'])}</span></div><button id='themeBtn' class='themebtn' onclick='toggleTheme()'>&#9728; Day</button></header>
<nav>{nav}</nav>
<main>{sections}</main>
<script>
function show(i){{document.querySelectorAll('.panel').forEach(p=>p.style.display='none');
 document.getElementById('sec-'+i).style.display='block';
 document.querySelectorAll('.tab').forEach(b=>b.classList.remove('active'));
 event.target.classList.add('active');}}
</script>
<script>{THEME_JS}</script>
<script>var PF_LIGHT_SWAP={json.dumps(LIGHT_SWAP)}, PF_LIGHT_TEXT={json.dumps(LIGHT_TEXT)}, PF_LIGHT_BG='{THEME_VARS["light"]["bg"]}';</script>
<script>{EXPORT_JS}</script>
<script>{READING_JS}</script>
</body></html>"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(doc)
    return path
