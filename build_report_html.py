"""
Self-contained Pathfinder HTML report.

Fully offline rendering: figures are generated as SVG on the Python side (no CDN,
no charting JS). Palette consistent with the pipeline dashboard.
"""
from __future__ import annotations
import html
import numpy as np
import pandas as pd

BG = "#0d1117"; PANEL = "#161b22"; GRID = "#30363d"; FG = "#c9d1d9"; MUT = "#8b949e"
BLUE = "#388bfd"; UP = "#f85149"; DOWN = "#388bfd"; GREEN = "#3fb950"; PURP = "#bc8cff"; AMBER = "#d29922"
COND_PALETTE = ["#388bfd", "#3fb950", "#bc8cff", "#d29922", "#f85149", "#39c5cf", "#db61a2", "#a5d6ff",
                "#7ee787", "#ffa657", "#d2a8ff", "#f2cc60", "#79c0ff", "#ff7b72", "#56d364", "#e3b341"]
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
    cw, ch = 46, 20
    left, top = 360, 150
    W = left + len(contrasts) * cw + 30
    H = top + len(m) * ch + 42
    s = [f"<svg viewBox='0 0 {W} {H}' width='100%' style='max-width:{W}px'>"]
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


def svg_scatter(df, xcol, ycol, cat_col, title_x, title_y, size=(520, 420)):
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

    s = [f"<svg viewBox='0 0 {W} {Hh}' width='100%' style='max-width:{W}px'>"]
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


def svg_scatter_value(df, xcol, ycol, valcol, diverging, size=(560, 440), r=3):
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
    s = [f"<svg viewBox='0 0 {W} {Hh}' width='100%' style='max-width:{W}px'>",
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
    s = [f"<svg viewBox='0 0 {W} {Hh}' width='100%' style='max-width:{W}px'>",
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
<canvas id='go_cv' width='1120' height='880' style='width:100%;max-width:560px;height:auto;background:#161b22;border-radius:6px'></canvas>
<div id='go_leg' class='golegend'></div>
<script>
(function(){
 var D=__DATA__, C=D.coords, M=D.memb, T=D.terms, PAL=D.palette;
 var cv=document.getElementById('go_cv'), ctx=cv.getContext('2d');
 var xs=C.map(function(p){return p[0];}), ys=C.map(function(p){return p[1];});
 var xmin=Math.min.apply(null,xs), xmax=Math.max.apply(null,xs);
 var ymin=Math.min.apply(null,ys), ymax=Math.max.apply(null,ys);
 var W=560,H=440,pad=30; ctx.scale(2,2);
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
    W, Hh = 520, 240; pad = 40
    n = len(ve); bw = (W - 2 * pad) / max(n, 1)
    vmax = max(ve) or 1
    s = [f"<svg viewBox='0 0 {W} {Hh}' width='100%' style='max-width:{W}px'>",
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


def svg_contrast_corr(corr: pd.DataFrame, shared: pd.DataFrame):
    """contrast x contrast similarity heatmap (diverging, hierarchical order)."""
    labels = list(corr.columns)
    n = len(labels)
    if n < 2:
        return ""
    cell = min(26, max(9, int(520 / n)))
    lab_w = min(190, 8 + max(len(str(l)) for l in labels) * 5)
    left, top = lab_w, lab_w
    W = left + n * cell + 24
    H = top + n * cell + 44
    s = [f"<svg viewBox='0 0 {W} {H}' width='100%' style='max-width:{W}px'>"]
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
            tip = (f"{_esc(li)} / {_esc(lj)} : "
                   + ("r=n/a" if pd.isna(v) else f"r={float(v):+.2f}")
                   + f" · {sh} shared terms")
            s.append(f"<rect x='{left+j*cell}' y='{top+i*cell}' width='{cell-1}' height='{cell-1}' "
                     f"fill='{fill}'><title>{tip}</title></rect>")
    s.append(_colorbar(left, H - 28, 1.0, width=150, uid="cbctr"))
    return "".join(s) + "</svg>"


def svg_corr(corr: pd.DataFrame):
    labels = list(corr.columns)
    n = len(labels); cell = min(28, max(10, int(560 / max(n, 1))))
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

    s = [f"<svg viewBox='0 0 {W} {H}' width='100%' style='max-width:{W}px'>"]
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
def _resume_html(resume):
    if resume is None or not resume.get("sections"):
        return "<p class='mut'>Summary unavailable.</p>"
    blocks = []
    for title, lines in resume["sections"]:
        items = []
        for l in lines:
            cls = "warn" if l.strip().startswith("⚠") else ""
            items.append(f"<li class='{cls}'>{_esc(l)}</li>")
        blocks.append(f"<h3>{_esc(title)}</h3><ul class='resume'>{''.join(items)}</ul>")
    return "".join(blocks)


def build(path, meta_go, explore, provenance: dict, kpis: dict, resume=None, branding=None, pubs=None):
    b = {"title": "Pathfinder",
         "subtitle": "Pathway And Term Harmonisation For INtegrated Design Exploration & Reporting",
         "affiliation": "Proteogen · Université de Caen",
         "icon_b64": None, "icon_svg": None}
    if branding:
        b.update({k: v for k, v in branding.items() if v is not None})
    tabs = ["Summary", "Overview", "Meta-GO", "Exploration", "Protein UMAP"]
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
            cc = svg_contrast_corr(ccorr, meta_go.get("contrast_shared"))
            if cc:
                off = ccorr.where(~np.eye(ccorr.shape[0], dtype=bool)).stack().dropna()
                extra = ""
                if len(off):
                    top_pair = off.idxmax(); top_v = off.max()
                    neg = off.min(); neg_pair = off.idxmin()
                    extra = (f" Most similar: <b>{_esc(top_pair[0])}</b> ↔ <b>{_esc(top_pair[1])}</b> "
                             f"(r={top_v:+.2f}).")
                    if neg < -0.3:
                        extra += (f" Most opposed: <b>{_esc(neg_pair[0])}</b> ↔ "
                                  f"<b>{_esc(neg_pair[1])}</b> (r={neg:+.2f}).")
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
        top = load.reindex(load["PC1"].abs().sort_values(ascending=False).index).head(15)
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

    panels = {"Summary": _resume_html(resume), "Overview": overview,
              "Meta-GO": metago, "Exploration": explore_html, "Protein UMAP": umap_p_html}

    # Publications (online axis) - tab added only when present
    if pubs is not None and any(pubs.get("results", {}).values()):
        pub_parts = ["<p class='mut'>Papers over-represented in each contrast signature (STRING "
                     "PMID enrichment, ranked by FDR). Known bias: heavily studied proteins surface "
                     "more readily — cross-check against the LFC.</p>"]
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
        f"<section id='sec-{i}' class='panel' style='display:{'block' if i==0 else 'none'}'>{panels[t]}</section>"
        for i, t in enumerate(tabs))

    doc = f"""<!doctype html><html lang='en'><head><meta charset='utf-8'>
<meta name='viewport' content='width=device-width, initial-scale=1'>
<title>Pathfinder — Exploration & Synthesis Report</title>
<style>
*{{box-sizing:border-box}} body{{margin:0;background:{BG};color:{FG};
 font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;font-size:14px}}
header{{padding:14px 22px;border-bottom:1px solid {GRID};display:flex;align-items:center;gap:14px}}
header .logo{{width:46px;height:46px;flex:none;display:block}}
header .logo svg{{width:46px;height:46px;display:block}}
header img.logo{{width:46px;height:auto;border-radius:8px;flex:none}}
header .brand{{display:flex;flex-direction:column;gap:1px;line-height:1.25}}
header .brand .t{{font-size:20px;font-weight:700;color:#fff}}
header .brand .s{{font-size:11px;color:{MUT};letter-spacing:.02em}}
header .brand .a{{font-size:11px;color:{BLUE}}}
nav{{display:flex;gap:6px;padding:10px 24px;border-bottom:1px solid {GRID};flex-wrap:wrap}}
.tab{{background:{PANEL};color:{FG};border:1px solid {GRID};border-radius:6px;padding:7px 14px;cursor:pointer;font-size:13px}}
.tab.active{{background:{BLUE};border-color:{BLUE};color:#fff}}
main{{padding:20px 24px;max-width:1200px}}
h3{{color:#fff;font-size:14px;margin:22px 0 8px;font-weight:600}}
.cards{{display:flex;gap:12px;flex-wrap:wrap;margin-bottom:8px}}
.kpi{{background:{PANEL};border:1px solid {GRID};border-radius:8px;padding:14px 18px;min-width:120px}}
.kv{{font-size:24px;font-weight:700;color:#fff}} .kl{{color:{MUT};font-size:11px;margin-top:2px}}
.grid2{{display:grid;grid-template-columns:1fr 1fr;gap:18px}}
@media(max-width:820px){{.grid2{{grid-template-columns:1fr}}}}
table{{border-collapse:collapse;width:100%;margin:6px 0 12px;font-size:12px}}
th,td{{border:1px solid {GRID};padding:5px 8px;text-align:left}} th{{background:{PANEL};color:{FG}}}
td{{color:{FG}}} tbody tr:nth-child(even){{background:#11161d}}
.mut{{color:{MUT};font-size:12px}} .warn{{color:{AMBER};font-size:12px;border-left:3px solid {AMBER};padding-left:8px}}
ul.resume{{margin:6px 0 4px;padding-left:20px}} ul.resume li{{margin:5px 0;line-height:1.5}}
ul.resume li.warn{{list-style:none;margin-left:-12px}}
code{{background:{PANEL};padding:1px 5px;border-radius:4px;color:{PURP}}}
.goctl{{display:flex;gap:16px;align-items:center;flex-wrap:wrap;margin:6px 0 10px;font-size:12px;color:{FG}}}
.goctl select{{background:{PANEL};color:{FG};border:1px solid {GRID};border-radius:5px;padding:3px 6px;font-size:12px}}
.golegend{{display:flex;flex-wrap:wrap;gap:6px 16px;margin:8px 0 4px;font-size:11px;color:{FG}}}
.golegend .legit{{display:flex;align-items:center;gap:5px}}
.golegend .sw{{width:10px;height:10px;border-radius:2px;display:inline-block}}
</style></head><body>
<header>{(f"<span class='logo'>{b['icon_svg']}</span>" if b.get('icon_svg') else (f"<img class='logo' src='data:image/png;base64,{b['icon_b64']}' alt=''/>" if b['icon_b64'] else ""))}<div class='brand'><span class='t'>{_esc(b['title'])}</span><span class='s'>{_esc(b['subtitle'])}</span><span class='a'>{_esc(b['affiliation'])}</span></div></header>
<nav>{nav}</nav>
<main>{sections}</main>
<script>
function show(i){{document.querySelectorAll('.panel').forEach(p=>p.style.display='none');
 document.getElementById('sec-'+i).style.display='block';
 document.querySelectorAll('.tab').forEach(b=>b.classList.remove('active'));
 event.target.classList.add('active');}}
</script></body></html>"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(doc)
    return path
