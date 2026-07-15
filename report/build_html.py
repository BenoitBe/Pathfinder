"""
Rapport HTML autonome Pathfinder.

Rendu 100% offline : figures générées en SVG côté Python (pas de CDN, pas de JS
de charting). Palette cohérente avec le dashboard du pipeline.
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
SRC_LABEL = {"GO:BP": "GO:BP — Processus biologique", "GO:MF": "GO:MF — Fonction moléculaire",
             "GO:CC": "GO:CC — Composant cellulaire", "KEGG": "KEGG", "REAC": "Reactome"}


def _sources_present(meta):
    present = [s for s in SRC_ORDER if (meta["source"] == s).any()]
    others = [s for s in meta["source"].dropna().unique() if s not in SRC_ORDER]
    return present + sorted(others)


def _esc(s):
    return html.escape(str(s))


def _diverging(v, vmax):
    """z signé -> couleur (rouge=up, bleu=down), NaN -> gris foncé."""
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
def svg_meta_heatmap(matrix: pd.DataFrame, meta: pd.DataFrame, max_terms=40):
    if matrix.empty:
        return "<p class='mut'>Aucun terme GO significatif.</p>"
    m = matrix.reindex(meta.sort_values("n_sig", ascending=False).index).head(max_terms)
    names = [meta.loc[k, "term_name"] for k in m.index]
    contrasts = list(m.columns)
    vmax = np.nanmax(np.abs(m.values)) or 1.0
    cw, ch = 46, 20
    left, top = 360, 150
    W = left + len(contrasts) * cw + 30
    H = top + len(m) * ch + 20
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
    # légende
    ly = H - 6
    s.append(f"<text x='{left}' y='{ly}' fill='{MUT}' font-size='10'>bleu = down · rouge = up · gris = non-sig  (|z|max={vmax:.1f})</text>")
    s.append("</svg>")
    return "".join(s)


def svg_scatter(df, xcol, ycol, cat_col, title_x, title_y, size=(520, 420)):
    W, Hh = size
    xs = df[xcol].astype(float).values
    ys = df[ycol].astype(float).values
    cats = list(dict.fromkeys(df[cat_col].tolist())) if cat_col else [None]
    cmap = {c: COND_PALETTE[i % len(COND_PALETTE)] for i, c in enumerate(cats)}

    # gouttière de légende : le nuage se rétracte, la légende ne le recouvre jamais
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
    # légende dans la gouttière
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
    # barre de légende
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
    # points gris d'abord, colorés au-dessus
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
    """Widget interactif : UMAP protéines recoloriée par terme GO, avec sélecteurs
    source et nombre de termes. Données embarquées, rendu canvas côté client."""
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

    # tokens par protéine + appartenance (indices de termes)
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
  <label>Source <select id='go_src'><option value='ALL'>Toutes</option>__SRC__</select></label>
  <label>Termes <select id='go_k'><option>4</option><option selected>8</option><option>12</option><option>16</option><option>24</option></select></label>
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
     d.innerHTML="<span class='sw' style='background:"+PAL[j%PAL.length]+"'></span>"+T[ti].n+" ["+T[ti].s+"] · "+T[ti].k+" contr."; leg.appendChild(d); });
   document.getElementById('go_info').textContent=ncol+' / '+C.length+' protéines coloriées';
 }
 document.getElementById('go_src').addEventListener('change',draw);
 document.getElementById('go_k').addEventListener('change',draw);
 draw();
})();
</script>"""
    return tmpl.replace("__SRC__", src_opts).replace("__DATA__", _json.dumps(data))


import re  # noqa: E402 (utilisé par _go_overlay)


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
    s.append(f"<text x='{W/2}' y='16' fill='{FG}' font-size='11' text-anchor='middle'>Variance expliquée (%)</text>")
    s.append("</svg>")
    return "".join(s)


def svg_corr(corr: pd.DataFrame):
    labels = list(corr.columns)
    n = len(labels); cell = min(28, max(10, int(560 / max(n, 1))))
    left, top = 130, 130
    W = left + n * cell + 20; H = top + n * cell + 20
    s = [f"<svg viewBox='0 0 {W} {H}' width='100%' style='max-width:{W}px'>"]
    for i, ri in enumerate(labels):
        s.append(f"<text x='{left-6}' y='{top+i*cell+cell*0.7}' fill='{FG}' font-size='9' text-anchor='end'>{_esc(ri)}</text>")
        x = left + i * cell + cell / 2
        s.append(f"<text x='{x}' y='{top-6}' fill='{MUT}' font-size='9' text-anchor='end' transform='rotate(-50 {x} {top-6})'>{_esc(ri)}</text>")
    for i, ri in enumerate(labels):
        for j, cj in enumerate(labels):
            v = float(corr.iloc[i, j])
            t = max(0.0, min(1.0, (v - 0.7) / 0.3)) if v >= 0 else 0
            g = int(40 + t * 180)
            col = f"rgb(30,{g},60)"
            s.append(f"<rect x='{left+j*cell}' y='{top+i*cell}' width='{cell-1}' height='{cell-1}' fill='{col}'>"
                     f"<title>{_esc(ri)} / {_esc(cj)} : r={v:.3f}</title></rect>")
    s.append(f"<text x='{left}' y='{H-4}' fill='{MUT}' font-size='9'>Pearson r (log₂ LFQ), ordonné par clustering</text>")
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
        return "<p class='mut'>Résumé non disponible.</p>"
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
         "subtitle": "Proteomic Exploration & Cross-Comparison Synthesis Suite",
         "affiliation": "Proteogen · Université de Caen",
         "icon_b64": None}
    if branding:
        b.update({k: v for k, v in branding.items() if v is not None})
    tabs = ["Résumé", "Vue d'ensemble", "Méta-GO", "Exploration", "UMAP protéines"]
    # KPI cards
    kcards = "".join(
        f"<div class='kpi'><div class='kv'>{_esc(v)}</div><div class='kl'>{_esc(k)}</div></div>"
        for k, v in kpis.items())

    # Vue d'ensemble
    prov_tbl = _table(pd.DataFrame(list(provenance.items()), columns=["clé", "valeur"]), max_rows=30)
    overview = f"<div class='cards'>{kcards}</div><h3>Provenance</h3>{prov_tbl}"

    # Méta-GO — ventilé par source
    matrix, meta, direction = meta_go["matrix"], meta_go["meta"], meta_go["directionality"]
    if matrix.empty:
        metago = "<p class='mut'>Enrichissement GO non fourni — synthèse fonctionnelle non disponible.</p>"
    else:
        n_clusters = meta["redundancy_cluster"].nunique() if (meta["redundancy_cluster"] >= 0).any() else 0
        parts = [f"<p class='mut'>{matrix.shape[0]} termes significatifs sur {matrix.shape[1]} "
                 f"contraste(s) · {n_clusters} groupes de redondance (Jaccard). "
                 f"Rouge = up, bleu = down, gris = non-sig.</p>"]
        for src in _sources_present(meta):
            keys = meta.index[meta["source"] == src]
            sub_m = matrix.loc[keys]
            sub_meta = meta.loc[keys]
            sub_dir = direction[direction["source"] == src].drop(columns=["source"], errors="ignore")
            hm = svg_meta_heatmap(sub_m, sub_meta)
            dtbl = _table(sub_dir, max_rows=25)
            parts.append(f"<h3>{_esc(SRC_LABEL.get(src, src))} "
                         f"<span class='mut'>· {len(keys)} termes</span></h3>{hm}{dtbl}")
        metago = "".join(parts)

    # Exploration
    pca = explore["pca"]; assoc = explore["pc_factor"]
    scatter = svg_scatter(pca["scores"], "PC1", "PC2", "condition",
                          f"PC1 ({pca['var_exp'][0]:.1f}%)",
                          f"PC2 ({pca['var_exp'][1]:.1f}%)" if len(pca['var_exp']) > 1 else "PC2")
    scree = svg_scree(pca["var_exp"])
    corr = svg_corr(explore["correlation"])
    assoc_tbl = _table(assoc, max_rows=10)
    # top loadings PC1
    load = pca["loadings"]
    load_tbl = ""
    if "PC1" in load.columns:
        top = load.reindex(load["PC1"].abs().sort_values(ascending=False).index).head(15)
        load_tbl = _table(top[["PC1"] + [c for c in ["PC2", "PC3"] if c in top.columns]].reset_index()
                          .rename(columns={"index": "protéine"}), max_rows=15)
    # association vide (design non résolu) -> note au lieu d'un tableau vide
    if assoc["p_value"].isna().all():
        assoc_block = ("<p class='warn'>⚠ Test PC↔condition non calculable : conditions non résolues "
                       "(chaque échantillon isolé) ou sans réplicats. Fournir un design "
                       "(<code>--design</code>) pour activer ce test.</p>")
    else:
        assoc_block = _table(assoc, max_rows=10)

    # UMAP échantillons (design)
    umap_s_html = ""
    us = explore.get("umap_samples")
    if us:
        if us.get("df") is not None:
            umap_s_html = ("<h3>UMAP échantillons</h3>"
                           + svg_scatter(us["df"], "UMAP1", "UMAP2", "condition", "UMAP1", "UMAP2"))
        elif us.get("warning"):
            umap_s_html = f"<p class='mut'>UMAP échantillons : {_esc(us['warning'])}</p>"

    # multi-facteurs : quel facteur porte quel axe
    multi = explore.get("pc_factor_multi")
    multi_block = ""
    if multi is not None and not multi.empty:
        sig = multi[(multi["p_adj_BH"] < 0.1) | (multi["effet"] >= 0.14)]
        show = (sig if len(sig) else multi).sort_values(["effet"], ascending=False)
        multi_block = ("<h3>Quel facteur porte quel axe ? <span class='mut'>· PC × facteurs du design</span></h3>"
                       + _table(show, max_rows=24)
                       + "<p class='mut'>Catégoriel : Kruskal-Wallis, effet = ε². Numérique (ex. ordre "
                         "d'injection) : Spearman, effet = |rho|. Un facteur technique associé à une PC "
                         "de forte variance ⇒ effet batch / dérive à corriger.</p>")

    explore_html = (
        f"<div class='grid2'><div><h3>PCA — PC1 vs PC2</h3>{scatter}</div>"
        f"<div><h3>Scree</h3>{scree}</div></div>"
        f"<h3>Association PC ↔ condition (Kruskal-Wallis, ε² = taille d'effet)</h3>{assoc_block}"
        f"<p class='mut'>ε² élevé + p faible sur une PC ⇒ cet axe capture le facteur biologique ; "
        f"sur une PC non désirée, suspecter un effet batch.</p>"
        f"{multi_block}"
        f"<div class='grid2'><div><h3>Corrélation inter-échantillons</h3>{corr}</div>"
        f"<div><h3>Top contributeurs PC1</h3>{load_tbl}</div></div>"
        f"{umap_s_html}")

    # UMAP protéines (onglet dédié)
    up = explore.get("umap_proteins")
    if up and up.get("df") is not None:
        udf = up["df"]
        blocks = ["<p class='mut'>Chaque point = une protéine, positionnée selon la similarité de "
                  "son profil d'intensité entre échantillons. Le coloriage révèle la structure : "
                  "zones de forte/faible abondance, et régions régulées par contraste.</p>"]
        blocks.append("<h3>Coloriage par abondance moyenne (log2)</h3>"
                      + svg_scatter_value(udf, "UMAP1", "UMAP2", "abondance", diverging=False))
        for c in up.get("contrasts_lfc", [])[:3]:
            blocks.append(f"<h3>Coloriage par LFC — {_esc(c)} <span class='mut'>(rouge=up, bleu=down)</span></h3>"
                          + svg_scatter_value(udf, "UMAP1", "UMAP2", f"LFC:{c}", diverging=True))
        # coloriage par terme GO — widget interactif (source + nombre de termes)
        if not meta_go["matrix"].empty and meta_go.get("term_genes"):
            widget = _go_umap_widget(udf, meta_go["meta"], meta_go["term_genes"])
            if widget:
                blocks.append("<h3>Coloriage par terme GO <span class='mut'>· interactif</span></h3>"
                              + widget
                              + "<p class='mut'>Choisis la source (ex. GO:CC pour les compartiments) et le "
                                "nombre de termes. Une fonction concentrée dans une région ⇒ co-régulation spatiale.</p>")
        umap_p_html = "".join(blocks)
    elif up and up.get("warning"):
        umap_p_html = f"<p class='mut'>{_esc(up['warning'])}</p>"
    else:
        umap_p_html = "<p class='mut'>UMAP protéines indisponible.</p>"

    panels = {"Résumé": _resume_html(resume), "Vue d'ensemble": overview,
              "Méta-GO": metago, "Exploration": explore_html, "UMAP protéines": umap_p_html}

    # Publications (axe A) — onglet ajouté seulement si présent
    if pubs is not None and any(pubs.get("results", {}).values()):
        pub_parts = ["<p class='mut'>Articles sur-représentés dans la signature de chaque contraste "
                     "(enrichissement PMID STRING, classé par FDR). Biais connu : les protéines très "
                     "étudiées ressortent plus — à croiser avec le LFC.</p>"]
        for c, rows in pubs["results"].items():
            if not rows:
                continue
            trs = []
            for r in rows:
                link = f"<a href='{_esc(r.get('url',''))}' target='_blank' style='color:{BLUE}'>PMID:{_esc(r.get('pmid',''))}</a>"
                fdr = r.get("fdr")
                fdr_s = f"{fdr:.1e}" if isinstance(fdr, float) and fdr == fdr else ""
                trs.append(f"<tr><td>{link}</td><td>{_esc(r.get('titre',''))}</td>"
                           f"<td>{_esc(r.get('n_genes',''))}</td><td>{fdr_s}</td></tr>")
            pub_parts.append(f"<h3>{_esc(c)}</h3><table><thead><tr><th>PMID</th><th>Titre</th>"
                             f"<th>#gènes</th><th>FDR</th></tr></thead><tbody>{''.join(trs)}</tbody></table>")
        panels["Publications"] = "".join(pub_parts)
        tabs = tabs + ["Publications"]
    nav = "".join(f"<button class='tab{' active' if i==0 else ''}' onclick=\"show({i})\">{_esc(t)}</button>"
                  for i, t in enumerate(tabs))
    sections = "".join(
        f"<section id='sec-{i}' class='panel' style='display:{'block' if i==0 else 'none'}'>{panels[t]}</section>"
        for i, t in enumerate(tabs))

    doc = f"""<!doctype html><html lang='fr'><head><meta charset='utf-8'>
<meta name='viewport' content='width=device-width, initial-scale=1'>
<title>Pathfinder — Synthèse exploratoire</title>
<style>
*{{box-sizing:border-box}} body{{margin:0;background:{BG};color:{FG};
 font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;font-size:14px}}
header{{padding:14px 22px;border-bottom:1px solid {GRID};display:flex;align-items:center;gap:14px}}
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
<header>{f"<img class='logo' src='data:image/png;base64,{b['icon_b64']}' alt=''/>" if b['icon_b64'] else ""}<div class='brand'><span class='t'>{_esc(b['title'])}</span><span class='s'>{_esc(b['subtitle'])}</span><span class='a'>{_esc(b['affiliation'])}</span></div></header>
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
