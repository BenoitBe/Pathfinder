# Pathfinder

*Proteomic Exploration & Cross-Comparison Synthesis Suite*
Couche d'exploration et de synthèse au-dessus des sorties de **Deimos** (module DEP/DIA) et du module d'enrichissement GO du pipeline **Proteogen**.

Pathfinder ne recalcule aucune statistique : il lit les classeurs Excel produits en amont, les réorganise et les visualise pour donner une **vision d'ensemble** multi-contrastes, un **résumé interprétatif** automatique, et un rapport HTML autonome.

---

## Ce que fait Pathfinder

Trois axes complémentaires :

1. **Méta-GO** — croise tous les contrastes sur leurs termes GO/pathway (matrice terme × contraste de z-score signé, directionnalité croisée up/down/mixte, regroupement de redondance par Jaccard). Ventilé par source (BP / MF / CC / KEGG / Reactome).
2. **Exploration** — PCA enrichie (scree, loadings, **association PC ↔ condition** par Kruskal-Wallis + ε²), corrélation inter-échantillons, UMAP échantillons (design) et **UMAP protéines** interactive coloriable par abondance, LFC, et **terme GO** (sélecteurs source / nombre de termes).
3. **Publications** *(online, opt-in)* — enrichissement **PMID STRING** par contraste : quels articles sont sur-représentés dans la signature, avec cache daté (snapshot).

Un **résumé interprétatif déterministe** (capture du design, suspicion de batch, échantillons aberrants, synthèse fonctionnelle par source) est généré en Markdown, injecté dans le rapport HTML et l'Excel.

## Sorties

- `pathfinder_report.html` — rapport autonome (onglets Résumé, Vue d'ensemble, Méta-GO, Exploration, UMAP protéines, Publications), icône et données embarquées, aucun accès réseau à l'ouverture.
- `pathfinder_synthese.xlsx` — matrices et tables (provenance, résumé, méta-GO, PCA, associations, corrélation, publications).
- `pathfinder_resume.md` — résumé interprétatif seul.

---

## Installation

```bash
python -m venv .venv && source .venv/bin/activate      # Windows : .venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

`umap-learn` (UMAP) et `requests` (axe Publications) sont optionnels : sans eux, les modules concernés sont ignorés proprement.

## Utilisation

```bash
# Offline (auto-détection des Excel par contenu dans le dossier courant)
python pathfinder.py

# Fichiers explicites
python pathfinder.py --stats rapportstatistique.xlsx --go rapportenrichissementGO.xlsx

# Design expérimental (si les colonnes échantillon ne codent pas la condition)
python pathfinder.py --design ExperimentalDesign.csv

# Axe Publications (online, opt-in) — nécessite un réseau et l'espèce
python pathfinder.py --online --org rnorvegicus        # ou hsapiens, mmusculus… ou un taxid
python pathfinder.py --online --org 10116 --refresh    # ignorer le cache

# Diagnostic (colonnes/sheets détectées, tailles de groupe, matrice PCA retenue)
python pathfinder.py --diagnose
```

### Principales options

| Option | Rôle |
|---|---|
| `--stats` | Classeur stats (auto-détecté par contenu si absent) |
| `--go` | Classeur d'enrichissement GO (peut être le même classeur) |
| `--design` | Mapping échantillon → condition (CSV/xlsx) |
| `--outdir` | Répertoire de sortie |
| `--pthr` | Seuil p pour masquer la matrice GO (défaut 0.05) |
| `--jaccard` | Distance de coupe pour la redondance GO (défaut 0.5) |
| `--online` / `--org` / `--refresh` | Axe Publications (STRING PMID) + cache |
| `--no-umap` | Désactive les UMAP |
| `--diagnose` | Inspecte la détection puis quitte |

### Entrées attendues

- **Stats** : un classeur dont une sheet porte des colonnes `<contraste>_p.val` (+ `_diff`, `_p.adj`, `Pi_Score_*`…) et une matrice d'intensités par échantillon (sheet dédiée type `Z-score_Heatmap` ou colonnes d'intensité).
- **GO** : des sheets `GO_<abrégé>` avec `term_name`, `source`, `p_value`, `intersection` (accessions) et un z-score.
- **Design** (optionnel) : colonnes identifiant échantillon (`sample`/`label`/`run`/`file`…) + condition (`condition`/`group`/`class`…).

La détection des fichiers et colonnes se fait **par contenu**, tolérante au nommage.

---

## Exemple reproductible

```bash
cd examples && python make_synthetic.py        # génère rapportstatistique.xlsx + rapportenrichissementGO.xlsx
python ../pathfinder.py --stats rapportstatistique.xlsx --go rapportenrichissementGO.xlsx
```

## Icône

L'icône (rover *Sojourner* sur sol martien, dans le style des lunes Deimos/Phobos) est générée par script :

```bash
python generate_icon.py        # écrit dans assets/
```

---

## Structure

```
pathfinder.py          orchestration + CLI
io_layer.py            lecture Excel + design → objets communs
axis_meta_go.py        axe Méta-GO
axis_explore.py        axe Exploration (PCA/UMAP)
axis_publications.py   axe Publications (STRING PMID, online)
summary.py             résumé interprétatif déterministe
species.py             résolution org gprofiler → taxid NCBI
cache.py               cache disque des appels online
report/build_html.py   rapport HTML autonome
report/build_xlsx.py   classeur de synthèse
generate_icon.py       génération de l'icône
```

## Notes méthodologiques

- Aucune statistique n'est recalculée ; tout dérive des Excel amont.
- La corrélation QC est calculée sur les intensités **brutes** (les z-scores sont décorrélés par construction).
- L'association PC ↔ condition requiert des réplicats ; sans design résolu ou sans réplicats, le test est signalé non calculable plutôt que faussement interprété.
- L'axe Publications est reproductible via cache + date de fetch, malgré l'évolution du corpus littérature. Biais connu : sur-représentation des protéines très étudiées (gènes d'entrée triés par |LFC| pour l'atténuer).

## Nommage

Pathfinder complète les briques **Deimos** et **Phobos** (les deux lunes de Mars) : la sonde *Mars Pathfinder* explore la surface — ici, ce que les modules amont ont produit.
