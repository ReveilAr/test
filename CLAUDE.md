# MS/MS → FBMN Molecular Networking Pipeline

## Objective
Process MS/MS data end-to-end: raw vendor files → aligned feature map → feature
identification via MS/MS spectrum matching → feature-based molecular networking
(FBMN) → unknown-molecule identity prediction using network topology →
statistical analysis. Delivered as a Snakemake-orchestrated pipeline with a
Panel-based interactive dashboard.

## Pipeline stages (Snakemake rule groups, in dependency order)
1. **Raw conversion** — ThermoRawFileParser → mzML.
2. **Preprocessing** — pyOpenMS (centroiding if not already done by the parser;
   lock-mass/internal recalibration if available).
3. **Alignment** — pyOpenMS (FeatureFinderMetabo → MetaboliteAdductDecharger for
   IIMN adduct grouping → MapAlignerPoseClustering → IDMapper) → `GNPSExport`
   TOPP tool to produce an MGF + feature quantification table in the standard
   GNPS/FBMN schema (this schema is what FBMN-STATS and pyMolNetEnhancer expect
   downstream — don't reinvent the export format).
4. **Re-quantification** — `FeatureFinderMetaboIdent` targeted re-extraction
   (gap-filling) so features triggered for MS2 in only some samples still get
   quantified everywhere.
5. **Network scoring** — matchms; primary edge score MS2DeepScore (ONNX-exported
   for inference speed), modified cosine as a fallback/complementary score.
   Filter spectra by independent QC (peak count, precursor intensity, blank
   subtraction) *before* scoring — never filter on the score itself.
6. **Network construction** — build via top-N pool → threshold + matching-peaks
   filter → top-K per node → iterative weakest-edge pruning for oversized
   components (the arteMIS/GNPS-style algorithm). **Don't hand-pick GNPS
   defaults (cutoff 0.7, top-K 10, max component 100)** — they're benchmarked
   as suboptimal; tune per-score via parameter search (Latin Hypercube Sampling
   + composite topology/chemistry score), or adapt arteMIS's own implementation
   (networkx + matchms based, MIT-compatible stack). Family assignment: Louvain
   (`nx.community.louvain_communities`, tune `resolution` upward from 1.0) or
   Leiden (`leidenalg`, stronger well-connectedness guarantee) — **not
   HDBSCAN**, which clusters embeddings directly rather than consuming graph
   edges. graph-tool is a viable CPU-only accelerator specifically for
   community detection (stochastic block model, avoids Louvain's resolution
   limit) — run it as an isolated Docker/conda process (flat edge-list in,
   community-label table out), both for GPL isolation and because it doesn't
   support Windows.
7. **Feature ID / annotation:**
   - Spectral library matching (matchms/pyOpenMS) — decide if pyOpenMS's
     matcher is redundant with MS2Query below.
   - MS2Query — exact + analog library search.
   - SIRIUS via the PySirius client (local Nightsky REST API, background
     service — not the CLI in a loop): formula (+ ZODIAC network-based
     refinement), CSI:FingerID structure candidates, CANOPUS compound class.
   - MIST-CF as an independent formula cross-check (post-SIRIUS agreement/
     disagreement flag, not a pre-SIRIUS input — no clean hook exists to
     inject external priors into SIRIUS's internal search).
   - **Network-aware rescoring** (runs *after* stage 6, not before — this was
     a dependency bug in an earlier draft): pull full CSI:FingerID candidate
     lists via the API; define high-confidence anchors per cluster (library
     hits, or SIRIUS confidence above a validated threshold; fall back to
     CANOPUS-class match where no anchor exists); combine Tanimoto/class
     consistency with CSI:FingerID's own score into a re-ranked list. Single
     pass, not iterative — guard against propagating a wrong anchor's error
     through its whole cluster. Feed the *rescored* output into
     pyMolNetEnhancer, not the raw pass.
   - SNAP-MS as an optional complementary annotator (MS1-mass-distribution
     based compound-family assignment against NPAtlas/COCONUT) — under
     consideration, not yet integrated.
8. **MS2LDA 2.0** — substructure/Mass2Motif discovery, parallel track over the
   same spectra, feeds pyMolNetEnhancer.
9. **Stats** — FBMN-STATS (has both R and Python notebooks; shape your
   feature/edge/node tables to the standard GNPS-FBMN schema to reuse its code
   with minimal adaptation). Can run in parallel with stage 10.
10. **Merge** — pyMolNetEnhancer: combine SIRIUS/CANOPUS classes + MS2LDA
    motifs onto the network's node/edge tables.
11. **Harmonize annotations** — one script reconciling spectral-library hits,
    SIRIUS/CSI:FingerID, MS2Query, CANOPUS into a single best-annotation field
    per node, tagged with a confidence tier (Schymanski 5-level scale:
    confirmed by reference standard > library-spectrum match > diagnostic
    evidence > tentative candidate > formula-only) plus annotation source.
12. **Visualization** — Panel app; pyOpenMS-viz for standard spectrum/
    chromatogram/peakmap plots (Bokeh backend); HoloViews + Datashader
    specifically for large rasterized views; `hv.Graph.from_networkx` +
    `bundle_graph`/`dynspread` for the network, standardized on Bokeh
    throughout (not Plotly) so click-driven cross-panel linking works cleanly.
13. **Optional Cytoscape export** — py4cytoscape, carrying all annotation
    layers as node/edge table columns.

## UI architecture
- A single served Panel app (`panel serve`), not a static page and not a
  notebook — notebooks are for prototyping individual visualizations only.
- Layout: sidebar (global Snakemake parameters, rule on/off profile, Run
  button, status/log) + main area `pn.Tabs` (Network / Chromatograms /
  Feature List). Annotation info is a persistent detail panel bound to a
  shared "selected node" state, not a fourth peer tab.
- Cross-tab linking via one shared reactive state object (`param.Parameterized`
  or `pn.bind`) — network tap stream, Tabulator row-selection, etc. all read/
  write the same "selected feature" value.
- Snakemake triggered via an `async def` button callback (subprocess or
  `snakemake.api.SnakemakeApi`, Snakemake 8+), streaming stdout to a log pane.
  `pn.state.add_periodic_callback` polls a lightweight status file/DuckDB
  table that Snakemake writes at the end of each major rule — don't poll the
  raw output directory. Enable `nthreads` so polling doesn't compete with the
  websocket.

## Licensing notes (checked; re-verify before any commercial/distributed use)
- **SIRIUS**: source is AGPL-3.0; CSI:FingerID/CANOPUS/MSNovelist web services
  are academic/non-commercial only (Bright Giant GmbH license needed
  otherwise). Fine for internal research use as planned.
- **graph-tool, leidenalg/python-igraph**: GPL-3. No issue for private/
  internal use; only matters if the combined codebase is redistributed —
  isolate as a separate process (already the plan above) to keep this clean.
- **matchms, MS2DeepScore**: Apache-2.0 (confirmed). Spec2Vec: same lab,
  presumed Apache-2.0. **pyMolNetEnhancer**: MIT (confirmed). **SNAP-MS**:
  code is MIT (confirmed) — but its default reference database, NPAtlas, is
  CC BY-NC 4.0 (non-commercial); COCONUT is the alternative reference DB,
  license not yet verified.
- Not yet verified — check each repo's LICENSE before depending on them:
  MS2LDA 2.0, MS2Query, DreaMS, MIST-CF, FBMN-STATS (no explicit license
  found), py4cytoscape, Cytoscape, and STEP/MT-GEM/arteMIS (very recent
  preprints — research code at that stage often has no license file yet).

## Tools under consideration, not yet committed
DreaMS (transformer embeddings, alternative/complementary edge score),
STEP/MT-GEM (very new, May 2026 preprint — discrete structural-transformation
prediction from spectrum pairs; treat as an experimental, clearly-labeled
optional module once the core pipeline is stable, not a day-one dependency),
SNAP-MS (see stage 7 above).
