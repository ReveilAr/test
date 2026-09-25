# MS/MS → FBMN Molecular Networking Pipeline

## Objective
Process LC-MS/MS data end to end: raw vendor files → aligned feature map →
feature-based molecular networking (FBMN) → annotation with Schymanski
confidence levels (lipid-focused) → statistics. The pipeline is orchestrated by
Snakemake, with a Panel app to set up, run and explore a project.

The full design is in [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md). This file
keeps the decisions and the rules every session must follow.

## Status
Architecture agreed on 2026-09-25 (draft under review). The dataset comes after
the architecture is validated. No code yet. `fbmnflow` is a working name for
the package.

## Decisions (from the design Q&A)
- **Platform:** Linux only for now (Windows maybe later, don't design for it
  yet). Conda environments. Runs locally on laptop/desktop, with or without a
  GPU. 10–100 raw files per project. One project open at a time.
- **Input:** instrument, polarity and adducts are parameters edited in the
  Panel app. Conversion uses ThermoRawFileParser for Thermo files and
  ProteoWizard msconvert (Docker) for the other vendors. Current data is
  positive mode only, with no QC-based drift correction.
- **Preprocessing:** UmetaFlow is the base, ported from its OpenMS command-line
  tools to pyOpenMS. Its step order is kept: FFM → align → decharge → IDMapper →
  link → FeatureFinderMetaboIdent gap filling → re-link → export.
- **Network:** the user chooses modified cosine *or* MS2DeepScore in the app.
  Parameters start from placeholders and are tuned **once, after the first real
  network**. Families are connected components, with Louvain communities
  inside them. GraphML is always written, and a py4cytoscape push to Cytoscape
  is in v1.
- **Annotation:** matchms library search (user MSP/MGF/JSON + public libraries),
  MS2Query, SIRIUS 6 through its **REST API** (PySirius), a lipid module,
  MS2LDA 2.0. Pretrained models only.
- **Confidence (Schymanski):** Level 1 only from libraries flagged as reference
  standards (with RT). CSI:FingerID is **capped at Level 3** whatever its COSMIC
  confidence. Formula disagreement (SIRIUS vs MIST-CF, in v2) is a **flag, not
  a demotion**. No manual curation in v1. In-silico library hits are at most
  Level 3. Lipids also carry a Liebisch structural level.
- **Stats:** a FBMN-STATS port (BSD-3). Cleanup (blank removal, imputation,
  normalization) runs in Snakemake; tests and plots are interactive in the
  app. Metadata is entered manually in the app.
- **App tabs:** Setup | Network (clickable network + pyOpenMS-viz plots +
  feature list with the top-1 annotation) | Annotation (deep dive on the
  selected node or cluster: candidates, motifs) | Statistics. A sidebar holds
  the parameters, Run and the log.
- **Deferred to v2:** MIST-CF, network-aware rescoring, manual curation,
  graph-tool/Leiden, DreaMS, SNAP-MS, STEP/MT-GEM, negative mode, and the
  PubChem / LOTUS / NPAtlas / COCONUT / GNPS2 plugins.
- **Distribution:** academic only, no commercial distribution.

## Rules to keep
- **Formats:** don't reinvent the export format. The GNPS FBMN schema
  (quant table + MGF + metadata) is the contract that FBMN-STATS-style code and
  Cytoscape expect. Internal Parquet tables mirror it.
- **Filtering:** filter spectra with independent QC (peak count, precursor
  intensity, blank ratio) *before* scoring. Never filter on the similarity
  score itself.
- **Rescoring (v2):** network-aware rescoring runs *after* network
  construction, in a single pass (not iterative), and its rescored output is
  what feeds the harmonization.
- **Confidence levels:** harmonization may lower a plugin's proposed level but
  never raise it. Disagreements become flags; candidates are never silently
  dropped.
- **Code placement:** all logic lives in the `fbmnflow` package. Snakemake
  scripts are thin. Tools with conflicting dependencies (MS2Query, MS2LDA,
  SIRIUS client) get their own conda env and exchange files only.
- **Parameters:** they are defined once, as `param` classes in
  `fbmnflow.config`. The app renders them and Snakemake validates with the same
  classes.
- **App:** the Panel app never runs pipeline tools itself. It writes the
  project config, launches Snakemake as a subprocess and reads `results/`.
- **Visualization:** Bokeh everywhere (not Plotly) so cross-panel linking works.
  Datashader for large rasterized views such as peak maps.
- **Resources:** never commit spectral libraries, models or raw data. They are
  downloaded to a cache or referenced by path.
- **Tests:** tests must not need network access, a SIRIUS account or GPUs.
  Mock the external tools.

## Verified facts (2026-09)
- **UmetaFlow** (`biosustain/snakemake_UmetaFlow`, Apache-2.0): Linux/macOS
  only, calls OpenMS 3.2 command-line tools (not pyOpenMS), ships SIRIUS 5.8.6,
  last commit January 2025.
- **Conflicting dependencies:** ms2query 1.5.4 pins matchms ≤ 0.26.4,
  ms2deepscore == 2.0.0 and torch < 2.6. ms2lda 2.0.1 needs matchms ≥ 0.27 and
  Python 3.11–3.12. pyopenms-viz 1.2 needs Python ≥ 3.12. Result: core env is
  Python 3.12, and MS2Query and MS2LDA get their own envs.
- **SIRIUS:** 6.5.4 is current. The PySirius client (`py-sirius-ms`,
  conda-forge) is versioned together with SIRIUS. The API exposes El Gordo
  lipid annotation (`LipidAnnotation`), lipid-class tagging of CSI:FingerID
  candidates, custom structure databases and spectral library search.
- **MIST-CF:** frozen since 2023 (Python 3.8, PyTorch 1.9, CUDA 11.1). The
  public model is trained on NPLIB1 only; the NIST20 model needs a NIST licence
  (not available to us).
- **pyMolNetEnhancer:** last release 0.1.9 (2019), built for GNPS1 outputs.

## Licensing (checked 2026-09; re-verify before any redistribution)
| Tool | License |
|---|---|
| pyOpenMS, pyopenms-viz | BSD-3-Clause |
| FBMN-STATS | BSD-3-Clause |
| matchms, MS2DeepScore, MS2Query | Apache-2.0 |
| UmetaFlow, PySirius client | Apache-2.0 |
| MS2LDA 2.0, MIST-CF, pyMolNetEnhancer, pygoslin, py4cytoscape | MIT |
| leidenalg | GPL-3.0-or-later (keep optional, in its own process) |
| graph-tool | GPL-3 per the original summary, but upstream may have moved to LGPL-3: to verify (v2 only, isolated process) |
| SIRIUS | AGPL-3.0 source. CSI:FingerID/CANOPUS/MSNovelist web services are academic/non-commercial only (Bright Giant license otherwise). |
| NPAtlas | CC BY-NC 4.0 (data) |
| MassBank | licenses set per record, some NC |
| LipidBlast / MS-DIAL lipid libraries, LIPID MAPS LMSD, COCONUT | to verify before bundling |

The code licenses above all allow commercial use. What restricts commercial
use is the SIRIUS web-service terms and the NC data licenses.

## Open items
- JSON library flavour (GNPS, MoNA or MassBank JSON).
- Chromatography (RP assumed for the lipid RT/ECN model) and the sample
  matrix, which decide the lipid classes to prioritize.
- Replacing pyMolNetEnhancer with a reimplementation of its
  family-consensus logic (proposed).
- Final package name and the license of our own code (BSD-3-Clause proposed).
- Docker on the user's machines (needed for msconvert on Linux).
