# ATLAS-MS: MS/MS → FBMN Molecular Networking Pipeline

## Objective
Process LC-MS/MS data end to end: raw vendor files → aligned feature map →
feature-based molecular networking (FBMN) → annotation with Schymanski
confidence levels (lipid-focused) → statistics. The pipeline is orchestrated by
Snakemake, with a Panel app to set up, run and explore a project.

The full design is in [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md). This file
keeps the decisions and the rules every session must follow.

## Status
- **Milestone 1 (preprocessing) implemented:** raw → aligned, gap-filled
  feature table + GNPS export. It is tested end to end on synthetic data.
- **Real-data runs** (user's bacterial lipids, 4 Thermo runs of *different
  species*, so the low overlap is expected): about 2,000–2,400 features per
  run, 4,122 linked, 1,134 found in all runs.
  - Run 1: 4,016 features after gap filling and 1,348 with MS2, but only
    1,256 spectra in the MGF.
  - Run 2, after the fixes: 4,120 features, 1,374 with MS2, 1,374 spectra
    written, 13% missing values.
  - Run 2 also showed that gap filling mixed two intensity scales, and the
    OpenMS peptide isotope model used for targets. Gap filling was reworked
    (see Preprocessing).
  - Run 3, with the reworked gap filling: about 930–1,110 detected /
    re-extracted pairs per run and 1,237–1,346 gaps filled per run. Scale
    factors (detected / re-extracted) were 1.43, 1.24, 1.10 and 1.21: the
    symmetric model under-integrates tailing peaks.
  - Intensity check (weakest vs. strongest quarter of pairs): 1.12 vs 1.14,
    1.19 vs 1.20, 1.20 vs 1.29, 1.29 vs 1.56. There is no real intensity
    dependence except in legio_19 (gaps scaled about 10% high there), so one
    factor per run is kept. The per-feature ratio IQR is wide (about 0.85–2.0):
    single gap-filled values are approximate. They are therefore flagged in
    `results/quant_gap_filled.parquet` (and `features.n_gap_filled`) for the
    statistics.
  - **Preprocessing is considered validated on real data.**
- **Milestone 2 done** (awaiting the user's first real network).
  - `atlas_ms.network`: spectrum QC → modified cosine or MS2DeepScore
    candidate pool → GNPS-style network, families, Louvain communities,
    layout, GraphML.
  - `atlas_ms.app` (`atlas-ms app`): sidebar open / create / run; Setup tab
    (samples, presets, adducts, parameters); Network tab (network ↔ table
    selection, MS2 spectrum, chromatograms from the mzML).
  - Next: tuning on the real network (milestone 5), then milestone 3
    (annotation).
- **Name:** ATLAS-MS is a placeholder (Python package `atlas_ms`, command
  `atlas-ms`). The repo is private and licensing is decided later.

## How to work in this repo
- **Code style (the user's top priority):** keep the code **simple and
  extensively commented**, and as efficient as the libraries allow.
  - Docstrings explain the science and the *why*, not only the *what*.
  - Prefer pyOpenMS/numpy/pandas built-ins over Python loops over peaks.
  - No clever abstractions.
- **Environment:** `conda env create -f environment.yml` installs the
  package in editable mode. Python dependencies are listed once, in
  `pyproject.toml`.
- **Tests:** `pytest` (about 1.5 min: matchms and PyTorch imports are slow). `tests/synthetic.py` generates the
  LC-MS runs. Every change to a processing step needs a test with an
  expected value on the synthetic study.
- **Running:** `atlas-ms init <project> <files>` then `atlas-ms run <project>`.
  Snakemake runs with `--directory <project>`, so rule paths are relative to
  the project folder. `--sdm conda` needs the `conda` command on PATH; the
  tests run without it (mzML input only).
- **After a change:** Snakemake does not re-run a step when only the
  `atlas_ms` library changed (rule, script, input and parameter changes
  do trigger it). When a fix needs a re-run, give the user the
  `atlas-ms run <project> -- --forcerun <rule>` to use. The user runs the
  code from their own clone: they must `git pull` it (merging on GitHub is
  not enough). Every rule log's first line gives the commit and folder of
  the code (`logs.code_version`).

## Decisions (from the design Q&A)
- **Platform:** Linux only for now (Windows maybe later, don't design for it
  yet). Conda environments. Runs locally on laptop/desktop, with or without a
  GPU. 10–100 raw files per project. One project open at a time.
- **Input:** instrument, polarity and adducts are parameters edited in the
  Panel app. **Thermo only for now:** ThermoRawFileParser converts `.raw`
  files and centroided `.mzML` files are linked as they are. msconvert and
  Docker were dropped. Current data is positive mode only, RP
  chromatography, with no QC-based drift correction.
- **Samples:** mostly wastewater, so anything can be present. The lipid rules
  must cover broad classes. Choosing a subset of classes is a v2 option. The
  user's own test data are **bacterial lipids**, so the rule set must include
  bacterial classes (e.g. PE, PG, cardiolipins, lyso forms, ornithine lipids,
  glycolipids), not only the mammalian ones.
- **Preprocessing:** UmetaFlow is the base, ported from its OpenMS command-line
  tools to pyOpenMS. Its step order is kept: FFM → align → decharge → IDMapper →
  link → FeatureFinderMetaboIdent gap filling → re-link → export. Two
  deviations, both fixing data loss:
  - MS2 spectra without peaks are not attached to features
    (`annotate.drop_empty_ms2`);
  - **Gap filling only fills gaps** (`gap_filling.merge_gap_filled`).
    Detected features are never replaced by re-extracted ones: re-extraction
    can pick a neighbouring isomer. Re-extracted values are added only where
    a run had no feature. They are converted to the detected scale: first
    to the monoisotopic share, then by the per-run median ratio
    detected/re-extracted over targets that have both. Targets carry their
    measured isotope pattern (FFM `masstrace_intensity`) instead of OpenMS's
    peptide model. UmetaFlow instead re-extracts incomplete features in
    every run, which loses values when extraction fails and mixes two
    intensity definitions (trace area vs. model area summed over M+M+1).
- **Network:** the user chooses modified cosine *or* MS2DeepScore in the app.
  Parameters start from placeholders and are tuned **once, after the first real
  network**. Families are connected components, with Louvain communities
  inside them. GraphML is always written, and a py4cytoscape push to Cytoscape
  is in v1.
- **Annotation:** matchms library search (user MSP/MGF/JSON, JSON = GNPS
  flavour, + public libraries),
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
- **Distribution:** academic only, no commercial distribution. The repo is
  private and licensing is decided later.

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
- **Code placement:** all logic lives in the `atlas_ms` package. Snakemake
  scripts are thin. Tools with conflicting dependencies (MS2Query, MS2LDA,
  SIRIUS client) get their own conda env and exchange files only.
- **Parameters:** they are defined once, as `param` classes in
  `atlas_ms.config`. The app renders them and Snakemake validates with the same
  classes.
- **App:** the Panel app never runs pipeline tools itself. It writes the
  project config, launches Snakemake as a subprocess and reads `results/`.
- **Visualization:** Bokeh everywhere (not Plotly) so cross-panel linking works.
  Datashader for large rasterized views such as peak maps.
- **Resources:** never commit spectral libraries, models or raw data. They are
  downloaded to a cache or referenced by path.
- **mzML files:** never write a modified copy of an mzML file. Corrections are
  small side files (`work/features/<s>.precursors.tsv` for precursor m/z,
  `work/alignment/<s>.trafoXML` for RT), re-applied in memory by
  `msdata.load_run`.
- **Feature ids:** features with MS2 come first, so
  `feature_id` = GNPS `row ID` = MGF `SCANS` for them. Features without MS2
  follow, in the internal tables only.
- **Rule parameters:** a rule's `params` must hold only values that rule uses.
  Snakemake re-runs a rule whenever one of its params changes (e.g. the
  preset names live in their own `presets` section for this reason).
- **Tests:** tests must not need network access, a SIRIUS account or GPUs.
  Mock the external tools.

## pyOpenMS 3.5 quirks (found while building milestone 1)
- **MassTraceDetection** loses most traces when some MS1 scans have no peak
  above `noise_threshold_int`. `msdata.ms1_experiment(min_intensity=...)`
  drops those scans first, which is harmless.
- **IDMapper** labels its identification run `UNKNOWN_SEARCH_RUN_IDENTIFIER`
  in every file, and the linked consensusXML then can't be saved. Pass a
  `ProteinIdentification` with a unique identifier (the sample name).
- **Getters that fill a list:** `FeatureMap.getPrimaryMSRunPath(list)` (and
  similar getters) fill a list argument instead of returning a value.
- **PrecursorCorrection** looks spectra up by RT: two MS2 scans with exactly
  the same RT get mixed up. Real data never has this; synthetic data must
  avoid it.
- **IIMN meta values** on consensus features are named `row ID`,
  `best ion`, `partners` and `annotation network number`.
- **writeSupplementaryPairTable** writes no file when there are no adduct
  pairs.
- **GNPSMGFFile** (source checked) numbers entries `SCANS = position + 1` in the
  consensus map. For each feature it takes the MS2 spectrum from the run where
  the feature is most intense, and skips the feature if that spectrum is
  empty, even when other runs have good spectra. Hence `drop_empty_ms2`.
- **FeatureFinderMetaboIdent:**
  - Without a formula or an isotope pattern, it logs "No sum formula given…
    using estimation method for peptides" once per target.
  - Its feature intensity is a fitted-model area **summed over the extracted
    isotope traces**. Each trace subordinate has `native_id` `…_i<k>` and
    `isotope_probability`.
  - "SignalToNoiseEstimatorMedian: 100% of all windows were sparse" warnings
    are harmless.
- **PyYAML** reads `1.0e4` as a string (YAML 1.1). Write `10000.0`.

## Other library notes (milestone 2)
- **matchms 0.33:** `ModifiedCosine` is now `ModifiedCosineGreedy` (and
  `ModifiedCosineHungarian`, exact but slower). `.matrix()` returns a
  structured array with fields `score` and `matches`.
- **matchms import time:** importing matchms takes about 10 s. Nothing the
  Snakefile imports may import matchms. `network/model_files.py` exists for
  this reason, and `blank_ratios` lives in `graph.py`.
- **MS2DeepScore 2.11:**
  - It imports `onnxruntime` without declaring it, so it is listed in
    `pyproject.toml`.
  - The pretrained model is `ms2deepscore_model.pt` from Zenodo record
    17826815, loaded with `load_model(path, allow_legacy=True)`.
  - Tests use a tiny untrained `SiameseSpectralModel(SettingsMS2Deepscore(base_dims=(32,), embedding_dim=8))`.
- **Panel 1.9:** widgets take `label=` (not `name=`) and buttons take `color=`
  (not `button_type=`); the old names give deprecation warnings.
- **pyopenms-viz 1.2:** `df.plot(kind="spectrum" | "chromatogram",
  backend="ms_bokeh", show_plot=False)` returns a Bokeh figure. Its
  `__version__` string still says 1.0.1.
- **Checking the app:** serve it (`atlas-ms app <project> --no-browser
  --port N`) and screenshot it with Playwright. Use the preinstalled
  Chromium (`executable_path="/opt/pw-browsers/chromium"`): the pip
  Playwright wants a newer build. Never `pkill -f "atlas-ms app"`: the
  pattern matches the shell running it.
- **Modified cosine with one matched fragment** links unrelated spectra by
  coincidence. Example in the synthetic data: TG vs CE at 0.78, through one
  shifted fragment. Keep `min_matched_peaks` above 1.
- **ThermoRawFileParser 1.4.5 (bioconda, runs on Mono):** `--input=`,
  `--output=` (a file) and `--format=2` (indexed mzML). Vendor peak picking
  is on by default.

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
- Replacing pyMolNetEnhancer with a reimplementation of its family-consensus
  logic (proposed, not yet confirmed).
- Final name and license (later).
- Instrument auto-detection from ThermoRawFileParser metadata (`--metadata`),
  to pre-select the preset in the app (idea for milestone 2).
