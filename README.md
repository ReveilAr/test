# ATLAS-MS

LC-MS/MS pipeline for feature-based molecular networking (FBMN) with
lipid-focused annotation. It is orchestrated by Snakemake, and a Panel app
is planned.

> Working name, private project. Design: [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md).
> Decisions and project rules: [`CLAUDE.md`](CLAUDE.md).

**Status: milestones 1–3.** The pipeline goes from raw files to an aligned,
gap-filled feature table, a GNPS FBMN export, a molecular network (modified
cosine or MS2DeepScore) and annotations with Schymanski confidence levels:
spectral library search, rule-based lipid annotation, and optionally
SIRIUS 6 and MS2Query, combined into one best annotation per feature. The
app has Setup, Network and Annotation tabs. Statistics come next.

## Install (Linux)

```bash
git clone https://github.com/ReveilAr/test.git atlas-ms   # a clone, not a ZIP: updates are a `git pull`
cd atlas-ms
conda env create -f environment.yml
conda activate atlas-ms
```

The package is installed in editable mode, so the `workflow/` folder is found
next to it. Tools with conflicting dependencies (e.g. ThermoRawFileParser)
get their own conda environments. Snakemake creates them on first use, in
`~/.cache/atlas-ms/conda`.

## Use

```bash
# 1. Create a project from Thermo .raw (or centroided .mzML) files
atlas-ms init my_study data/*.raw --instrument orbitrap --adducts positive_lipids

# 2. Optional: edit my_study/samples.tsv (sample_type = sample/blank/standard,
#    add ATTRIBUTE_* metadata columns) and my_study/project.yaml (parameters)

# 3. Process
atlas-ms run my_study --cores 4

# Or do all of it in the app (opens in the browser)
atlas-ms app my_study   # port 5006, or the next free one if it is taken
```

The app's sidebar opens or creates a project and runs the pipeline. The
**Setup** tab edits the sample table (sample types, metadata columns),
presets, adducts and all parameters. The **Network** tab shows the
molecular network (scroll to zoom; nodes never overlap), coloured by family,
community, retention time, intensity, gap-filled values or spectrum QC, and
sized by mean intensity or by the intensity in one sample. Clicking a node
or a row shows its MS2 spectrum and its chromatogram in every run, to the
right of the network. Two selected features (shift-click, or ctrl-click in
the table) give a mirror plot. The MS2 search circles the features whose
spectrum holds given fragments or neutral losses. The feature table has a
search box above every column, and shows each feature's best annotation
and its confidence level (the network can be coloured by level or by lipid
class). The **Annotation** tab lists every candidate of the selected
feature, from every source, with its evidence: a mirror plot against the
library spectrum, or the diagnostic lipid ions marked on the spectrum.

**Spectral libraries** (MSP, MGF or GNPS JSON) are added in the Setup tab.
Mark a library as *reference standards* only if its spectra were measured
on your own method, with retention times: only those can give level 1.
In-silico libraries (e.g. LipidBlast) give at most level 3.

**SIRIUS 6** (formulas, El Gordo lipids, CSI:FingerID structures, CANOPUS
classes) and **MS2Query** (analog search) are switched on in the Setup
tab's parameters (the SIRIUS and MS2Query cards). The pipeline installs both in their
own conda environments on first use.

- SIRIUS needs a free academic account: sign in once in the SIRIUS 6
  GUI (or `sirius login`). The pipeline uses a SIRIUS that is already
  running, or starts its own. The SIRIUS project is kept in
  `work/sirius/project.sirius` for the GUI. Expect hours for thousands of
  features.
- MS2Query downloads its library and models (a few GB) once, to
  `~/.cache/atlas-ms/ms2query/`.

The first MS2DeepScore run downloads the pretrained model (about 130 MB)
to `~/.cache/atlas-ms/models/`.

Results in `my_study/results/`:

| File | Content |
|---|---|
| `features.parquet` | one row per feature: m/z, RT, charge, number of MS2 spectra, adduct groups |
| `quant.parquet` | feature intensities, one column per sample |
| `quant_gap_filled.parquet` | same shape, `True` where the value was re-extracted by gap filling (less precise than a detected value) |
| `gnps/` | GNPS FBMN input in "OpenMS" format: `ms2_spectra.mgf`, `quantification_table.txt`, `metadata.tsv`, `iimn_supplementary_pairs.csv` |
| `network/` | `nodes.parquet` (family, community, layout, spectrum QC), `edges.parquet` (spectral and adduct edges), `network.graphml` (open in Cytoscape; with the best annotations) |
| `annotations/` | `library.parquet`, `lipid_rules.parquet`, `sirius.parquet`, `ms2query.parquet` (candidates of each source), `candidates.parquet` (all, with final levels), `best.parquet` (one best annotation per feature, flags, family class consensus) |

Feature ids are the same everywhere. Features with MS2 are numbered first,
so `feature_id` = GNPS `row ID` = MGF `SCANS`.

## Update

```bash
cd <the clone the atlas-ms environment was created from>
git pull
python -c "import atlas_ms; print(atlas_ms.__file__)"   # check which clone is used
```

The package is installed in editable mode, so it always runs the code of
that clone. Every log file in `<project>/logs/` starts with the version,
commit and folder of the code that wrote it.

Snakemake re-runs a step when its rule, script, inputs or parameters
change, but not when only the `atlas_ms` library changes. To redo a step
and everything after it: `atlas-ms run my_study -- --forcerun <rule>`,
e.g. `--forcerun plan_gap_filling`.

## Tests

```bash
pytest
```

The tests generate small synthetic LC-MS/MS runs (`tests/synthetic.py`) and
run the real workflow on them (about 1.5 min: matchms and PyTorch load slowly). They need neither raw data nor
network access.
