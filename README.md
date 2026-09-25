# ATLAS-MS

LC-MS/MS pipeline for feature-based molecular networking (FBMN) with
lipid-focused annotation. It is orchestrated by Snakemake, and a Panel app
is planned.

> Working name, private project. Design: [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md).
> Decisions and project rules: [`CLAUDE.md`](CLAUDE.md).

**Status: milestones 1–2.** The pipeline goes from raw files to an aligned,
gap-filled feature table, a GNPS FBMN export and a molecular network
(modified cosine or MS2DeepScore). The app has Setup and Network tabs.
Annotation and statistics come next.

## Install (Linux)

```bash
conda env create -f environment.yml   # from the repository root
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
atlas-ms app my_study
```

The app's sidebar opens or creates a project and runs the pipeline. The
**Setup** tab edits the sample table (sample types, metadata columns),
presets, adducts and all parameters. The **Network** tab shows the
molecular network, colourable by family, community, intensity, gap-filled
values or spectrum QC. Clicking a node or a row shows its MS2 spectrum and
its chromatogram in every run.

The first MS2DeepScore run downloads the pretrained model (about 130 MB)
to `~/.cache/atlas-ms/models/`.

Results in `my_study/results/`:

| File | Content |
|---|---|
| `features.parquet` | one row per feature: m/z, RT, charge, number of MS2 spectra, adduct groups |
| `quant.parquet` | feature intensities, one column per sample |
| `quant_gap_filled.parquet` | same shape, `True` where the value was re-extracted by gap filling (less precise than a detected value) |
| `gnps/` | GNPS FBMN input in "OpenMS" format: `ms2_spectra.mgf`, `quantification_table.txt`, `metadata.tsv`, `iimn_supplementary_pairs.csv` |
| `network/` | `nodes.parquet` (family, community, layout, spectrum QC), `edges.parquet` (spectral and adduct edges), `network.graphml` (open in Cytoscape) |

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
