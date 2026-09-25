# ATLAS-MS

LC-MS/MS pipeline for feature-based molecular networking (FBMN) with
lipid-focused annotation. It is orchestrated by Snakemake, and a Panel app
is planned.

> Working name, private project. Design: [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md).
> Decisions and project rules: [`CLAUDE.md`](CLAUDE.md).

**Status: milestone 1 (preprocessing).** The pipeline goes from raw files to
an aligned, gap-filled feature table and a GNPS FBMN export. Networking,
annotation, statistics and the app come next.

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
```

Results in `my_study/results/`:

| File | Content |
|---|---|
| `features.parquet` | one row per feature: m/z, RT, charge, number of MS2 spectra, adduct groups |
| `quant.parquet` | feature intensities, one column per sample |
| `gnps/` | GNPS FBMN input in "OpenMS" format: `ms2_spectra.mgf`, `quantification_table.txt`, `metadata.tsv`, `iimn_supplementary_pairs.csv` |

Feature ids are the same everywhere. Features with MS2 are numbered first,
so `feature_id` = GNPS `row ID` = MGF `SCANS`.

## Tests

```bash
pytest
```

The tests generate small synthetic LC-MS/MS runs (`tests/synthetic.py`) and
run the real workflow on them (about 15 s). They need neither raw data nor
network access.
