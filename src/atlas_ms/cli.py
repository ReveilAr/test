"""
Command line interface.

    atlas-ms init  <project> <raw files...> [--instrument orbitrap] [--adducts positive_lipids]
    atlas-ms run   <project> [--cores N] [--dry-run] [-- extra snakemake options]
    atlas-ms app   [project] [--port 5006] [--no-browser]
    atlas-ms cache                      what ATLAS-MS keeps outside the projects, and its size
"""

import argparse
import sys
from pathlib import Path

from atlas_ms.config import load_presets
from atlas_ms.credentials import ACCOUNT_FILE
from atlas_ms.project import Project
from atlas_ms.runner import CONDA_PREFIX, available_cores, run_workflow

CACHE = Path.home() / ".cache" / "atlas-ms"


def folder_size_mb(folder: Path) -> float:
    """Total size of the files in a folder (symbolic links not followed)."""
    return sum(f.stat().st_size for f in folder.rglob("*") if f.is_file() and not f.is_symlink()) / 1e6


def show_cache() -> None:
    """
    Everything ATLAS-MS stores outside the project folders. All of it is
    shared by every project and made once (not per analysis), so nothing
    accumulates with each run.
    """
    print(f"ATLAS-MS files outside the projects ({CACHE}):")
    folders = {
        CONDA_PREFIX: "conda environments of the tools (ThermoRawFileParser, SIRIUS, MS2Query), made on first use",
        CACHE / "models": "MS2DeepScore model",
        CACHE / "ms2query": "MS2Query library and models",
    }
    for folder, what in folders.items():
        size = f"{folder_size_mb(folder):9.0f} MB" if folder.exists() else "   (none)   "
        print(f"  {size}  {folder}: {what}")
    print(f"  SIRIUS account: {'saved in ' + str(ACCOUNT_FILE) if ACCOUNT_FILE.exists() else 'none saved'}")
    # (Not Snakemake's --conda-cleanup-envs: it deletes the environments the
    # workflow uses, not the old ones.)
    print(f"\nTo free space, delete {CONDA_PREFIX} (the environments are made again when a step\n"
          f"needs them, which takes a few minutes) or all of {CACHE} (models re-downloaded too).")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="atlas-ms", description="ATLAS-MS LC-MS/MS pipeline")
    commands = parser.add_subparsers(dest="command", required=True)

    init = commands.add_parser("init", help="create a project folder for a set of raw files")
    init.add_argument("project", help="folder to create")
    init.add_argument("files", nargs="+", help="Thermo .raw or centroided .mzML files")
    init.add_argument("--instrument", default="orbitrap", choices=sorted(load_presets("instruments")))
    init.add_argument("--adducts", default="positive_lipids", choices=sorted(load_presets("adducts")))

    run = commands.add_parser("run", help="process a project")
    run.add_argument("project", help="project folder")
    run.add_argument("--cores", type=int, default=None,
                     help=f"maximum CPU cores used at once (default: all available, {available_cores()} here)")
    run.add_argument("--dry-run", action="store_true", help="only show what would be done")
    run.add_argument("snakemake_args", nargs=argparse.REMAINDER,
                     help="extra Snakemake options, after '--'")

    app = commands.add_parser("app", help="open the app in the browser")
    app.add_argument("project", nargs="?", help="project folder to open")
    app.add_argument("--port", type=int, default=5006)
    app.add_argument("--no-browser", action="store_true", help="do not open a browser tab")

    commands.add_parser("cache", help="show what ATLAS-MS keeps outside the projects")

    args = parser.parse_args(argv)
    if args.command == "cache":
        show_cache()
        return 0
    if args.command == "app":
        from atlas_ms.app.main import serve  # imported here: loading Panel takes a moment

        serve(args.project, port=args.port, show=not args.no_browser)
        return 0
    if args.command == "init":
        project = Project.create(args.project, args.files, instrument=args.instrument, adducts=args.adducts)
        print(f"Created {project.root}\n"
              f"  Edit samples.tsv (sample_type, ATTRIBUTE_ columns) and project.yaml if needed,\n"
              f"  then run: atlas-ms run {args.project}")
        return 0
    extra = [a for a in args.snakemake_args if a != "--"]
    return run_workflow(args.project, cores=args.cores, dry_run=args.dry_run, extra_args=extra)


if __name__ == "__main__":
    sys.exit(main())
