"""
Command line interface.

    atlas-ms init <project> <raw files...> [--instrument orbitrap] [--adducts positive_lipids]
    atlas-ms run  <project> [--cores 4] [--dry-run] [-- extra snakemake options]
"""

import argparse
import sys

from atlas_ms.config import load_presets
from atlas_ms.project import Project
from atlas_ms.runner import run_workflow


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
    run.add_argument("--cores", type=int, default=4, help="CPU cores Snakemake may use (default 4)")
    run.add_argument("--dry-run", action="store_true", help="only show what would be done")
    run.add_argument("snakemake_args", nargs=argparse.REMAINDER,
                     help="extra Snakemake options, after '--'")

    args = parser.parse_args(argv)
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
