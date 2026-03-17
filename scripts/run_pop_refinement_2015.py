import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = REPO_ROOT / "src"
for path in [REPO_ROOT, SRC_DIR]:
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from semcog_popsim.forecast_refinement.refinement_runner import run_population_refinement

DEFAULT_OUTPUT_DIR = Path("/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2015(2017)/092324_run_2015")
DEFAULT_REFINEMENT_EXCEL = Path("/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/semcog_estimates/PopHHEstimates715.xlsx")
DEFAULT_HDF_PATH = Path("/mnt/hgfs/urbansim/RDF2045/data/base_year/all_semcog_data_02-02-18-final-forecast.h5")


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--refinement-excel", default=str(DEFAULT_REFINEMENT_EXCEL))
    parser.add_argument("--hdf-path", default=str(DEFAULT_HDF_PATH))
    parser.add_argument("--household-csv", default=str(DEFAULT_OUTPUT_DIR / "households.csv"))
    parser.add_argument("--person-csv", default=str(DEFAULT_OUTPUT_DIR / "persons.csv"))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--refine-geo", default="semmcd")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    run_population_refinement(
        refinement_excel=args.refinement_excel,
        hdf_path=args.hdf_path,
        household_csv=args.household_csv,
        person_csv=args.person_csv,
        output_dir=args.output_dir,
        refine_geo=args.refine_geo,
    )


if __name__ == "__main__":
    main()
    print("Done.")
