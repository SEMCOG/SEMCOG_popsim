import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = REPO_ROOT / "src"
for path in [REPO_ROOT, SRC_DIR]:
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from semcog_popsim.forecast_refinement.placement_runner import (
    load_buildings_from_hdf,
    run_household_placement,
)

RUN_NUMBER = "092324_run_2020"
DEFAULT_HDF_PATH = Path("/mnt/hgfs/urbansim/RDF2045/data/base_year/all_semcog_data_02-02-18-final-forecast.h5")
DEFAULT_HOUSEHOLD_CSV = Path("/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2015(2017)/synthetic_households.csv")
DEFAULT_PERSON_CSV = Path("/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2015(2017)/synthetic_persons.csv")
DEFAULT_VOTER_REGISTRATION_CSV = Path("/mnt/hgfs/da/Staff/Nutting/RDF2050/Qualified Voter File/placement_2020.csv")
DEFAULT_OUTPUT_DIR = REPO_ROOT / "output" / RUN_NUMBER


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-number", default=RUN_NUMBER)
    parser.add_argument("--hdf-path", default=str(DEFAULT_HDF_PATH))
    parser.add_argument("--household-csv", default=str(DEFAULT_HOUSEHOLD_CSV))
    parser.add_argument("--person-csv", default=str(DEFAULT_PERSON_CSV))
    parser.add_argument("--voter-registration-csv", default=str(DEFAULT_VOTER_REGISTRATION_CSV))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    hdf_path = Path(args.hdf_path)
    run_household_placement(
        run_number=args.run_number,
        hdf_path=hdf_path,
        household_csv=args.household_csv,
        person_csv=args.person_csv,
        voter_registration_csv=args.voter_registration_csv,
        output_dir=args.output_dir,
        buildings_loader=lambda: load_buildings_from_hdf(hdf_path),
        pretransform_persons_source=lambda hdf: hdf["persons"],
    )


if __name__ == "__main__":
    main()
