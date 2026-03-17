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
HDF_PATH = "/mnt/hgfs/urbansim/RDF2045/data/base_year/all_semcog_data_02-02-18-final-forecast.h5"
HOUSEHOLD_CSV = "/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2015(2017)/synthetic_households.csv"
PERSON_CSV = "/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2015(2017)/synthetic_persons.csv"
VOTER_REGISTRATION_CSV = "/mnt/hgfs/da/Staff/Nutting/RDF2050/Qualified Voter File/placement_2020.csv"


def main():
    run_household_placement(
        run_number=RUN_NUMBER,
        hdf_path=HDF_PATH,
        household_csv=HOUSEHOLD_CSV,
        person_csv=PERSON_CSV,
        voter_registration_csv=VOTER_REGISTRATION_CSV,
        output_dir=REPO_ROOT / "output" / RUN_NUMBER,
        buildings_loader=lambda: load_buildings_from_hdf(HDF_PATH),
        pretransform_persons_source=lambda hdf: hdf["persons"],
    )


if __name__ == "__main__":
    main()
