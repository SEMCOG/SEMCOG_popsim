import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = REPO_ROOT / "src"
for path in [REPO_ROOT, SRC_DIR]:
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from semcog_popsim.forecast_refinement.refinement_runner import run_population_refinement


def main():
    output_dir = Path("/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2015(2017)/092324_run_2015")
    run_population_refinement(
        refinement_excel="/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/semcog_estimates/PopHHEstimates715.xlsx",
        hdf_path="/mnt/hgfs/urbansim/RDF2045/data/base_year/all_semcog_data_02-02-18-final-forecast.h5",
        household_csv=output_dir / "households.csv",
        person_csv=output_dir / "persons.csv",
        output_dir=output_dir,
    )


if __name__ == "__main__":
    main()
    print("Done.")
