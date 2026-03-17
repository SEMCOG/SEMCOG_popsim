import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = REPO_ROOT / "src"
for path in [REPO_ROOT, SRC_DIR]:
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from semcog_popsim.forecast_refinement.refinement_runner import run_population_refinement

RUN_NUMBER = "100124_run_2020"


def main():
    output_dir = REPO_ROOT / "output" / RUN_NUMBER
    run_population_refinement(
        refinement_excel="/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/semcog_estimates/PopHHEstimate720.xlsx",
        hdf_path="/home/da/share/urbansim/RDF2050/model_inputs/base_hdf/forecast_data_input_031523.h5",
        household_csv="/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2020(2022)/100124_run_2020/households.csv",
        person_csv="/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2020(2022)/100124_run_2020/persons.csv",
        output_dir=output_dir,
    )


if __name__ == "__main__":
    main()
    print("Done.")
