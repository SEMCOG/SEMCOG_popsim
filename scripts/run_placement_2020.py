import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = REPO_ROOT / "src"
for path in [REPO_ROOT, SRC_DIR]:
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from semcog_popsim.forecast_refinement.placement_runner import (
    load_buildings_from_hdf,
    load_buildings_from_sql,
    run_household_placement,
)

RUN_NUMBER = "100124_run_2020"
HDF_PATH = "/home/da/share/urbansim/RDF2050/model_inputs/base_hdf/forecast_data_input_031523.h5"
HOUSEHOLD_CSV = "/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2020(2022)/synthetic_households_ybl.csv"
PERSON_CSV = "/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2020(2022)/synthetic_persons.csv"
VOTER_REGISTRATION_CSV = "/mnt/hgfs/da/Staff/Nutting/RDF2050/Qualified Voter File/placement_2020.csv"
LOAD_FROM_HDF = False
SQL = """
      SELECT 
      urbansim_buildings.building_id,
      urbansim_buildings.parcel_id, 
      urbansim_buildings.nonres_sqft AS non_residential_sqft,
      urbansim_buildings.year_built,
      urbansim_buildings.residential_units,
      urbansim_buildings.owner_units,
      urbansim_buildings.building_type_id,
      urbansim_buildings.sqft_per_unit,
      urbansim_buildings.city_id as mcd,
      urbansim_buildings.stories,
      urbansim_buildings.market_value,
      urbansim_buildings.land_area,
      p.county_id AS county,
      substring(p.census_block_id, 1, 6)::INT tract,
      substring(p.census_block_id, 8, 3)::INT block,
      substring(p.census_block_id, 7, 1)::INT bg
  FROM urbansim_buildings
      LEFT JOIN urbansim_parcels as p ON 
          urbansim_buildings.parcel_id=p.parcel_id;
"""
CONNECTION_STRING = "postgresql://USER:PSWD@SERVER:PORT/DBNAME"


def main():
    if LOAD_FROM_HDF:
        loader = lambda: load_buildings_from_hdf(HDF_PATH)
    else:
        loader = lambda: load_buildings_from_sql(SQL, CONNECTION_STRING, HDF_PATH)

    run_household_placement(
        run_number=RUN_NUMBER,
        hdf_path=HDF_PATH,
        household_csv=HOUSEHOLD_CSV,
        person_csv=PERSON_CSV,
        voter_registration_csv=VOTER_REGISTRATION_CSV,
        output_dir=REPO_ROOT / "output" / RUN_NUMBER,
        buildings_loader=loader,
    )


if __name__ == "__main__":
    main()
