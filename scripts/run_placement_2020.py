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
    load_buildings_from_sql,
    run_household_placement,
)

RUN_NUMBER = "100124_run_2020"
DEFAULT_HDF_PATH = Path("/home/da/share/urbansim/RDF2050/model_inputs/base_hdf/forecast_data_input_031523.h5")
DEFAULT_HOUSEHOLD_CSV = Path("/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2020(2022)/synthetic_households_ybl.csv")
DEFAULT_PERSON_CSV = Path("/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2020(2022)/synthetic_persons.csv")
DEFAULT_VOTER_REGISTRATION_CSV = Path("/mnt/hgfs/da/Staff/Nutting/RDF2050/Qualified Voter File/placement_2020.csv")
DEFAULT_OUTPUT_DIR = REPO_ROOT / "output" / RUN_NUMBER
DEFAULT_SQL = """
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
DEFAULT_CONNECTION_STRING = "postgresql://USER:PSWD@SERVER:PORT/DBNAME"


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-number", default=RUN_NUMBER)
    parser.add_argument("--hdf-path", default=str(DEFAULT_HDF_PATH))
    parser.add_argument("--household-csv", default=str(DEFAULT_HOUSEHOLD_CSV))
    parser.add_argument("--person-csv", default=str(DEFAULT_PERSON_CSV))
    parser.add_argument("--voter-registration-csv", default=str(DEFAULT_VOTER_REGISTRATION_CSV))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--load-from-hdf", action="store_true")
    parser.add_argument("--sql", default=DEFAULT_SQL)
    parser.add_argument("--connection-string", default=DEFAULT_CONNECTION_STRING)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    hdf_path = Path(args.hdf_path)
    if args.load_from_hdf:
        loader = lambda: load_buildings_from_hdf(hdf_path)
    else:
        loader = lambda: load_buildings_from_sql(args.sql, args.connection_string, hdf_path)

    run_household_placement(
        run_number=args.run_number,
        hdf_path=hdf_path,
        household_csv=args.household_csv,
        person_csv=args.person_csv,
        voter_registration_csv=args.voter_registration_csv,
        output_dir=args.output_dir,
        buildings_loader=loader,
    )


if __name__ == "__main__":
    main()
