from forecast_input.input import load_from_files, load_from_sql
from forecast_input.transform import calculate_improvement_values, transform_buildings, transform_hh, transform_persons, households_add_n_18plus
from forecast_input.input import load_from_files, load_from_sql, load_pop_synthetic_csv
from forecast_input.placement import run_placement
import pandas as pd
import yaml
import os

run_number = '100124_run_2020'

# create output folder
if not os.path.exists(os.path.join('output', run_number)):
    os.makedirs(os.path.join('output', run_number))

#########################################
###          Set up Phase             ###
#########################################
# 2020 HDF
hdf_path = '/home/da/share/urbansim/RDF2050/model_inputs/base_hdf/forecast_data_input_031523.h5'
# read HDF cache
hdf = pd.HDFStore(hdf_path, 'r')

# sync hh and persons
sim_hh = pd.read_csv('/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2020(2022)/synthetic_households_ybl.csv')
sim_p = pd.read_csv('/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2020(2022)/synthetic_persons.csv')

voter_registration = pd.read_csv('/mnt/hgfs/da/Staff/Nutting/RDF2050/Qualified Voter File/placement_2020.csv')

#########################################
###           Input Phase             ###
#########################################
sim_hh = load_pop_synthetic_csv(sim_hh)
# 2020 hh synthesis doesn't have ybl
# if "ybl" not in sim_hh.columns:
#     sim_hh["ybl"] = 0
sim_hh = households_add_n_18plus(sim_hh, sim_p)

sim_p = transform_persons(sim_p)

# syn_hdf = pd.HDFStore(os.path.join(directory, "starter6_20170609-1702.h5"), mode="r")

LOAD_FROM_HDF = False
if LOAD_FROM_HDF:
    b = hdf['buildings']
    # fill in non-exist columns
    b['owner_units'] = 0
    parcels = hdf['parcels']
    parcels['tract'] = parcels['census_bg_id'].astype(str).str.slice(0, 6)
    parcels['bg'] = parcels['census_bg_id'].astype(str).str.slice(6, 7)
    parcels['mcd'] = parcels['semmcd']
    parcels['county'] = parcels['county_id']
    b = b.join(parcels[['tract', 'bg', 'county', 'mcd']], on='parcel_id')
    b = b.fillna(-1)
    b = transform_buildings(b)
else:
    sql = """
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
    b = pd.read_sql(sql, 'postgresql://gisad:forecast20@plannerprojection:5432/land', index_col='building_id')
    b = calculate_improvement_values(b, hdf["parcels"])
    b = b.fillna(-1)
    b = transform_buildings(b)

placement_households = run_placement(
	sim_hh,
	b,  # buildings
	voter_registration,
    run_number
)

placement_households = placement_households.loc[
    placement_households.matched_household_id != -1,
    ["matched_household_id", "building_id"],
]

placement_households = placement_households.rename(
    columns={"matched_household_id": "household_id"}
)

households = sim_hh

households.insert(1, "building_id", placement_households["building_id"].values)

households = households[
    [
        col
        for col in households.columns
        if col.lower()
        not in [
            "puma",
            "tract",
            "hincp",
            "r18",
            "hhisp",
            "adjinc",
            "ybl",
            "bld",
            "grntp",
            "adjhsg",
            "type",
            "valp",
        ]
    ]
]

households.loc[households["income"] < 0, "income"] = 0

persons = sim_p
households = transform_hh(households, persons)

households.to_csv(os.path.join('output', run_number, 'households.csv'))
persons.to_csv(os.path.join('output', run_number, 'persons.csv'))