from forecast_input.input import load_from_files, load_from_sql
from forecast_input.transform import calculate_improvement_values, transform_buildings, transform_hh, transform_persons, households_add_n_18plus
from forecast_input.input import load_from_files, load_from_sql, load_pop_synthetic_csv
from forecast_input.placement import run_placement
import pandas as pd
import yaml
import os

#########################################
###          Set up Phase             ###
#########################################
# 2015 HDF
hdf_path = '/mnt/hgfs/urbansim/RDF2045/data/base_year/all_semcog_data_02-02-18-final-forecast.h5'
# read HDF cache
hdf = pd.HDFStore(hdf_path, 'r')

# sync hh and persons
sim_hh = pd.read_csv('/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2015(2017)/synthetic_households.csv')
sim_p = pd.read_csv('/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2015(2017)/synthetic_persons.csv')

voter_registration = pd.read_csv('/mnt/hgfs/da/Staff/Nutting/RDF2050/Qualified Voter File/placement_2020.csv')

#########################################
###           Input Phase             ###
#########################################
sim_hh = load_pop_synthetic_csv(sim_hh)
sim_hh = households_add_n_18plus(sim_hh, sim_p)

sim_p = transform_persons(sim_p)

# syn_hdf = pd.HDFStore(os.path.join(directory, "starter6_20170609-1702.h5"), mode="r")

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

placement_households = run_placement(
	sim_hh,
	b,  # buildings
	voter_registration,
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

households = transform_hh(households.set_index("household_id"), hdf["persons"])

persons = sim_p
persons = transform_persons(persons)
households = transform_hh(households.set_index("household_id"), persons)

households.to_csv("test_run/092324_run_2015/households.csv")
persons.to_csv("test_run/092324_run_2015/persons.csv")