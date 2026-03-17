from pathlib import Path

import pandas as pd

from semcog_popsim.forecast_input.input import load_pop_synthetic_csv
from semcog_popsim.forecast_input.placement import run_placement
from semcog_popsim.forecast_input.transform import (
    calculate_improvement_values,
    households_add_n_18plus,
    transform_buildings,
    transform_hh,
    transform_persons,
)


DROP_HOUSEHOLD_COLUMNS = [
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


def load_buildings_from_hdf(hdf_path):
    hdf = pd.HDFStore(hdf_path, "r")
    buildings = hdf["buildings"]
    buildings["owner_units"] = 0
    parcels = hdf["parcels"]
    parcels["tract"] = parcels["census_bg_id"].astype(str).str.slice(0, 6)
    parcels["bg"] = parcels["census_bg_id"].astype(str).str.slice(6, 7)
    parcels["mcd"] = parcels["semmcd"]
    parcels["county"] = parcels["county_id"]
    buildings = buildings.join(parcels[["tract", "bg", "county", "mcd"]], on="parcel_id")
    buildings = buildings.fillna(-1)
    buildings = transform_buildings(buildings)
    return hdf, buildings


def load_buildings_from_sql(sql, connection_string, hdf_path):
    hdf = pd.HDFStore(hdf_path, "r")
    buildings = pd.read_sql(sql, connection_string, index_col="building_id")
    buildings = calculate_improvement_values(buildings, hdf["parcels"])
    buildings = buildings.fillna(-1)
    buildings = transform_buildings(buildings)
    return hdf, buildings


def run_household_placement(
    run_number,
    hdf_path,
    household_csv,
    person_csv,
    voter_registration_csv,
    output_dir,
    buildings_loader,
    pretransform_persons_source=None,
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    hdf, buildings = buildings_loader()
    households = pd.read_csv(household_csv)
    persons = pd.read_csv(person_csv)
    voter_registration = pd.read_csv(voter_registration_csv)

    households = load_pop_synthetic_csv(households)
    households = households_add_n_18plus(households, persons)
    persons = transform_persons(persons)

    placed = run_placement(households, buildings, voter_registration, run_number)
    placed = placed.loc[
        placed.matched_household_id != -1,
        ["matched_household_id", "building_id"],
    ].rename(columns={"matched_household_id": "household_id"})

    households.insert(1, "building_id", placed["building_id"].values)
    households = households[
        [col for col in households.columns if col.lower() not in DROP_HOUSEHOLD_COLUMNS]
    ]
    households.loc[households["income"] < 0, "income"] = 0

    if pretransform_persons_source is not None:
        households = transform_hh(households, pretransform_persons_source(hdf))
    households = transform_hh(households, persons)

    households.to_csv(output_dir / "households.csv")
    persons.to_csv(output_dir / "persons.csv")
    hdf.close()
    return output_dir
