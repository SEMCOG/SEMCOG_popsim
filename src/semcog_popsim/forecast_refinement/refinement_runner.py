from pathlib import Path

import pandas as pd

from semcog_popsim.forecast_input.pop_refinement import refine_pop_single_year


def run_population_refinement(
    refinement_excel,
    hdf_path,
    household_csv,
    person_csv,
    output_dir,
    refine_geo="semmcd",
):
    refinement_excel = Path(refinement_excel)
    hdf_path = Path(hdf_path)
    household_csv = Path(household_csv)
    person_csv = Path(person_csv)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    refinement = pd.read_excel(refinement_excel, sheet_name=0)
    refinement.columns = refinement.columns.str.lower()
    refinement = refinement[refinement[refine_geo] < 8000]
    refinement = refinement.set_index(refine_geo)
    refinement.index = refinement.index.astype(int)
    refinement = refinement.fillna(0).astype(int)
    refinement = refinement["totalpop"] - refinement["gqpop"]
    refinement = refinement[refinement > 0]

    hdf = pd.HDFStore(hdf_path, "r")
    buildings = hdf["buildings"]
    parcels = hdf["parcels"]
    buildings = buildings.join(parcels[[refine_geo]], on="parcel_id")
    buildings = buildings[~buildings[refine_geo].isna()]
    buildings = buildings.astype({refine_geo: int})

    households = pd.read_csv(household_csv, index_col=0)
    persons = pd.read_csv(person_csv, index_col=0)

    households = households.join(buildings[[refine_geo]], on="building_id")
    persons = persons.join(households[[refine_geo]], on="household_id")
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
    households = households.rename(
        columns={
            "AGEHOH": "age_of_head",
            "VEN": "cars",
            "NP": "persons",
            "HRACE": "race_id",
            "HWORKERS": "workers",
            "HHT": "hht",
            "hh_id": "seed_id",
        }
    )

    new_households, new_persons = refine_pop_single_year(households, persons, buildings, refinement, refine_geo)
    new_households.to_csv(output_dir / "households_after_refinement.csv")
    new_persons.to_csv(output_dir / "persons_after_refinement.csv")
    review = pd.DataFrame(
        {
            "pop_refinement": refinement,
            "before_refinement": households.groupby(refine_geo).sum()["persons"],
            "after_refinement": new_households.groupby(refine_geo).sum()["persons"],
        }
    ).fillna(0).astype(int)
    review["diff"] = review["after_refinement"] - review["pop_refinement"]
    review.to_csv(output_dir / "pop_refinement_review.csv")
    return output_dir
