# %% [markdown]
# UrbanSim refinement input preparation
#
# This script converts full UrbanSim model outputs into a PopulationSim-ready
# refinement package: geo crosswalk, controls, control totals, and seed tables.
#
# It is a script version of `urbansim_refine_input.ipynb` with added cleanup:
# - output paths are centralized
# - CSVs are written without pandas index columns
# - chained assignment is avoided where practical
# - overlapping-sample list construction no longer aliases household/person lists
# - command-line arguments make it easier to rerun for different years/scenarios

# %%
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


# %% [markdown]
# ## Defaults

# %%
DEFAULT_SAMPLING_METHOD = "ALL"
DEFAULT_TARGET_GEO = "TAZ"
DEFAULT_SAMPLE_GEO = "CITY"
DEFAULT_YEAR = "2045"
DEFAULT_GEO_UNITS = ["COUNTY", "CITY", "TAZ", "TRACTID", "BLKGRPID", "REGION"]
DEFAULT_WEIGHT_COL = "residential_units"
DEFAULT_PROJECT_NAME = "refine"
DEFAULT_NEW_ID_BASE = 2 * 10**10
DEFAULT_COUNTY = 125
DEFAULT_HDF_INPUT = "run4032_45.h5"
DEFAULT_HDF_TARGET = "run4032_taz_draft_ypsi.h5"
DEFAULT_PERSON_CONTROL_CATEGORIES = {
    "age": [[0, 17], [18, 24], [25, 64], [65, -1]],
    "sex": [[1, 1], [2, 2]],
    "race_id": [[1, 1], [2, 2], [3, 3], [4, 4]],
}

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT_DIR = SCRIPT_DIR


# %% [markdown]
# ## Helper functions

# %%
def output_path(output_dir: Path, filename: str) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir / filename


def building_geos(store: pd.HDFStore, year: str, weight_col: str | None = None) -> pd.DataFrame:
    """
    Add a standard set of geographies to buildings based on UrbanSim HDF data.
    """
    rename_map = {
        "b_city_id": "CITY",
        "county_id": "COUNTY",
        "census_bg_id": "BLKGRPID",
        "b_zone_id": "TAZ",
    }
    cols = ["parcel_id", "b_city_id", "b_zone_id"]
    if weight_col:
        cols.append(weight_col)

    buildings = store[f"{year}/buildings"][cols]
    parcels = store[f"{year}/parcels"][["census_bg_id", "county_id"]]
    dfgeo = buildings.merge(parcels, left_on="parcel_id", right_index=True, how="left")
    dfgeo = dfgeo.rename(columns=rename_map)
    dfgeo["BLKGRPID"] = dfgeo["COUNTY"] * (10**7) + dfgeo["BLKGRPID"]
    dfgeo["TRACTID"] = dfgeo["BLKGRPID"] // 10
    dfgeo["REGION"] = 2
    return dfgeo


def create_geo_crosswalk(
    dfgeo: pd.DataFrame,
    target_geo: str,
    sample_geo: str,
    weight_col: str,
    sampling_method: str,
    new_id_base: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Build the primary geo crosswalk and keep the full overlap table for
    multi-area sampling strategies.
    """
    grouped = (
        dfgeo.groupby([target_geo, sample_geo], as_index=False)[weight_col]
        .sum()
        .sort_values(by=[target_geo, weight_col], ascending=False)
    )
    geocross = grouped.drop_duplicates(target_geo).copy()
    geocross["SAMPLEGEO"] = geocross[sample_geo]

    geocross_dup = grouped.loc[grouped.duplicated(target_geo, keep=False)].copy()
    if sampling_method != "LARGEST" and not geocross_dup.empty:
        overlap_targets = geocross_dup[target_geo].unique()
        geocross.loc[geocross[target_geo].isin(overlap_targets), "SAMPLEGEO"] = (
            geocross.loc[geocross[target_geo].isin(overlap_targets), target_geo] + new_id_base
        )

    geocross["REGION"] = 2
    geocross = geocross.astype("int64")
    return geocross, grouped, geocross_dup


def cats_to_ctrl(
    dict_cats: dict[str, Iterable[Iterable[int]]],
    target_geo: str,
    seed_tbl: str,
    dfctrl: pd.DataFrame | None = None,
    importance: int = 500,
) -> pd.DataFrame:
    """
    Convert UrbanSim control categories to a PopulationSim control spec.
    """
    if dfctrl is None:
        dfctrl = pd.DataFrame(
            columns=["target", "geography", "seed_table", "importance", "control_field", "expression"]
        )
        indv = 0
    else:
        indv = len(dfctrl)

    for category, bounds in dict_cats.items():
        target_name = f"hh{category}" if seed_tbl == "households" else category
        ccount = 0
        for vmin, vmax in bounds:
            indv += 1
            ccount += 1
            if vmin == vmax:
                expression = f"({seed_tbl}.{category} == {vmin})"
            elif vmax == -1:
                expression = f"({seed_tbl}.{category} >= {vmin})"
            else:
                expression = f"({seed_tbl}.{category} >= {vmin}) & ({seed_tbl}.{category} <= {vmax})"
            dfctrl.loc[indv] = [
                f"{target_name}{vmin}",
                target_geo,
                seed_tbl,
                importance,
                f"{target_name.upper()}{ccount}",
                expression,
            ]

    return dfctrl


def extract_control_categories(ctotals: pd.DataFrame) -> dict[str, np.ndarray]:
    """
    Extract distinct min/max category boundaries from UrbanSim household controls.
    """
    ctotals = ctotals.copy()
    if "race_id" in ctotals.columns:
        ctotals["race_id_min"] = ctotals["race_id"]
        ctotals["race_id_max"] = ctotals["race_id"]

    categories = [col[:-4] for col in ctotals.columns if col.endswith("_min")]
    return {
        category: ctotals[[f"{category}_min", f"{category}_max"]]
        .drop_duplicates()
        .sort_values(by=f"{category}_min")
        .values
        for category in categories
    }


def attach_households_and_persons(
    store: pd.HDFStore,
    year: str,
    geo_units: list[str],
    filter_county: int | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Attach geography columns to UrbanSim households and persons.
    """
    bldg_geos = building_geos(store, year)
    households = store[f"{year}/households"].merge(
        bldg_geos[geo_units], left_on="building_id", right_index=True, how="left"
    )
    if filter_county is not None:
        households = households.loc[households["COUNTY"] == filter_county].copy()
    households = households.dropna(axis=0).astype("int64")

    persons = store[f"{year}/persons"].merge(
        households[geo_units], left_on="household_id", right_index=True, how="left"
    )
    persons = persons.dropna(axis=0).astype("int64")
    return households, persons


def summarize_geo_controls(
    simctrl: pd.DataFrame,
    seed_frames: dict[str, pd.DataFrame],
    target_geo: str,
) -> pd.DataFrame:
    """
    Summarize control totals from households/persons using the generated control spec.
    """
    control_series = []
    for _, row in simctrl.iterrows():
        seed_df = seed_frames[row.seed_table]
        mask = eval(row.expression, {"np": np}, seed_frames)
        series = seed_df.loc[mask].groupby(row.geography).size()
        series.name = row.control_field
        control_series.append(series)

    geo_ctrl = pd.concat(control_series, axis=1)
    geo_ctrl.index.name = target_geo
    geo_ctrl["REGION"] = 2
    geo_ctrl = geo_ctrl.fillna(0).astype("int64")
    return geo_ctrl


def sample_largest(df_total: pd.DataFrame, geocross: pd.DataFrame, sample_geo: str) -> pd.DataFrame:
    sampled = df_total.loc[df_total[sample_geo].isin(geocross[sample_geo].unique())].copy()
    sampled["SAMPLEGEO"] = sampled[sample_geo]
    return sampled.astype("int64")


def build_seed_tables(
    hhs_total: pd.DataFrame,
    pps_total: pd.DataFrame,
    geocross: pd.DataFrame,
    geocross_dup: pd.DataFrame,
    sample_geo: str,
    target_geo: str,
    sampling_method: str,
    new_id_base: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Create seed household and person tables for the requested sampling method.
    """
    base_geocross = geocross.loc[geocross["SAMPLEGEO"] <= new_id_base].copy()

    hhs = sample_largest(hhs_total, base_geocross, sample_geo)
    hhs["WGTP"] = 1
    pps = sample_largest(pps_total, base_geocross, sample_geo)
    pps["PWGTP"] = 1

    if sampling_method == "LARGEST":
        hhs.index.name = "household_id"
        pps.index.name = "person_id"
        return hhs, pps

    hhlst = []
    pplst = []
    hhind = new_id_base
    ppind = new_id_base

    for geo, dfg in geocross_dup.groupby(target_geo):
        if sampling_method == "OVERLAP":
            hh_lc = pd.concat(
                [
                    hhs_total.loc[hhs_total[sample_geo] == dfg[sample_geo].values[0]],
                    hhs_total.loc[
                        (hhs_total[sample_geo].isin(dfg[sample_geo].values[1:]))
                        & (hhs_total[target_geo] == geo)
                    ],
                ]
            )
        elif sampling_method == "ALL":
            hh_lc = hhs_total.loc[hhs_total[sample_geo].isin(dfg[sample_geo].values)].copy()
        else:
            raise ValueError(f"Unknown sampling method: {sampling_method}")

        hh_lc = hh_lc.copy()
        hh_lc["SAMPLEGEO"] = geo + new_id_base
        hh_lc["new_hhid"] = range(hhind, hhind + len(hh_lc))
        hhlst.append(hh_lc)
        hhind += len(hh_lc)

        pp_lc = pps_total.loc[pps_total.household_id.isin(hh_lc.index.values)].copy()
        pp_lc = pp_lc.merge(
            hh_lc[["new_hhid", "SAMPLEGEO"]],
            left_on="household_id",
            right_index=True,
            how="left",
        )
        pp_lc.index = range(ppind, ppind + len(pp_lc))
        pplst.append(pp_lc)
        ppind += len(pp_lc)

    if hhlst:
        dfhhs = pd.concat(hhlst)
        dfhhs["WGTP"] = 1
        dfhhs.index = dfhhs["new_hhid"]
        dfhhs.index.name = "household_id"
        dfhhs = dfhhs.drop(columns=["new_hhid"]).astype("int64")
        hhs = pd.concat([hhs, dfhhs]).astype("int64")

    if pplst:
        dfpps = pd.concat(pplst)
        dfpps["PWGTP"] = 1
        dfpps["household_id"] = dfpps["new_hhid"]
        dfpps = dfpps.drop(columns=["new_hhid"])
        pps = pd.concat([pps, dfpps]).astype("int64")

    hhs.index.name = "household_id"
    pps.index.name = "person_id"
    return hhs, pps


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare PopulationSim refinement inputs from UrbanSim outputs.")
    parser.add_argument("--year", default=DEFAULT_YEAR, help="UrbanSim year key in the HDF stores.")
    parser.add_argument("--input-hdf", default=DEFAULT_HDF_INPUT, help="Base UrbanSim HDF used for sampling.")
    parser.add_argument("--target-hdf", default=DEFAULT_HDF_TARGET, help="Target UrbanSim HDF used for refinement controls.")
    parser.add_argument("--project-name", default=DEFAULT_PROJECT_NAME, help="Output file prefix.")
    parser.add_argument("--sampling-method", default=DEFAULT_SAMPLING_METHOD, choices=["LARGEST", "OVERLAP", "ALL"], help="Seed sampling strategy.")
    parser.add_argument("--target-geo", default=DEFAULT_TARGET_GEO, help="Target geography for refinement controls.")
    parser.add_argument("--sample-geo", default=DEFAULT_SAMPLE_GEO, help="Sampling geography used to source seed households/persons.")
    parser.add_argument("--weight-col", default=DEFAULT_WEIGHT_COL, help="Building weight column used to rank overlaps.")
    parser.add_argument("--county", type=int, default=DEFAULT_COUNTY, help="Optional county filter. Use a negative value to disable filtering.")
    parser.add_argument("--new-id-base", type=int, default=DEFAULT_NEW_ID_BASE, help="Large offset used for synthetic SAMPLEGEO and copied household ids.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR), help="Directory for generated refinement input CSVs.")
    parser.add_argument("--add-person-controls", action="store_true", help="Include the default person control categories in the control file.")
    return parser.parse_args()


def build_refinement_inputs(args: argparse.Namespace) -> dict[str, Path]:
    output_dir = Path(args.output_dir).resolve()
    filter_county = None if args.county is not None and args.county < 0 else args.county

    sampling_method = args.sampling_method
    if args.target_geo == args.sample_geo and sampling_method != "LARGEST":
        print("target_geo equals sample_geo; forcing sampling method to LARGEST")
        sampling_method = "LARGEST"

    with pd.HDFStore(args.input_hdf, "r") as store_input, pd.HDFStore(args.target_hdf, "r") as store_target:
        # %% [markdown]
        # ## Step 1. Build geo crosswalk
        # %%
        dfgeo = building_geos(store_input, args.year, args.weight_col)
        if filter_county is not None:
            dfgeo = dfgeo.loc[dfgeo["COUNTY"] == filter_county].copy()

        geocross, geocross_all, geocross_dup = create_geo_crosswalk(
            dfgeo,
            target_geo=args.target_geo,
            sample_geo=args.sample_geo,
            weight_col=args.weight_col,
            sampling_method=sampling_method,
            new_id_base=args.new_id_base,
        )

        geo_crosswalk_file = output_path(output_dir, f"{args.project_name}_geo_cross_walk.csv")
        geocross.to_csv(geo_crosswalk_file, index=False)

        print("sampling_method:", sampling_method)
        print("geo_cross_walk:", geo_crosswalk_file)
        print("geocross:", len(geocross), "   geocross_all:", len(geocross_all))

        # %% [markdown]
        # ## Step 2. Build PopulationSim control spec
        # %%
        ctotals = store_target["/base/annual_household_control_totals"]
        dict_cats = extract_control_categories(ctotals)

        simctrl = cats_to_ctrl(dict_cats, args.target_geo, "households")
        simctrl.loc[len(simctrl) + 1] = [
            "num_hh",
            args.target_geo,
            "households",
            10000000,
            "HHBASE",
            "(households.persons > 0)",
        ]

        if args.add_person_controls:
            simctrl = cats_to_ctrl(
                DEFAULT_PERSON_CONTROL_CATEGORIES,
                args.target_geo,
                "persons",
                dfctrl=simctrl,
            )

        controls_file = output_path(output_dir, f"{args.project_name}_controls.csv")
        simctrl.to_csv(controls_file, index=False)
        print("config control table:", controls_file)

        # %% [markdown]
        # ## Step 3. Generate target control totals
        # %%
        households_target, persons_target = attach_households_and_persons(
            store_target,
            args.year,
            DEFAULT_GEO_UNITS,
            filter_county=filter_county,
        )

        geo_ctrl = summarize_geo_controls(
            simctrl,
            seed_frames={"households": households_target, "persons": persons_target},
            target_geo=args.target_geo,
        )

        control_totals_file = output_path(
            output_dir, f"{args.project_name}_control_totals_{args.target_geo.lower()}.csv"
        )
        geo_ctrl.to_csv(control_totals_file, index=False)
        print("geography control table:", control_totals_file)

        # %% [markdown]
        # ## Step 4. Build seed households and persons
        # %%
        bldg_geos = building_geos(store_input, args.year)
        hhs_total = store_input[f"{args.year}/households"].merge(
            bldg_geos[DEFAULT_GEO_UNITS], left_on="building_id", right_index=True, how="left"
        )
        pps_total = store_input[f"{args.year}/persons"].merge(
            hhs_total[DEFAULT_GEO_UNITS], left_on="household_id", right_index=True, how="left"
        )

        if filter_county is not None:
            hhs_total = hhs_total.loc[hhs_total["COUNTY"] == filter_county].copy()
            pps_total = pps_total.loc[pps_total["COUNTY"] == filter_county].copy()

        hhs, pps = build_seed_tables(
            hhs_total=hhs_total,
            pps_total=pps_total,
            geocross=geocross,
            geocross_dup=geocross_dup,
            sample_geo=args.sample_geo,
            target_geo=args.target_geo,
            sampling_method=sampling_method,
            new_id_base=args.new_id_base,
        )

        seed_households_file = output_path(output_dir, f"{args.project_name}_seed_households.csv")
        seed_persons_file = output_path(output_dir, f"{args.project_name}_seed_persons.csv")
        hhs.to_csv(seed_households_file, index=False)
        pps.to_csv(seed_persons_file, index=False)

        print("seed households table:", seed_households_file)
        print("seed persons table:", seed_persons_file)

    return {
        "geo_crosswalk": geo_crosswalk_file,
        "controls": controls_file,
        "control_totals": control_totals_file,
        "seed_households": seed_households_file,
        "seed_persons": seed_persons_file,
    }


# %% [markdown]
# ## Run

# %%
def main() -> None:
    args = parse_args()
    outputs = build_refinement_inputs(args)
    print("created files:")
    for name, path in outputs.items():
        print(f"  {name}: {path}")


if __name__ == "__main__":
    main()
