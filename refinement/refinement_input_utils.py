from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

SUPPORTED_GEO_ALIASES = {
    "COUNTY": "COUNTY",
    "MCD": "MCD",
    "CITY": "MCD",
    "TAZ": "TAZ",
    "TRACT": "TRACT",
    "TRACTID": "TRACT",
    "BLKGRP": "BLKGRP",
    "BLKGRPID": "BLKGRP",
}
STANDARD_GEO_UNITS = ["COUNTY", "MCD", "TAZ", "TRACT", "BLKGRP", "REGION"]
DEFAULT_PERSON_CONTROL_CATEGORIES = {
    "age": [[0, 17], [18, 24], [25, 64], [65, -1]],
    "sex": [[1, 1], [2, 2]],
    "race_id": [[1, 1], [2, 2], [3, 3], [4, 4]],
}
DEFAULT_NEW_ID_BASE = 2 * 10**10


def normalize_geo_name(name: str) -> str:
    normalized = str(name).strip().upper().replace(" ", "_")
    if normalized not in SUPPORTED_GEO_ALIASES:
        supported = ", ".join(sorted(SUPPORTED_GEO_ALIASES))
        raise ValueError(f"Unsupported geography '{name}'. Supported values: {supported}")
    return SUPPORTED_GEO_ALIASES[normalized]


def ensure_columns(df: pd.DataFrame, required: Iterable[str], label: str) -> None:
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def ensure_hdf_table(store: pd.HDFStore, key: str, required_columns: Iterable[str], label: str) -> pd.DataFrame:
    if key not in store:
        raise KeyError(f"Missing HDF table '{key}' for {label}")
    df = store[key]
    ensure_columns(df, required_columns, label)
    return df.copy()


def cast_int_columns(df: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    for column in columns:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="raise").astype("int64")
    return df


def building_geos(store: pd.HDFStore, year: str, weight_col: str | None = None) -> pd.DataFrame:
    building_cols = ["parcel_id", "b_city_id", "b_zone_id"]
    if weight_col:
        building_cols.append(weight_col)
    buildings = ensure_hdf_table(store, f"{year}/buildings", building_cols, f"{year}/buildings")
    parcels = ensure_hdf_table(store, f"{year}/parcels", ["census_bg_id", "county_id"], f"{year}/parcels")

    dfgeo = buildings.merge(parcels[["census_bg_id", "county_id"]], left_on="parcel_id", right_index=True, how="left")
    if dfgeo[["b_city_id", "b_zone_id", "census_bg_id", "county_id"]].isnull().any().any():
        raise ValueError("Geography join produced null values. Check buildings/parcels linkage before refinement generation.")

    dfgeo = dfgeo.rename(
        columns={
            "b_city_id": "MCD",
            "county_id": "COUNTY",
            "census_bg_id": "BLKGRP",
            "b_zone_id": "TAZ",
        }
    )
    dfgeo["BLKGRP"] = dfgeo["COUNTY"] * (10**7) + dfgeo["BLKGRP"]
    dfgeo["TRACT"] = dfgeo["BLKGRP"] // 10
    dfgeo["REGION"] = 2
    return cast_int_columns(dfgeo, ["COUNTY", "MCD", "TAZ", "TRACT", "BLKGRP", "REGION"])


def create_geo_crosswalk(
    dfgeo: pd.DataFrame,
    target_geo: str,
    sample_geo: str,
    weight_col: str,
    sampling_method: str,
    new_id_base: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ensure_columns(dfgeo, [target_geo, sample_geo, weight_col], "refinement geography frame")
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
    geocross = cast_int_columns(geocross, [target_geo, sample_geo, "SAMPLEGEO", "REGION"])
    return geocross, grouped, geocross_dup


def cats_to_ctrl(
    dict_cats: dict[str, Iterable[Iterable[int]]],
    target_geo: str,
    seed_tbl: str,
    dfctrl: pd.DataFrame | None = None,
    importance: int = 500,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = [] if dfctrl is None else dfctrl.to_dict("records")
    for category, bounds in dict_cats.items():
        target_name = f"hh{category}" if seed_tbl == "households" else category
        for index, (vmin, vmax) in enumerate(bounds, start=1):
            if vmin == vmax:
                expression = f"({seed_tbl}.{category} == {vmin})"
            elif vmax == -1:
                expression = f"({seed_tbl}.{category} >= {vmin})"
            else:
                expression = f"({seed_tbl}.{category} >= {vmin}) & ({seed_tbl}.{category} <= {vmax})"
            rows.append(
                {
                    "target": f"{target_name}{vmin}",
                    "geography": target_geo,
                    "seed_table": seed_tbl,
                    "importance": importance,
                    "control_field": f"{target_name.upper()}{index}",
                    "expression": expression,
                }
            )
    return pd.DataFrame(rows)


def extract_control_categories(ctotals: pd.DataFrame) -> dict[str, np.ndarray]:
    ctotals = ctotals.copy()
    if "race_id" in ctotals.columns:
        ctotals["race_id_min"] = ctotals["race_id"]
        ctotals["race_id_max"] = ctotals["race_id"]

    categories = [column[:-4] for column in ctotals.columns if column.endswith("_min")]
    if not categories:
        raise ValueError("No *_min control category columns found in annual_household_control_totals")
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
    bldg_geos = building_geos(store, year)

    households = ensure_hdf_table(store, f"{year}/households", ["building_id"], f"{year}/households")
    households["household_id"] = households.index.astype("int64")
    households = households.merge(bldg_geos[geo_units], left_on="building_id", right_index=True, how="left")
    if filter_county is not None:
        households = households.loc[households["COUNTY"] == filter_county].copy()
    households = households.dropna(subset=geo_units)
    if households.empty:
        raise ValueError("No household records remain after applying geography join and county filter.")
    households = cast_int_columns(households, ["household_id", *geo_units])

    persons = ensure_hdf_table(store, f"{year}/persons", ["household_id"], f"{year}/persons")
    persons["person_id"] = persons.index.astype("int64")
    persons = persons.merge(households[["household_id", *geo_units]], on="household_id", how="left")
    persons = persons.dropna(subset=geo_units)
    if persons.empty:
        raise ValueError("No person records remain after applying geography join and county filter.")
    persons = cast_int_columns(persons, ["person_id", "household_id", *geo_units])
    return households, persons


def evaluate_control_mask(seed_df: pd.DataFrame, seed_table: str, expression: str) -> pd.Series:
    clean_expression = expression.replace(f"{seed_table}.", "")
    try:
        mask = seed_df.eval(clean_expression)
    except Exception as exc:  # pragma: no cover - detailed message path
        raise ValueError(f"Failed to evaluate control expression '{expression}': {exc}") from exc
    if mask.dtype != bool:
        raise ValueError(f"Control expression '{expression}' did not evaluate to a boolean mask")
    return mask


def summarize_geo_controls(
    simctrl: pd.DataFrame,
    seed_frames: dict[str, pd.DataFrame],
    target_geo: str,
) -> pd.DataFrame:
    control_series = []
    for _, row in simctrl.iterrows():
        seed_df = seed_frames[row.seed_table]
        ensure_columns(seed_df, [row.geography], f"{row.seed_table} seed frame")
        mask = evaluate_control_mask(seed_df, row.seed_table, row.expression)
        series = seed_df.loc[mask].groupby(row.geography).size()
        series.name = row.control_field
        control_series.append(series)

    if not control_series:
        raise ValueError("No control totals were generated from the control specification.")

    geo_ctrl = pd.concat(control_series, axis=1)
    geo_ctrl.index.name = target_geo
    geo_ctrl["REGION"] = 2
    geo_ctrl = geo_ctrl.fillna(0)
    geo_ctrl = cast_int_columns(geo_ctrl.reset_index(), [target_geo, "REGION", *[c for c in geo_ctrl.columns if c != "REGION"]])
    return geo_ctrl


def sample_largest(
    df_total: pd.DataFrame,
    geocross: pd.DataFrame,
    sample_geo: str,
    allow_empty: bool = False,
) -> pd.DataFrame:
    sampled = df_total.loc[df_total[sample_geo].isin(geocross[sample_geo].unique())].copy()
    if sampled.empty and not allow_empty:
        raise ValueError(f"No seed records matched sample geography '{sample_geo}'.")
    sampled["SAMPLEGEO"] = sampled[sample_geo]
    return cast_int_columns(sampled, [sample_geo, "SAMPLEGEO"])


def finalize_seed_households(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["household_id"] = pd.to_numeric(df["household_id"], errors="raise").astype("int64")
    df["hh_id"] = df["household_id"]
    df["WGTP"] = 1
    return df


def finalize_seed_persons(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["person_id"] = pd.to_numeric(df["person_id"], errors="raise").astype("int64")
    df["household_id"] = pd.to_numeric(df["household_id"], errors="raise").astype("int64")
    df["hh_id"] = df["household_id"]
    df["PWGTP"] = 1
    return df


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
    base_geocross = geocross.loc[geocross["SAMPLEGEO"] <= new_id_base].copy()

    allow_empty_base = sampling_method != "LARGEST"
    hhs = finalize_seed_households(sample_largest(hhs_total, base_geocross, sample_geo, allow_empty=allow_empty_base))
    pps = finalize_seed_persons(sample_largest(pps_total, base_geocross, sample_geo, allow_empty=allow_empty_base))

    if sampling_method == "LARGEST":
        return hhs, pps

    hh_frames: list[pd.DataFrame] = []
    person_frames: list[pd.DataFrame] = []
    next_household_id = new_id_base if hhs.empty else max(new_id_base, int(hhs["household_id"].max()) + 1)
    next_person_id = new_id_base if pps.empty else max(new_id_base, int(pps["person_id"].max()) + 1)

    for geo, dfg in geocross_dup.groupby(target_geo):
        if sampling_method == "OVERLAP":
            hh_local = pd.concat(
                [
                    hhs_total.loc[hhs_total[sample_geo] == dfg[sample_geo].values[0]],
                    hhs_total.loc[
                        (hhs_total[sample_geo].isin(dfg[sample_geo].values[1:]))
                        & (hhs_total[target_geo] == geo)
                    ],
                ]
            ).copy()
        elif sampling_method == "ALL":
            hh_local = hhs_total.loc[hhs_total[sample_geo].isin(dfg[sample_geo].values)].copy()
        else:
            raise ValueError(f"Unknown sampling method: {sampling_method}")

        if hh_local.empty:
            continue

        hh_local["SAMPLEGEO"] = geo + new_id_base
        hh_local["household_id"] = np.arange(next_household_id, next_household_id + len(hh_local), dtype="int64")
        next_household_id += len(hh_local)
        hh_frames.append(finalize_seed_households(hh_local))

        person_local = pps_total.loc[pps_total["household_id"].isin(hh_local.index.values)].copy()
        if person_local.empty:
            continue
        household_lookup = hh_local[["household_id", "SAMPLEGEO"]].copy()
        household_lookup.index.name = "source_household_id"
        person_local = person_local.merge(
            household_lookup,
            left_on="household_id",
            right_index=True,
            how="left",
        )
        person_local["household_id"] = person_local["household_id_y"]
        person_local = person_local.drop(columns=["household_id_x", "household_id_y"])
        person_local["person_id"] = np.arange(next_person_id, next_person_id + len(person_local), dtype="int64")
        next_person_id += len(person_local)
        person_frames.append(finalize_seed_persons(person_local))

    if hh_frames:
        hhs = pd.concat([hhs, *hh_frames], ignore_index=True)
    if person_frames:
        pps = pd.concat([pps, *person_frames], ignore_index=True)

    return hhs, pps


def runtime_control_columns(df: pd.DataFrame) -> pd.DataFrame:
    runtime_cols = [
        column
        for column in ["target", "geography", "seed_table", "importance", "control_field", "expression"]
        if column in df.columns
    ]
    return df.loc[:, runtime_cols].copy()


def build_output_synthetic_population(households: pd.DataFrame, persons: pd.DataFrame) -> dict[str, Any]:
    return {
        "household_id": "household_id",
        "households": {
            "filename": "synthetic_households.csv",
            "columns": list(households.columns),
        },
        "persons": {
            "filename": "synthetic_persons.csv",
            "columns": list(persons.columns),
        },
    }


def build_runtime_settings(
    template: dict[str, Any],
    target_geo: str,
    geo_crosswalk_file: Path,
    controls_file: Path,
    control_totals_file: Path,
    seed_households_file: Path,
    seed_persons_file: Path,
    households: pd.DataFrame,
    persons: pd.DataFrame,
) -> dict[str, Any]:
    settings = dict(template)
    settings["geographies"] = ["REGION", "SAMPLEGEO", target_geo]
    settings["seed_geography"] = "SAMPLEGEO"
    settings["data_dir"] = "data"
    settings["input_table_list"] = [
        {"tablename": "households", "filename": seed_households_file.name, "index_col": "hh_id"},
        {"tablename": "persons", "filename": seed_persons_file.name},
        {"tablename": "geo_cross_walk", "filename": geo_crosswalk_file.name},
        {"tablename": f"{target_geo}_control_data", "filename": control_totals_file.name},
    ]
    settings["household_weight_col"] = "WGTP"
    settings["household_id_col"] = "hh_id"
    settings["total_hh_control"] = "num_hh"
    settings["control_file_name"] = controls_file.name
    settings["output_tables"] = {"action": "include", "tables": [f"summary_{target_geo}"]}
    settings["output_synthetic_population"] = build_output_synthetic_population(households, persons)
    settings["models"] = [
        "input_pre_processor",
        "setup_data_structures",
        "initial_seed_balancing",
        "meta_control_factoring",
        "final_seed_balancing",
        "integerize_final_seed_weights",
        f"sub_balancing.geography={target_geo}",
        "expand_households",
        "summarize",
        "write_tables",
        "write_synthetic_population",
    ]
    settings.pop("run_list", None)
    return settings
