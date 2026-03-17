import argparse
import re
import shutil
import time
from pathlib import Path

import pandas as pd
from census import Census

from semcog_popsim.input_prep.utils import (
    CensusDownloader,
    combine_puma_data,
    group_pums_data,
    marginal_summary,
    preprocess_pums,
    pums_update,
    read_tract_puma_crosswalk,
)

try:
    import oyaml as yaml
except ImportError:
    import yaml


SORT_ORDER = ["REGION", "PUMA", "MCD", "SAMPLEGEO", "TRACT", "BLKGRP", "TAZ", "BLK", "BUILDING"]
CROSSWALK_MAP = {
    ("TRACT10", "PUMA00"): "2010_Census_Tract_to_2000_PUMA_SEMCOG.csv",
    ("TRACT10", "PUMA10"): "2010_Census_Tract_to_2010_PUMA_SEMCOG.csv",
    ("TRACT20", "PUMA10"): "2020_Census_Tract_to_2010_PUMA_SEMCOG.csv",
    ("TRACT20", "PUMA20"): "2020_Census_Tract_to_2020_PUMA_SEMCOG.csv",
}
ACS_PUMA_MAP = {
    (2010, 2011): [("TRACT10", "PUMA00")],
    (2012, 2013, 2014, 2015, 2016): [("TRACT10", "PUMA00"), ("TRACT10", "PUMA10")],
    (2017, 2018, 2019): [("TRACT10", "PUMA10")],
    (2020, 2021): [("TRACT20", "PUMA10")],
    (2022, 2023, 2024, 2025, 2026): [("TRACT20", "PUMA10"), ("TRACT20", "PUMA20")],
}


def _load_yaml(path):
    with path.open("r") as stream:
        return yaml.load(stream, Loader=yaml.FullLoader)


def _resolve_base_dir(config_path):
    if config_path.parent.name.isdigit():
        return config_path.parent.parent
    return config_path.parent


def _resolve_path(path_value, base_dir, repo_root):
    candidate = Path(path_value)
    if candidate.is_absolute():
        return candidate
    for root in [base_dir, repo_root, Path.cwd()]:
        resolved = (root / candidate).resolve()
        if resolved.exists():
            return resolved
    return (base_dir / candidate).resolve()


def _resolve_geo_dir(base_dir, repo_root):
    candidates = [
        base_dir / "geo",
        repo_root / "projects" / "geo",
        repo_root / "input_prep" / "geo",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def _find_puma_definition(acs_year):
    for years, definitions in ACS_PUMA_MAP.items():
        if acs_year in years:
            return definitions
    raise ValueError(f"No PUMA crosswalk mapping defined for ACS year {acs_year}")


def run_input_prep(api_key, config_path):
    config_path = Path(config_path).resolve()
    base_dir = _resolve_base_dir(config_path)
    repo_root = base_dir.parent
    geo_dir = _resolve_geo_dir(base_dir, repo_root)
    conf = _load_yaml(config_path)
    start = time.time()

    project = conf["project"]
    project_name = project["name"]
    target = project["target"]
    acs_year = project["acs_year"]
    acs_sample = project["acs_sample"]
    project_folder = base_dir / str(acs_year)
    settings_file = _resolve_path(project["settings"].format(str(acs_year)), project_folder, repo_root)
    pre_control = _resolve_path(project["pre_control"].format(str(acs_year)), project_folder, repo_root)
    household_pums_csv = _resolve_path(project["h_pums_csv"], base_dir, repo_root)
    person_pums_csv = _resolve_path(project["p_pums_csv"], base_dir, repo_root)

    geography = conf["geography"]
    state = geography["state"][0]
    counties = geography["counties"]

    output_folder = project_folder / "data"
    output_folder.mkdir(parents=True, exist_ok=True)
    output_geo_cross = f"{project_name}_{acs_year}_geo_cross_walk.csv"
    output_control = f"{project_name}_{acs_year}_control_totals_.csv"
    output_seed_households = f"{project_name}_{acs_year}_seed_households.csv"
    output_seed_persons = f"{project_name}_{acs_year}_seed_persons.csv"

    print(f"\n *** preparing {target} synthesisdata for year {acs_year} ***")
    if acs_year < 2010 or acs_year >= 2026:
        raise ValueError("synthesis year should be between 2010 and 2026")

    census_client = Census(api_key, year=acs_year)
    print(f"\nPreparing Census geographies crosswalk: \n\tstate: {state}  \n\tcounty: {counties}")
    geography_reader = CensusDownloader(census_client.acs5, state, counties, "*", "*")
    df_geo = pd.DataFrame.from_dict(geography_reader.download("NAME")).drop("NAME", axis=1)
    df_geo["tractid"] = df_geo["state"] + df_geo["county"] + df_geo["tract"]
    df_geo["blkgrpid"] = df_geo["tractid"] + df_geo["block group"]
    df_geo.columns = [column.upper() for column in df_geo.columns]

    puma_definitions = _find_puma_definition(acs_year)
    if len(puma_definitions) == 1:
        df_tract_puma = read_tract_puma_crosswalk(puma_definitions[0], CROSSWALK_MAP, geo_dir)
    else:
        df_tract_puma = read_tract_puma_crosswalk(puma_definitions[0], CROSSWALK_MAP, geo_dir)
        df_tract_puma_extra = read_tract_puma_crosswalk(puma_definitions[1], CROSSWALK_MAP, geo_dir)
        df_tract_puma["PUMA"] = df_tract_puma["PUMA"] + df_tract_puma_extra["PUMA"]
    df_tract_puma = df_tract_puma.reset_index()

    df_geo_cross = pd.merge(df_geo, df_tract_puma, on="TRACTID", how="left")
    df_geo_cross["REGION"] = 2
    df_geo_cross = df_geo_cross[["TRACTID", "BLKGRPID", "PUMA", "COUNTYID", "REGION"]]
    print(f"  saving geo cross walk to: {output_folder / output_geo_cross}")
    df_geo_cross.to_csv(output_folder / output_geo_cross)

    census_sample = getattr(census_client, acs_sample)
    print(f"\nCreating popsim marginal controls: {acs_year}")
    print("  downloading Census variables ...")
    df_controls = pd.read_csv(pre_control)
    marginals = {}
    for geo_name, geo_frame in df_controls.groupby("geography"):
        full_vars = list(set(re.findall(r"[B-C][0-9]{5}[A-Z]{0,1}_[0-9]{3}E", str(list(geo_frame.acs_variables)))))
        if geo_name == "BLKGRP":
            downloader = CensusDownloader(census_sample, state, counties, "*", "*")
            geo_cols = ["state", "county", "tract", "block group"]
        elif geo_name == "TRACT":
            downloader = CensusDownloader(census_sample, state, counties, "*")
            geo_cols = ["state", "county", "tract"]
        elif geo_name == "COUNTY":
            downloader = CensusDownloader(census_sample, state, counties)
            geo_cols = ["state", "county"]
        else:
            raise ValueError(f"Unsupported geography {geo_name}")
        print(f"\t{geo_name} marginals ")
        marginals[geo_name] = downloader.download(full_vars).set_index(geo_cols)
        if "GEO_ID" in marginals[geo_name].columns:
            marginals[geo_name].drop("GEO_ID", axis=1, inplace=True)

    print("  compiling popsim control fields ...")
    for geo_name, geo_frame in df_controls.groupby("geography"):
        marginals[geo_name] = marginals[geo_name].astype(float).fillna(0)
        for _, row in geo_frame.iterrows():
            marginals[geo_name][row["control_field"]] = marginals[geo_name].eval(row["acs_variables"].replace('"', ""))
        marginals[geo_name] = marginals[geo_name][list(geo_frame.control_field)]

    control_outputs = {}
    mapping_dict = {4: "BLKGRPID", 3: "TRACTID", 2: "COUNTYID"}
    for geo_name, df_margin in marginals.items():
        df_margin[mapping_dict[df_margin.index.nlevels]] = df_margin.index.map("".join)
        df_margin.reset_index(drop=True, inplace=True)
        df_margin.fillna(0, inplace=True)
        df_margin.columns = [column.upper() for column in df_margin.columns]
        output_control_name = output_control.replace(".csv", geo_name.lower() + ".csv")
        control_outputs[geo_name] = output_control_name
        print(f"  saving controls to: {output_folder / output_control_name}")
        df_margin.to_csv(output_folder / output_control_name)
        marginal_summary(df_margin)

    print("\nExtrating PUMS seed households and persons from state samples")
    puma_list = df_geo_cross.PUMA.unique()
    household_pums = pd.read_csv(household_pums_csv, dtype={"SERIALNO": str, "PUMA": str}).set_index("SERIALNO")
    person_pums = pd.read_csv(person_pums_csv, dtype={"SERIALNO": str, "PUMA": str})
    if acs_year in conf["pums_var_updates"]:
        household_pums = pums_update(household_pums, conf["pums_var_updates"][acs_year])
        person_pums = pums_update(person_pums, conf["pums_var_updates"][acs_year])

    if len(puma_definitions) == 1:
        household_pums = household_pums.loc[household_pums.PUMA.isin(puma_list)]
        person_pums = person_pums.loc[person_pums.PUMA.isin(puma_list)]
    else:
        grouped = group_pums_data(household_pums, person_pums, puma_definitions[0][1], puma_definitions[1][1])
        household_samples, person_samples = combine_puma_data(puma_list, puma_definitions[0][1], puma_definitions[1][1], grouped)
        household_pums = pd.concat(household_samples)
        person_pums = pd.concat(person_samples)

    if target != "housing_units":
        household_pums = household_pums.loc[(household_pums.TYPE == 1) & (household_pums.NP > 0)]
    person_pums = person_pums.loc[person_pums["SERIALNO"].isin(household_pums.index)]
    household_pums, person_pums = preprocess_pums(household_pums, person_pums)
    household_pums["hh_id"] = range(len(household_pums))
    person_pums = pd.merge(person_pums, household_pums[["hh_id"]], left_on="SERIALNO", right_index=True, how="left")

    print(f"- saving seed households: {output_folder / output_seed_households}.| total {len(household_pums)} records")
    household_pums.to_csv(output_folder / output_seed_households)
    print(f"- saving seed persons: {output_folder / output_seed_persons}.| total {len(person_pums)} records")
    person_pums.to_csv(output_folder / output_seed_persons)

    print("\nupdate popsim settings")
    sorted_geos = list(control_outputs.keys())
    sorted_geos.sort(key=lambda value: SORT_ORDER.index(value))
    project_settings = _load_yaml(settings_file)
    geographies = sorted_geos + ["REGION", "PUMA"]
    geographies.sort(key=lambda value: SORT_ORDER.index(value))
    project_settings["geographies"] = geographies
    project_settings["seed_geography"] = "PUMA"
    project_settings["data_dir"] = f"data/{acs_year}"

    remaining_controls = control_outputs.copy()
    for item in list(project_settings["input_table_list"]):
        if item["tablename"] == "households":
            item["filename"] = output_seed_households
        if item["tablename"] == "persons":
            item["filename"] = output_seed_persons
        if item["tablename"] == "geo_cross_walk":
            item["filename"] = output_geo_cross
        if "_control_data" in item["tablename"]:
            geo_name = item["tablename"].replace("_control_data", "")
            if geo_name not in remaining_controls:
                project_settings["input_table_list"].remove(item)
            else:
                item["filename"] = remaining_controls[geo_name]
                del remaining_controls[geo_name]
    for geo_name, filename in remaining_controls.items():
        project_settings["input_table_list"].append({"tablename": geo_name + "_control_data", "filename": filename})

    project_settings["control_file_name"] = f"{project_name}_{acs_year}_controls.csv"
    project_settings["output_tables"] = {"action": "include", "tables": ["summary_" + geo for geo in sorted_geos]}
    sub_balance_steps = ["sub_balancing.geography=" + geo for geo in sorted_geos]
    project_settings["run_list"]["steps"] = [
        "input_pre_processor",
        "setup_data_structures",
        "initial_seed_balancing",
        "meta_control_factoring",
        "final_seed_balancing",
        "integerize_final_seed_weights",
    ] + sub_balance_steps + ["expand_households", "summarize", "write_tables", "write_synthetic_population"]

    settings_output = output_folder / f"{project_name}_{acs_year}_settings.yaml"
    with settings_output.open("w") as yaml_file:
        yaml.dump(project_settings, yaml_file, default_style=None, default_flow_style=False, sort_keys=False)

    print("copy popsim master control file")
    shutil.copy(pre_control, output_folder / f"{project_name}_{acs_year}_controls.csv")
    print(
        f"\ntotal time: {round(time.time() - start, 1)} seconds",
        f"\nDone. All files are saved to {output_folder}",
        '\nTo run Populationsim:',
        '\n\t copy new settings and controls to configs folder and rename "xxx_settings.yaml" to "settings.yaml"',
        f"\n\t copy other files in {acs_year}/data folder to data/{acs_year}/",
    )



def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("key", help="Census API key")
    parser.add_argument("yaml", help="yaml configuration file name")
    return parser



def main(argv=None):
    args = build_parser().parse_args(argv)
    run_input_prep(args.key, args.yaml)


if __name__ == "__main__":
    main()
