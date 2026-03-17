import argparse
import re
import time
from pathlib import Path

import pandas as pd
from census import Census

from semcog_popsim.input_prep.maker import _load_yaml, _resolve_base_dir, _resolve_path
from semcog_popsim.input_prep.utils import CensusDownloader, marginal_summary


def run_control_adjustment(api_key, config_path):
    config_path = Path(config_path).resolve()
    base_dir = _resolve_base_dir(config_path)
    repo_root = base_dir.parent
    conf = _load_yaml(config_path)
    start = time.time()

    project = conf["project"]
    project_name = project["name"]
    target = project["target"]
    acs_year = project["acs_year"]
    acs_sample = project["acs_sample"]
    project_folder = base_dir / str(acs_year)
    pre_control = _resolve_path(project["pre_control"].format(str(acs_year)), project_folder, repo_root)

    geography = conf["geography"]
    state = geography["state"][0]
    counties = geography["counties"]

    output_folder = project_folder / "data"
    output_folder.mkdir(parents=True, exist_ok=True)
    output_control = f"{project_name}_{acs_year}_control_totals_.csv"
    print(f"\n *** download {target} for year {acs_year} ***")

    census_client = Census(api_key, year=acs_year)
    census_sample = getattr(census_client, acs_sample)
    print(f"\nMaking popsim adjustment controls: {acs_year}")
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

    mapping_dict = {4: "BLKGRPID", 3: "TRACTID", 2: "COUNTYID"}
    for geo_name, df_margin in marginals.items():
        df_margin[mapping_dict[df_margin.index.nlevels]] = df_margin.index.map("".join)
        df_margin.reset_index(drop=True, inplace=True)
        df_margin.fillna(0, inplace=True)
        df_margin.columns = [column.upper() for column in df_margin.columns]
        if "HHBASE" in df_margin.columns:
            df_margin = df_margin.loc[df_margin.HHBASE > 0]
        marginal_summary(df_margin)
        output_name = output_control.replace(".csv", geo_name.lower() + "_adj.csv")
        print(f"  saving control file: {output_folder / output_name}")
        df_margin.to_csv(output_folder / output_name)

    print(f"\ntotal time: {round(time.time() - start, 1)} seconds")


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("key", help="Census API key")
    parser.add_argument("yaml", help="yaml configuration file name")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    run_control_adjustment(args.key, args.yaml)


if __name__ == "__main__":
    main()
