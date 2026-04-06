from __future__ import annotations

import argparse
import shutil
from copy import deepcopy
from pathlib import Path
from typing import Any

import pandas as pd

try:
    import oyaml as yaml
except ModuleNotFoundError:
    import yaml

try:
    from .refinement_input_utils import (
        DEFAULT_NEW_ID_BASE,
        DEFAULT_PERSON_CONTROL_CATEGORIES,
        STANDARD_GEO_UNITS,
        attach_households_and_persons,
        build_runtime_settings,
        build_seed_tables,
        building_geos,
        cast_int_columns,
        cats_to_ctrl,
        create_geo_crosswalk,
        extract_control_categories,
        normalize_geo_name,
        runtime_control_columns,
        summarize_geo_controls,
    )
except ImportError:
    from refinement_input_utils import (  # type: ignore
        DEFAULT_NEW_ID_BASE,
        DEFAULT_PERSON_CONTROL_CATEGORIES,
        STANDARD_GEO_UNITS,
        attach_households_and_persons,
        build_runtime_settings,
        build_seed_tables,
        building_geos,
        cast_int_columns,
        cats_to_ctrl,
        create_geo_crosswalk,
        extract_control_categories,
        normalize_geo_name,
        runtime_control_columns,
        summarize_geo_controls,
    )

DEFAULT_SAMPLING_METHOD = "ALL"
DEFAULT_TARGET_GEO = "TAZ"
DEFAULT_SAMPLE_GEO = "MCD"
DEFAULT_YEAR = "2045"
DEFAULT_WEIGHT_COL = "residential_units"
DEFAULT_PROJECT_NAME = "refine"
DEFAULT_COUNTY = None

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_SETTINGS_TEMPLATE = REPO_ROOT / "configs" / "templates" / "refinement_settings_template.yaml"
DEFAULT_RUN_ROOT = REPO_ROOT.parent / "d_drive" / "popsim" / "runs"
ALLOWED_MODEL_GEOS = {"COUNTY", "MCD", "TAZ", "TRACT", "BLKGRP"}


def load_yaml(path: Path) -> dict[str, Any]:
    with open(path, "r") as stream:
        loaded = yaml.load(stream, Loader=yaml.FullLoader)
    return loaded or {}


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a full PopulationSim refinement input package from UrbanSim model outputs. "
            "This script creates configs/, data/, and output/ under d_drive/popsim/runs/<project>_<target_geo>/; "
            "it does not run refinement synthesis itself."
        )
    )
    parser.add_argument("config", help="YAML config file for refinement input generation")
    parser.add_argument("--project-name", help="Override the base project/package name from the YAML config")
    parser.add_argument("--year", help="Override the UrbanSim year key from the YAML config")
    parser.add_argument("--target-geo", help="Override the target geography from the YAML config")
    parser.add_argument("--sample-geo", help="Override the sample geography from the YAML config")
    parser.add_argument("--input-hdf", help="Override the base UrbanSim HDF path from the YAML config")
    parser.add_argument("--target-hdf", help="Override the target UrbanSim HDF path from the YAML config")
    parser.add_argument("--output-root", help="Override the run package root directory")
    parser.add_argument(
        "--sampling-method",
        choices=["LARGEST", "OVERLAP", "ALL"],
        help="Override the seed sampling method",
    )
    parser.add_argument("--weight-col", help="Override the building weight column")
    parser.add_argument(
        "--county",
        type=int,
        help="Optionally filter to a single county code. Omit this to keep all counties.",
    )
    parser.add_argument(
        "--new-id-base",
        type=int,
        help="Override the synthetic id offset for overlap/all sampling",
    )
    parser.add_argument(
        "--add-person-controls",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Override whether default person control categories are included",
    )
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def resolve_existing_path(path_str: str, config_dir: Path, label: str) -> Path:
    if not path_str:
        raise ValueError(f"Missing required config value for {label}")

    candidate = Path(path_str)
    if candidate.is_absolute():
        resolved = candidate
    else:
        for root in [config_dir, REPO_ROOT, SCRIPT_DIR]:
            probe = (root / candidate).resolve()
            if probe.exists():
                resolved = probe
                break
        else:
            resolved = (config_dir / candidate).resolve()

    if not resolved.exists():
        raise FileNotFoundError(f"Could not find {label}: {resolved}")
    if resolved.is_dir():
        raise ValueError(f"Expected a file path for {label}, but found a directory: {resolved}")
    return resolved


def resolve_output_root(path_str: str, config_dir: Path) -> Path:
    if not path_str:
        return DEFAULT_RUN_ROOT.resolve()
    candidate = Path(path_str)
    if candidate.is_absolute():
        return candidate.resolve()
    return (config_dir / candidate).resolve()


def validate_model_geo(name: str, label: str) -> str:
    geo = normalize_geo_name(name)
    if geo not in ALLOWED_MODEL_GEOS:
        allowed = ", ".join(sorted(ALLOWED_MODEL_GEOS))
        raise ValueError(f"{label} must be one of: {allowed}")
    return geo


def load_runtime_config(args: argparse.Namespace) -> dict[str, Any]:
    config_path = Path(args.config).resolve()
    if not config_path.exists():
        raise FileNotFoundError(f"Refinement config not found: {config_path}")

    conf = load_yaml(config_path)
    project = conf.get("project", {})
    paths = conf.get("paths", {})
    options = conf.get("options", {})

    project_name = args.project_name or project.get("name", DEFAULT_PROJECT_NAME)
    year = str(args.year or project.get("year", DEFAULT_YEAR))
    target_geo = validate_model_geo(
        args.target_geo or project.get("target_geography", DEFAULT_TARGET_GEO),
        "target geography",
    )
    sample_geo = validate_model_geo(
        args.sample_geo or project.get("sample_geography", DEFAULT_SAMPLE_GEO),
        "sample geography",
    )

    input_hdf = resolve_existing_path(args.input_hdf or paths.get("input_hdf", ""), config_path.parent, "input_hdf")
    target_hdf = resolve_existing_path(args.target_hdf or paths.get("target_hdf", ""), config_path.parent, "target_hdf")
    settings_template = resolve_existing_path(
        paths.get("settings_template", str(DEFAULT_SETTINGS_TEMPLATE)),
        config_path.parent,
        "settings_template",
    )
    output_root = resolve_output_root(
        args.output_root or paths.get("output_root", str(DEFAULT_RUN_ROOT)),
        config_path.parent,
    )

    sampling_method = args.sampling_method or options.get("sampling_method", DEFAULT_SAMPLING_METHOD)
    weight_col = args.weight_col or options.get("weight_col", DEFAULT_WEIGHT_COL)
    county = args.county if args.county is not None else options.get("county", DEFAULT_COUNTY)
    new_id_base = args.new_id_base if args.new_id_base is not None else options.get("new_id_base", DEFAULT_NEW_ID_BASE)
    add_person_controls = (
        args.add_person_controls
        if args.add_person_controls is not None
        else options.get("add_person_controls", False)
    )

    filter_county = None if county is not None and county < 0 else county
    if target_geo == sample_geo and sampling_method != "LARGEST":
        sampling_method = "LARGEST"

    run_name = f"{project_name}_{target_geo.lower()}"
    run_folder = output_root / run_name
    config_output_dir = run_folder / "configs"
    data_output_dir = run_folder / "data"
    runtime_output_dir = run_folder / "output"

    config_output_dir.mkdir(parents=True, exist_ok=True)
    data_output_dir.mkdir(parents=True, exist_ok=True)
    runtime_output_dir.mkdir(parents=True, exist_ok=True)

    return {
        "config_path": config_path,
        "project_name": project_name,
        "year": year,
        "target_geo": target_geo,
        "sample_geo": sample_geo,
        "input_hdf": input_hdf,
        "target_hdf": target_hdf,
        "settings_template": settings_template,
        "output_root": output_root,
        "run_name": run_name,
        "run_folder": run_folder,
        "configs_dir": config_output_dir,
        "data_dir": data_output_dir,
        "output_dir": runtime_output_dir,
        "sampling_method": sampling_method,
        "weight_col": weight_col,
        "filter_county": filter_county,
        "new_id_base": int(new_id_base),
        "add_person_controls": bool(add_person_controls),
    }


def build_refinement_inputs(config: dict[str, Any]) -> dict[str, Path]:
    project_name = config["project_name"]
    target_geo = config["target_geo"]

    geo_crosswalk_file = config["data_dir"] / f"{project_name}_geo_cross_walk.csv"
    control_totals_file = config["data_dir"] / f"{project_name}_control_totals_{target_geo}.csv"
    seed_households_file = config["data_dir"] / f"{project_name}_seed_households.csv"
    seed_persons_file = config["data_dir"] / f"{project_name}_seed_persons.csv"
    controls_file = config["configs_dir"] / "controls.csv"
    settings_file = config["configs_dir"] / "settings.yaml"
    config_copy = config["configs_dir"] / "refinement_input_config.yaml"

    with pd.HDFStore(config["input_hdf"], "r") as store_input, pd.HDFStore(config["target_hdf"], "r") as store_target:
        dfgeo = building_geos(store_input, config["year"], config["weight_col"])
        if config["filter_county"] is not None:
            dfgeo = dfgeo.loc[dfgeo["COUNTY"] == config["filter_county"]].copy()
        if dfgeo.empty:
            raise ValueError("No building geography rows remain after applying the county filter.")

        geocross, geocross_all, geocross_dup = create_geo_crosswalk(
            dfgeo,
            target_geo=config["target_geo"],
            sample_geo=config["sample_geo"],
            weight_col=config["weight_col"],
            sampling_method=config["sampling_method"],
            new_id_base=config["new_id_base"],
        )
        geocross.to_csv(geo_crosswalk_file, index=False)

        if "/base/annual_household_control_totals" not in store_target:
            raise KeyError("Missing '/base/annual_household_control_totals' in the target HDF")
        ctotals = store_target["/base/annual_household_control_totals"]
        dict_cats = extract_control_categories(ctotals)

        simctrl = cats_to_ctrl(dict_cats, config["target_geo"], "households")
        simctrl = pd.concat(
            [
                simctrl,
                pd.DataFrame(
                    [
                        {
                            "target": "num_hh",
                            "geography": config["target_geo"],
                            "seed_table": "households",
                            "importance": 10000000,
                            "control_field": "HHBASE",
                            "expression": "(households.persons > 0)",
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )
        if config["add_person_controls"]:
            simctrl = cats_to_ctrl(
                DEFAULT_PERSON_CONTROL_CATEGORIES,
                config["target_geo"],
                "persons",
                dfctrl=simctrl,
            )
        runtime_control_columns(simctrl).to_csv(controls_file, index=False)

        households_target, persons_target = attach_households_and_persons(
            store_target,
            config["year"],
            STANDARD_GEO_UNITS,
            filter_county=config["filter_county"],
        )
        geo_ctrl = summarize_geo_controls(
            simctrl,
            seed_frames={"households": households_target, "persons": persons_target},
            target_geo=config["target_geo"],
        )
        geo_ctrl.to_csv(control_totals_file, index=False)

        households_seed, persons_seed = attach_households_and_persons(
            store_input,
            config["year"],
            STANDARD_GEO_UNITS,
            filter_county=config["filter_county"],
        )
        hhs, pps = build_seed_tables(
            hhs_total=households_seed,
            pps_total=persons_seed,
            geocross=geocross,
            geocross_dup=geocross_dup,
            sample_geo=config["sample_geo"],
            target_geo=config["target_geo"],
            sampling_method=config["sampling_method"],
            new_id_base=config["new_id_base"],
        )
        cast_int_columns(hhs, ["household_id", "hh_id", "WGTP", "SAMPLEGEO"])
        cast_int_columns(pps, ["person_id", "household_id", "hh_id", "PWGTP", "SAMPLEGEO"])
        hhs.to_csv(seed_households_file, index=False)
        pps.to_csv(seed_persons_file, index=False)

        settings_template = deepcopy(load_yaml(config["settings_template"]))
        settings = build_runtime_settings(
            settings_template,
            target_geo=config["target_geo"],
            geo_crosswalk_file=geo_crosswalk_file,
            controls_file=controls_file,
            control_totals_file=control_totals_file,
            seed_households_file=seed_households_file,
            seed_persons_file=seed_persons_file,
            households=hhs,
            persons=pps,
        )

    with open(settings_file, "w") as stream:
        yaml.dump(settings, stream, default_flow_style=False, sort_keys=False)
    shutil.copy2(config["config_path"], config_copy)

    print(f"sampling_method: {config['sampling_method']}")
    print(f"geo_cross_walk: {geo_crosswalk_file}")
    print(f"geocross: {len(geocross)}   geocross_all: {len(geocross_all)}")
    print(f"config control table: {controls_file}")
    print(f"geography control table: {control_totals_file}")
    print(f"seed households table: {seed_households_file}")
    print(f"seed persons table: {seed_persons_file}")

    return {
        "run_folder": config["run_folder"],
        "configs_dir": config["configs_dir"],
        "data_dir": config["data_dir"],
        "output_dir": config["output_dir"],
        "geo_crosswalk": geo_crosswalk_file,
        "controls": controls_file,
        "control_totals": control_totals_file,
        "seed_households": seed_households_file,
        "seed_persons": seed_persons_file,
        "settings": settings_file,
    }


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    config = load_runtime_config(args)
    outputs = build_refinement_inputs(config)
    print("created refinement input package:")
    print(f"  run folder: {outputs['run_folder']}")
    print(f"  configs: {outputs['configs_dir']}")
    print(f"  data: {outputs['data_dir']}")
    print(f"  output: {outputs['output_dir']}")
    print(
        "review handoff:\n"
        f"  edit and replace {outputs['control_totals']} in-place after review,\n"
        "  then use this package for the later refinement/synthesis run."
    )


if __name__ == "__main__":
    main()
