import argparse
from pathlib import Path

import pandas as pd

from input_utils import marginal_summary


def csv_reader(fname, key_col):
    df = pd.read_csv(fname, index_col=key_col)
    df.drop(df.columns[0], axis=1, inplace=True)
    return df


def make_countyid(df):
    df["COUNTYID"] = df.index.astype(str).str[:5].astype(int)
    df = df.sort_index()
    return df


def update_by_category(df_controls, county_adj):
    index_name = df_controls.index.name
    df_sum = df_controls.groupby("COUNTYID").sum()
    df_diff = county_adj[df_sum.columns] / df_sum
    df_controls = df_controls.reset_index().set_index("COUNTYID")
    df_controls[df_diff.columns] = df_controls[df_diff.columns] * df_diff
    df_controls = df_controls.reset_index().set_index(index_name)
    return df_controls


def update_by_total(df_controls, county_totals):
    index_name = df_controls.index.name
    df_sum = df_controls.groupby("COUNTYID").sum()
    household_cols = [col for col in df_sum.columns if "HH" in col]
    person_cols = [col for col in df_sum.columns if col not in household_cols]
    df_controls = df_controls.reset_index().set_index("COUNTYID")

    if household_cols:
        if "HHBASE" not in household_cols:
            attr = "".join(char for char in household_cols[0] if not char.isdigit())
            cols = [col for col in household_cols if col.startswith(attr)]
            df_sum["HHBASE"] = df_sum[cols].sum(axis=1)
        diff = county_totals["HHBASE"] / df_sum["HHBASE"]
        df_controls[household_cols] = df_controls[household_cols].apply(lambda x: x * diff, axis=0)

    if person_cols:
        if "POPBASE" not in person_cols:
            attr = "".join(char for char in household_cols[0] if not char.isdigit())
            cols = [col for col in person_cols if col.startswith(attr)]
            df_sum["POPBASE"] = df_sum[cols].sum(axis=1)
        diff = county_totals["POPBASE"] / df_sum["POPBASE"]
        df_controls[person_cols] = df_controls[person_cols].apply(lambda x: x * diff, axis=0)

    df_controls = df_controls.reset_index().set_index(index_name)
    return df_controls


def integerize(series):
    add_count = int(series.sum().round() - (series // 1).sum())
    add_index = (series % 1).nlargest(add_count).index
    series = series // 1
    series.loc[add_index] = series // 1 + 1
    return series


def integerize_df(df_controls, int_col):
    for _, df_group in df_controls.groupby(int_col):
        for col in df_group.columns:
            df_controls.loc[df_group.index, col] = integerize(df_group[col])
    return df_controls


def prepare(csv_file, geoid):
    df = csv_reader(csv_file, geoid)
    return make_countyid(df)


def adjust_by_county_cat(geo_control, county_control):
    geo_control = update_by_category(geo_control, county_control)
    geo_control = integerize_df(geo_control, "COUNTYID")
    return geo_control


def adjust_by_county_total(geo_control, county_control):
    geo_control = update_by_total(geo_control, county_control)
    geo_control = integerize_df(geo_control, "COUNTYID")
    return geo_control


def build_parser():
    parser = argparse.ArgumentParser(description="Adjust blockgroup and tract controls to ACS1 county controls.")
    parser.add_argument("--year", type=int, default=2019)
    parser.add_argument("--mode", choices=["category", "total"], default="category")
    parser.add_argument("--base-dir", default=".", help="Directory containing the year folder.")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    base_dir = Path(args.base_dir).resolve()
    year_dir = base_dir / str(args.year) / "data"

    bg_file = year_dir / f"SEMCOG_{args.year}_control_totals_blkgrp.csv"
    trt_file = year_dir / f"SEMCOG_{args.year}_control_totals_tract.csv"
    county_file = year_dir / f"SEMCOG_{args.year}_control_totals_county_adj.csv"
    output_suffix = "adj"

    df_county = prepare(county_file, "COUNTYID")
    df_bg = prepare(bg_file, "BLKGRPID")
    df_trt = prepare(trt_file, "TRACTID")

    print("\nBefore blockgroup adjustment")
    marginal_summary(df_bg)
    print("\nBefore tract adjustment")
    marginal_summary(df_trt)

    if args.mode == "category":
        df_bg = adjust_by_county_cat(df_bg, df_county)
        df_trt = adjust_by_county_cat(df_trt, df_county)
    else:
        df_bg = adjust_by_county_total(df_bg, df_county)
        df_trt = adjust_by_county_total(df_trt, df_county)
        output_suffix = "total_adj"

    print("\nAfter blockgroup adjustment")
    marginal_summary(df_bg)
    print("\nAfter tract adjustment")
    marginal_summary(df_trt)

    bg_output = year_dir / f"SEMCOG_{args.year}_control_totals_blkgrp_{output_suffix}.csv"
    trt_output = year_dir / f"SEMCOG_{args.year}_control_totals_tract_{output_suffix}.csv"
    df_bg.to_csv(bg_output)
    df_trt.to_csv(trt_output)

    print(f"Saved {bg_output}")
    print(f"Saved {trt_output}")


if __name__ == "__main__":
    main()
