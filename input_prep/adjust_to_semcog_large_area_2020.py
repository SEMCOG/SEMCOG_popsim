import argparse
from pathlib import Path

import pandas as pd


def csv_reader(fname, key_col):
    df = pd.read_csv(fname, index_col=key_col)
    df.drop(df.columns[0], axis=1, inplace=True)
    return df


def ctrl_ratio_by_geo(df, col_val, col_geo, df_ctrl):
    ds_sum = df.groupby([col_geo])[col_val].sum()
    return df_ctrl[col_val] / ds_sum


def adjust_to_new_totals(df, ds_total):
    return df.div(df.sum(axis=1), axis=0).mul(ds_total, axis=0)


def integerize(series):
    add_count = int(series.sum().round() - (series // 1).sum())
    add_index = (series % 1).nlargest(add_count).index
    series = series // 1
    series.loc[add_index] = series // 1 + 1
    return series


def total_by_dict(df, groups):
    for key, values in groups.items():
        if values[0] in df.columns:
            print(key, df[values].sum().sum())


def build_parser():
    parser = argparse.ArgumentParser(description="Adjust 2020 controls to SEMCOG large-area totals.")
    parser.add_argument("--base-dir", default=".", help="Directory containing the year folder.")
    parser.add_argument("--year", type=int, default=2020)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    base_dir = Path(args.base_dir).resolve()
    year_dir = base_dir / str(args.year) / "data"

    input_files = {
        "TRACT": [year_dir / f"SEMCOG_{args.year}_control_totals_tract.csv", "TRACTID"],
        "BLOCKGROUP": [year_dir / f"SEMCOG_{args.year}_control_totals_blkgrp.csv", "BLKGRPID"],
        "BLOCK": [year_dir / f"SEMCOG_{args.year}_control_totals_blk.csv", "BLKID"],
    }
    large_area_file = year_dir / "final_hhpop2020_largearea.csv"
    large_area_key = "LARGE_AREA_ID"
    geography_file = year_dir / "Census2020_Tract_BG_LARGEAREA.csv"

    df_geo = pd.read_csv(geography_file)
    df_large_area = pd.read_csv(large_area_file, index_col=large_area_key) * 1.0

    data = {}
    for geo_name, (filename, key) in input_files.items():
        print(f"\n{geo_name}")
        data[geo_name] = csv_reader(filename, key)
        df_geo_type = df_geo[df_geo.GEOTYPE == geo_name].set_index("GEOID20")
        data[geo_name]["LARGE_AREA_ID"] = df_geo_type["LARGE_AREA_ID"]
        print(data[geo_name].head(2))

    data["BLOCK"]["BLOCKGROUP"] = data["BLOCK"].index.astype(str).str[:-3].astype(int)
    df_block_to_group = data["BLOCK"].groupby("BLOCKGROUP")[["HHBASE", "POPBASE"]].sum()
    data["BLOCKGROUP"][["HHBASE_old", "POPBASE_old"]] = data["BLOCKGROUP"][["HHBASE", "POPBASE"]]
    data["BLOCKGROUP"][["HHBASE", "POPBASE"]] = df_block_to_group[["HHBASE", "POPBASE"]]

    data["BLOCKGROUP"]["TRACTID"] = data["BLOCKGROUP"].index.astype(str).str[:-1].astype(int)
    df_group_to_tract = data["BLOCKGROUP"].groupby("TRACTID")[["HHBASE_old", "POPBASE_old", "HHBASE", "POPBASE"]].sum()
    data["TRACT"][["HHBASE_old", "POPBASE_old", "HHBASE", "POPBASE"]] = df_group_to_tract[["HHBASE_old", "POPBASE_old", "HHBASE", "POPBASE"]]

    dict_adj_total = {
        "HHAGE": ["HHAGE1", "HHAGE2", "HHAGE3", "HHAGE4"],
        "HHRACE": ["HHRACE1", "HHRACE2", "HHRACE3", "HHRACE4"],
        "HHHISP": ["HHHISP1", "HHHISP2"],
        "HHCHD": ["HHCHD1", "HHCHD2"],
        "HHINC": ["HHINC1", "HHINC2", "HHINC3", "HHINC4"],
        "HHCAR": ["HHCAR0", "HHCAR1", "HHCAR2"],
        "HHPERSONS": ["HHPERSONS1", "HHPERSONS2", "HHPERSONS3", "HHPERSONS4", "HHPERSONS5", "HHPERSONS6", "HHPERSONS7"],
        "HHTENURE": ["HHTENURE1", "HHTENURE0"],
        "HHWORKER": ["HHWORKER0", "HHWORKER1", "HHWORKER2"],
        "AGEP": ["AGEP1", "AGEP2", "AGEP3", "AGEP4", "AGEP5", "AGEP6"],
        "RACE": ["RACE1", "RACE2", "RACE3", "RACE4", "RACE5"],
        "SEX": ["SEX1", "SEX2"],
    }
    dict_adj_ratio = {
        "IND": [
            "INDUSTRY1", "INDUSTRY2", "INDUSTRY3", "INDUSTRY4", "INDUSTRY5", "INDUSTRY6", "INDUSTRY7",
            "INDUSTRY8", "INDUSTRY9", "INDUSTRY10", "INDUSTRY11", "INDUSTRY12", "INDUSTRY13", "INDUSTRY14",
        ],
        "EMP": ["EMPWORKER"],
    }

    for geo_name in ["BLOCKGROUP", "TRACT"]:
        for groups in [dict_adj_total, dict_adj_ratio]:
            total_by_dict(data[geo_name], groups)

    for key, values in dict_adj_total.items():
        for geo_name in ["BLOCKGROUP", "TRACT"]:
            if values[0] in data[geo_name].columns:
                print(key)
                total_col = "HHBASE" if key[:2] == "HH" else "POPBASE"
                df_new = adjust_to_new_totals(data[geo_name][values], data[geo_name][total_col])
                df_new = df_new.apply(lambda x: integerize(x))
                data[geo_name][values] = df_new[values]

    for key, values in dict_adj_ratio.items():
        for geo_name in ["BLOCKGROUP", "TRACT"]:
            if values[0] in data[geo_name].columns:
                print(key)
                total_col = "HHBASE" if key[:2] == "HH" else "POPBASE"
                new_ratio = data[geo_name][total_col] / data[geo_name][total_col + "_old"]
                data[geo_name][values] = data[geo_name][values].mul(new_ratio, axis=0)

    data["BLOCKGROUP"].fillna(0, inplace=True)
    data["TRACT"].fillna(0, inplace=True)
    data["BLOCKGROUP"].drop(["HHBASE_old", "POPBASE_old", "TRACTID"], axis=1, inplace=True)
    data["TRACT"].drop(["HHBASE_old", "POPBASE_old", "HHBASE", "POPBASE"], axis=1, inplace=True)

    data["BLOCKGROUP"].to_csv(year_dir / f"SEMCOG_{args.year}_control_totals_blkgrp_census.csv")
    data["TRACT"].to_csv(year_dir / f"SEMCOG_{args.year}_control_totals_tract_census.csv")

    for geo_name in ["BLOCKGROUP", "TRACT"]:
        for groups in [dict_adj_total, dict_adj_ratio]:
            total_by_dict(data[geo_name], groups)

    large_area_pop_ratios = df_large_area.POPBASE / data["BLOCKGROUP"].groupby("LARGE_AREA_ID").POPBASE.sum()
    dict_update_cat = {
        "BLOCKGROUP": ["POPBASE", "AGEP1", "AGEP2", "AGEP3", "AGEP4", "AGEP5", "AGEP6", "RACE1", "RACE2", "RACE3", "RACE4", "RACE5"]
    }
    for geo_name, columns in dict_update_cat.items():
        for col in columns:
            print(geo_name, col)
            ratios = ctrl_ratio_by_geo(data[geo_name], col, "LARGE_AREA_ID", df_large_area)
            index_name = data[geo_name].index.name
            data[geo_name].set_index("LARGE_AREA_ID", append=True, inplace=True)
            data[geo_name][col] = ratios.mul(data[geo_name][col], level=1, axis=0)
            data[geo_name] = data[geo_name].reset_index().set_index(index_name)

    dict_update_tot = {
        "BLOCKGROUP": ["SEX1", "SEX2"],
        "TRACT": [
            "INDUSTRY1", "INDUSTRY2", "INDUSTRY3", "INDUSTRY4", "INDUSTRY5", "INDUSTRY6", "INDUSTRY7",
            "INDUSTRY8", "INDUSTRY9", "INDUSTRY10", "INDUSTRY11", "INDUSTRY12", "INDUSTRY13", "INDUSTRY14", "EMPWORKER",
        ],
    }
    for geo_name, columns in dict_update_tot.items():
        print(geo_name)
        data[geo_name].set_index("LARGE_AREA_ID", append=True, inplace=True)
        for col in columns:
            print(col)
            data[geo_name][col] = large_area_pop_ratios.mul(data[geo_name][col], level=1, axis=0)
        data[geo_name] = data[geo_name].droplevel("LARGE_AREA_ID")
        print(data[geo_name].head(2))

    data["BLOCKGROUP"].to_csv(year_dir / f"SEMCOG_{args.year}_control_totals_blkgrp_update.csv")
    data["TRACT"].to_csv(year_dir / f"SEMCOG_{args.year}_control_totals_tract_update.csv")

    for geo_name in ["BLOCKGROUP", "TRACT"]:
        for groups in [dict_adj_total, dict_adj_ratio]:
            total_by_dict(data[geo_name], groups)

    print(data["TRACT"].isnull().sum().sum())
    print(data["BLOCKGROUP"].isnull().sum().sum())


if __name__ == "__main__":
    main()
