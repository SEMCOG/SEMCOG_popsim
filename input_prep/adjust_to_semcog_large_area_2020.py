# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
import pandas as pd
import os
from input_utils import *

#%%

####### === #####
## this step is to correct wrong SEMMCD and assign correct large area id. it is optional.
##create new large area ids, when SEMMCD is incorrect
## all three levels, BG, BLOCK and TRACT are consistent and cover same area, even a couple of blocks might belong to other MCD

# df_geo = pd.read_csv("2020/data/Census2020_Tract_BG_BLK_000.csv") #input table
# df_geo.rename(columns= {"SEMMCD":"SEMMCD_incorrect"}, inplace=True)
# df_geo["LARGE_AREA_ID"] = df_geo.GEOID20.astype(str).str[2:5].astype(int)
# df_geo.loc[df_geo.LARGE_AREA_ID == 163, "LARGE_AREA_ID"] = 3
# df_geo["GEOID20_BG"] = -1
# df_geo.loc[df_geo.GEOTYPE == "BLOCK", "GEOID20_BG"] = df_geo.GEOID20//1000

# df_det = pd.read_csv("2020/data/detroit_2020_blockgroups.csv")
# df_det["TRACTID"] = df_det.GEOID20//10


# df_geo.loc[(df_geo.GEOID20.isin(df_det.GEOID20)) & (df_geo.GEOTYPE == "BLOCKGROUP"), "LARGE_AREA_ID"] = 5
# df_geo.loc[(df_geo.GEOID20.isin(df_det.TRACTID)) & (df_geo.GEOTYPE == "TRACT"), "LARGE_AREA_ID"] = 5
# df_geo.loc[(df_geo.GEOID20_BG.isin(df_det.GEOID20)) & (df_geo.GEOTYPE == "BLOCK"), "LARGE_AREA_ID"] = 5

# print(df_geo.loc[df_geo.GEOTYPE == "TRACT"].LARGE_AREA_ID.value_counts())
# print(df_geo.loc[df_geo.GEOTYPE == "BLOCKGROUP"].LARGE_AREA_ID.value_counts())
# print(df_geo.loc[df_geo.GEOTYPE == "BLOCK"].LARGE_AREA_ID.value_counts())

# df_geo.to_csv("2020/data/Census2020_Tract_BG_LARGEAREA.csv") #output table

####### === #####


### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ###
### start here 2020 control updates
# %%
def csv_reader(fname, key_col):
    df = pd.read_csv(fname, index_col=key_col)
    df.drop(df.columns[0], axis=1, inplace=True)
    return df


def make_countyid(df):
    df["COUNTYID"] = df.index.astype(str).str[:5].astype(int)
    df = df.sort_index()  # important for later multiplication
    return df


def make_large_area_id(df, df_geo):
    # df contains a equivalant geo id as index, such as track, block group etc
    # df_geo has all equivalant geo id, SEMMCD and County ID
    df[["SEMMCD", "COUNTY"]] = df_geo[["SEMMCD", "COUNTY"]]
    df["LARGE_AREA_ID"] = df["COUNTY"]
    df.loc[df.COUNTY == 163, "LARGE_AREA_ID"] = 3
    df.loc[df.SEMMCD == 5, "LARGE_AREA_ID"] = 5
    df.drop(["SEMMCD", "COUNTY"], axis=1, inplace=True)
    return df


def ctrl_ratio_by_geo(df, col_val, col_geo, df_ctrl):
    ds_sum = df.groupby([col_geo])[col_val].sum()
    return df_ctrl[col_val] / ds_sum


def col_update_subgroup_ratio(ds, subgrp_ratio, lev_name):
    # subgrp_ratio: ratio of new/old subgroup totals
    # ds: series with index/multiindex, one level has subgroup ID
    # lev_name: level name
    # return pd series with decimal values
    return subgrp_ratio.mul(ds, level=lev_name, axis=0)


def col_update_subgroup_totals(df, df_la, col_name):
    # use ratios between large area controls to existing sum to update cell values
    # df_la: large area controls/new subgroup totals
    # col_name: same name for df value column and df_la control column
    # df: dataframe with value column and LARGE_AREA_ID/subgroup ID, index is the detailed geo id
    # return pd series with decimal values
    ds_sum = df.groupby(["LARGE_AREA_ID"])[col_name].sum()
    ds_ratio = df_la[col_name] / ds_sum
    df.set_index("LARGE_AREA_ID", append=True, inplace=True)
    return ds_ratio.mul(df[col_name], level=1, axis=0)


def adjust_to_new_totals(df, ds_total):
    return df.div(df.sum(axis=1), axis=0).mul(ds_total, axis=0)


def intergerize(ds):
    add_count = int(ds.sum().round() - (ds // 1).sum())
    add_ind = (ds % 1).nlargest(add_count).index
    ds = ds // 1
    ds.loc[add_ind] = ds // 1 + 1
    return ds


def total_by_dict(df, dict):
    for k, v in dict.items():
        if v[0] in df.columns:
            print(k, df[v].sum().sum())


# step 1. prepare all inputs
# %%
# BLKGRP and TRACT marginals are from 2020 5-year ACS
# BLK HHs and HHpop  Census 2020 PL94 (all we have right now)
dic_files = {
    "TRACT": ["2020/data/SEMCOG_2020_control_totals_tract.csv", "TRACTID"],
    "BLOCKGROUP": ["2020/data/SEMCOG_2020_control_totals_blkgrp.csv", "BLKGRPID"],
    "BLOCK": ["2020/data/SEMCOG_2020_control_totals_blk.csv", "BLKID"],
}
la_file = "2020/data/final_hhpop2020_largearea.csv"
la_key = "LARGE_AREA_ID"
geo = "2020/data/Census2020_Tract_BG_LARGEAREA.csv"

# read and process inputs
df_geo = pd.read_csv(geo)
df_la = pd.read_csv(la_file, index_col=la_key) * 1.0

#%%
df = {}
for k, v in dic_files.items():
    print("\n" + k)
    df[k] = csv_reader(v[0], v[1])
    dfg = df_geo[df_geo.GEOTYPE == k].set_index("GEOID20")
    df[k]["LARGE_AREA_ID"] = dfg["LARGE_AREA_ID"]
    print(df[k].head(2))
    # df[k].drop([])


# step 2. adjust to Census 2020 BLK HH and HHpop
# %%
# replace blkgrp HH and HHpop with BLK numbers
df["BLOCK"]["BLOCKGROUP"] = df["BLOCK"].index.astype(str).str[:-3].astype(int)
df_blk2grp = df["BLOCK"].groupby("BLOCKGROUP")[["HHBASE", "POPBASE"]].sum()
df["BLOCKGROUP"][["HHBASE_old", "POPBASE_old"]] = df["BLOCKGROUP"][
    ["HHBASE", "POPBASE"]
]
df["BLOCKGROUP"][["HHBASE", "POPBASE"]] = df_blk2grp[["HHBASE", "POPBASE"]]

df["BLOCKGROUP"]["TRACTID"] = df["BLOCKGROUP"].index.astype(str).str[:-1].astype(int)
df_grp2trt = (
    df["BLOCKGROUP"]
    .groupby("TRACTID")[["HHBASE_old", "POPBASE_old", "HHBASE", "POPBASE"]]
    .sum()
)
df["TRACT"][["HHBASE_old", "POPBASE_old", "HHBASE", "POPBASE"]] = df_grp2trt[
    ["HHBASE_old", "POPBASE_old", "HHBASE", "POPBASE"]
]

# %%
dict_adj_total = {
    "HHAGE": ["HHAGE1", "HHAGE2", "HHAGE3", "HHAGE4"],
    "HHRACE": ["HHRACE1", "HHRACE2", "HHRACE3", "HHRACE4"],
    "HHHISP": ["HHHISP1", "HHHISP2"],
    "HHCHD": ["HHCHD1", "HHCHD2"],
    "HHINC": ["HHINC1", "HHINC2", "HHINC3", "HHINC4"],
    "HHCAR": ["HHCAR0", "HHCAR1", "HHCAR2"],
    "HHPERSONS": [
        "HHPERSONS1",
        "HHPERSONS2",
        "HHPERSONS3",
        "HHPERSONS4",
        "HHPERSONS5",
        "HHPERSONS6",
        "HHPERSONS7",
    ],
    "HHTENURE": ["HHTENURE1", "HHTENURE0"],
    "HHWORKER": ["HHWORKER0", "HHWORKER1", "HHWORKER2"],
    "AGEP": ["AGEP1", "AGEP2", "AGEP3", "AGEP4", "AGEP5", "AGEP6"],
    "RACE": ["RACE1", "RACE2", "RACE3", "RACE4", "RACE5"],
    "SEX": ["SEX1", "SEX2"],
}
dict_adj_ratio = {
    "IND": [
        "INDUSTRY1",
        "INDUSTRY2",
        "INDUSTRY3",
        "INDUSTRY4",
        "INDUSTRY5",
        "INDUSTRY6",
        "INDUSTRY7",
        "INDUSTRY8",
        "INDUSTRY9",
        "INDUSTRY10",
        "INDUSTRY11",
        "INDUSTRY12",
        "INDUSTRY13",
        "INDUSTRY14",
    ],
    "EMP": ["EMPWORKER"],
}

#%%
# proportionally adjust subgroup numbers by each Census geo
for geo in ["BLOCKGROUP", "TRACT"]:
    for dict in [dict_adj_total, dict_adj_ratio]:
        total_by_dict(df[geo], dict)

for k, v in dict_adj_total.items():
    for cen_geo in ["BLOCKGROUP", "TRACT"]:
        if v[0] in df[cen_geo].columns:
            print(k)
            total = "HHBASE" if k[:2] == "HH" else "POPBASE"
            df_new = adjust_to_new_totals(
                df[cen_geo][v], df[cen_geo][total]
            )  # proportion adjust to new total "HHBASE" or "POPBASE"
            df_new = df_new.apply(
                lambda x: intergerize(x)
            )  # integerize all rows and meet row totals
            df[cen_geo][v] = df_new[v]


for k, v in dict_adj_ratio.items():
    for cen_geo in ["BLOCKGROUP", "TRACT"]:
        if v[0] in df[cen_geo].columns:
            print(k)
            total = "HHBASE" if k[:2] == "HH" else "POPBASE"
            new_ratio = df[cen_geo][total] / df[cen_geo][total + "_old"]
            df[cen_geo][v] = df[cen_geo][v].mul(new_ratio, axis=0)

df["BLOCKGROUP"].fillna(0, inplace=True)
df["TRACT"].fillna(0, inplace=True)

#%%
df["BLOCKGROUP"].drop(["HHBASE_old", "POPBASE_old", "TRACTID"], axis=1, inplace=True)
df["TRACT"].drop(
    ["HHBASE_old", "POPBASE_old", "HHBASE", "POPBASE"], axis=1, inplace=True
)

df["BLOCKGROUP"].to_csv("2020/data/SEMCOG_2020_control_totals_blkgrp_census.csv")
df["TRACT"].to_csv("2020/data/SEMCOG_2020_control_totals_tract_census.csv")

for geo in ["BLOCKGROUP", "TRACT"]:
    for dict in [dict_adj_total, dict_adj_ratio]:
        total_by_dict(df[geo], dict)

# step 3. adjust to SEMCOG large area controls
#%%
la_pop_ratios = df_la.POPBASE / df["BLOCKGROUP"].groupby("LARGE_AREA_ID").POPBASE.sum()

# %%
dict_update_cat = {
    "BLOCKGROUP": [
        "POPBASE",
        "AGEP1",
        "AGEP2",
        "AGEP3",
        "AGEP4",
        "AGEP5",
        "AGEP6",
        "RACE1",
        "RACE2",
        "RACE3",
        "RACE4",
        "RACE5",
    ]
}

# adjust by each attribute and large area totals
for geo, cols in dict_update_cat.items():
    for col in cols:
        print(geo, col)
        ratios = ctrl_ratio_by_geo(df[geo], col, "LARGE_AREA_ID", df_la)
        ind_name = df[geo].index.name
        df[geo].set_index("LARGE_AREA_ID", append=True, inplace=True)
        newds = col_update_subgroup_ratio(df[geo][col], ratios, "LARGE_AREA_ID")
        df[geo][col] = newds
        df[geo] = df[geo].reset_index().set_index(ind_name)


# %%
dict_update_tot = {
    "BLOCKGROUP": ["SEX1", "SEX2"],
    "TRACT": [
        "INDUSTRY1",
        "INDUSTRY2",
        "INDUSTRY3",
        "INDUSTRY4",
        "INDUSTRY5",
        "INDUSTRY6",
        "INDUSTRY7",
        "INDUSTRY8",
        "INDUSTRY9",
        "INDUSTRY10",
        "INDUSTRY11",
        "INDUSTRY12",
        "INDUSTRY13",
        "INDUSTRY14",
        "EMPWORKER",
    ],
}
for geo, cols in dict_update_tot.items():
    print(geo)
    df[geo].set_index("LARGE_AREA_ID", append=True, inplace=True)
    for col in cols:
        print(col)
        df[geo][col] = la_pop_ratios.mul(df[geo][col], level=1, axis=0)
    df[geo] = df[geo].droplevel("LARGE_AREA_ID")
    print(df[geo].head(2))


# df["BLOCK"].to_csv("2020/data/SEMCOG_2020_control_totals_blk_update.csv")
df["BLOCKGROUP"].to_csv("2020/data/SEMCOG_2020_control_totals_blkgrp_update.csv")
df["TRACT"].to_csv("2020/data/SEMCOG_2020_control_totals_tract_update.csv")

#%%
for geo in ["BLOCKGROUP", "TRACT"]:
    for dict in [dict_adj_total, dict_adj_ratio]:
        total_by_dict(df[geo], dict)

#######


#%%
df["TRACT"].isnull().sum().sum()

#%%
df["BLOCKGROUP"].isnull().sum().sum()


# %%
# {def update_by_category(df_controls, county_adj):
#     ind_name = df_controls.index.name
#     df_sum = df_controls.groupby('COUNTYID').sum()
#     df_diff = county_adj[df_sum.columns]/df_sum

#     # change index to COUNTYID, prepare for calculation
#     df_controls = df_controls.reset_index().set_index('COUNTYID')
#     #df_controls, df_diff = df_controls.align(df_diff)
#     # print(df_controls)
#     df_controls[df_diff.columns] = df_controls[df_diff.columns] * df_diff
#     df_controls = df_controls.reset_index().set_index(ind_name)

#     return df_controls


# def update_by_total(df_controls, county_totals):
#     # county_totals should have county_id as index, and HHBASE, POPBASE

#     ind_name = df_controls.index.name
#     df_sum = df_controls.groupby('COUNTYID').sum()
#     hcols = [c for c in df_sum.columns if "HH" in c]
#     pcols = [c for c in df_sum.columns if c not in hcols]
#     df_controls = df_controls.reset_index().set_index('COUNTYID')

#     if hcols != []:
#         if not 'HHBASE' in hcols:
#             at = ''.join(i for i in hcols[0] if not i.isdigit())
#             cols = [c for c in hcols if c.startswith(at)]
#             df_sum['HHBASE'] = df_sum[cols].sum(axis=1)
#         #df_controls, tt = df_controls.align(tt, axis=1)
#         diff = county_totals['HHBASE']/df_sum['HHBASE']
#         df_controls[hcols] = df_controls[hcols].apply(
#             lambda x: x * diff, axis=0)

#     if pcols != []:
#         if not 'POPBASE' in pcols:
#             at = ''.join(i for i in hcols[0] if not i.isdigit())
#             cols = [c for c in pcols if c.startswith(at)]
#             df_sum['POPBASE'] = df_sum[cols].sum(axis=1)
#         diff = county_totals['POPBASE']/df_sum['POPBASE']
#         df_controls[pcols] = df_controls[pcols].apply(
#             lambda x: x * diff, axis=0)

#     df_controls = df_controls.reset_index().set_index(ind_name)

#     return df_controls


# def integerize_df(df_controls, int_col):
#     # int_col is the columns for aggregation, integrized sum should add up to the original total
#     # of this column
#     for grp, dfg in df_controls.groupby(int_col):
#         for col in dfg.columns:
#             df_controls.loc[dfg.index, col] = intergerize(dfg[col])
#     return df_controls


# def prepare(csv_file, geoid):
#     df = csv_reader(csv_file, geoid)
#     return make_countyid(df)


# def adjust_by_county_cat(geo_control, cnty_control):
#     geo_control = update_by_category(geo_control, cnty_control)
#     geo_control = integerize_df(geo_control, 'COUNTYID')
#     return geo_control


# def adjust_by_county_total(geo_control, cnty_control):
#     geo_control = update_by_total(geo_control, cnty_control)
#     geo_control = integerize_df(geo_control, 'COUNTYID')
#     return geo_control}
