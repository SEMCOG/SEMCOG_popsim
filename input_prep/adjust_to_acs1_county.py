# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
import pandas as pd
import os
from input_utils import *


# %%
def csv_reader(fname, key_col):
    df = pd.read_csv(fname, index_col=key_col)
    df.drop(df.columns[0], axis=1, inplace=True)
    return df


def make_countyid(df):
    df['COUNTYID'] = df.index.astype(str).str[:5].astype(int)
    df = df.sort_index()  # important for later multiplication
    return df


# %%
def update_by_category(df_controls, county_adj):
    ind_name = df_controls.index.name
    df_sum = df_controls.groupby('COUNTYID').sum()
    df_diff = county_adj[df_sum.columns]/df_sum

    # change index to COUNTYID, prepare for calculation
    df_controls = df_controls.reset_index().set_index('COUNTYID')
    #df_controls, df_diff = df_controls.align(df_diff)
    # print(df_controls)
    df_controls[df_diff.columns] = df_controls[df_diff.columns] * df_diff
    df_controls = df_controls.reset_index().set_index(ind_name)

    return df_controls


def update_by_total(df_controls, county_totals):
    # county_totals should have county_id as index, and HHBASE, POPBASE

    ind_name = df_controls.index.name
    df_sum = df_controls.groupby('COUNTYID').sum()
    hcols = [c for c in df_sum.columns if "HH" in c]
    pcols = [c for c in df_sum.columns if c not in hcols]
    df_controls = df_controls.reset_index().set_index('COUNTYID')

    if hcols != []:
        if not 'HHBASE' in hcols:
            at = ''.join(i for i in hcols[0] if not i.isdigit())
            cols = [c for c in hcols if c.startswith(at)]
            df_sum['HHBASE'] = df_sum[cols].sum(axis=1)
        #df_controls, tt = df_controls.align(tt, axis=1)
        diff = county_totals['HHBASE']/df_sum['HHBASE']
        df_controls[hcols] = df_controls[hcols].apply(
            lambda x: x * diff, axis=0)

    if pcols != []:
        if not 'POPBASE' in pcols:
            at = ''.join(i for i in hcols[0] if not i.isdigit())
            cols = [c for c in pcols if c.startswith(at)]
            df_sum['POPBASE'] = df_sum[cols].sum(axis=1)
        diff = county_totals['POPBASE']/df_sum['POPBASE']
        df_controls[pcols] = df_controls[pcols].apply(
            lambda x: x * diff, axis=0)

    df_controls = df_controls.reset_index().set_index(ind_name)

    return df_controls


def intergerize(ds):
    add_count = int(ds.sum().round() - (ds//1).sum())
    add_ind = (ds % 1).nlargest(add_count).index
    ds = ds//1
    ds.loc[add_ind] = ds//1 + 1

    return ds


def integerize_df(df_controls, int_col):
    # int_col is the columns for aggregation, integrized sum should add up to the original total of this column
    for grp, dfg in df_controls.groupby(int_col):
        for col in dfg.columns:
            df_controls.loc[dfg.index, col] = intergerize(dfg[col])
    return df_controls


def prepare(csv_file, geoid):
    df = csv_reader(csv_file, geoid)
    return make_countyid(df)


def adjust_by_county_cat(geo_control, cnty_control):
    geo_control = update_by_category(geo_control, cnty_control)
    geo_control = integerize_df(geo_control, 'COUNTYID')
    return geo_control


def adjust_by_county_total(geo_control, cnty_control):
    geo_control = update_by_total(geo_control, cnty_control)
    geo_control = integerize_df(geo_control, 'COUNTYID')
    return geo_control


# %%
bg_file = "2019/data/SEMCOG_2019_control_totals_blkgrp.csv"
bg_key = "BLKGRPID"
trt_file = "2019/data/SEMCOG_2019_control_totals_tract.csv"
trt_key = "TRACTID"
cnty_file = "2019/data/SEMCOG_2019_control_totals_county_adj.csv"
cnty_key = "COUNTYID"


# %%
# adjust by all categories
df_cnty = prepare(cnty_file, cnty_key)
df_bg = prepare(bg_file, bg_key)
marginal_summary(df_bg)
df_bg = adjust_by_county_cat(df_bg, df_cnty)
marginal_summary(df_bg)


# %%
df_trt = prepare(trt_file, trt_key)
marginal_summary(df_trt)
df_trt = adjust_by_county_cat(df_trt, df_cnty)
marginal_summary(df_trt)

# %%
df_bg.to_csv("2019/data/SEMCOG_2019_control_totals_blkgrp_adj.csv")
df_trt.to_csv("2019/data/SEMCOG_2019_control_totals_tract_adj.csv")


# %%
%cd "d: \projects\populationsim\SEMCOG_popsim\input_prep\"

# %%
# adjust by county totals

# %%
df_cnty = prepare(cnty_file, cnty_key)
df_bg = prepare(bg_file, bg_key)
marginal_summary(df_bg)
df_bg = adjust_by_county_total(df_bg, df_cnty)
print('\n after ___________')
marginal_summary(df_bg)
# %%
df_trt = prepare(trt_file, trt_key)
marginal_summary(df_trt)
df_trt = adjust_by_county_cat(df_trt, df_cnty)
marginal_summary(df_trt)

# %%
