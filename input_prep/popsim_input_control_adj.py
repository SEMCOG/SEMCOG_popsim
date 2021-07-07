# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
# This program download ACS 1-year control file for adjustment purpose
# Syntax:  > python popsim_input_control_adj.py key adj.yaml (key: Census API key; yaml: a yaml configuration file such as region.yaml)

# Inputs:
# [year]/region_[year]_adj.yaml  (yaml config input for input_maker.py)
# [year]/control_pre_[year]_adj.csv   (control file with all columns the same as 5-year but geography values are COUNTYID)

# outputs:
# [year]/data/SEMCOG_[year]_control_totals_county_adj.csv   (marginal controls by geography)


# %%
import os
import re
import time
import pandas as pd
from census import Census
import oyaml as yaml
from collections import defaultdict
from input_utils import *
import argparse
import shutil

# %%
parser = argparse.ArgumentParser()
parser.add_argument("key", help="Census API key")
parser.add_argument("yaml", help="yaml configuration file name")
args = parser.parse_args()
t0 = time.time()

# %%
conf = yaml.load(open("./" + args.yaml, "r"), Loader=yaml.Loader)
#conf = yaml.load(open("./region_2019.yaml", "r"), Loader=yaml.Loader)

# %%
prj = conf["project"]
prj_name = prj["name"]
target = prj["target"]
acs_year = prj["acs_year"]
acs_sample = prj["acs_sample"]
prj_folder = f"{acs_year}/"
pre_control = prj_folder + prj["pre_control"].format(str(acs_year))

geo = conf["geography"]
state = geo["state"][0]
counties = geo["counties"]

output_folder = prj_folder + "data/"
output_control = "{}_{}_control_totals_.csv".format(prj_name, str(acs_year))
print(f"\n *** download {target} for year {acs_year} ***")


# %% [markdown]
# step 2. download and compile maginal control files
c = Census(args.key, year=acs_year)
cs = eval(f'c.{acs_sample}')  # set API type

# %%
# download Census marginal controls
# Census controls variables are defined in "controls_pre" table(popsim "controls" table + "acs_variables" field )
# "acs_variables" contains Census variables and expressions
#
print("\nMaking popsim adjustment controls: ", acs_year)
print("  downloading Census variables ...")
dfc = pd.read_csv(pre_control)
dic_margs = {}
for geo, dfgeo in dfc.groupby("geography"):
    full_vars = list(
        set(
            re.findall(r"[B-C][0-9]{5}[A-Z]{0,1}_[0-9]{3}E", str(list(dfgeo.acs_variables)))
        )
    )
    if geo == "BLKGRP":
        ac5 = Census_Downloader(cs, state, counties, "*", "*")
        geo_cols = ["state", "county", "tract", "block group"]
    elif geo == "TRACT":
        ac5 = Census_Downloader(cs, state, counties, "*")
        geo_cols = ["state", "county", "tract"]
    elif geo == "COUNTY":
        ac5 = Census_Downloader(cs, state, counties)
        geo_cols = ["state", "county"]

    print("\t" + geo + " marginals ")
    dic_margs[geo] = ac5.download(full_vars).set_index(geo_cols)
    if "GEO_ID" in dic_margs[geo].columns:
        # new Census API downloads an extra "GEO_ID" column
        dic_margs[geo].drop('GEO_ID', axis=1, inplace=True)

# %%
# Compute popsim control variables from Census marginals
print("  compiling popsim control fields ...")
for geo, dfg in dfc.groupby("geography"):
    dic_margs[geo] = dic_margs[geo].astype(float).fillna(0)
    for ind, r in dfg.iterrows():
        # print(r['control_field'] + ":  " + r['acs_variables'])
        dic_margs[geo][r["control_field"]] = dic_margs[geo].eval(
            r["acs_variables"].replace('"', "")
        )
    dic_margs[geo] = dic_margs[geo][list(dfg.control_field)]  # keep only control fields


# %%
# add unique geoids
ctr_geos = {}
mapping_dict = {4: 'BLKGRPID', 3: "TRACTID", 2: "COUNTYID"}
for geo, dfm in dic_margs.items():
    dfm[mapping_dict[dfm.index.nlevels]] = dfm.index.map(''.join)
    dfm.reset_index(drop=True, inplace=True)

    dfm.fillna(0, inplace=True)
    dfm.columns = [col.upper() for col in dfm.columns]
    if "HHBASE" in dfm.columns:
        dfm = dfm.loc[dfm.HHBASE > 0]
    marginal_summary(dfm)

    f_output_control = output_control.replace(".csv", geo.lower() + "_adj.csv")
    ctr_geos[geo] = f_output_control
    print("  saving control file: " + output_folder + f_output_control)
    dfm.to_csv(output_folder + f_output_control)
