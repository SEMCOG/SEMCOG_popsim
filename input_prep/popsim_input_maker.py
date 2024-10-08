# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'


# %%
# This program prepares all inputs files needed by PopulationSim(RSG), including configuration, geo_crosswalk, Census marginals and PUMS sample households and persons
# Syntax:  > python popsim_input_maker.py key yaml (key: Census API key; yaml: a yaml configuration file such as region.yaml)

# Inputs:
# [year]/region_[year].yaml  (yaml config input for input_maker.py)
# [year]/settings_[year].yaml   (a generic setting file for populationsim)
# [year]/control_pre_[year].csv   (control file same as populationsim, with additional "acs_variables" column)
# PUMS/xxxxhxx.csv   (PUMS household samples)
# PUMS/xxxxpxx.csv   (PUMS person samples)
# geo/tract0010_to_puma.csv   (tract00 and tract10 to PUMA crosswalk file)
# geo/tract10_detroit_city.csv    (tract10 to city_id, including Detroit neighborhood, equiv file)
# geo/2010_Census_Tract_to_2010_PUMA.txt    (tract2010 to PUMA)

# outputs:
# [year]/data/SEMCOG_[year]_control_totals_[geo].csv   (marginal controls by geography)
# [year]/data/SEMCOG_[year]_controls.csv   (same as control_pre)
# [year]/data/SEMCOG_[year]_geo_cross_walk.csv   (geo cross walk)
# [year]/data/SEMCOG_[year]_seed_households.csv   (PUMS HH samples selected by synthesizing region)
# [year]/data/SEMCOG_[year]_seed_persons.csv   (PUMS person samples selected by synthesizing region)
# [year]/data/SEMCOG_[year]_settings.yaml   (modified setting file)
# %%
import os
#change working directory
os.chdir("/home/da/populationsim/SEMCOG_popsim/input_prep")

#%%
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
# set up project information from yaml
conf = yaml.load(open("./2022/prepare_2022.yaml", "r"), Loader=yaml.Loader)

prj = conf["project"]
prj_name = prj["name"]
target = prj["target"]
acs_year = prj["acs_year"]
acs_sample = prj["acs_sample"]  # acs5 or acs1
prj_folder = f"{acs_year}/"
settings_file = prj_folder + prj["settings"].format(str(acs_year))
pre_control = prj_folder + prj["pre_control"].format(str(acs_year))
h_pums_csv = prj["h_pums_csv"]
p_pums_csv = prj["p_pums_csv"]

geo = conf["geography"]
state = geo["state"][0]
counties = geo["counties"]

output_folder = prj_folder + "data/"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
output_geo_cross = "{}_{}_geo_cross_walk.csv".format(prj_name, str(acs_year))
output_control = "{}_{}_control_totals_.csv".format(prj_name, str(acs_year))
output_seed_hhs = "{}_{}_seed_households.csv".format(prj_name, str(acs_year))
output_seed_persons = "{}_{}_seed_persons.csv".format(prj_name, str(acs_year))

print(f"\n *** preparing {target} synthesisdata for year {acs_year} ***")

# %% [markdown]
# step 1. make geographic cross work file
if (acs_year<2010) or (acs_year>=2026):
    print("synthesis year should be between 2010 and 2026")
    exit()

c = Census(args.key, year=acs_year)

# %%
print(
    f"\nPreparing Census geographies crosswalk: \n\tstate: {state}  \n\tcounty: {counties}"
)
# download Census BGs for this region
acgeo = Census_Downloader(c.acs5, state, counties, "*", "*")

df_geo = pd.DataFrame.from_dict(acgeo.download("NAME")).drop("NAME", axis=1)
df_geo["tractid"] = df_geo["state"] + df_geo["county"] + df_geo["tract"]
df_geo["blkgrpid"] = df_geo["tractid"] + df_geo["block group"]
df_geo.columns = [col.upper() for col in df_geo.columns]

# %%
# depends on synthesizing year, switch to different tract-PUMA files
# %%
dict_cross = {
    ("TRACT10", "PUMA00"): "2010_Census_Tract_to_2000_PUMA_SEMCOG.csv",
    ("TRACT10", "PUMA10"): "2010_Census_Tract_to_2010_PUMA_SEMCOG.csv",
    ("TRACT20", "PUMA10"): "2020_Census_Tract_to_2010_PUMA_SEMCOG.csv",
    ("TRACT20", "PUMA20"): "2020_Census_Tract_to_2020_PUMA_SEMCOG.csv",
}

dict_acs_pums = {
    (2010,2011):[("TRACT10", "PUMA00")],
    (2012,2013,2014,2015,2016):[("TRACT10", "PUMA00"),("TRACT10", "PUMA10")],
    (2017,2018,2019):[("TRACT10", "PUMA10")],
    (2020,2021):[("TRACT20", "PUMA10")],
    (2022,2023,2024,2025,2026):[("TRACT20", "PUMA10"),("TRACT20", "PUMA20")]
}

# %%
# compile a tract-to-PUMA crosswalk file
for yrs, vals in dict_acs_pums.items():
    if acs_year in yrs:
        if len(vals) == 1:
            df_tract_puma = read_tract_puma_crosswalk(vals[0], dict_cross)
        else:
            df_tract_puma = read_tract_puma_crosswalk(vals[0], dict_cross)
            df_tract_puma1 = read_tract_puma_crosswalk(vals[1], dict_cross)
            df_tract_puma['PUMA'] = df_tract_puma['PUMA'] + df_tract_puma1['PUMA']
        df_tract_puma = df_tract_puma.reset_index()

# %%
# join tract-PUMA and creat geo cross file
df_geo_cross = pd.merge(df_geo, df_tract_puma, on="TRACTID", how="left")
df_geo_cross["REGION"] = 2
df_geo_cross = df_geo_cross[["TRACTID", "BLKGRPID", "PUMA", "COUNTYID", "REGION"]]
print("  saving geo cross walk to: " + output_folder + output_geo_cross)
df_geo_cross.to_csv(output_folder + output_geo_cross)


# %% [markdown]
# step 2. download and compile maginal control files

# %%
# download Census marginal controls
# Census controls variables are defined in "controls_pre" table(popsim "controls" table + "acs_variables" field )
# "acs_variables" contains Census variables and expressions
#

cs = eval(f"c.{acs_sample}")  # set API type

print("\nCreating popsim marginal controls: ", acs_year)
print("  downloading Census variables ...")
dfc = pd.read_csv(pre_control)
dic_margs = {}
for geo, dfgeo in dfc.groupby("geography"):
    full_vars = list(
        set(
            re.findall(
                r"[B-C][0-9]{5}[A-Z]{0,1}_[0-9]{3}E", str(list(dfgeo.acs_variables))
            )
        )
    )
    if geo == "BLKGRP":
        acd = Census_Downloader(cs, state, counties, "*", "*")
        geo_cols = ["state", "county", "tract", "block group"]
    elif geo == "TRACT":
        acd = Census_Downloader(cs, state, counties, "*")
        geo_cols = ["state", "county", "tract"]
    elif geo == "COUNTY":
        acd = Census_Downloader(cs, state, counties)
        geo_cols = ["state", "county"]

    print("\t" + geo + " marginals ")
    dic_margs[geo] = acd.download(full_vars).set_index(geo_cols)

    if "GEO_ID" in dic_margs[geo].columns:
        # new Census API downloads an extra "GEO_ID" column
        dic_margs[geo].drop("GEO_ID", axis=1, inplace=True)

# %%
# Compute popsim control variables from Census marginals
print("  compiling popsim control fields ...")
for geo, dfg in dfc.groupby("geography"):
    dic_margs[geo] = dic_margs[geo].astype(float).fillna(0)
    for ind, r in dfg.iterrows():
        dic_margs[geo][r["control_field"]] = dic_margs[geo].eval(
            r["acs_variables"].replace('"', "")
        )
    dic_margs[geo] = dic_margs[geo][list(dfg.control_field)]  # keep only control fields


# %%
ctr_geos = {}
mapping_dict = {4: "BLKGRPID", 3: "TRACTID", 2: "COUNTYID"}
for geo, dfm in dic_margs.items():
    dfm[mapping_dict[dfm.index.nlevels]] = dfm.index.map("".join)
    dfm.reset_index(drop=True, inplace=True)
    dfm.fillna(0, inplace=True)
    dfm.columns = [col.upper() for col in dfm.columns]

    # if "HHBASE" in dfm.columns:
    #    dfm = dfm.loc[dfm.HHBASE >= 0]

    f_output_control = output_control.replace(".csv", geo.lower() + ".csv")
    ctr_geos[geo] = f_output_control

    print("  saving controls to: " + output_folder + f_output_control)
    dfm.to_csv(output_folder + f_output_control)
    marginal_summary(dfm)


# %% [markdown]
# step 3. extract PUMS seed households and persons
print("\nExtrating PUMS seed households and persons from state samples")

puma_lst = df_geo_cross.PUMA.unique()

h_pums = pd.read_csv(h_pums_csv, dtype={"SERIALNO": str, "PUMA": str})
h_pums = h_pums.set_index("SERIALNO")  # keep index as string, one-liner doesn't work
p_pums = pd.read_csv(p_pums_csv, dtype={"SERIALNO": str, "PUMA": str})

# Census might change variable names by year, changed variables are in region config file
# https://www2.census.gov/programs-surveys/acs/tech_docs/pums/ACS2019_PUMS_README.pdf?
if acs_year in conf["pums_var_updates"]:
    h_pums = pums_update(h_pums, conf["pums_var_updates"][acs_year])
    p_pums = pums_update(p_pums, conf["pums_var_updates"][acs_year])

emp_df = pd.DataFrame()
h_samples, p_samples = [], []
count = 0

#%%
#between pums 2012 and 2016 5-year ACS, PUMS samples are avaiable at two different PUMA areas(00 and 10), need to combine both
#start from pums 2022 5-year ACS, PUMS samples are avaiable at two different PUMA areas(10 and 20), need to combine both

for yrs, vals in dict_acs_pums.items():
    if acs_year in yrs:
        if len(vals) == 1:
            h_pums = h_pums.loc[h_pums.PUMA.isin(puma_lst)]
            p_pums = p_pums.loc[p_pums.PUMA.isin(puma_lst)]
        else:
            pums_grp = group_pums_data(h_pums, p_pums, vals[0][1], vals[1][1])
            h_samples, p_samples = combine_puma_data(puma_lst, vals[0][1], vals[1][1], pums_grp)
            h_pums = pd.concat(h_samples)
            p_pums = pd.concat(p_samples)


if target != "housing_units":
    h_pums = h_pums.loc[
        (h_pums.TYPE == 1) & (h_pums.NP > 0)
    ]  # remove group quarters and empty units
p_pums = p_pums.loc[p_pums["SERIALNO"].isin(h_pums.index)]

h_pums, p_pums = preprocess_pums(h_pums, p_pums)

# h_pums["hh_id"] = h_pums.index.values
h_pums["hh_id"] = range(len(h_pums))
p_pums = pd.merge(
    p_pums, h_pums[["hh_id"]], left_on="SERIALNO", right_index=True, how="left"
)

print(
    f"- saving seed households: {output_folder+output_seed_hhs}.| total {str(len(h_pums))} records"
)
h_pums.to_csv(output_folder + output_seed_hhs)
print(
    f"- saving seed persons: {output_folder+output_seed_persons}.| total {str(len(p_pums))} records"
)
p_pums.to_csv(output_folder + output_seed_persons)


# %% [markdown]
# step 4. modify popsim config
print("\nupdate popsim settings")

# %%
SORT_ORDER = [
    "REGION",
    "PUMA",
    "MCD",
    "SAMPLEGEO",
    "TRACT",
    "BLKGRP",
    "TAZ",
    "BLK",
    "BUILDING",
]
sorted_geos = list(ctr_geos.keys())
sorted_geos.sort(key=lambda val: SORT_ORDER.index(val))


# %%
prj_settings = yaml.load(open(settings_file, "r"), Loader=yaml.FullLoader)

# %%
geos = sorted_geos + ["REGION", "PUMA"]
geos.sort(key=lambda val: SORT_ORDER.index(val))
prj_settings["geographies"] = geos  # sort by predefined order)
prj_settings["seed_geography"] = "PUMA"

prj_settings["data_dir"] = "data/" + str(acs_year)

# %%
for litem in prj_settings["input_table_list"]:
    if litem["tablename"] == "households":
        litem["filename"] = output_seed_hhs
    if litem["tablename"] == "persons":
        litem["filename"] = output_seed_persons
    if litem["tablename"] == "geo_cross_walk":
        litem["filename"] = output_geo_cross
    if "_control_data" in litem["tablename"]:
        geo = litem["tablename"].replace("_control_data", "")
        if geo not in ctr_geos.keys():
            prj_settings["input_table_list"].remove(litem)
        else:
            litem["filename"] = ctr_geos[geo]
            del ctr_geos[geo]
for k in ctr_geos:
    prj_settings["input_table_list"].append(
        {"tablename": k + "_control_data", "filename": ctr_geos[k]}
    )

# %%
prj_settings["control_file_name"] = "{}_{}_controls.csv".format(prj_name, str(acs_year))
prj_settings["output_tables"] = {
    "action": "include",
    "tables": ["summary_" + x for x in sorted_geos],
}


sub_bal_lst = ["sub_balancing.geography=" + x for x in sorted_geos]
prj_settings["run_list"]["steps"] = (
    [
        "input_pre_processor",
        "setup_data_structures",
        "initial_seed_balancing",
        "meta_control_factoring",
        "final_seed_balancing",
        "integerize_final_seed_weights",
    ]
    + sub_bal_lst
    + ["expand_households", "summarize", "write_tables", "write_synthetic_population"]
)

with open(
    "{}{}_{}_settings.yaml".format(output_folder, prj_name, acs_year), "w"
) as yaml_file:
    yaml.dump(prj_settings, yaml_file, default_style=None, default_flow_style=False)


# %% [markdown]
# # copy pre control file
print("copy popsim master control file")
shutil.copy(
    pre_control, output_folder + "{}_{}_controls.csv".format(prj_name, str(acs_year)),
)

# %%
print(
    "\ntotal time: {} seconds".format(round(time.time() - t0, 1)),
    "\nDone. All files are saved to " + output_folder,
    "\nTo run Populationsim:",
    '\n\t copy new settings and controls to configs folder and rename "xxx_settings.yaml" to "settings.yaml"',
    f"\n\t copy other files in {acs_year}/data folder to data/{acs_year}/",
)

