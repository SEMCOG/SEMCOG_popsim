# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'


# %%
## This program prepares all inputs files needed by PopulationSim(RSG), including configuration, geo_crosswalk, Census marginals and PUMS sample households and persons
## Syntax:  > python popsim_input_maker.py key yaml (key: Census API key; yaml: a yaml configuration file such as region.yaml)
## Inputs (defined in yaml config):
##      preprocess/semcog_res_setting.yaml   (a generic setting file for populationsim)
##      preprocess/control_pre.csv   (control file same as populationsim, with additional "acs_variables" column)
##      preprocess/ACS/xxxxhxx.csv   (PUMS household samples)
##      preprocess/ACS/xxxxpxx.csv   (PUMS person samples)
##      preprocess/ACS/tract0010_to_puma.csv   (tract00 and tract10 to PUMA crosswalk file)


# %%
import os, re, time
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


# %%
acs_year = conf["year"]
print("\nsynthersizing population for year ", acs_year)

c = Census(args.key, year=acs_year).acs5  # define dataset
prj_name = conf["region"]["name"]
data_folder = conf["data_folder"]
state = conf["region"]["state"][0]
counties = conf["region"]["counties"]
pre_folder = conf["preprocess"]["folder"]
settings_file = "../" + pre_folder + conf["preprocess"]["settings_file"]
h_pums_csv = (
    "../" + pre_folder + conf["preprocess"]["h_pums_csv"].format(str(acs_year - 2000))
)
p_pums_csv = (
    "../" + pre_folder + conf["preprocess"]["p_pums_csv"].format(str(acs_year - 2000))
)

output_folder = "../{}/{}/".format(data_folder, str(acs_year))
output_geo_cross = "{}_{}_geo_cross_walk.csv".format(prj_name, str(acs_year))
output_control = "{}_{}_control_totals_.csv".format(prj_name, str(acs_year))
output_seed_hhs = "{}_{}_seed_households.csv".format(prj_name, str(acs_year))
output_seed_persons = "{}_{}_seed_persons.csv".format(prj_name, str(acs_year))


# %% [markdown]
# step 1. make geographic cross work file

# %%
print(
    "\nPreparing Census geographies: \n\tstate: {}  \n\tcounty: {}".format(
        state, counties
    )
)
acgeo = Census_Downloader(
    c, state, counties, "*", "*"
)  # download Census BGs for this region
df_geo = pd.DataFrame.from_dict(acgeo.download("NAME")).drop("NAME", axis=1)
df_geo["tractid"] = df_geo["state"] + df_geo["county"] + df_geo["tract"]
df_geo["blkgrpid"] = df_geo["tractid"] + df_geo["block group"]
df_geo.columns = [col.upper() for col in df_geo.columns]

# %%
# depends on synthesizing year, switch to different tract-PUMA files
if acs_year >= 2017:
    df_tract_puma = pd.read_csv(conf["geographies"]["tract_puma_file"], dtype=str)
    df_tract_puma.rename(columns={"PUMA5CE": "PUMA"}, inplace=True)
elif (acs_year >= 2010) & (acs_year <= 2016):
    df_tract_puma = pd.read_csv(
        "../" + pre_folder + conf["geographies"]["tract_puma0010_file"], dtype=str
    )
    df_tract_puma.columns = [col.upper() for col in df_tract_puma.columns]
    if acs_year >= 2012:
        df_tract_puma.fillna("00000", inplace=True)
        df_tract_puma["PUMA"] = df_tract_puma["PUMA10_ID"] + df_tract_puma["PUMA00_ID"]
    else:
        df_tract_puma.rename(columns={"PUMA00_ID": "PUMA"}, inplace=True)
else:
    print("synthesis year should be 2010 or later")
    exit()

df_tract_puma = df_tract_puma.loc[df_tract_puma.STATEFP == acgeo.states]
df_tract_puma["TRACTID"] = (
    df_tract_puma["STATEFP"] + df_tract_puma["COUNTYFP"] + df_tract_puma["TRACTCE"]
)  # remove  for now to reduce ID length as very long IDs cause calculation errors in populationsim


# %%
# join tract-PUMA and creat geo cross file
df_geo_cross = pd.merge(df_geo, df_tract_puma, on="TRACTID", how="left")
df_geo_cross["REGION"] = 1
df_geo_cross = df_geo_cross[["TRACTID", "BLKGRPID", "PUMA", "REGION"]]
print("  saving geo cross walk file: " + output_folder + output_geo_cross)
df_geo_cross.to_csv(output_folder + output_geo_cross)


# %% [markdown]
# step 2. download and compile maginal control files

# %%
# download Census marginal controls
# Census controls variables are defined in "controls_pre" table(popsim "controls" table + "acs_variables" field )
# "acs_variables" contains Census variables and expressions
#
print("\nMaking popsim controls")
print("  downloading Census variables ...")
dfc = pd.read_csv("../" + pre_folder + conf["preprocess"]["pre_control"])
dic_margs = {}
for geo, dfgeo in dfc.groupby("geography"):
    full_vars = list(
        set(
            re.findall(r"B[0-9]{5}[A-Z]{0,1}_[0-9]{3}E", str(list(dfgeo.acs_variables)))
        )
    )
    if geo == "BLKGRP":
        ac5 = Census_Downloader(c, state, counties, "*", "*")
        geo_cols = ["state", "county", "tract", "block group"]
    elif geo == "TRACT":
        ac5 = Census_Downloader(c, state, counties, "*")
        geo_cols = ["state", "county", "tract"]
    print("\t" + geo + " marginals ")
    dic_margs[geo] = ac5.download(full_vars).set_index(geo_cols)


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
for geo, dfm in dic_margs.items():

    # dfcross = pd.read_csv(geo_cross_csv, dtype = str)
    if dfm.index.nlevels == 3:
        dfm["TRACTID"] = [l1 + l2 + l3 for l1, l2, l3 in dfm.index]  # remove l1 for now
    elif dfm.index.nlevels == 4:
        dfm["BLKGRPID"] = [
            l1 + l2 + l3 + l4 for l1, l2, l3, l4 in dfm.index
        ]  # remove l1 for now
    dfm.reset_index(drop=True, inplace=True)

    dfm.fillna(0, inplace=True)
    dfm.columns = [col.upper() for col in dfm.columns]
    if "HHBASE" in dfm.columns:
        dfm = dfm.loc[dfm.HHBASE > 0]

    f_output_control = output_control.replace(".csv", geo.lower() + ".csv")
    ctr_geos[geo] = f_output_control
    print("  saving control file: " + output_folder + f_output_control)
    dfm.to_csv(output_folder + f_output_control)

# %% [markdown]
# step 3. extract PUMS seed households and persons

print("\nExtrating PUMS seed households and persons")

# %%
# puma_lst = df_tract_puma.loc[
#     (df_tract_puma.STATEFP == acgeo.states)
#     & (df_tract_puma.COUNTYFP.isin(acgeo.counties.split(",")))
# ].PUMA.unique()
puma_lst = df_geo_cross.PUMA.unique()

# %%
h_pums = pd.read_csv(h_pums_csv, index_col="SERIALNO", dtype={"SERIALNO": str})
p_pums = pd.read_csv(p_pums_csv, dtype={"SERIALNO": str})
emp_df = pd.DataFrame()
h_samples, p_samples = [], []
count = 0

if (acs_year <= 2011) or (acs_year >= 2017):
    h_pums = h_pums.loc[h_pums.PUMA.isin(puma_lst)]
    p_pums = p_pums.loc[p_pums.PUMA.isin(puma_lst)]
else:
    pums_grp = {}
    for pma in ["PUMA10", "PUMA00"]:
        pums_grp[pma] = defaultdict(dict)
        for indx, grp in h_pums.loc[h_pums[pma] != -9].groupby(pma):
            pums_grp[pma]["households"][indx] = grp
        for indx, grp in p_pums.loc[p_pums[pma] != -9].groupby(pma):
            pums_grp[pma]["persons"][indx] = grp
        pums_grp[pma]["households"][0] = emp_df
        pums_grp[pma]["persons"][0] = emp_df

    for puma in puma_lst:
        count += 1
        h_puma = pd.concat(
            [
                pums_grp["PUMA10"]["households"][int(puma[:5])],
                pums_grp["PUMA00"]["households"][int(puma[5:])],
            ]
        )
        h_puma["PUMA"] = puma
        h_puma.index = count * (10 ** 13) + h_puma.index
        h_samples.append(h_puma)
        p_puma = pd.concat(
            [
                pums_grp["PUMA10"]["persons"][int(puma[:5])],
                pums_grp["PUMA00"]["persons"][int(puma[5:])],
            ]
        )
        p_puma["PUMA"] = puma
        p_puma["SERIALNO"] = count * (10 ** 13) + p_puma["SERIALNO"]
        p_samples.append(p_puma)

    h_pums = pd.concat(h_samples)
    p_pums = pd.concat(p_samples)

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
    "  saving seed households: {}   . {} records".format(
        output_folder + output_seed_hhs, str(len(h_pums))
    )
)
h_pums.to_csv(output_folder + output_seed_hhs)
print(
    "  saving seed persons: {}   . {} records".format(
        output_folder + output_seed_persons, str(len(p_pums))
    )
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
    "BUILDING",
]
sorted_geos = list(ctr_geos.keys())
sorted_geos.sort(key=lambda val: SORT_ORDER.index(val))


# %%
prj_config = yaml.load(open(settings_file, "r"), Loader=yaml.FullLoader)

# %%
geos = sorted_geos + ["REGION", "PUMA"]
geos.sort(key=lambda val: SORT_ORDER.index(val))
prj_config["geographies"] = geos  # sort by predefined order)
prj_config["seed_geography"] = "PUMA"

prj_config["data_dir"] = "data/" + str(acs_year)

# %%
for litem in prj_config["input_table_list"]:
    if litem["tablename"] == "households":
        litem["filename"] = output_seed_hhs
    if litem["tablename"] == "persons":
        litem["filename"] = output_seed_persons
    if litem["tablename"] == "geo_cross_walk":
        litem["filename"] = output_geo_cross
    if "_control_data" in litem["tablename"]:
        geo = litem["tablename"].replace("_control_data", "")
        if geo not in ctr_geos.keys():
            prj_config["input_table_list"].remove(litem)
        else:
            litem["filename"] = ctr_geos[geo]
            del ctr_geos[geo]
for k in ctr_geos:
    prj_config["input_table_list"].append(
        {"tablename": k + "_control_data", "filename": ctr_geos[k]}
    )

# %%
prj_config["control_file_name"] = "{}_{}_controls.csv".format(prj_name, str(acs_year))
prj_config["output_tables"] = {
    "action": "include",
    "tables": ["summary_" + x for x in sorted_geos],
}


sub_bal_lst = ["sub_balancing.geography=" + x for x in sorted_geos]
prj_config["run_list"]["steps"] = (
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
    yaml.dump(prj_config, yaml_file, default_style=None, default_flow_style=False)


# %% [markdown]
# # copy pre control file

print("copy popsim master control file")
shutil.copy(
    "../" + pre_folder + conf["preprocess"]["pre_control"],
    output_folder + "{}_{}_controls.csv".format(prj_name, str(acs_year)),
)


# %% [markdown]

print(
    "\ntotal time: {} seconds".format(round(time.time() - t0, 1)),
    "\nDone. All files are saved to " + output_folder,
    '\nTo run Populationsim, copy new settings and controls to configs folder and rename "xxx_settings.yaml" to "settings.yaml"',
)

