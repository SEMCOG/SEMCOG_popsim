# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
from IPython import get_ipython

# %%
# generate control totals needed by populationsim(RSG)


# %%
import os
import pandas as pd
import re
from census import Census

######### Census Downloader  ######
# %%
class Census_Downloader:
    def __init__(
        self, census_reader, states, counties=None, tract_ids=None, blockgroup_ids=None
    ):
        self.states = states
        self.counties = counties
        self.tracts = tract_ids
        self.blockgroups = blockgroup_ids
        self.cread = census_reader
        self.udpate_states_counties([self.states, self.counties])

    def state_download(self, vars):
        return self.cread.get(vars, geo={"for": "state:{}".format(self.states)})

    def county_download(self, vars):
        print(self.states, self.counties)
        return self.cread.get(
            vars,
            geo={
                "for": "county:{}".format(self.counties),
                "in": "state:{}".format(self.states),
            },
        )

    def tract_download(self, vars):
        return self.cread.get(
            vars,
            geo={
                "for": "tract:{}".format(self.tracts),
                "in": "state:{} county:{}".format(self.states, self.counties),
            },
        )

    def blockgroup_download(self, vars):
        clst = self.counties.split(",")
        cm = []
        for cn in clst:
            cm += self.cread.get(
                vars,
                geo={
                    "for": "block group:{}".format(self.blockgroups),
                    "in": "state:{} county:{} tract:{}".format(
                        self.states, cn, self.tracts
                    ),
                },
            )
        return cm

    def fips_lookup(self, states, counties=None):
        if counties == "*":
            counties = None
        fips_table = pd.read_csv(
            "https://www2.census.gov/geo/docs/reference/codes/files/national_county.txt",
            header=None,
            names=["state", "state_fips", "county_fips", "county", "type"],
            dtype=str,
        )
        qstr = "(state in {})".format(states)
        if counties:
            qstr += " & (county in {})".format(counties)
            dfq = fips_table.query(qstr)
            return list(dfq.state_fips.unique()), list(dfq.county_fips.unique())
        dfq = fips_table.query(qstr)
        print(list(dfq.state_fips.unique()))
        return list(dfq.state_fips.unique()), None

    def udpate_states_counties(self, geos):
        for i in [0, 1]:
            if geos[i] != None:
                if (type(geos[i]) != list) & (geos[i] != "*"):
                    geos[i] = [geos[i]]
                geos[i] = [str(x) for x in geos[i]]
        if (type(geos[0]) == list) & (geos[0][0].isdigit() == False):
            geos[0], geos[1] = self.fips_lookup(geos[0], geos[1])
        for i in [0, 1]:
            if geos[i]:
                geos[i] = ",".join(
                    [str(x).zfill(i + 2) for x in geos[i]]
                )  # i+2 cause state and county need 2 and 3 0s in lead
        self.states = geos[0]
        self.counties = geos[1]

    def download(self, variables):
        dfm = pd.DataFrame()
        if not (self.counties):
            downv = self.state_download(variables)
        elif not (self.tracts):
            downv = self.county_download(variables)
        elif not (self.blockgroups):
            downv = self.tract_download(variables)
        else:
            downv = self.blockgroup_download(variables)
        dfm = pd.DataFrame.from_dict(downv)

        return dfm

def read_tract_puma_crosswalk(tup, dict_cross):
    """
    Read a tract to PUMA crosswalk file and return a dataframe with the appropriate columns renamed and set as the index.

    Inputs:
    tup : A tuple containing the column names for the tract and puma variables.
    dict_cross : A dictionary mapping tuples of tract and puma column names to crosswalk file paths.

    Returns
    A dataframe with the tractid column set as the index and the puma column renamed to "PUMA".
    """

    df = pd.read_csv("./geo/" + dict_cross[tup], dtype=str)
    df["COUNTYID"] = df["STATEFP"] + df["COUNTYFP"]
    df["TRACTID"] = df["STATEFP"] + df["COUNTYFP"] + df[tup[0]]
    df = df.set_index("TRACTID")
    df.rename(columns={tup[1]: "PUMA"}, inplace=True)

    return df[["COUNTYID","PUMA"]]

def preprocess_pums(h_pums, p_pums):
    # h_pums and p_pums must have the same index column 'SERIALNO'
    p_pums = p_pums.set_index("SERIALNO")

    # add AGEHOH to PUMS sample
    h_pums["AGEHOH"] = p_pums.loc[p_pums.RELP == 0].AGEP

    # add HRACE to PUMS sample
    # RAC1P Character 1
    # Recoded detailed race code
    # 1 .White alone
    # 2 .Black or African American alone
    # 3 .American Indian alone
    # 4 .Alaska Native alone
    # 5 .American Indian and Alaska Native tribes specified; or
    # .American Indian or Alaska Native, not specified and no other
    # .races
    # 6 .Asian alone
    # 7 .Native Hawaiian and Other Pacific Islander alone
    # 8 .Some Other Race alone
    # 9 .Two or More Races
    rac_map = {1: 1, 2: 2, 6: 3}
    h_pums["HRACE"] = p_pums.loc[p_pums.RELP == 0].RAC1P
    h_pums["HRACE"] = h_pums["HRACE"].map(rac_map).fillna(4)

    # add HHISP to PUMS sample
    # HISP Character 2
    # Recoded detailed Hispanic origin
    # 01 .Not Spanish/Hispanic/Latino
    # 02 - 24, hisp
    hisp_map = {1: 0}
    h_pums["HHISP"] = p_pums.loc[p_pums.RELP == 0].HISP
    h_pums["HHISP"] = h_pums["HHISP"].map(hisp_map).fillna(1)

 
    # adjust hh income to current ACS year (release year)
    h_pums['income'] = (h_pums.HINCP * h_pums.ADJINC / 1000000).astype(int)

    # add HWORKERS to PUMS sample
    # ESR Character 1
    # Employment status recode
    # b .N/A (less than 16 years old)
    # 1 .Civilian employed, at work
    # 2 .Civilian employed, with a job but not at work
    # 3 .Unemployed
    # 4 .Armed forces, at work
    # 5 .Armed forces, with a job but not at work
    # 6 .Not in labor force
    p_pums = p_pums.reset_index()
    h_pums["HWORKERS"] = (
        p_pums.loc[p_pums.ESR.isin([1, 2, 4, 5])].groupby("SERIALNO").ESR.size()
    )

    h_pums.fillna(0, inplace=True)
    # for v in ["AGEHOH","HRACE", "HHISP","HWORKERS" ]:
    #    print (v, h_pums.groupby(v).size())

    # recode NAICSP to industry
    # https://www2.census.gov/programs-surveys/acs/tech_docs/code_lists/2018_ACS_Code_Lists.pdf
    # summary table example https://www.socialexplorer.com/data/ACS2017_5yr/metadata/?ds=ACS17_5yr&table=B08126
    dict_naics2ind = {
        "11": 1,
        "21": 1,
        "23": 2,
        "31": 3,
        "32": 3,
        "33": 3,
        "3M": 3,
        "42": 4,
        "44": 5,
        "45": 5,
        "4M": 5,
        "48": 6,
        "49": 6,
        "22": 6,
        "51": 7,
        "52": 8,
        "53": 8,
        "54": 9,
        "55": 9,
        "56": 9,
        "61": 10,
        "62": 10,
        "71": 11,
        "72": 11,
        "81": 12,
        "92": 13,
        "99": 0,
    }
    p_pums["industry"] = p_pums.NAICSP.str[:2]
    p_pums.industry.replace(dict_naics2ind, inplace=True)
    p_pums.loc[p_pums.NAICSP.str[:6] == "928110", "industry"] = 14
    p_pums.loc[p_pums.NAICSP.isnull(), "industry"] = 0
    p_pums.industry = p_pums.industry.astype(int)
 
    # adjust person income to current ACS year (release year)
    p_pums['pincome'] = p_pums['PINCP'] * p_pums['ADJINC'] / 1000000
    p_pums['pincome'] = p_pums.loc[~p_pums['pincome'].isnull(), 'pincome'].astype(int)

    return h_pums, p_pums


def pums_update(df, dic_var):
    for k in dic_var.keys():
        if k in df.columns:
            conv = dic_var[k]
            df[conv['std_variable']] = df[k]
            if 'std_codes' in conv:
                df[conv['std_variable']].replace(conv['std_codes'], inplace=True)

    return df

def group_pums_data(h_pums, p_pums, pma0_col, pma1_col):
    """
    Group households and persons data by PUMA10 and PUMA00 columns.

    Args:
        h_pums (DataFrame): Households data.
        p_pums (DataFrame): Persons data.
        pma0_col (str): Name of first PUMA column.
        pma1_col (str): Name of second PUMA column.

    Returns:
        dict: Nested dictionary with grouped data.
    """
    pums_grp = {}
    emp_df = pd.DataFrame()  # Create an empty DataFrame
    for pma in [pma0_col, pma1_col]:
        pums_grp[pma] = {"households": {}, "persons": {}}
        for indx, grp in h_pums.loc[h_pums[pma] != -9].groupby(pma):
            pums_grp[pma]["households"][indx] = grp
        for indx, grp in p_pums.loc[p_pums[pma] != -9].groupby(pma):
            pums_grp[pma]["persons"][indx] = grp
        pums_grp[pma]["households"][0] = emp_df
        pums_grp[pma]["persons"][0] = emp_df
    return pums_grp

def combine_puma_data(puma_lst, pma0_col, pma1_col, pums_grp):
    """
    Combine household and person data from two different PUMA areas.

    Parameters:
    puma_lst (list): List of PUMA values to combine.
    pma0_col (str): Key for first PUMA  in pums_grp.
    pma1_col (str): Key for second PUMA in pums_grp.
    pums_grp (dict): Dictionary containing PUMA data.

    Returns:
    h_samples (list): List of combined household DataFrames.
    p_samples (list): List of combined person DataFrames.
    """
    h_samples = []
    p_samples = []
    count = 0

    for puma in puma_lst:
        count += 1
        h_puma = pd.concat(
            [
                pums_grp[pma0_col]["households"][int(puma[:5])],
                pums_grp[pma1_col]["households"][int(puma[5:])],
            ]
        )
        h_puma["PUMA"] = puma
        h_puma.index = str(count) + h_puma.index
        h_samples.append(h_puma)

        p_puma = pd.concat(
            [
                pums_grp[pma0_col]["persons"][int(puma[:5])],
                pums_grp[pma1_col]["persons"][int(puma[5:])],
            ]
        )
        p_puma["PUMA"] = puma
        p_puma["SERIALNO"] = str(count) + p_puma["SERIALNO"]
        p_samples.append(p_puma)

    return h_samples, p_samples
# def control_summary(df_controls, df_margins):
#     for _, dft in df_controls[['geography', 'attr', 'control_field']].groupby('geography'):
#         dft = dft[['attr', 'control_field']].groupby('attr').agg(list)
#         att_dict = dft.to_dict()
#     for k in att_dict.keys():
#         print(k, df_margins[att_dict[k]].sum().sum())


def marginal_summary(df_margin):
    print('\n * * * verify maringal sums:')
    at_lst = []
    for col in df_margin.columns:
        c = ''.join(i for i in col if not i.isdigit())
        if 'ID' not in c:
            at_lst.append(c)
    for at in sorted(set(at_lst)):
        cols = [c for c in df_margin.columns if c.startswith(at)]
        print(at, df_margin[cols].sum().sum(), cols)
