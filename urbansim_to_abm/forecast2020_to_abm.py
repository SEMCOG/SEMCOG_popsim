
# Convert 2020 SEMCOG Urbansim forecast demographic data to ABM test inputs

## ABM HHs: HHID, TAZ, *TYPE, HINCP, *ADJINC, NP, *HHT, VEH (* new variables)
## ABM Persons: HHID, PERID, AGEP, SEX, *ESR, *WKHP, *WKW, *SCHG, *MIL, *PINCP(for GQ)
## note. both ABM HH and Person tables include GQ HH and Person records

## see document for the details of task purpose and methodology

# 04-22-2020
# fix member id: no duplicates in same HH
# zone_id: remove invalid zone_id
# household id: consider shorter int number(in the future)
# NAs: Set all NAs to -9
# HHT: for GQ records is set to 0
# HINCP: Use PINCP (personal income) for GQ records instead of  (household income).
# HINCP adjustment: this dataset uses original HINCP and PINCP from PUMS, so not adjustment was implemented. The "income" variable in official 2015 model base year was adjusted.
# document: provide detailed data dictionary


# %%
import pandas as pd
import numpy as np
import math
import time
import os
#from pandas_profiling import ProfileReport
from datetime import date
import random
import concurrent.futures as conc
from itertools import repeat

random.seed(1)

# %%
def find_matching_index(df1, df2, keycols):

    rstart, rend = df1.index[0], df1.index[-1]
    print(f"process records {rstart} -- {rend}")

    v1, v2 = df1.index.name, df2.index.name
    df1 = df1.reset_index().set_index([v1] + keycols).sort_index()

    qrylst = df1.index  # assign searching sequence
    dfsample = []
    df_remain = pd.DataFrame(columns=["start", "end", "keycols", "matched", "remain"])
    while (len(qrylst) > 0) and (len(keycols) > 0):
        df2 = df2.reset_index()
        df2 = df2.set_index(keycols).sort_index()  # set index for query
        qrylst_remain = ([])  # if no matches found, remove end keycol and save the record
        for k in qrylst:
            try:
                dfsample.append([k[0], df2.loc[k[1:]].sample(1)[v2].values[0]])
            except:
                qrylst_remain.append(k[:-1])  # store with 1 less key
        qrylst = qrylst_remain
        df_remain.loc[len(df_remain)] = [rstart,rend,keycols,len(dfsample), len(qrylst)]
        keycols = keycols[:-1]  # drop one key

    if len(qrylst_remain) > 0:
        print("warning! still have unmatched records! ")
        df_remain.loc[len(df_remain)] = [rstart, rend, "", len(dfsample), list(qrylst_remain)]

    df1 = pd.merge(df1.reset_index(), pd.DataFrame(dfsample, columns=[v1, v2]),
                                    left_on=v1,right_on=v1, how="left").set_index(v1)

    return [df1, df_remain]


# %%
def concur_match(df1, df2, keycols):
    """ for every record in df1, identify and sample a matching record from df2
        df1: table to search (iterate rows)
        df2: table to sample from (query and sample)
        keycols: list of columns to be used as matching key, both df1 and df2 must have same keycols 
        df1 ands df2 must have valid index name. final table will have df2 index column joined to df1
    """

    v1, v2 = df1.index.name, df2.index.name
    assert (type(v1) == str) & (type(v2) == str) == True
    print(f"matching indices:  df1|{v1}  --  df|{v2}")

    csize = len(df1)//8 + 1 #8 is the # of cores
    while csize > 100_000:  #arbitrary size
        csize = csize//2 + 1
    print('chunk size:', csize)
    list_df = [df1[i : i + csize] for i in range(0, df1.shape[0], csize)]

    with conc.ProcessPoolExecutor() as executor:
        results = executor.map(
            find_matching_index, list_df, repeat(df2), repeat(keycols)
        )

    dfout, dfremain = pd.DataFrame(), pd.DataFrame()
    for result in results:
        dfout = pd.concat([dfout, result[0]], axis = 0)
        dfremain = pd.concat([dfremain, result[1]], axis = 0)


    return dfout, dfremain


# %%
def building_geos(st, year, df_tract_puma):
    ### associate building with other geographic units: county, tract, PUMA, etc ###
    bldgeos = pd.merge(
        st[year + "/buildings"][["parcel_id", "b_zone_id"]],
        st[year + "/parcels"][["census_bg_id", "county_id","large_area_id"]],
        left_on="parcel_id",
        right_index=True,
        how="left",
    )
    bldgeos["tract"] = bldgeos.census_bg_id // 10
    bldgeos = pd.merge(
        bldgeos.reset_index(),
        df_tract_puma[["COUNTYFP", "TRACTCE", "PUMA5CE"]],
        left_on=["county_id", "tract"],
        right_on=["COUNTYFP", "TRACTCE"],
        how="left",
    ).set_index("building_id")
    bldgeos = bldgeos.rename(columns={"b_zone_id": "zone_id", "PUMA5CE": "PUMA"})
    bldgeos = bldgeos.fillna(na_value).astype(int)
    bldgeos.columns = [x.lower() for x in bldgeos.columns]

    return bldgeos


# %%
def pums_cleaner(df):
    df.columns = [x.lower() for x in df.columns]
    df.index.names = [(None if x is None else x.lower()) for x in df.index.names]
    df.columns = [str(c) for c in df.columns]
    return df

def recode_values(df, ctrl_bins):
    ### recode/digitize df colum values to control bins """"
    for k in ctrl_bins.keys():
        df["match_" + k] = pd.np.digitize(df[k], ctrl_bins[k])
    return df

def update_PUMS_attrs(pums_hhs, pums_pps):
    # process age, race, worker, sex
    pums_pps = pums_cleaner(pums_pps)
    pums_pps.rename(
        columns={
            "agep": "age",
            "rac1p": "race_id",
            "relp": "relate",
            "sporder": "member_id",
        },
        inplace=True,
    )

    pums_pps.loc[pums_pps.race_id > 2, "race_id"] = 4
    pums_pps.loc[pums_pps.hisp > 1, "race_id"] = 3

    pums_pps["worker"] = 0
    pums_pps.loc[
        pums_pps.esr.isin(range(1, 6)), "worker"
    ] = 1  # consistent with forecast definition

    pums_pps["child"] = 0
    pums_pps.loc[(pums_pps.age < 18), "child"] = 1

    # process hhs
    pums_hhs = pums_cleaner(pums_hhs)
    pums_hhs.rename(columns={"np": "persons", "veh": "cars"}, inplace=True)

    # race_id and age
    sp = pums_pps.loc[pums_pps.relate == 0][["serialno", "race_id", "age"]].set_index("serialno")
    pums_hhs["race_id"] = sp["race_id"]
    pums_hhs["age_of_head"] = sp["age"]

    # income
    pums_hhs["income"] = pums_hhs["hincp"]
    pums_hhs["income"] *= pums_hhs["adjinc"] / 1000000.0

    # workers: update number of workers from persons table since only family workers counted in HH table
    pums_hhs["workers"] = pums_pps.groupby("serialno").worker.sum()

    # children: update number of childrens(AGEP<18) from persons table since NOC has only own children
    pums_hhs["children"] = pums_pps.groupby("serialno").child.sum()

    gq_pps = pums_pps.loc[pums_pps.relate > 15]
    pums_pps = pums_pps.loc[pums_pps.relate <= 15]  # 16 and 17 are GQ
    gq_hhs = pums_hhs.loc[(pums_hhs.type > 1) & (pums_hhs.persons > 0)]
    pums_hhs = pums_hhs.loc[(pums_hhs.type == 1) & (pums_hhs.persons > 0)]  # type 1=res, 2,3=gq

    assert len(gq_pps) == len(gq_hhs)

    return pums_hhs, pums_pps, gq_hhs, gq_pps

def fix_member_relate(dfp):

    #fix duplicated member_id
    dfp = dfp.reset_index()
    dfp = dfp.sort_values(by=['household_id', 'member_id', 'person_id']) #sorting by existing ids
    dfp['member_id'] = dfp.groupby('household_id').cumcount() + 1
    dfp = dfp.set_index('person_id')

    #fix duplicated partner records, update relate 1 to 10
    dfp_relp = dfp.loc[(dfp.relate==1) & (dfp.duplicated(['household_id', 'relate']))].index 
    dfp.loc[dfp_relp, 'relate'] = 10

    return dfp

def remove_invalid_zone(df):
    n = len(df)
    df = df.loc[df.zone_id.notnull()] #remove records with invalid zone id, in case
    print(f"removed {n - len(df)} records with invalid zone id")

    return df

if __name__ == "__main__":

    # %%
    # step 1. input setup
    # all data under ABM/test_data_012020/
    year = "2020"
    infolder = "inputs"
    outfolder = "outputs"
    hdf_file = "run4032_taz_draft_ypsi.h5"  # model base year, final version
    tract_puma_url = "https://www2.census.gov/geo/docs/maps-data/data/rel/2010_Census_Tract_to_2010_PUMA.txt"
    na_value = -9 #fill na with this value

    pums_hhs_file = "psam_h26.csv"
    pums_pps_file = "psam_p26.csv"
    new_hh_vars = ["TYPE", "HINCP", "ADJINC", "HHT"]
    new_person_vars = ["ESR", "WKHP", "WKW", "SCHG", "MIL"]

    #raw output files  #outfolder, outhhs_raw,outpersons_raw,,outgq_raw,outgq_hhs_raw
    out_hhs = "abm2020_{}_hhs.csv".format(str(date.today()))
    out_pps = "abm2020_{}_persons.csv".format(str(date.today()))

    # ctrl bins
    hhs_bins = {
        "age_of_head": [0, 16, 18, 25, 35, 45, 55, 65, 75, 85],
        "persons": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "children": [0, 1, 2, 3, 4, 5],
        "cars": [0, 1, 2, 3, 4, 5],
        "workers": [0, 1, 2, 3, 4, 5],
        "income": list(range(0, 110_000, 10_000)) + [120_000, 150_000],
        "race_id": [1, 2, 3, 4],
    }

    pps_bins = {
        "age": [0, 16, 18, 25, 35, 45, 55, 65, 75, 85],
        "race_id": [1, 2, 3, 4],
        "sex": [1, 2],
        "worker": [0, 1],
        "relate": list(range(16)),
    }

    gq_bins = {
        "relate": [16, 17],
        "age": [0, 16, 18, 25, 35, 45, 55, 65, 75, 85],
        "race_id": [1, 2, 3, 4],
    } # no gender in model data

    # %%
    stm = pd.HDFStore(os.path.join(infolder, hdf_file), "r")  # model inputs
    final_hvars = [c.lower() for c in list(stm["base/households"].columns) + new_hh_vars]
    final_pvars = [c.lower() for c in list(stm["base/persons"].columns) + new_person_vars]

    #  tract to puma cross-walk for MI
    df_tract_puma = pd.read_csv(tract_puma_url)
    df_tract_puma = df_tract_puma.loc[df_tract_puma.STATEFP == 26]

    bldgeos = building_geos(stm, year, df_tract_puma)
    puma_lst = list(bldgeos.puma.unique())



    ## step 2. --  prepare PUMS samples --
    # need new variables -> persons ['ESR', 'WKHP', 'WKW', 'SCHG', 'MIL'], hhs  ['TYPE', 'HINCP', 'ADJINC', 'HHT']
    pums_hhs = pd.read_csv(
        os.path.join(infolder, pums_hhs_file),
        usecols=["SERIALNO", "NP", "VEH", "PUMA"] + new_hh_vars,
        dtype={"SERIALNO": object},
        index_col="SERIALNO",
    )
    pums_hhs = pums_hhs.loc[pums_hhs.PUMA.isin(puma_lst)]

    pums_pps = pd.read_csv(
        os.path.join(infolder, pums_pps_file),
        usecols=[
            "SERIALNO",
            "SPORDER",
            "AGEP",
            "RELP",
            "SEX",
            "HISP",
            "RAC1P",
            "PUMA",
            "PINCP"
        ]
        + new_person_vars,
        dtype={"SERIALNO": object}
    )

    pums_pps = pums_pps.loc[pums_pps.SERIALNO.isin(pums_hhs.index.values)]
    pums_pps.index = range(len(pums_pps))
    pums_pps.index.name = "pid"

    # update and recode
    pums_hhs, pums_pps, pums_gq_hhs, pums_gq_pps = update_PUMS_attrs(pums_hhs, pums_pps)
    pums_hhs = recode_values(pums_hhs, hhs_bins)
    pums_pps = recode_values(pums_pps, pps_bins)
    pums_gq_pps = recode_values(pums_gq_pps, gq_bins)
    print(pums_pps.columns)



    ## step 3. --  prepare 2020 forecast data  --
    fcast_hhs = pd.merge(
        stm[year + "/households"],
        bldgeos,
        left_on="building_id",
        right_index=True,
        how="left",
        suffixes=['','_y']
    )
    fcast_hhs = remove_invalid_zone(fcast_hhs)
    fcast_hhs = recode_values(fcast_hhs, hhs_bins)
    fcast_hhs.index.name = "household_id"
    

    # %% add geo ids to forecast persons and recode it for matches
    fcast_pps = pd.merge(
        stm[year + "/persons"],
        fcast_hhs[bldgeos.columns],
        left_on="household_id",
        right_index=True,
        how="left",
    )
    fcast_pps = remove_invalid_zone(fcast_pps)
    fcast_pps = recode_values(fcast_pps, pps_bins)
    fcast_pps.index.name = "person_id"
    



    # step 4.  ---- match PUMS samples to forecast HHs and persons

    t0 = time.time()
    print("matching hhs starts ... \n")
    keycols = ["match_" + c for c in ["persons","race_id",
            "income", "age_of_head","children", "workers","cars"]]
    keycols.insert(3, "puma")
    fcast_hhs, hhs_remain = concur_match(fcast_hhs, pums_hhs, keycols)

    print("match hhs time", time.time() - t0)
    fcast_hhs = pd.merge(
        fcast_hhs, pums_hhs[list(set(pums_hhs.columns)-set(fcast_hhs.columns))], left_on="serialno", right_index=True, how="left"
    )
    fcast_hhs.to_csv(outfolder + "/" + out_hhs.replace('.csv','_raw.csv'))
    hhs_remain.to_csv(outfolder + "/" + out_hhs.replace('.csv', '_raw_remain.csv'))

    # %%  ---- match PUMS samples to forecast persons
    t0 = time.time()
    print("matching persons starts...\n")
    keycols = ["match_" + c for c in ["sex", "race_id", "relate", "worker", "age"]]
    fcast_pps, pps_remain = concur_match(fcast_pps, pums_pps, keycols)

    print("match persons time", time.time() - t0)
    fcast_pps = pd.merge(
        fcast_pps, pums_pps[list(set(pums_pps.columns)-set(fcast_pps.columns))], left_on="pid", right_index=True, how="left"
    )
    fcast_pps.to_csv(outfolder + "/" + out_pps.replace('.csv','_raw.csv'))
    pps_remain.to_csv(outfolder + "/" + out_pps.replace('.csv', '_raw_remain.csv'))



    ## step 5. --  Match PUMS to group quarter population --

    fcast_gq = pd.merge(stm[year + "/group_quarters"], bldgeos, left_on = 'building_id', right_index=True, how='left')
    fcast_gq = remove_invalid_zone(fcast_gq)
    fcast_gq = fcast_gq.loc[fcast_gq.zone_id.notnull()] # some SEMCOG GQ buildings have no valid zone_id/
    fcast_gq.index.name = 'person_id'

    # GQ types
    # 101-401 Institutional
    # # 601 MILITARY, NOT PRESENT IN THIS DATASET
    # # 101	cfa	Correctional Facilities for Adults
    # # 201	juv	Juvenile Facilities
    # # 301	nur	Nursing Homes
    # # 401	oif	Other Institutional
    # # 501	csh	College/Student Housing
    # # 701	onf	Other NonInstitutional
    # # 789	HL	Homeless Population

    fcast_gq["relate"] = 16
    fcast_gq.loc[fcast_gq.gq_code >= 501,'relate'] = 17  #for GQ HH type
    #fcast_gq['puma']= fcast_gq['puma'].fillna(na_value).astype(int)
    fcast_gq= recode_values(fcast_gq, gq_bins)

    print("matching GQ starts ... \n")
    t0 = time.time()
    keycols = ["match_" + c for c in ["relate", "age", "race_id"]]
    fcast_gq, gq_remain = concur_match(fcast_gq, pums_gq_pps, keycols)

    print("match GQ time", time.time() - t0)
    fcast_gq = pd.merge(
        fcast_gq, pums_gq_pps[list(set(pums_gq_pps.columns)-set(fcast_gq.columns))], left_on="pid", right_index=True, how="left"
    )
    fcast_gq['household_id'] = range(fcast_hhs.index.max() + 1, fcast_hhs.index.max() + len(fcast_gq) + 1)
    fcast_gq.to_csv(outfolder + '/' + out_pps.replace('.csv','_gq_raw.csv'))
    gq_remain.to_csv(outfolder + '/' + out_pps.replace('.csv','_gq_raw._remain.csv'))


    fcast_gq_hhs = pd.merge(fcast_gq, pums_gq_pps[['serialno'] + list(set(pums_gq_pps.columns) - set(fcast_gq.columns))], 
                                    left_on = 'serialno', right_on = 'serialno', how = 'left')
    fcast_gq_hhs.rename(columns = {'age':'age_of_head',
                                'veh': 'cars',
                                'noc': 'children',
                                'np': 'persons',
                                'worker': 'workers',
                                }, inplace = True)

    fcast_gq_hhs['hht'] = 0 #based on RSG suggestion
    fcast_gq_hhs['hincp'] = fcast_gq_hhs['pincp']
    fcast_gq_hhs.drop('pincp', axis=1, inplace=True)
    
    fcast_gq_hhs = pd.merge(fcast_gq_hhs, pums_gq_hhs[list(set(pums_gq_hhs.columns) - set(fcast_gq_hhs.columns))], 
                            left_on = 'serialno', right_index = True, how = 'left')
    fcast_gq_hhs = fcast_gq_hhs.set_index('household_id')
    fcast_gq_hhs.to_csv(outfolder + '/' + out_hhs.replace('.csv','_gq_hhs_raw.csv'))


    ## step 6. combine resi and gq outputs
    fcast_hhs = fcast_hhs[[c.lower() for c in list(stm["base/households"].columns) + new_hh_vars]]
    fcast_pps = fcast_pps[[c.lower() for c in list(stm["base/persons"].columns) + new_person_vars]]
    fcast_gq_hhs = fcast_gq_hhs[[c.lower() for c in list(stm["base/households"].columns) + new_hh_vars]]
    fcast_gq = fcast_gq[[c.lower() for c in list(stm["base/persons"].columns) + new_person_vars]]
    final_hhs = pd.concat([fcast_hhs, fcast_gq_hhs], axis = 0)
    final_pps = pd.concat([fcast_pps, fcast_gq], axis = 0)

    # additional fixes
    na_value = -9 #RSG suggestion to fill na with -9
    final_hhs.fillna(na_value, inplace =True)
    final_pps.fillna(na_value, inplace =True)
    final_pps =  fix_member_relate(final_pps)

    final_hhs.to_csv(os.path.join(outfolder, out_hhs))
    final_pps.to_csv(os.path.join(outfolder, out_pps))
