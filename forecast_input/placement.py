# placement.py
# Household placement process
import time
import os
import subprocess as sbp
import math

import numpy as np
import pandas as pd
import numba
import yaml
from numba.typed import List

from forecast_input.transform import calculate_housing_units_by_bg

with open('configs/placement.yaml', 'r') as f:
    placement_config = yaml.load(f, Loader=yaml.CLoader)
with open('configs/mcd.yaml', 'r') as f:
    mcd_mapping = yaml.load(f, Loader=yaml.CLoader)

# get around the numba type checking
recode_building_type = List([82, 83, 84])
recode_hh_bld = List([2, 4, 5])
recode_nu_hh = List([3, 4, 5, 6, 7, 8, 9])
recode_nu_building = List([1, 2, 3, 5, 10, 20, 50])
recode_biv = List([ 15051, 22263, 27672, 32979, 38158, 44372, 52168,
    61582, 71645, 82921, 94136, 108549, 135510, 179252])
recode_hpv = List([ 30000, 47000, 60000, 75000, 90000, 100000, 120000,
    135000, 150000, 170000, 189000, 210000, 260000, 330000])
recode_rent = List([ 396, 589, 650, 690, 730, 772, 823, 890, 955, 1021,
    1216, 1387, 1605, 1823])
recode_income = List([ 11816, 19456, 26249, 33883, 40961, 48853, 57344,
    65694, 75276, 86724, 100842, 116168, 138270, 177482])
recode_yb = List([1940, 1950, 1960, 1970, 1980, 1990, 2000, 2005, 2006,
    2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017,
    2018, 2019, 2020, 2021, 2022])
recode_p_over18 = List([ 2, 3, 4, 5, 6, 7])

def preprocess_households(households, mcd_by_bg, hu_by_mcd_bg):
    # households["tenure"] = 1
    # households.loc[households['hh_tenure_agehead'].str.startswith("r"), "tenure"] = 0
    # process tract and block group
    households['county'] = households['tract'] // 1000000 % 1000
    households['tract'] = households['tract'].abs() % 1000000 # last 6 digits
    households['block group'] = households['block group'] % 10 # last 1 digit

    # adjusting housing value by ADJHSG
    households["valp"] = households["valp"] * households['ADJHSG'] / 1000000
    # replace 0 with nan
    households["valp"] = households["valp"].replace(0, np.nan)
    # get mcd
    def bg_to_mcd(row):
        # helper func
        geo = (row.county, row.tract, row['block group'])
        # if the geo in the hashmap
        if geo in mcd_by_bg:
            # calculate weight by divide HU by the BG HU total
            options = mcd_by_bg[geo]
            hu = hu_by_mcd_bg[geo]
            w = [x/sum(hu) for x in hu]
            # make the choice based on weight
            return np.random.choice(options, 1, p=w)[0]
        else:
        # the geo not found in the hashmap, prob means that the bg not in the building table
        # return 9999 as dummy mcd code
            return 9999
    t1 = time.time()
    households['mcd'] = households[['county', 'tract', 'block group']].apply(bg_to_mcd, axis=1)
    t2 = time.time()
    print('HH MCD generating run finshed, total run time %s' % str(t2-t1))

    # tenure 1-2 owner(1), 3-4 renter(0)
    households.loc[households['tenure'].lt(3), 'tenure'] = 1
    households.loc[households['tenure'].ge(3), 'tenure'] = 0
    households = households.reset_index()[['household_id', 'county', 'mcd', 'tract', 'block group', 'tenure', 'bld',
                                  'valp', 'rent', 'ybl', 'income', 'p_over18']]
    return households.fillna(-1).astype(int)

@numba.jit(nopython=True)
def recode_households_np(households, recode_hh_bld, recode_nu_hh, recode_hpv, recode_rent, recode_yb, recode_income, recode_p_over18):
    recoded = np.zeros((households.shape[0], 12))
    # hh_id, county, mcd, tract, bg, tenure, buildingtype: keep
    recoded[:, :6] = households[:, :6]
    # recode # of unit
    recoded[:, 6] = np.digitize(
        households[:, 6],
        recode_hh_bld 
    )
    # recode nu
    recoded[:, 7] = np.digitize(
        households[:, 6],
        recode_nu_hh 
    )
    # recode hpv
    hpv = np.digitize(
        households[:, 7],
        recode_hpv
    )
    # recode rent
    rent = np.digitize(
        households[:, 8],
        recode_rent 
    )

    # calculated housing cost for both owner and renter
    recoded[:, 8] = hpv * recoded[:, 5] + rent * (1 - recoded[:, 5])

    # recode yb
    recoded[:, 9] = np.digitize(
        households[:, 9],
        recode_yb 
    )
    # recode income
    recoded[:, 10] = np.digitize(
        households[:, 10],
        recode_income 
    )
    # recode number of people over 18
    recoded[:, 11] = np.digitize(
        households[:, 11],
        recode_p_over18
    )
    return recoded

def preprocess_buildings(buildings):
    # Read selected fields
    data = buildings[
        ['parcel_id', 'non_residential_sqft', 'year_built', 
        'residential_units', 'owner_units', 'building_type_id',
         'improvement_value', 'county', 'tract', 'bg', 'mcd']]
    # Remove records with duplicated records
    data = data[~data.index.duplicated()].copy()
    # Calculate Residential sqft
    data['res_sqft'] = data.residential_units * buildings[~buildings.index.duplicated()].sqft_per_unit

    # Calculate tax per sqft
    data['tax_per_sqft'] = 0.5 * data.improvement_value / (1.0 + data.res_sqft + data.residential_units)
    # Set tax per sqft to NaN if the values < 0.2
    data.loc[~(data.tax_per_sqft > 0.2), "tax_per_sqft"] = np.nan

    # Filter out buildings without residential units
    data = data[data['residential_units'] > 0]
    # Calculate Tax per unit
    data["tax_per_units"] = data["tax_per_sqft"] * data["res_sqft"] / data["residential_units"]

    # TODO: Ask Jeff about this check
        # there are 2 building records with owner_units > residential units
    # # Hope we don't need in the future
    # print "buildings with more owner_units than residential_units", (data.owner_units > data.residential_units).sum()
    data["owner_units"] = data[["owner_units", "residential_units"]].min(axis=1)

    # export hu table
    data.reset_index()[[
       'building_id', 'county', 'tract', 'bg', 
       'year_built', 'building_type_id', 'tax_per_units', 
       'res_sqft', 'residential_units', 'owner_units'
    ]].to_csv('temp_for_placement/rust_test/HU.csv', index=False)

    data = data.reset_index()[[
        'building_id', 'county', 'mcd', 'tract', 'bg',
        'building_type_id','tax_per_units',
        'year_built', 'residential_units', 'owner_units']]
    return data.fillna(-1).astype(int)

# @numba.jit(nopython=True)
def recode_buildings_table_np(buildings, recode_building_type, recode_nu_building, recode_biv, recode_yb):
    # for each building row 
    # number of building copy = #of unit - owner_unit
    # [
        # 'building_id', 'county', 'mcd', 'tract', 
        # ''buildingtype', blockgroup', 'tenure',  
        # 'recode_nu', 'recode_biv', 'recode_yb', 
        # ]
    recoded = np.zeros((buildings.shape[0], 11))
    # hh_id, county, mcd, tract, bg, tenure, buildingtype: keep
    recoded[:, :5] = buildings[:, :5]
    # recode tenure default 0:renter
    recoded[:, 5] = np.zeros(buildings.shape[0])
    # recode building type
    recoded[:, 6] = np.digitize( 
        buildings[:, 5],
        recode_building_type
    )
    # recode nu
    recoded[:, 7] = np.digitize(
        buildings[:, 5],
        recode_nu_building 
    )
    # recode biv
    recoded[:, 8] = np.digitize(
        buildings[:, 6],
        recode_biv 
    )
    # recode yb
    recoded[:, 9] = np.digitize(
        buildings[:, 7],
        recode_yb 
    )
    # biv as last col verse hh_income
    recoded[:, 10] = recoded[:, 8]
    # residential unit(ind8) and owner units(ind9)
    # recoded[:, 10:12] = buildings[:, 8:10]
    # add renter buildings
    renters_arr = buildings[:, 8] - buildings[:, 9]
    renters = np.repeat(recoded, renters_arr, axis=0)
    owners_arr = buildings[:, 9]
    owners = np.repeat(recoded, owners_arr, axis=0)
    # set tenure to 1 as owner
    owners[:, 5] = 1
    return np.vstack((renters, owners)).astype(int)

def generate_n_over18_for_hu(hu, voter):
    # filter out voters w/o building_id
    voter = voter[voter['BUILDING_2020'] != 0]
    # filter out invalid building_id
    voter = voter[voter['BUILDING_2020'] < 100000000]
    voter_by_bid_and_hh_id = voter[
        ['PERSON_ID', 'BUILDING_2020', 'HOUSEHOLD_ID']
    ].groupby(['BUILDING_2020', 'HOUSEHOLD_ID']).count().reset_index()
    voter_by_bid_and_hh_id.columns = ['building_id', 'household_id', 'p_over18']
    # Array[Array[building_id, hh_code, persons_count]]
    hu_n = len(hu)
    # initialize p_over18 and default as 0
    out_hu = np.hstack((hu, np.zeros(hu_n).reshape(-1, 1))).astype(int)
    # create df from array for merging 
    hu_df = pd.DataFrame(
        out_hu, 
        columns=[
            'building_id', 'county', 'mcd', 'tract', 
            'blockgroup', 'tenure', 'buildingtype', 
            'recode_nu', 'recode_biv', 'recode_yb', 'recode_bc', 
            'recode_p_over18'
            ],
        dtype=int
    )
    # calculate hu count per building_id
    hu_by_building = hu_df[['county', 'building_id']].groupby('building_id').count().county
    hu_by_building.name = 'hu_count'
    # join back to the table
    hu_df = hu_df.merge(hu_by_building, how='left', on='building_id')
    # populate index for later usage
    hu_df = hu_df.reset_index()

    # get HU whose building has only 1 HU in it
    single_hu = hu_df[hu_df['hu_count'] == 1]
    # merge with the voters' table
    single_hu_merged = single_hu.merge(voter_by_bid_and_hh_id, how='left', on='building_id')
    # ** Important ** May worth discussion
    # some building which has only 1 HU in the building table, contains more than 1 voters family
    # r = single_hu_merged.groupby('index').count()['p_over18']
    # single_hu_merged.loc[single_hu_merged['index'].isin(r.sort_values().tail(10).index)]['building_id'].unique()
    # pick the first one in this case
    first_pick_voter_hh = single_hu_merged.groupby('index').first()['p_over18']
    # cleaning
    single_hu.loc[:, 'recode_p_over18'] = first_pick_voter_hh.fillna(0).astype(int)
    # pass the result back to the main table
    hu_df.loc[hu_df['hu_count'] == 1, :] = single_hu

    # get HU whose building has more than 1 HU in it
    multi_hu = hu_df[hu_df['hu_count'] > 1]
    # convert to numpy array
    multi_hu_np = multi_hu.iloc[:, 1:-1].to_numpy()
    # set up iteration variable for preventing reusing voters family
    visited = set()
    for hu_ind in range(len(multi_hu_np)):
        unit = multi_hu_np[hu_ind]
        unit_bid = unit[0]
        # for a unit, find matching voters
        # if building_id cannot be found, skip
        options = voter_by_bid_and_hh_id.loc[voter_by_bid_and_hh_id['building_id'] == unit_bid]
        if len(options) == 0:
            continue
        # increase ind until found the unused voter family
        ind = 0
        while ind < len(options) and options.index[ind] in visited:
            ind += 1
        # continue if exceed the length
        if ind > len(options) - 1:
            continue
        option = options.iloc[ind]
        # assign p_over18 with the person count
        unit[11] = option['p_over18']
        # drop the chosen row
        # print(options.index[ind], hu_ind)
        visited.add(options.index[ind])
    # pass the result back the main table
    hu_df.loc[hu_df['hu_count'] > 1, 'recode_p_over18'] = multi_hu_np[:, -1]
    ### TODO: need recode? default value?
    # export as numpy array
    # remove the first(index) and the last column(hu_count)
    return hu_df.iloc[:, 1:-1].to_numpy()

# TODO: convert to numba func
# @numba.jit(nopython=True)
def row_rank( match_matrix, row):
    '''
    Row ranking

    args:
    match_matrix:   np.ndarray matrix to ranking against
    row:            np.array row

    return:
    prefer_ranking  np.array
    '''
    weight = np.array([4000, 3000, 2000, 1000, 1000, 500, 100, 100, 100, 100, 50], dtype=np.int64)
    # bool match geo
    scores_mat = np.zeros((match_matrix.shape[0], 8), dtype=np.int64)
    scores_mat[:, 0]    = ((match_matrix[:, 1:5] != row[1:5])       * weight[0:4]).sum(axis=1)
    # bool tanure, building_type
    scores_mat[:, 1: 3] = (match_matrix[:, 5:7] != row[5:7])         * weight[4:6]
    # numeric num_of_unit, property_values/rent, year_built, and income/biv
    scores_mat[:, 3:8]  = np.abs(match_matrix[:, 7:12] - row[7:12]) * weight[6:11] 
    return scores_mat.sum(axis=1).argsort()

#credits goes to Geekforgeek
@numba.jit(nopython=True)
def wPrefersM1OverM(prefer, w, m, m1, N):
    for i in range(N):
        if (prefer[w][i] == m1):
            return True
        if (prefer[w][i] == m):
            return False
            
@numba.jit(nopython=True)
def stable_matching(man_prefer, wmn_prefer):
    # return list(ind=wmn, value=man)
    N = len(man_prefer)
    wPartner = [-1 for i in range(N)]
    mFree = [False for i in range(N)]
 
    freeCount = N
    while (freeCount > 0):
        m = 0
        while (m < N):
            if (mFree[m] == False):
                break
            m += 1
 
        i = 0
        while i < N and mFree[m] == False:
            w = man_prefer[m][i]
            if wPartner[w] == -1:
                wPartner[w] = m
                mFree[m] = True
                freeCount -= 1
            else:
                m1 = wPartner[w]
                if (wPrefersM1OverM(wmn_prefer, w, m, m1, N) == False):
                    wPartner[w] = m
                    mFree[m] = True
                    mFree[m1] = False
            i += 1
    return wPartner

# @numba.jit(nopython=True)
def placement(encoded_buildings, encoded_households): 
    '''
    solving the stable matching program given encoded buildings and hh
    Input: 
        - Encoded Building record matrix: np.ndarray
        - Encoded HH record matrix: np.ndarray
    Return:
        - placement result
    '''
    # # of HH rows
    n = encoded_households.shape[0]
    # # of B rows
    m = encoded_buildings.shape[0]
    # calculate Building prefer matrix and HH prefer matrix 
    hh_prefer_b_rank = np.zeros((n, m))
    i = 0
    for hh in encoded_households:
        hh_prefer_b_rank[i, :] = row_rank(encoded_buildings, hh)
        i += 1
    building_prefer_hh_rank = np.zeros((m, n))
    i = 0
    for b in encoded_buildings:
        building_prefer_hh_rank[i, :] = row_rank(encoded_households, b)
        i += 1
    if n < m:
        # Adding household dummies
        dummies_for_b = np.repeat(np.arange(n, m), np.array([m])).reshape(
            (m-n, m)).transpose()
        # append dummies
        building_prefer_hh_rank = np.hstack((building_prefer_hh_rank, dummies_for_b))
        dummies_for_hh = np.repeat(
            np.arange(m), np.array([m-n])).reshape((m, m-n)).transpose()
        # append dummies
        hh_prefer_b_rank = np.vstack((hh_prefer_b_rank, dummies_for_hh))
    elif n > m:
        # Adding building dummies
        dummies_for_hh = np.repeat(np.arange(m, n), np.array([n])).reshape(
            (n-m, n)).transpose()
        # append dummies
        hh_prefer_b_rank = np.hstack((hh_prefer_b_rank, dummies_for_hh))
        dummies_for_b = np.repeat(
            np.arange(n), np.array([n-m])).reshape((n, n-m)).transpose()
        # append dummies
        building_prefer_hh_rank = np.vstack((building_prefer_hh_rank, dummies_for_b))
    # list: ind=building_id value=hh_id
    # take min(n, m) count as final result
    # the remains HH&building used for higher level matching
    result = stable_matching(hh_prefer_b_rank.astype(np.int64), building_prefer_hh_rank.astype(np.int64))
    # only keep l result excluding dummies
    l = min(n, m)
    if n > m:
        # if more hh, take the first m record, all buildings have been assigned
        # match_building_ind = np.arange(m)
        match_building_ind = [i for i in range(m)] 
        result_wo_dum = np.array(result[0: l], dtype=np.int64)
    else:
        # if more or equal b, take the first n record, all hh have been assigned
        # get the buildings which have valid HH assigned
        # init the list with dtype int
        match_building_ind = [1]
        # remove 1 from list 
        match_building_ind.pop()
        for i, hh_ind in enumerate(result):
            if hh_ind in range(n):
                match_building_ind.append(i)
        result_wo_dum = np.array([result[i]
                                  for i in match_building_ind], dtype=np.int64)
    return (
        result_wo_dum,
        np.array(match_building_ind, dtype=np.int64)
    )

def chunked_match(out, cant_match, con_l, syn_bg_count):
    '''
    chunked match function
    Arguments:
        out: buildings dataframe 
        cant_match: households to match 
        con_l: function to get geo id from Building
        syn_bg_count: HH synthesis in bg summary
    return:
        cant_match: hh records which can't be matched
    '''
    # check arguments
    if cant_match.shape[0] == 0:
        print("Empty cant_match skip, chunked_match")
        return out, cant_match
    out = out.sort_values(by=con_l)
    # group by building by the geo_id con_l
    building_groupby = out.groupby(by=con_l)
    building_groups = building_groupby.groups
        
    # group by HH by the geo_id con_l 
    hh_groupby = cant_match.groupby(by=con_l)
    hh_groups = hh_groupby.groups

    # init new_cant_match_ind with HH groups which cannot match
    new_cant_match_ind = []
    # debug
    # total_assigned_count = 0
    for geo_id, _ in building_groups.items():
        geo_slice = building_groupby.get_group(geo_id)
        ##
        # use only the non-occupied building units, not unplacing hh which 
        # have been assigned in the previous step
        geo_slice = geo_slice.loc[geo_slice["matched_household_id"] == -1]
        # if this geo don't have any HH in it, skip
        if geo_id not in hh_groups:
            continue
        hh_insert = hh_groupby.get_group(geo_id)
        # convert both building and hh slice to matrix
        geo_slice_np = geo_slice.to_numpy().astype(int)
        hh_insert_np = hh_insert.to_numpy().astype(int)

        ##
        # max matrix size set to reduce the matrix size for faster speed
        max_matrix_size = 1000
        n = hh_insert_np.shape[0]
        m = geo_slice_np.shape[0]

        # initialized HU index for the geo
        hu_index = np.arange(m)
        out_index = geo_slice.index
        # initialize start index for HH
        start_ind = 0
        # debug
        # init_assigned_hh = (out.loc[geo_slice.index, 'matched_household_id'] != -1).sum()
        while hu_index.shape[0] > 0:
            # if #of hh exceed max size, create chunk with size min(n, max_matrix_size)
            num_hu = min(m, max_matrix_size)
            num_hh = min(n, max_matrix_size)
            # if index exceed range, use n
            end_ind = min(start_ind + num_hh, n)
            hh_sample_ind = np.arange(start_ind, end_ind, dtype=np.int64)
            # get the HH split by index 
            hh_split_i = hh_insert_np[start_ind:end_ind, :]
            # get group by bg
            hu_bg_unique, hu_bg_unique_count = np.unique(geo_slice_np[hu_index][:, [1, 3, 4]], axis=0, return_counts=True)
            if len(hu_bg_unique_count) > 1:
                # if more than one bg in this geo
                # generate importance score based on syn HH summary table by BG
                bg = np.array([syn_bg_count[np.all(np.isin(syn_bg_count[:, :3], ar), axis=1)][0]
                            # if the bg is found in the synthesis count
                        if np.all(np.isin(syn_bg_count[:, :3], ar), axis=1).sum()!=0
                            # if not found, use count 1
                        else np.append(ar, 1)
                        for ar in hu_bg_unique])
                # bg = syn_bg_count[np.all(np.isin(syn_bg_count[:, :3], hu_bg_unique), axis=1)]
                bg_count = np.repeat(bg[:, -1], hu_bg_unique_count) if len(bg) > 0 else np.array([])
                p = np.divide(bg_count, np.sum(bg_count)) if len(bg_count) > 0 else None
            else:
                # if only one bg in this geo, skip this step
                p = None

            # get sampled HU split by random sampling 
            hu_sample_ind = np.random.choice(m, num_hu, replace=False, p=p)

            # not including the last col(matched_hh_id)
            hu_split_i = geo_slice_np[hu_index][hu_sample_ind, :-1]

            # run placement algo on the set of Building and HH
            # building_result list: ind=building_id value=hh_id

            match_result, match_building_ind = placement(hu_split_i, hh_split_i)
            # update buildings with match hh_id
            hu_ind_to_update = out_index[hu_index[match_building_ind]]
            out.loc[hu_ind_to_update, 'matched_household_id'] = [
                hh_insert_np[hh_sample_ind[hh_ind], 0] for hh_ind in match_result]

            ## important
            # remove assigned hu from the pool
            hu_index = np.delete(hu_index, match_building_ind)

            # current_geo_slice = geo_slice_np[[hu_index], :]
            m = hu_index.shape[0]

            # update start_ind
            start_ind += num_hh
            # if next start ind > hh length, stop while loop
            if start_ind >= n:
                # debug
                # hh_assigned = (out.loc[geo_slice.index, 'matched_household_id'] != -1).sum() - init_assigned_hh
                # total_assigned_count += hh_assigned
                # print("break while loop", "hh_assigned", hh_assigned, "out of", n, "hu", num_hu, "total_assigned_count", total_assigned_count)
                break
    # add all unmatched HH to cant_match
    matched_hh = out[out.matched_household_id != -1].matched_household_id.values
    # get all unmatched HH
    new_cant_match_ind += cant_match[ 
                    # didn't get assigned to HU
                    ~cant_match.household_id.isin(matched_hh) &
                    # avoid duplicates
                    ~cant_match.household_id.isin(new_cant_match_ind)].index.to_list()
    # avoid duplicates
    new_cant_match_ind = list(set(new_cant_match_ind))
    return (out, cant_match.loc[new_cant_match_ind])

def run_placement(households, buildings, voters_registration):
    ###
    skip_preprocess = False 
    if not skip_preprocess:
        print("Processing buildings...")
        buildings = preprocess_buildings(buildings)
        print("Calculating housing units by bg...")
        mcd_by_bg, hu_by_mcd_bg = calculate_housing_units_by_bg(buildings)
        print("Processing households...")
        households = preprocess_households(households, mcd_by_bg, hu_by_mcd_bg)
        print("Saving buildings and households table...")
        households.to_csv(placement_config['hh_path'], float_format='%.0f', index=False)
        buildings.to_csv(placement_config["buildings_path"], index=False)
        syn_bg_count = households[['county', 'tract', 'block group', 'household_id']
                                ].groupby(by=['county', 'tract', 'block group']
                                            ).count().reset_index().to_numpy()
        print("Recoding households table...")
        recoded_hh_np = recode_households_np(
            households.to_numpy(),
            recode_hh_bld, recode_nu_hh, recode_hpv, recode_rent, recode_yb, recode_income, recode_p_over18
        )
        recoded_households = pd.DataFrame(
            recoded_hh_np, 
            columns=[
                'household_id', 'county', 'mcd', 'tract', 
                'blockgroup', 'tenure', 'buildingtype', 
                'recode_nu', 'recode_hc', 
                'recode_yb', 'recode_income', 'recode_p_over18'
                ],
            dtype=int
        )

        print("Recoding buildings table...")
        # generate HU from building table
        recoded_buildings_np = recode_buildings_table_np(
            buildings.to_numpy(),
            recode_building_type, recode_nu_building, recode_biv, recode_yb
        )
        # generate p_over18 for HU table using voters' dataset
        t0 = time.time()
        print("Generating n_over18 for HU table...")
        recoded_buildings_np = generate_n_over18_for_hu(recoded_buildings_np, voters_registration)
        t1 = time.time()
        hu_columns = [
                'building_id', 'county', 'mcd', 'tract', 
                'blockgroup', 'tenure', 'buildingtype', 
                'recode_nu', 'recode_biv', 'recode_yb', 'recode_bc',
                'recode_p_over18'
                ]
        recoded_buildings = pd.DataFrame(
            recoded_buildings_np, 
            columns=hu_columns,
            dtype=int
        )
        print("Saving recoded households and housing_units tables...")
        recoded_buildings.to_csv(placement_config['recoded_building_path'], index=False)
        recoded_households.to_csv(placement_config['recoded_hh_path'], index=False)
        print("saved recoded done.")
    else:
        print("Skipping preprocessing households and buildings table...")
        households = pd.read_csv(placement_config['hh_path'])
        buildings = pd.read_csv(placement_config['buildings_path'])
        print("Loading recoded buildings table from %s..." % placement_config['recoded_building_path'])
        recoded_buildings = pd.read_csv(placement_config['recoded_building_path'])
        print("Loading recoded households table from %s..." % placement_config['recoded_hh_path'])
        recoded_households = pd.read_csv(placement_config['recoded_hh_path'])
    ### use buildings earlier than 2020 year_built != 2020
    # buildings = buildings[buildings.year_built != 2020]
    ### use only couty 125 as testing data
    # households = households.loc[households.county == 125]
    # buildings = buildings.loc[buildings.county == 125]
    # generate syn_HH bg summary


    print("reading from %s" % placement_config["buildings_path"])
    print("read %s lines" % buildings.shape[0])

    print("reading from %s" % placement_config["hh_path"])
    print("read %s lines" % households.shape[0])

    print("matching %s HH with %s HU" % (recoded_households.shape[0], recoded_buildings.shape[0]))

    # Rust implementation
    #############
    # binary_path = 'temp_for_placement/bildings'
    # with sbp.Popen([binary_path], stdout=sbp.PIPE) as proc:
    #     print('Start running placement executable...')
    #     while proc.poll() is None: # Check if still running
    #         print(proc.stdout.readline().decode('utf-8'))
    #     print('end')

    # match = read_match()
    # match.index.names = households.index.names
    #############
    # Python impl
    #############
    # Add the two new col
    out = recoded_buildings.copy(deep=True)
    out['matched_household_id'] = -1
    out['placed_at'] = None
    # out['compared_result'] = None
    # shuffle out 
    out = out.sample(random_state=42, frac=1)
    # move HH data to cant_match
    cant_match = recoded_households.copy(deep=True)
    # sort by hh id
    cant_match = cant_match.loc[
        cant_match['household_id'].sort_values().index
    ]
    geo_to_iter = [
        ['county', 'mcd', 'tract','blockgroup'],
        ['county', 'mcd', 'tract'], 
        # geos below take long time to run
        ['county', 'mcd'], 
        ['county']
    ]
    # iter through geo types
    t_i = time.time()
    for i in range(len(geo_to_iter)):
        # by all geo_id
        match_result, cant_match = chunked_match(
            out[
                hu_columns + ['matched_household_id']
            ],
            cant_match,
            geo_to_iter[i],
            syn_bg_count 
        )
        non_placed_hu = out[out['placed_at'].isnull()]
        placed_hu = match_result[
            match_result['matched_household_id'] != -1]
        matched_index = np.intersect1d(non_placed_hu.index, placed_hu.index)
        out.loc[matched_index, 'placed_at'] = ','.join(geo_to_iter[i])
        out['matched_household_id'] = match_result['matched_household_id']
        print(
            "matching ", geo_to_iter[i],
            "took %.2fs with %s left to match, %s HH got assigned, %s placed_at got editted" %
            (time.time() - t_i,
             cant_match.shape[0],
             out[out['matched_household_id'] != -1].shape[0],
             len(matched_index)
             )
        )
        t_i = time.time()
    cant_match.to_csv('test_run/092324_run_2015/cant_match_region.csv', index=False)
    out.to_csv('test_run/092324_run_2015/match_result_region_temp.csv', index=False)
    return out 

