""" This script will forcefully ensure hh + gq_pop = target_pop by
adding/dropping household
** This step will add or drop hh, be careful running this script **
"""

import pandas as pd
import numpy as np

BASED_HDF_PATH = '~/semcog_urbansim/runs/run2120_TAZ_draft_final.h5'
BASE_SIM_RUN_PATH = '/mnt/hgfs/RDF2050/run2120_may_refinement_v3.h5'

REFINE_PATH = '/home/da/share/urbansim/RDF2050/model_outputs/demographic/PopByCity_allAges_Final_gqadj.xlsx'
GQ_PATH = '/mnt/hgfs/RDF2050/gq/gq_%s_match_gq_age.csv'

def main(based_hdf_path, base_sim_run_path, refine_path, gq_path):
    hh_hdf = pd.HDFStore(based_hdf_path, 'a')
    hdf = pd.HDFStore(base_sim_run_path, 'r')
    # loop through baseyear and forecast years
    for year in ['base'] + list(range(2025,2051,5)):
        hh = hh_hdf['/%s/households'%year]
        hh_col = hh.columns
        p = hh_hdf['/%s/persons'%year]
        p_col = p.columns
        gq = pd.read_csv(gq_path % year)
        b = hdf['/%s/buildings'%year]
        hh = hh.join(b[['city_id']], on='building_id')
        p = p.join(hh[['city_id']], on='household_id')
        gq = gq.join(b[['city_id']], on='building_id')
        # load total_pop refinement 
        refinement = pd.read_excel(refine_path, sheet_name='TotalPop')
        # filter out large area(<500) and total(8888)
        refinement = refinement[refinement.city_id>500]
        refinement = refinement[refinement.city_id<8000]
        # set index to city_id
        refinement = refinement.set_index('city_id')
        refinement.index = refinement.index.astype(int)
        refinement = refinement[2020 if year=='base' else year]
        refinement = refinement[refinement>0]
        # calculate adj
        p_by_city = p.groupby('city_id').size()
        gq_by_city = gq.groupby('city_id').size()
        adj_yr = -refinement.sub(p_by_city, fill_value=0).sub(gq_by_city, fill_value=0).astype(int)
        adj_yr.name = 'to_drop'
        adj_yr = adj_yr.reset_index()
        
        # adj_yr = adj[adj.year==str(year)]
        # adj_yr = adj_yr.groupby(['city_id']).sum().to_drop.reset_index()
        hh_id_to_drop = []
        hh_id_to_add = []
        remains = []
        for ind, row in adj_yr.iterrows():
            # only adjusting hhsize <7
            to_drop = row.to_drop
            # loop until meet the pop target
            while to_drop != 0:
                # the hh pool for the city
                small_hh_pool = hh[(hh.city_id==row.city_id)&(hh.persons<7)&(~hh.index.isin(hh_id_to_drop))&(~hh.index.isin(hh_id_to_add))]
                if small_hh_pool.persons.sum() < to_drop:
                    # if not enough small hh to drop, use the whole hh pool
                    pool = hh[(hh[refine_geo]==row[refine_geo])&(~hh.index.isin(hh_id_to_drop))&(~hh.index.isin(hh_id_to_add))]
                else:
                    # less than 7 persons hh, not yet been touched
                    pool = small_hh_pool
                # if need to add and number is less than 5 and there is no hh in pool that could satisfy that
                if to_drop < 0 and abs(to_drop) < 5 and -to_drop not in pool.persons.values:
                    # add to remains to solve it later
                    remains.append(ind)
                    adj_yr.loc[ind, 'to_drop'] = to_drop
                    print('cannot solve', row)
                    break
                # if need to drop
                if to_drop > 0:
                    if to_drop < 7 and pool[pool.persons==to_drop].shape[0]>0:
                        # if to drop number less that 7, sample it from pool
                        local_pool = pool[pool.persons==to_drop]
                    else:
                        local_pool = pool
                    selected = local_pool.sample(1)
                    hh_id_to_drop += list(selected.index)
                    # update to_drop
                    to_drop -= selected.persons.sum()
                # if need to add
                else:
                    if to_drop > -7 and pool[pool.persons==-to_drop].shape[0]>0:
                        # if to add number less that 7, sample it from pool
                        local_pool = pool[pool.persons==-to_drop]
                    else:
                        local_pool = pool
                    if local_pool.shape[0] == 0:
                        # if empty pool, add to remains to solve it later
                        remains.append(ind)
                        print('empty city', row)
                        break
                    selected = local_pool.sample(1)
                    hh_id_to_add += list(selected.index)
                    # update to_drop
                    to_drop += selected.persons.sum()

        # drop hh and members
        print('dropping %s hh and %s persosn' % (len(hh_id_to_drop), hh.loc[hh_id_to_drop].persons.sum()))
        new_hh = hh[~hh.index.isin(hh_id_to_drop)]
        new_p  = p[p.household_id.isin(new_hh.index)]
        # add hh (**index**)
        print('adding %s hh and %s persosn' % (len(hh_id_to_add), hh.loc[hh_id_to_add].persons.sum()))
        hh_id_to_add.sort()
        add_hh = hh.loc[hh_id_to_add].copy()
        add_per = p[p.household_id.isin(add_hh.index)].sort_values('household_id').copy()
        hh_id_start = new_hh.index.max()+1
        per_id_start = add_per.index.max()+1
        # update index
        add_hh.index = np.arange(hh_id_start, hh_id_start+len(add_hh), 1)
        add_per.index = np.arange(per_id_start, per_id_start+len(add_per), 1)
        # update ppl household_id
        add_per.loc[:, 'household_id'] = np.repeat(add_hh.index, add_hh.persons)
        new_hh = pd.concat((new_hh, add_hh), axis=0)
        new_p = pd.concat((new_p, add_per), axis=0)
        end_hh_ind = new_hh.index.max()+1
        end_ppl_ind = new_p.index.max()+1
        # fix remains in the end
        for ind in remains:
            row = adj_yr.loc[ind]
            to_drop = row.to_drop
            # all hh as pool
            pool = hh
            bpool = b[(b.city_id==row.city_id)&(b.residential_units>0)]
            # make sure to_drop < 0, all process below are adding hh
            assert to_drop < 0
            while to_drop < 0:
                if to_drop > -7:
                    local_pool = pool[pool.persons==-to_drop]
                else:
                    local_pool = pool
                selected_hh = local_pool.sample(1).copy()
                selected_ppl = p[p.household_id.isin(selected_hh.index)]
                selected_hh.index = [end_hh_ind]
                selected_ppl.index = np.arange(end_ppl_ind, end_ppl_ind+selected_hh.persons.sum(), 1)
                hu = bpool.sample(1)
                selected_hh['building_id'] = hu.index[0]
                selected_ppl['household_id'] = end_hh_ind
                new_hh = pd.concat((new_hh, selected_hh), axis=0)
                new_p = pd.concat((new_p, selected_ppl), axis=0)
                end_hh_ind += 1
                end_ppl_ind += selected_hh.persons.sum()
                to_drop += selected_hh.persons.sum()


        hh_hdf['/%s/households_adj'%year] = new_hh[hh_col]
        hh_hdf['/%s/persons_adj'%year] = new_p[p_col]
    hdf.close()
    hh_hdf.close()
    return

def refine_pop_single_year(hh, p, b, refinement, refine_geo='city_id'):
    hh_col = hh.columns
    p_col = p.columns
    # hh = hh.join(b[[refine_geo]], on='building_id')
    # p = p.join(hh[[refine_geo]], on='household_id')
    # load total_pop refinement 
    # calculate adj
    p_by_city = p.groupby(refine_geo).size()
    adj_yr = -refinement.sub(p_by_city, fill_value=0).astype(int)
    adj_yr.name = 'to_drop'
    adj_yr = adj_yr.reset_index()

    hh_id_to_drop = []
    hh_id_to_add = []
    remains = []
    for ind, row in adj_yr.iterrows():
        # only adjusting hhsize <7
        to_drop = row.to_drop
        # loop until meet the pop target
        while to_drop != 0:
            # the hh pool for the city
            # less than 7 persons hh, not yet been touched
            pool = hh[(hh[refine_geo]==row[refine_geo])&(~hh.index.isin(hh_id_to_drop))&(~hh.index.isin(hh_id_to_add))]
            # if need to add and number is less than 5 and there is no hh in pool that could satisfy that
            if to_drop < 0 and abs(to_drop) < 5 and -to_drop not in pool.persons.values:
                # add to remains to solve it later
                remains.append(ind)
                adj_yr.loc[ind, 'to_drop'] = to_drop
                print('cannot solve', row)
                break
            # if need to drop
            if to_drop > 0:
                if to_drop < 7 and pool[pool.persons==to_drop].shape[0]>0:
                    # if to drop number less that 7, sample it from pool
                    local_pool = pool[pool.persons==to_drop]
                else:
                    local_pool = pool
                selected = local_pool.sample(1)
                hh_id_to_drop += list(selected.index)
                # update to_drop
                to_drop -= selected.persons.sum()
            # if need to add
            else:
                if to_drop > -7 and pool[pool.persons==-to_drop].shape[0]>0:
                    # if to add number less that 7, sample it from pool
                    local_pool = pool[pool.persons==-to_drop]
                else:
                    local_pool = pool
                if local_pool.shape[0] == 0:
                    # if empty pool, add to remains to solve it later
                    remains.append(ind)
                    print('empty city', row)
                    break
                selected = local_pool.sample(1)
                hh_id_to_add += list(selected.index)
                # update to_drop
                to_drop += selected.persons.sum()

    # drop hh and members
    print('dropping %s hh and %s persosn' % (len(hh_id_to_drop), hh.loc[hh_id_to_drop].persons.sum()))
    new_hh = hh[~hh.index.isin(hh_id_to_drop)]
    new_p  = p[p.household_id.isin(new_hh.index)]
    # add hh (**index**)
    print('adding %s hh and %s persosn' % (len(hh_id_to_add), hh.loc[hh_id_to_add].persons.sum()))
    hh_id_to_add.sort()
    add_hh = hh.loc[hh_id_to_add].copy()
    add_per = p[p.household_id.isin(add_hh.index)].sort_values('household_id').copy()
    hh_id_start = new_hh.index.max()+1
    per_id_start = add_per.index.max()+1
    # update index
    add_hh.index = np.arange(hh_id_start, hh_id_start+len(add_hh), 1)
    add_per.index = np.arange(per_id_start, per_id_start+len(add_per), 1)
    # update ppl household_id
    add_per.loc[:, 'household_id'] = np.repeat(add_hh.index, add_hh.persons)
    new_hh = pd.concat((new_hh, add_hh), axis=0)
    new_p = pd.concat((new_p, add_per), axis=0)
    end_hh_ind = new_hh.index.max()+1
    end_ppl_ind = new_p.index.max()+1
    # fix remains in the end
    for ind in remains:
        row = adj_yr.loc[ind]
        to_drop = row.to_drop
        # all hh as pool
        pool = hh
        bpool = b[(b[refine_geo]==row[refine_geo])&(b.residential_units>0)]
        if bpool.shape[0] == 0:
            print(refine_geo, row[refine_geo], 'has 0 buildings, skipped')
            continue
        # make sure to_drop < 0, all process below are adding hh
        assert to_drop < 0
        while to_drop < 0:
            if to_drop > -7:
                local_pool = pool[pool.persons==-to_drop]
            else:
                local_pool = pool
            selected_hh = local_pool.sample(1).copy()
            selected_ppl = p[p.household_id.isin(selected_hh.index)]
            selected_hh.index = [end_hh_ind]
            selected_ppl.index = np.arange(end_ppl_ind, end_ppl_ind+selected_hh.persons.sum(), 1)
            hu = bpool.sample(1)
            selected_hh['building_id'] = hu.index[0]
            selected_ppl['household_id'] = end_hh_ind
            new_hh = pd.concat((new_hh, selected_hh), axis=0)
            new_p = pd.concat((new_p, selected_ppl), axis=0)
            end_hh_ind += 1
            end_ppl_ind += selected_hh.persons.sum()
            to_drop += selected_hh.persons.sum()
    return new_hh, new_p

def compile_final_hdf(based_hdf_path, gq_path):
    # add households, persons and gq to final hdf
    final_hdf = pd.HDFStore(based_hdf_path, 'a')
    reviewed_hdf = pd.HDFStore(BASE_SIM_RUN_PATH, 'r')

    for year in ['base'] + list(range(2025, 2051, 5)):
        final_hdf['/%s/households' % year] = final_hdf['/%s/households' % year].sort_index()
        final_hdf['/%s/persons' % year] = final_hdf['/%s/persons' % year].sort_index()
        USE_GQ_REFINEMENT = True
        if USE_GQ_REFINEMENT:
            final_hdf['/%s/group_quarters' % year] = pd.read_csv(gq_path % year)
        else:
            final_hdf['/%s/group_quarters' % year] = reviewed_hdf['/%s/group_quarters' % year].sort_index()
    final_hdf.close()
    reviewed_hdf.close()
    return

if __name__ == '__main__':
    # main(BASED_HDF_PATH, BASE_SIM_RUN_PATH, REFINE_PATH, GQ_PATH)
    compile_final_hdf(BASED_HDF_PATH, GQ_PATH)