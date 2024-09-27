"""
Group quarters placement
"""
import pandas as pd
import numpy as np

def run_gq_placement(gqs, ps):
    """
    Parameters:
    gqs: GroupQuarters with
        - building_id, county_id, city_id, gq_code, licensed_beds
    ps: Synthesis gq persons with
        - person_id, large_area_id, city_id, gq_code
    Return
    ps_with_bid: Synthesis gq persons with building_id
    """
    # group by cty_id, city_id, gq_code and b_id
    total_beds_by_cat = gqs.groupby([ 'city_id', 'gq_code', 'building_id']).sum()['licensed_beds']

    # initialize b_id with -1
    ps['building_id'] = -1
    for ind, row in ps.iterrows():
        # category index pairs
        cat = (row['city_id'], row['gq_code'])
        beds_by_b = total_beds_by_cat.loc[cat]
        if beds_by_b.sum() == 0:
            # not beds available
            print('persons ', cat, ' couldnot find a bed')
            continue
        # randomly pick one bed
        picked_bid = np.random.choice(beds_by_b.index.repeat(beds_by_b.values), 1)[0]
        # assign bid
        row['building_id'] = picked_bid
        # reduce beds count
        beds_by_b[picked_bid] -= 1
    return ps