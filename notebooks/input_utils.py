# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
from IPython import get_ipython

# %%
### generate control totals needed by populationsim(RSG)


# %%
import os
import pandas as pd
import re
from census import Census



# %%
class ACS5_downloader:

    def __init__(self, census_reader, states, counties = None, tract_ids = None, blockgroup_ids = None):
        self.states = states
        self.counties = counties
        self.tracts = tract_ids
        self.blockgroups = blockgroup_ids 
        self.cread = census_reader
        self.udpate_states_counties([self.states, self.counties])
            
    def state_download(self, vars):
        return self.cread.acs5.get(vars, geo={'for': 'state:{}'.format(self.states)})

    def county_download(self, vars):
        print (self.states,self.counties)
        print(vars)
        return self.cread.acs5.get(vars, geo={'for': 'county:{}'.format(self.counties), 
                                        'in': 'state:{}'.format(self.states)})

    def tract_download(self, vars):
        return self.cread.acs5.get(vars, geo={'for': 'tract:{}'.format(self.tracts), 
                                    'in': 'state:{} county:{}'.format(self.states, self.counties)})

    def blockgroup_download(self, vars):
        clst = self.counties.split(',')
        cm = []
        for cn in clst:
            cm += self.cread.acs5.get(vars, geo={'for': 'block group:{}'.format(self.blockgroups), 
                                'in': 'state:{} county:{} tract:{}'.format(self.states, cn,                                                                                     self.tracts)}) 
        return cm
    
    def fips_lookup(self, states, counties = None):
        if counties == "*": counties = None
        fips_table = pd.read_csv(
                "https://www2.census.gov/geo/docs/reference/codes/files/national_county.txt",
                header=None, names=['state','state_fips', 'county_fips', 'county' ,'type'],                          dtype=str)  
        qstr = '(state in {})'.format(states)
        if counties:
            qstr += ' & (county in {})'.format(counties)
            dfq = fips_table.query(qstr)
            return list(dfq.state_fips.unique()), list(dfq.county_fips.unique())
        dfq = fips_table.query(qstr)
        print(list(dfq.state_fips.unique()))
        return list(dfq.state_fips.unique()), None

    def udpate_states_counties(self, geos):
        for i in [0,1]:
            if geos[i] != None:
                if (type(geos[i]) != list) & (geos[i] != '*'):
                    geos[i] = [geos[i]]
                geos[i] = [str(x) for x in geos[i]]
        if (type(geos[0]) == list) & (geos[0][0].isdigit() == False):
            geos[0], geos[1] = self.fips_lookup(geos[0], geos[1])
        for i in [0,1]:
            if geos[i]: geos[i] = ','.join([str(x).zfill(i+2) for x in geos[i]]) #i+2 cause state and county need 2 and 3 0s in lead
        self.states = geos[0]
        self.counties = geos[1]


    def download(self, variables):
        dfm = pd.DataFrame()
        if not(self.counties):
            downv = self.state_download(variables)
        elif not(self.tracts):
            downv = self.county_download(variables)
        elif not(self.blockgroups):
            downv = self.tract_download(variables)
        else:
            downv = self.blockgroup_download(variables)
        dfm = pd.DataFrame.from_dict(downv)

        return dfm
 


def preprocess_pums(h_pums, p_pums):
    ### h_pums and p_pums must have the same index column 'SERIALNO'
    p_pums = p_pums.set_index('SERIALNO')

    # add AGEHOH to PUMS sample
    h_pums['AGEHOH'] = p_pums.loc[p_pums.RELP == 0].AGEP

    # add HRACE to PUMS sample
    rac_map = {1:1, 2:2, 6:3}
    h_pums['HRACE'] = p_pums.loc[p_pums.RELP == 0].RAC1P
    h_pums['HRACE'] = h_pums['HRACE'].map(rac_map).fillna(4)

    # add HHISP to PUMS sample
    hisp_map = {1:0}
    h_pums['HHISP'] = p_pums.loc[p_pums.RELP == 0].HISP
    h_pums['HHISP'] = h_pums['HHISP'].map(hisp_map).fillna(1)

    # add HWORKERS to PUMS sample
    p_pums = p_pums.reset_index()
    h_pums["HWORKERS"] = p_pums.loc[p_pums.ESR.isin([1,2,4,5])].groupby('SERIALNO').ESR.size()
    h_pums.loc[h_pums.HWORKERS >= 2, "HWORKERS"] = 2
    
    h_pums.fillna(0, inplace = True)
    #for v in ["AGEHOH","HRACE", "HHISP","HWORKERS" ]:
    #    print (v, h_pums.groupby(v).size())

    return h_pums, p_pums
