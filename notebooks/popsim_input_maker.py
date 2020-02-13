# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'


# %%
### prepare all inputs needed by populationsim(RSG)
## popsim inputs: geo_crosswalk, Census marginals and PUMS samples
## this programs runs with region.yml(as config), and input_utils.py
## this programs also needs controls_pre.csv and ACS sample HHs and persons


# %%
import os, re, time
import pandas as pd
from census import Census
import yaml
from input_utils  import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("key", help="Census API key")
parser.add_argument("yml", help="yaml configuration file name")
args = parser.parse_args()
t0 = time.time()

# %%
conf = yaml.load(open('./' + args.yml, 'r'), Loader = yaml.Loader)


# %%
c = Census(args.key, year=2017)
prj_name = conf['region']['name']
state = conf['region']['state'][0]
counties = conf['region']['counties']
pre_folder = conf['preprocess']['folder']
h_pums_csv = '../' + pre_folder + conf['preprocess']['h_pums_csv']
p_pums_csv = '../' + pre_folder + conf['preprocess']['p_pums_csv']

ouptut_geo_cross = "../{}{}_geo_cross_walk.csv".format(pre_folder,prj_name) 
output_control = '../{}{}_control_totals_.csv'.format(pre_folder, prj_name)
output_seed_hhs = '../{}{}_seed_households.csv'.format(pre_folder, prj_name)
output_seed_persons = '../{}{}_seed_persons.csv'.format(pre_folder, prj_name)

# %% [markdown]
# # make geographic cross work file

# %%
print('\nPreparing Census geographies: \n\tstate: {}  \n\tcounty: {}'.format(state, counties))
acgeo = ACS5_downloader(c, state, counties, '*', '*') #download Census BGs for this region
df_geo = pd.DataFrame.from_dict(acgeo.download('NAME')).drop('NAME', axis = 1)
df_geo['tractid'] = df_geo['state'] + df_geo['county'] + df_geo['tract']
df_geo['blkgrpid'] = df_geo['tractid'] + df_geo['block group']
df_geo.columns = [col.upper() for col in df_geo.columns]


# %%
df_tract_puma = pd.read_csv(conf['geographies']['tract_puma_file'], dtype = str)
df_tract_puma = df_tract_puma.loc[df_tract_puma.STATEFP == acgeo.states]
df_tract_puma['TRACTID'] = df_tract_puma['STATEFP'] + df_tract_puma['COUNTYFP'] + df_tract_puma['TRACTCE']
df_tract_puma.rename(columns={"PUMA5CE": "PUMA"}, inplace = True)


# %%
df_geo_cross = pd.merge(df_geo, df_tract_puma, on = 'TRACTID', how = 'left')
df_geo_cross['REGION'] = 1
print("  saving geo cross walk file: " + ouptut_geo_cross)
df_geo_cross.to_csv(ouptut_geo_cross)


# %% [markdown]
# # make maginal control files

# %%
# download all needed Census marginal controls
# Census marginal controls variables are from controls_pre table (popsim "controls" table + "acs_variables" field )
# "acs_variables" contains Census API variables and expressions
# 
print('\nMaking popsim controls')
print('  downloading Census variables ...')
dfc = pd.read_csv('../' + pre_folder + conf['preprocess']['pre_control'])
dic_margs = {}
for geo, dfgeo in dfc.groupby('geography'):  
    full_vars= list(set(re.findall(r'B[0-9]{5}[A-Z]{0,1}_[0-9]{3}E', 
                    str(list(dfgeo.acs_variables)))))
    if geo == 'BLKGRP':
        ac5 = ACS5_downloader(c, state, counties, "*", "*")
        geo_cols = ['state', 'county', 'tract', 'block group']
    elif geo == 'TRACT':
        ac5 = ACS5_downloader(c, state, counties, "*")
        geo_cols = ['state', 'county', 'tract']
    print('\t' + geo + ' marginals ')
    dic_margs[geo] = ac5.download(full_vars).set_index(geo_cols)


# %%
# Compute popsim control variables from Census marginals 
print('  compiling popsim control fields ...')
for geo, dfg in dfc.groupby('geography'):
    dic_margs[geo] = dic_margs[geo].astype(float).fillna(0)
    for ind, r in dfg.iterrows():
        #print(r['control_field'] + ":  " + r['acs_variables'])
        dic_margs[geo][r['control_field']] = dic_margs[geo].eval(r['acs_variables'].replace('"', ''))
    dic_margs[geo] = dic_margs[geo][list(dfg.control_field)] #keep only control fields


# %%
# add unique geoids and PUMA 
for geo, dfm in dic_margs.items():

    #dfcross = pd.read_csv(geo_cross_csv, dtype = str)
    if dfm.index.nlevels == 3:
        dfm['TRACTID'] = [l1 + l2 + l3 for l1, l2, l3 in dfm.index]
        df_geo_tract = df_geo_cross.drop_duplicates('TRACTID')
        dfm = pd.merge(dfm.reset_index(), df_geo_tract[['TRACTID','PUMA']], on = 'TRACTID', how = 'left')
    elif dfm.index.nlevels == 4:
        dfm['BLKGRPID'] = [l1 + l2 + l3 + l4 for l1, l2, l3, l4 in dfm.index]
        dfm = pd.merge(dfm.reset_index(), df_geo_cross[['BLKGRPID','PUMA']], on = 'BLKGRPID',how = 'left')
    dfm.columns = [col.upper() for col in dfm.columns]
    output_control = output_control.replace('.csv', geo.lower() + '.csv')
    print('  saving control file: ' + output_control)
    dfm.to_csv(output_control)

# %% [markdown]
# # extract PUMS seed households and persons

print('\nExtrating PUMS seed households and persons')

# %%
puma_lst=df_tract_puma.loc[(df_tract_puma.STATEFP == acgeo.states)
                        &(df_tract_puma.COUNTYFP.isin(acgeo.counties.split(',')))].PUMA.unique()


# %%
h_pums =  pd.read_csv(h_pums_csv, index_col="SERIALNO")
h_pums = h_pums.loc[h_pums.PUMA.isin(puma_lst)]
h_pums = h_pums.loc[(h_pums.TYPE == 1) & (h_pums.NP > 0 )] #remove group quarters and empty units

p_pums = pd.read_csv(p_pums_csv)
p_pums = p_pums.loc[p_pums.PUMA.isin(puma_lst)]
p_pums = p_pums.loc[p_pums["SERIALNO"].isin(h_pums.index)]


# %%
h_pums, p_pums = preprocess_pums(h_pums, p_pums)

h_pums['hh_id'] = h_pums.index.values
p_pums['hh_id'] = p_pums.SERIALNO

print('  saving seed households: {}   . {} records'.format(output_seed_hhs, str(len(h_pums))))
h_pums.to_csv(output_seed_hhs)
print('  saving seed persons: {}   . {} records'.format(output_seed_persons,str(len(p_pums))))
p_pums.to_csv(output_seed_persons)

print('Done.', '\ntotal time: {} seconds'.format(round(time.time()-t0, 1)))

