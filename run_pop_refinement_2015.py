from forecast_input.pop_refinement import refine_pop_single_year
import pandas as pd

refinement_excel = "/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/semcog_estimates/PopHHEstimates715.xlsx"
hdf_path = '/mnt/hgfs/urbansim/RDF2045/data/base_year/all_semcog_data_02-02-18-final-forecast.h5'
hh_path = '/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2015(2017)/092324_run_2015/households.csv'
p_path = '/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2015(2017)/092324_run_2015/persons.csv'

def main():
    refine_geo = 'semmcd'
    year = 2015
    # process refinemnt excel to get adjustment
    refinement = pd.read_excel(refinement_excel, sheet_name=0)
    refinement.columns = refinement.columns.str.lower()
    refinement = refinement[refinement[refine_geo]<8000]
    refinement = refinement.set_index(refine_geo)
    refinement.index = refinement.index.astype(int)
    refinement = refinement.fillna(0).astype(int)
    # refinement = refinement[2020 if year=='base' else year]
    refinement = refinement['totalpop'] - refinement['gqpop']
    refinement = refinement[refinement>0]
    # refine_pop_single_year(year, hh, p, b, refinement):

    # 2015 HDF
    # read HDF cache
    hdf = pd.HDFStore(hdf_path, 'r')
    # get building
    b = hdf['buildings']
    parcels = hdf['parcels']
    b = b.join(parcels[[refine_geo]], on='parcel_id')
    b = b[~b[refine_geo].isna()]
    b = b.astype({refine_geo: int})
    # hh and p
    hh = pd.read_csv(hh_path, index_col=0)
    p = pd.read_csv(p_path, index_col=0)

    hh = hh.join(b[[refine_geo]], on='building_id')
    p = p.join(hh[[refine_geo]], on='household_id')
    # hh.columns = hh.columns.str.lower()
    hh = hh[
        [
            col
            for col in hh.columns
            if col.lower()
            not in [
                "puma",
                "tract",
                "hincp",
                "r18",
                "hhisp",
                "adjinc",
                "ybl",
                "bld",
                "grntp",
                "adjhsg",
                "type",
                "valp",
            ]
        ]
    ]
    hh = hh.rename(columns={
        'AGEHOH': 'age_of_head',
        'VEN': 'cars',
        'NP': 'persons',
        'HRACE': 'race_id',
        'HWORKERS': 'workers',
        'HHT': 'hht',
        'hh_id': 'seed_id'
    })
    # semmcd 7027 have 0 buildings, 0 parcels
    new_hh, new_p = refine_pop_single_year(hh, p, b, refinement, refine_geo)
    new_hh.to_csv('/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2015(2017)/092324_run_2015/households_after_refinement.csv')
    new_p.to_csv('/mnt/hgfs/urbansim/RDF2050/population_synthesis/historical/2015(2017)/092324_run_2015/persons_after_refinement.csv')
    return
    
if __name__ == "__main__":
    main()
    print('Done.')

# empty city semmcd     7027
# to_drop     -10
# Name: 221, dtype: int64
# dropping 23248 hh and 54864 persosn
# adding 20890 hh and 49930 persosn
# semmcd 7027 has 0 buildings, skipped