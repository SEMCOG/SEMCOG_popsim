import pandas as pd
import numpy as np
# File includes checking functions 

def default_checker(df):
  check_null_values(df)
  return

def check_null_values(df):
  if df.isnull().any().any():
    print("Null value found in the table ", df.columns)
    print(df.head)

def find_dtype(s):
    infer = pd.api.types.infer_dtype(s)
    if infer == 'unicode' or infer == 'string':
        return s.astype(str)

    if infer == 'mixed':
        if (s.astype(unicode) == s).all():
            return s.astype(str)
        print(s.name)
        print(s.apply(type).drop_duplicates())

    def match(s, s_as):
        return abs(s - s_as).max() == 0

    try:
        s_as_int = s.astype(np.int64)
        if match(s, s_as_int):
            # we get overflows with np.int16
            for dtype in [np.int8, np.int32]:
                s_as = s.astype(dtype)
                if match(s, s_as):
                    return s_as
            return s_as_int
    except ValueError:
        pass

    try:
        s_as_float = s.astype(np.float64)
        if match(s, s_as_float):
            for dtype in [np.float16, np.float32]:
                s_as = s.astype(dtype)
                if match(s, s_as):
                    return s_as
            return s_as_float
    except ValueError:
        pass

    return s

def default_cleaner(df):
    # type: (pd.DataFrame) -> pd.DataFrame
    df.columns = [str(x).lower().strip() for x in df.columns]
    df.index.names = [(None if x is None else str(x).lower().strip()) for x in df.index.names]
    df = df.apply(find_dtype)

    if pd.api.types.infer_dtype(df.index.values) == 'unicode':
        df.index = df.index.astype(str)
    if hasattr(df.index, "levels"):
        for i in range(len(df.index.levels)):
            df.index = df.index.set_levels(find_dtype(df.index.levels[i]), level=i)
    else:
        df.index = find_dtype(df.index)

    return df

def checking_buildings(buildings):
    assert ~buildings.building_id.duplicated()
    assert (buildings.sqft_per_unit < 0).sum() == 0
    assert (buildings.res_sqft < 0).sum() == 0
    assert (buildings.market_value < 0).sum() == 0

def checking_households(households):
    assert ~households.household_id.duplicated()

def checking_persons(persons, hdf):
    assert ~persons.household_id.duplicated()
    # household_id in the households table index
    assert (~(persons.household_id.isin(hdf.households.index))).sum() == 0
    # check to make sure persons = household persons count 
    p_size = persons.groupby('household_id').size()
    assert len(hdf.households[hdf.households.persons - p_size != 0]) == 0