import numpy as np
import pandas as pd

# Transform.py
# Data manipulation and processing steps


def transform_buildings(buildings):
    # buildings.rename(columns={'city_id': 'b_city_id', 'zone_id': 'b_zone_id'}, inplace=True)
    buildings["owner_units"] = buildings["owner_units"].astype("int")
    buildings["tract"] = buildings["tract"].astype(np.int64)
    if "block" in buildings.columns:
        buildings["block"] = buildings["block"].astype(np.int64)
    buildings["bg"] = buildings["bg"].astype(np.int64)
    return buildings


def calculate_improvement_values(raw_buildings, parcels):
    """ Calculate improvement values from building and parcels table
  (pd.DataFrame, pd.DataFrame) -> pd.DataFrame
  """

    raw_buildings["res_sqft"] = (
        raw_buildings.residential_units * raw_buildings.sqft_per_unit
    )
    a = (
        raw_buildings[["parcel_id", "non_residential_sqft", "res_sqft"]]
        .groupby("parcel_id")
        .sum()
    )
    a["landvalue"] = parcels.landvalue
    a["sev_value"] = parcels.sev_value
    # set missing landvalue to 0
    a.loc[a.landvalue.isnull(), "landvalue"] = 0
    # parcel_sqft
    a["parcel_sqft"] = a.res_sqft + a.non_residential_sqft
    # join parcel_sqft, landvalue, sev_value back to raw_buildings
    raw_buildings = raw_buildings.join(a["parcel_sqft"], "parcel_id")
    raw_buildings = raw_buildings.join(a["landvalue"], "parcel_id")
    raw_buildings = raw_buildings.join(a["sev_value"], "parcel_id")
    # market_value = (2 * sev_value) * ( building_sqft / parcel_sqft)
    raw_buildings["market_value"] = (
        2
        * raw_buildings["sev_value"]
        * (
            (raw_buildings["res_sqft"] + raw_buildings["non_residential_sqft"])
            / (raw_buildings["parcel_sqft"] + 1)
        )
    )
    # improvement_value = (2 * sev_value - parcel_landvalue) * ( building_sqft / parcel_sqft)
    improvement_value = raw_buildings["market_value"] - raw_buildings["landvalue"] * (
        (raw_buildings["res_sqft"] + raw_buildings["non_residential_sqft"])
        / (raw_buildings["parcel_sqft"] + 1)
    )
    # set negative value to 0
    improvement_value[improvement_value < 0] = 0
    raw_buildings.loc[:, "improvement_value"] = improvement_value
    return raw_buildings


def calculate_housing_units_by_block(buildings):
    """Calculate housing unit count by block group summary table

  Args:
      buildings (pd.DataFrame): buildings table with residential units and block id
  """
    buildings["block_id"] = (
        buildings[["county", "tract", "bg", "block"]]
        .apply(
            lambda x: 26 * 1e13
            + x["county"] * 1e10
            + x["tract"] * 1e4
            + x["bg"] * 1e3
            + x["block"],
            axis=1,
        )
        .astype("int")
        .astype("string")
    )
    groupby = buildings[["block_id", "mcd", "residential_units"]].groupby(
        ["block_id", "mcd"]
    )
    units_by_block = groupby.sum()
    del buildings["block_id"]
    return units_by_block


def calculate_housing_units_by_bg(buildings):
    """Calculate housing unit count by block group summary table

  Args:
      buildings (pd.DataFrame): buildings table with residential units and block id
  """
    groupby = buildings[["county", "tract", "bg", "mcd", "residential_units"]].groupby(
        ["county", "tract", "bg", "mcd"]
    )
    units_by_bg = groupby.sum().reset_index()
    out_ru = {}
    out_mcd = {}
    for _, row in units_by_bg.iterrows():
        geo = (row.county, row.tract, row.bg)
        if geo not in out_ru:
            out_ru[geo] = []
            out_mcd[geo] = []
        out_ru[geo] += [row.residential_units]
        out_mcd[geo] += [row.mcd]
    return out_mcd, out_ru


def households_add_n_18plus(households, persons):
    # add 18plus variable to households table, generated from persons table
    persons = persons.reset_index()
    persons_over_18 = persons[persons["AGEP"] >= 18]
    over_18_by_hh = (
        persons_over_18[["index", "household_id"]].groupby("household_id").count()
    )
    households["p_over18"] = over_18_by_hh
    households["p_over18"] = households["p_over18"].fillna(0)
    households["p_over18"] = households["p_over18"].astype(int)
    return households


def transform_refiner_events(*args):
    """args are event dataframes with event_id as index"""
    max_event_id = 0
    max_transaction_id = 0
    for df in args:
        df.index += max_event_id
        df.transaction_id += max_transaction_id
        max_event_id += df.index.max() + 10
        max_transaction_id += df.transaction_id.max() + 10

    return pd.concat(args, axis=0)


def transform_persons(p):
    """
  Transform persons table to produce required forecast input
  Return: DataFrame
  """
    # age
    p.loc[:, "age"] = p["AGEP"]
    # member_id
    p.loc[:, "member_id"] = p.per_num

    # race_id
    # expect output for race_id: 1 non-hispanic white, 2. non-hispanic black, 3, hispanic, 4, non-hispanic others
    # RAC1P -> HRACE or race_id, 1, white, 2 black, 3, asian, 4 others
    # HISP 1 non-hisp; 2 and above: hispanic

    p.loc[:, "race_id"] = 4  # other races
    p.loc[p["RAC1P"] == 1, "race_id"] = 1  # white alone
    p.loc[p["RAC1P"] == 2, "race_id"] = 2  # black alone
    p.loc[p["RAC1P"] == 6, "race_id"] = 3  # asian alone

    p.loc[p.race_id > 2, "race_id"] = 4  # other non-hispanic races including asian
    p.loc[p.HISP > 1, "race_id"] = 3  # hispanic

    # relate
    p.loc[:, "relate"] = p[
        "RELP"
    ]  # use PUMS 2018 RELP code. Not RELSHIPP in 2019 and later
    # sex
    p.loc[:, "sex"] = p["SEX"]
    # worker
    p.loc[:, "worker"] = 0
    p.loc[p.ESR.isin([1, 2, 4, 5]), "worker"] = 1  # same as Census /travel model worker
    # school grade, carry-on person attribute for ABM
    # p.loc[:, "school_grade"] = p["SCHG"]
    p.loc[:, "naicsp"] = p["NAICSP"]
    # industry: add for ABM
    naics_18ind = {
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
        "22": 7,
        "51": 8,
        "52": 9,
        "53": 9,
        "54": 10,
        "55": 11,
        "56": 12,
        "61": 13,
        "6S": 14,  # ['621', '623', '624'] -> '6S'
        "62": 15,  # <- '622'
        "71": 16,
        "72": 16,
        "81": 17,
        "92": 18,
        "99": 0,
    }
    p.loc[:, "industry"] = p.loc[:, "naicsp"].str[:3]
    p.loc[:, "industry"].replace(["621", "623", "624"], "6S", inplace=True)
    p.loc[:, "industry"] = p.loc[:, "industry"].str[:2]
    p.loc[:, "industry"].replace(naics_18ind, inplace=True)
    p.loc[p["naicsp"].isnull(), "industry"] = 0
    # set index name to person_id
    p.index.rename("person_id", inplace=True)

    return p[
        [
            "age",
            "household_id",
            "member_id",
            "race_id",
            "relate",
            "sex",
            "worker",
            "industry",
        ]
    ]


def transform_hh(hh, p):
    ## transform hh's HRACE base on hh head race_id from persons table
    hh_head = p.loc[p["relate"] == 0, ["household_id", "race_id"]]
    hh_head_race_id = hh_head.set_index("household_id")["race_id"]
    hh["HRACE"] = hh_head_race_id
    hh_children = p.loc[p.age < 18].groupby("household_id").size()
    hh["children"] = 0
    hh.loc[hh_children.index, "children"] = hh_children
    return hh


def transform_gq_controls(df):
    """ transform 2045 gq controls table """
    df.columns = df.columns.str.lower().str.replace("gq", "").astype(int)
    df.columns.name = "year"
    df = df.unstack().to_frame("count").reset_index()
    return df


def transform_gq(df):
    """ transform 2045 gq pop table """
    df["age"] = df["AGEP"]
    df["race"] = df["Race"]
    df["type"] = df["Type"]
    df["building_id"] = df["BUILDING_I"]
    df["gq_code"] = df["GQ_CODE"]

    df.loc[df.race == "White", "race_id"] = 1
    df.loc[df.race == "Black", "race_id"] = 2
    df.loc[df.race == "Hispanic", "race_id"] = 3
    df.loc[df.race == "Other", "race_id"] = 4
    df.drop(["race"], axis=1, inplace=True)
    df = df[["age", "building_id", "gq_code", "race_id", "type"]]
    return df


def transform_gq_persons(p):
    """
    Transform persons table to produce required forecast input
    Return: DataFrame
    """
    # age
    p.loc[:, "age"] = p["AGEP"]
    p.loc[:, "member_id"] = p.per_num
    # race_id
    # expect output for race_id: 1 non-hispanic white, 2. non-hispanic black, 3, hispanic, 4, non-hispanic others
    # RAC1P -> HRACE or race_id, 1, white, 2 black, 3, asian, 4 others
    # HISP 1 non-hisp; 2 and above: hispanic

    p.loc[:, "race_id"] = 4  # other races
    p.loc[p["RAC1P"] == 1, "race_id"] = 1  # white alone
    p.loc[p["RAC1P"] == 2, "race_id"] = 2  # black alone
    p.loc[p["RAC1P"] == 6, "race_id"] = 3  # asian alone

    p.loc[p.race_id > 2, "race_id"] = 4  # other non-hispanic races including asian
    p.loc[p.HISP > 1, "race_id"] = 3  # hispanic
    p["household_id"] = p.index
    p["relate"] = p["RELSHIPP"] - 21
    p["sex"] = p["SEX"]
    p.loc[p.ESR.isin([1, 2, 4, 5]), "worker"] = 1
    p.loc[:, "naicsp"] = p["NAICSP"]
    # industry: add for ABM
    naics_18ind = {
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
        "22": 7,
        "51": 8,
        "52": 9,
        "53": 9,
        "54": 10,
        "55": 11,
        "56": 12,
        "61": 13,
        "6S": 14,  # ['621', '623', '624'] -> '6S'
        "62": 15,  # <- '622'
        "71": 16,
        "72": 16,
        "81": 17,
        "92": 18,
        "99": 0,
    }
    p.loc[:, "industry"] = p.loc[:, "naicsp"].str[:3]
    p.loc[:, "industry"].replace(["621", "623", "624"], "6S", inplace=True)
    p.loc[:, "industry"] = p.loc[:, "industry"].str[:2]
    p.loc[:, "industry"].replace(naics_18ind, inplace=True)
    p.loc[p["naicsp"].isnull(), "industry"] = 0
    p["industry"] = p["industry"].astype(int)
    # set index name to person_id
    p.index.rename("person_id", inplace=True)
    p = p.fillna(-9)

    return p[
        [
            "age",
            "household_id",
            "building_id",
            "hh_id",
            "gq_code",
            "race_id",
            "member_id",
            "relate",
            "sex",
            "worker",
            "industry",
        ]
    ]


def transform_gq_hh(hh, p):
    # remove extra columns
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
    # assign gq_pop race_id to hh HRACE
    # given same index
    hh.loc[:, "HRACE"] = p["race_id"]
    hh.hht = 0
    hh.index.name = "household_id"
    hh = hh.fillna(-9)
    return hh
