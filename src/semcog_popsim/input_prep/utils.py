import pandas as pd


class CensusDownloader:
    def __init__(self, census_reader, states, counties=None, tract_ids=None, blockgroup_ids=None):
        self.states = states
        self.counties = counties
        self.tracts = tract_ids
        self.blockgroups = blockgroup_ids
        self.cread = census_reader
        self.update_states_counties([self.states, self.counties])

    def state_download(self, variables):
        return self.cread.get(variables, geo={"for": f"state:{self.states}"})

    def county_download(self, variables):
        return self.cread.get(
            variables,
            geo={
                "for": f"county:{self.counties}",
                "in": f"state:{self.states}",
            },
        )

    def tract_download(self, variables):
        return self.cread.get(
            variables,
            geo={
                "for": f"tract:{self.tracts}",
                "in": f"state:{self.states} county:{self.counties}",
            },
        )

    def blockgroup_download(self, variables):
        county_list = self.counties.split(",")
        rows = []
        for county in county_list:
            rows += self.cread.get(
                variables,
                geo={
                    "for": f"block group:{self.blockgroups}",
                    "in": f"state:{self.states} county:{county} tract:{self.tracts}",
                },
            )
        return rows

    def fips_lookup(self, states, counties=None):
        if counties == "*":
            counties = None
        fips_table = pd.read_csv(
            "https://www2.census.gov/geo/docs/reference/codes/files/national_county.txt",
            header=None,
            names=["state", "state_fips", "county_fips", "county", "type"],
            dtype=str,
        )
        query = f"(state in {states})"
        if counties:
            query += f" & (county in {counties})"
            filtered = fips_table.query(query)
            return list(filtered.state_fips.unique()), list(filtered.county_fips.unique())
        filtered = fips_table.query(query)
        return list(filtered.state_fips.unique()), None

    def update_states_counties(self, geographies):
        for i in [0, 1]:
            if geographies[i] is not None:
                if not isinstance(geographies[i], list) and geographies[i] != "*":
                    geographies[i] = [geographies[i]]
                geographies[i] = [str(value) for value in geographies[i]]
        if isinstance(geographies[0], list) and not geographies[0][0].isdigit():
            geographies[0], geographies[1] = self.fips_lookup(geographies[0], geographies[1])
        for i in [0, 1]:
            if geographies[i]:
                geographies[i] = ",".join(str(value).zfill(i + 2) for value in geographies[i])
        self.states = geographies[0]
        self.counties = geographies[1]

    def download(self, variables):
        if not self.counties:
            data = self.state_download(variables)
        elif not self.tracts:
            data = self.county_download(variables)
        elif not self.blockgroups:
            data = self.tract_download(variables)
        else:
            data = self.blockgroup_download(variables)
        return pd.DataFrame.from_dict(data)


def read_tract_puma_crosswalk(crosswalk_key, crosswalk_map, geo_dir):
    df = pd.read_csv(geo_dir / crosswalk_map[crosswalk_key], dtype=str)
    df["COUNTYID"] = df["STATEFP"] + df["COUNTYFP"]
    df["TRACTID"] = df["STATEFP"] + df["COUNTYFP"] + df[crosswalk_key[0]]
    df = df.set_index("TRACTID")
    df.rename(columns={crosswalk_key[1]: "PUMA"}, inplace=True)
    return df[["COUNTYID", "PUMA"]]


def preprocess_pums(households, persons):
    persons = persons.set_index("SERIALNO")
    households["AGEHOH"] = persons.loc[persons.RELP == 0].AGEP
    race_map = {1: 1, 2: 2, 6: 3}
    households["HRACE"] = persons.loc[persons.RELP == 0].RAC1P
    households["HRACE"] = households["HRACE"].map(race_map).fillna(4)
    hisp_map = {1: 0}
    households["HHISP"] = persons.loc[persons.RELP == 0].HISP
    households["HHISP"] = households["HHISP"].map(hisp_map).fillna(1)
    households["income"] = (households.HINCP * households.ADJINC / 1000000).astype(int)
    persons = persons.reset_index()
    households["HWORKERS"] = persons.loc[persons.ESR.isin([1, 2, 4, 5])].groupby("SERIALNO").ESR.size()
    households.fillna(0, inplace=True)
    industry_map = {
        "11": 1, "21": 1, "23": 2, "31": 3, "32": 3, "33": 3, "3M": 3,
        "42": 4, "44": 5, "45": 5, "4M": 5, "48": 6, "49": 6, "22": 6,
        "51": 7, "52": 8, "53": 8, "54": 9, "55": 9, "56": 9, "61": 10,
        "62": 10, "71": 11, "72": 11, "81": 12, "92": 13, "99": 0,
    }
    persons["industry"] = persons.NAICSP.str[:2]
    persons.industry.replace(industry_map, inplace=True)
    persons.loc[persons.NAICSP.str[:6] == "928110", "industry"] = 14
    persons.loc[persons.NAICSP.isnull(), "industry"] = 0
    persons.industry = persons.industry.astype(int)
    persons["pincome"] = persons["PINCP"] * persons["ADJINC"] / 1000000
    persons["pincome"] = persons.loc[~persons["pincome"].isnull(), "pincome"].astype(int)
    return households, persons


def pums_update(df, variable_map):
    for column, conversion in variable_map.items():
        if column in df.columns:
            df[conversion["std_variable"]] = df[column]
            if "std_codes" in conversion:
                df[conversion["std_variable"]].replace(conversion["std_codes"], inplace=True)
    return df


def group_pums_data(households, persons, puma0_col, puma1_col):
    grouped = {}
    empty_df = pd.DataFrame()
    for puma_col in [puma0_col, puma1_col]:
        grouped[puma_col] = {"households": {}, "persons": {}}
        for key, group in households.loc[households[puma_col] != -9].groupby(puma_col):
            grouped[puma_col]["households"][key] = group
        for key, group in persons.loc[persons[puma_col] != -9].groupby(puma_col):
            grouped[puma_col]["persons"][key] = group
        grouped[puma_col]["households"][0] = empty_df
        grouped[puma_col]["persons"][0] = empty_df
    return grouped


def combine_puma_data(puma_list, puma0_col, puma1_col, grouped):
    household_samples = []
    person_samples = []
    count = 0
    for puma in puma_list:
        count += 1
        household_puma = pd.concat(
            [
                grouped[puma0_col]["households"][int(puma[:5])],
                grouped[puma1_col]["households"][int(puma[5:])],
            ]
        )
        household_puma["PUMA"] = puma
        household_puma.index = str(count) + household_puma.index
        household_samples.append(household_puma)

        person_puma = pd.concat(
            [
                grouped[puma0_col]["persons"][int(puma[:5])],
                grouped[puma1_col]["persons"][int(puma[5:])],
            ]
        )
        person_puma["PUMA"] = puma
        person_puma["SERIALNO"] = str(count) + person_puma["SERIALNO"]
        person_samples.append(person_puma)
    return household_samples, person_samples


def marginal_summary(df_margin):
    print("\n * * * verify maringal sums:")
    attribute_list = []
    for column in df_margin.columns:
        cleaned = "".join(char for char in column if not char.isdigit())
        if "ID" not in cleaned:
            attribute_list.append(cleaned)
    for attribute in sorted(set(attribute_list)):
        columns = [column for column in df_margin.columns if column.startswith(attribute)]
        print(attribute, df_margin[columns].sum().sum(), columns)
