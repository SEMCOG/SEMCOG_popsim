# Download ACS 5-year tables by block group for the SEMCOG counties:
#   B25002 (occupancy status)  -> adjust_to_mcd_2025.py, block-group vacancy rates
#   B25007 (tenure by age of householder) -> build_controls_mcd_2025.py, the split of
#          householders aged 65+ into 65-74 and 75+ (B19037 has only 65+)
#
# The Census API key is read from the CENSUS_API_KEY environment variable; it is
# never written to disk.
#
# Usage:
#     CENSUS_API_KEY=... python input_prep/scripts/fetch_acs_b25002_bg.py [acs_year]

import json
import os
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = REPO_ROOT.parent / "d_drive" / "popsim" / "inputs" / "acs_2024_bg"
COUNTIES = ["093", "099", "115", "125", "147", "161", "163"]
# table -> {ACS variable: output column}; columns with the same name are added up
TABLES = {
    "B25002": {"B25002_001E": "acs_hu", "B25002_002E": "acs_occ", "B25002_003E": "acs_vac"},
    "B25007": {"B25007_009E": "acs_hoh_65_74", "B25007_019E": "acs_hoh_65_74",       # owner, renter
               "B25007_010E": "acs_hoh_75_plus", "B25007_020E": "acs_hoh_75_plus",   # 75-84
               "B25007_011E": "acs_hoh_75_plus", "B25007_021E": "acs_hoh_75_plus"},  # 85+
}


def main():
    year = int(sys.argv[1]) if len(sys.argv) > 1 else 2024
    key = os.environ.get("CENSUS_API_KEY")
    if not key:
        sys.exit("set CENSUS_API_KEY")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for table, variables in TABLES.items():
        rows = []
        for county in COUNTIES:
            query = urllib.parse.urlencode({
                "get": ",".join(variables),
                "for": "block group:*",
                "in": f"state:26 county:{county} tract:*",
                "key": key,
            })
            url = f"https://api.census.gov/data/{year}/acs/acs5?{query}"
            with urllib.request.urlopen(url, timeout=120) as resp:
                data = json.load(resp)
            rows += [dict(zip(data[0], r)) for r in data[1:]]

        df = pd.DataFrame(rows)
        df["BLKGRPID"] = df["state"] + df["county"] + df["tract"] + df["block group"]
        vals = df[list(variables)].astype(int)
        vals.columns = list(variables.values())
        out_df = vals.T.groupby(level=0).sum().T  # add up columns with the same name
        out_df.insert(0, "BLKGRPID", df["BLKGRPID"])
        out = OUT_DIR / f"acs{year}_5yr_{table}_bg.csv"
        out_df.to_csv(out, index=False)
        print(f"{table}: {len(out_df)} block groups, totals {out_df.drop(columns='BLKGRPID').sum().to_dict()}")
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
