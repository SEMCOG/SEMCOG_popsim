# Download ACS 5-year B25002 (occupancy status) by block group for the SEMCOG
# counties. Used by adjust_to_mcd_2025.py for block-group occupancy rates.
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
VARS = {"B25002_001E": "acs_hu", "B25002_002E": "acs_occ", "B25002_003E": "acs_vac"}


def main():
    year = int(sys.argv[1]) if len(sys.argv) > 1 else 2024
    key = os.environ.get("CENSUS_API_KEY")
    if not key:
        sys.exit("set CENSUS_API_KEY")

    rows = []
    for county in COUNTIES:
        query = urllib.parse.urlencode({
            "get": ",".join(VARS),
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
    df = df.rename(columns=VARS)[["BLKGRPID"] + list(VARS.values())]
    df[list(VARS.values())] = df[list(VARS.values())].astype(int)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / f"acs{year}_5yr_B25002_bg.csv"
    df.to_csv(out, index=False)
    print(f"{len(df)} block groups, HU {df.acs_hu.sum():,}, occupied {df.acs_occ.sum():,}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
