# Household-person age targets for the 2025 synthesis, by large area x age band.
#
# The seed (PUMS 2020-24) and the ACS block-group age controls describe about 2022;
# the region aged fast in 2022-25 (baby boomers turning 75), so the base was too
# young (65+ -5.3%, 75+ -7.4% against RSQE 2025). Target per large area and band:
#     RSQE 2025 total population  -  GQ persons (GQ synthesis, by age)
# scaled so each large area adds up to the synthesis household population (POPBASE
# of the block-group targets, which follows the SEMCOG MCD estimate).
#
# Also writes the current synthesis (persons by large area x band) for comparison.
#
# Usage:
#     python input_prep/scripts/age_targets_2025.py --gq <group_quarters.csv> \
#         --final <synthesis final/> --out <targets folder>

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
D = REPO_ROOT.parent / "d_drive"
RSQE = (D / "forecast_inputs/forecast_controls/data/SEMCOG 2055 final/01 Baseline/"
        "RSQE_Final_Baseline_Forecast_for_SEMCOG.xlsx")
GQ_BUILDINGS = D / "forecast_inputs/group_quarters/outputs/GQ_BuildingLevel_2025.csv"
BG_TARGETS = D / "popsim/runs/2025_synthesis_mcd/data/targets/bg_targets_2025.csv"

REGION_LA = {"Rest of Wayne County": 3, "City of Detroit": 5, "Livingston County": 93,
             "Macomb County": 99, "Monroe County": 115, "Oakland County": 125,
             "St. Clair County": 147, "Washtenaw County": 161}
# first digit of a city_id -> large area (Detroit neighborhoods 5xx -> 5)
CITY_DIGIT_LA = {"1": 3, "2": 125, "3": 99, "4": 161, "5": 115, "6": 147, "7": 93}
# controlled bands (BG broad bands AGEP1-4 split further at PUMA) and diagnostics
BANDS = [(0, 17), (18, 24), (25, 34), (35, 54), (55, 64), (65, 74), (75, 200)]
DIAG = [(65, 200), (75, 84), (85, 200)]
# persons at the ages of the head-age control bands HHAGE1-3 (HHAGE4-5 = 65_74, 75_plus);
# the builder scales the head-age controls by target / current persons (constant headship)
HEAD = [(15, 24), (25, 44), (45, 64)]
YEAR = 2025


def band_name(lo, hi):
    return "%d_%s" % (lo, "plus" if hi >= 200 else hi)


def city_to_la(city_id):
    c = int(city_id)
    return 5 if 500 <= c < 600 else CITY_DIGIT_LA[str(c)[0]]


def bg_large_area(path=BG_TARGETS):
    """large_area_id per block group; a block group split between areas has none in the
    targets file, so it takes its county's large area (not possible in Wayne: 2 areas)."""
    t = pd.read_csv(path, usecols=["BLKGRPID", "POPBASE", "large_area_id"], dtype={"BLKGRPID": str})
    county_la = {"093": 93, "099": 99, "115": 115, "125": 125, "147": 147, "161": 161}
    miss = t.large_area_id.isna() & (t.POPBASE > 0)
    t.loc[miss, "large_area_id"] = t.loc[miss, "BLKGRPID"].str[2:5].map(county_la)
    assert t.loc[t.POPBASE > 0, "large_area_id"].notna().all(), "populated BG without large area"
    return t.set_index("BLKGRPID")


def by_band(ages, counts, la, bands):
    """Sum counts by large area x band (inclusive age bounds)."""
    out = {}
    for lo, hi in bands:
        m = (ages >= lo) & (ages <= hi)
        out[band_name(lo, hi)] = pd.Series(np.where(m, counts, 0)).groupby(np.asarray(la)).sum()
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gq", required=True, help="GQ persons (UrbanSim group_quarters.csv)")
    ap.add_argument("--final", required=True, help="synthesis final/ folder (current run)")
    ap.add_argument("--out", required=True, help="output folder")
    a = ap.parse_args()
    out_dir = Path(a.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    r = pd.read_excel(RSQE, sheet_name="1-Pop by Age, Gender, Race")
    r = r[(r.Race == "All Races") & (r.Gender == "Total") & r.Region.isin(REGION_LA)
          & (r.Age != "All Ages (0-100)")]
    r_age = r.Age.str.extract(r"(\d+)")[0].astype(int).to_numpy()
    r_la = r.Region.map(REGION_LA).to_numpy()

    gq = pd.read_csv(a.gq, usecols=["age", "building_id"])
    bld = pd.read_csv(GQ_BUILDINGS, usecols=["building_id", "city_id"]).drop_duplicates()
    assert not bld.building_id.duplicated().any(), "GQ building in two cities"
    gq_city = gq.building_id.map(bld.set_index("building_id").city_id)
    assert gq_city.notna().all(), "GQ building without city_id"
    gq_la = gq_city.map(city_to_la).to_numpy()

    bgs = bg_large_area()
    pop_la = bgs.groupby("large_area_id").POPBASE.sum()
    pop_la.index = pop_la.index.astype(int)

    bands = BANDS + DIAG + HEAD
    rsqe = by_band(r_age, r[YEAR].to_numpy(), r_la, bands)
    gqb = by_band(gq.age.to_numpy(), np.ones(len(gq)), gq_la, bands).reindex(rsqe.index).fillna(0)
    hh_raw = (rsqe - gqb).clip(lower=0)
    ctl = [band_name(*b) for b in BANDS]
    scale = pop_la.reindex(hh_raw.index) / hh_raw[ctl].sum(axis=1)
    target = hh_raw.mul(scale, axis=0)

    fin = pd.read_csv(Path(a.final) / "synthetic_persons.csv", usecols=["AGEP", "BLKGRP"], dtype={"BLKGRP": str})
    cur = by_band(fin.AGEP.to_numpy(), np.ones(len(fin)), fin.BLKGRP.map(bgs.large_area_id).astype(int).to_numpy(), bands)

    rows = []
    for la in target.index:
        for b in target.columns:
            rows.append({"large_area_id": la, "band": b, "controlled": b in ctl,
                         "rsqe_total": rsqe.at[la, b], "gq": gqb.at[la, b],
                         "target_hh": target.at[la, b], "current_hh": cur.at[la, b]})
    res = pd.DataFrame(rows)
    res["current_vs_target_pct"] = (res.current_hh / res.target_hh - 1) * 100
    res.to_csv(out_dir / "age_targets_2025.csv", index=False)
    print("scale to POPBASE per large area: %s" % scale.round(4).to_dict())
    reg = res.groupby("band", sort=False)[["rsqe_total", "gq", "target_hh", "current_hh"]].sum()
    reg["current_vs_target_pct"] = (reg.current_hh / reg.target_hh - 1) * 100
    print(reg.round(1).to_string())


if __name__ == "__main__":
    main()
