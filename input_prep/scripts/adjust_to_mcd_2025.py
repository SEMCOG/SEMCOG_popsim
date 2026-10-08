# Build 2025 block-group household (HHBASE) and household-population (POPBASE)
# targets from the SEMCOG MCD household estimate.
#
# Method (option B, October 2026):
#   - Household counts: SEMCOG MCD estimate, exact. Detroit is split into its 55
#     neighborhoods (city_id 501-555) by their 2025 HH shares, scaled to the
#     Detroit total.
#   - Household population: MCD (or neighborhood) POP - GQ, then one scalar per
#     large area so each large area equals RSQE 2025 total population - GQ. This
#     keeps the base year consistent with the RSQE-minus-GQ forecast controls.
#   - Area -> block group: through BG x area pieces of the base-year building
#     stock. The area's vacant units (units - HH) are spread over its pieces by
#     weight = units x ACS BG vacancy rate, and piece HH = units - vacant. All
#     vacancy rates in an area scale by one factor, so a piece is full only when
#     its ACS vacancy is 0. (Allocating HH by occupancy instead pushed ~60% of BGs
#     to 100% occupancy, because the estimate needs 95.5% occupancy outside
#     Detroit vs 93.0% from ACS.) Vacancy rate source: ACS B25002 (--vacancy-source
#     b25002) or 1 - ACS HH / BG units (acs_hh, no Census key needed).
#   - Piece HH pop = 1 person per HH, plus the area's persons above 1 per HH spread
#     by HH x (ACS BG average size - 1); ACS size capped at SIZE_CAP. Size >= 1.
#
# Attribute controls are built in a later step from these targets.
#
# Usage:
#     python input_prep/scripts/adjust_to_mcd_2025.py \
#         --estimate .../July1_2025_Population_revised.xlsx \
#         --base-hdf .../main_100226.h5

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- paths
REPO_ROOT = Path(__file__).resolve().parents[2]
PROJECT_ROOT = REPO_ROOT.parent  # d_drive is a sibling of the SEMCOG_popsim repo
D_DRIVE = PROJECT_ROOT / "d_drive"
RUN_DIR = D_DRIVE / "popsim" / "runs" / "2025_synthesis_mcd"
OUT_DIR = RUN_DIR / "data" / "targets"

ACS_CONTROLS = D_DRIVE / "popsim" / "runs" / "2024_synthesis" / "data" / "SEMCOG_2024_control_totals_blkgrp.csv"
GEO_XWALK = D_DRIVE / "popsim" / "runs" / "2024_synthesis" / "data" / "SEMCOG_2024_geo_cross_walk.csv"
ACS_B25002 = D_DRIVE / "popsim" / "inputs" / "acs_2024_bg" / "acs2024_5yr_B25002_bg.csv"
CITY_TABLE = D_DRIVE / "forecast_inputs" / "group_quarters" / "data" / "POP_GQ_HH_HU_by_CityID_2025.xlsx"
RSQE_WORKBOOK = (D_DRIVE / "forecast_inputs" / "forecast_controls" / "data" / "SEMCOG 2055 final"
                 / "01 Baseline" / "RSQE_Final_Baseline_Forecast_for_SEMCOG.xlsx")

YEAR = 2025
DETROIT = 5
DETROIT_NBHDS = range(501, 556)
# MCD code prefix -> large_area_id (the estimate's county summary codes differ)
PREFIX_TO_LA = {1: 3, 2: 125, 3: 99, 4: 161, 5: 115, 6: 147, 7: 93}
RSQE_TO_LA = {
    "Rest of Wayne County": 3, "City of Detroit": 5, "Livingston County": 93,
    "Macomb County": 99, "Monroe County": 115, "Oakland County": 125,
    "St. Clair County": 147, "Washtenaw County": 161,
}
SIZE_CAP = 6.0  # ACS BG average size used as a weight is capped here
SIZE_FLAG_MAX = 6.0  # BG average size above this is hard for the 7+ size bin (mean <= 10)


# ---------------------------------------------------------------- helpers
def largest_remainder(values, total, cap=None):
    """Round non-negative floats to integers summing to `total`, preserving
    proportions; with `cap`, no element ends above its cap."""
    values = np.asarray(values, dtype=float)
    out = np.floor(values).astype(np.int64)
    diff = int(round(total)) - int(out.sum())
    frac = values - out
    if diff > 0:
        room = np.ones(len(values), bool) if cap is None else out < np.asarray(cap)
        idx = [i for i in np.argsort(-frac, kind="stable") if room[i]][:diff]
        if len(idx) < diff:
            raise ValueError("no room to place rounding remainder under the cap")
        out[idx] += 1
    elif diff < 0:
        idx = [i for i in np.argsort(frac, kind="stable") if out[i] > 0][:-diff]
        out[idx] -= 1
    return out


def allocate_capped(total, weights, cap):
    """Split `total` in proportion to `weights` with x <= cap; the excess over a
    cap goes to the uncapped elements, again by weight."""
    w = np.asarray(weights, dtype=float)
    cap = np.asarray(cap, dtype=float)
    if total > cap.sum() + 1e-9:
        raise ValueError(f"total {total} exceeds capacity {cap.sum()}")
    x = np.zeros(len(w))
    free = cap > 0
    remaining = float(total)
    while remaining > 1e-9:
        wf = np.where(free, w, 0.0)
        if wf.sum() <= 0:  # no weight left: fall back to remaining capacity
            wf = np.where(free, cap - x, 0.0)
        add = remaining * wf / wf.sum()
        over = x + add > cap
        if not over.any():
            x += add
            break
        remaining -= (cap[over] - x[over]).sum()
        x[over] = cap[over]
        free &= ~over
    return x


# ---------------------------------------------------------------- inputs
def load_estimate(path):
    est = pd.read_excel(path)
    est.columns = ["mcd", "name", "pop", "hh"]
    # keep the 233 MCD rows only: drop region/county rows and the (total) rows
    # of split cities (8005-8020), which repeat their parts
    est = est[(est.mcd == DETROIT) | est.mcd.between(1000, 7999)].copy()
    est["large_area_id"] = np.where(est.mcd == DETROIT, DETROIT, (est.mcd // 1000).map(PREFIX_TO_LA))
    return est


def load_rsqe_hhpop(gq_by_la):
    pop = pd.read_excel(RSQE_WORKBOOK, sheet_name="1-Pop by Age, Gender, Race")
    tot = pop[(pop.Race == "All Races") & (pop.Gender == "Total") & (pop.Age == "All Ages (0-100)")]
    tot = tot[tot.Region.isin(RSQE_TO_LA)]
    tot = tot.set_index(tot.Region.map(RSQE_TO_LA))[YEAR]
    return tot - gq_by_la.reindex(tot.index)


def build_areas(est, city):
    """One row per allocation area: each non-Detroit MCD, and Detroit's
    neighborhoods. Columns: area_id, mcd, large_area_id, hh, pop, gq."""
    gq = city.set_index("CITYID")["GQ2025"]
    mcd = est[est.mcd != DETROIT].rename(columns={"mcd": "area_id"})
    mcd["mcd"] = mcd.area_id
    mcd["gq"] = mcd.area_id.map(gq).fillna(0)

    det = est.loc[est.mcd == DETROIT].iloc[0]
    nb = city[city.CITYID.isin(DETROIT_NBHDS)].copy()
    nb_hh = largest_remainder(nb.HH2025 * det.hh / nb.HH2025.sum(), det.hh)
    # neighborhood POP - GQ shares carry the Detroit POP - GQ
    det_gq = nb.GQ2025.sum()
    nb_hhpop_share = (nb.POP2025 - nb.GQ2025) / (nb.POP2025 - nb.GQ2025).sum()
    nbhd = pd.DataFrame({
        "area_id": nb.CITYID.values, "name": nb.NAME.values, "mcd": DETROIT,
        "large_area_id": DETROIT, "hh": nb_hh, "gq": nb.GQ2025.values,
    })
    nbhd["pop"] = nb_hhpop_share.values * (det["pop"] - det_gq) + nbhd.gq
    print("Detroit neighborhoods: pre-scale HH %d vs estimate %d (shortfall %+d, %.2f%%)"
          % (nb.HH2025.sum(), det.hh, det.hh - nb.HH2025.sum(),
             100 * (det.hh / nb.HH2025.sum() - 1)))

    cols = ["area_id", "name", "mcd", "large_area_id", "hh", "pop", "gq"]
    return pd.concat([mcd[cols], nbhd[cols]], ignore_index=True)


def load_pieces(base_hdf):
    """Residential units by block group x allocation area."""
    with pd.HDFStore(str(base_hdf), "r") as store:
        parcels = store["parcels"][["census_bg_id", "county_id", "semmcd", "city_id"]]
        bldg = store["buildings"][["parcel_id", "residential_units"]]
    b = bldg.join(parcels, on="parcel_id")
    b = b[b.residential_units > 0]
    b["BLKGRPID"] = ("26" + b.county_id.astype(str).str.zfill(3)
                     + b.census_bg_id.astype(str).str.zfill(7))
    b["area_id"] = np.where(b.semmcd == DETROIT, b.city_id, b.semmcd)
    return b.groupby(["BLKGRPID", "area_id"]).residential_units.sum().rename("units").reset_index()


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description="2025 BG HH / HH-pop targets from the SEMCOG MCD estimate")
    ap.add_argument("--estimate", required=True, type=Path, help="SEMCOG MCD estimate workbook")
    ap.add_argument("--base-hdf", required=True, type=Path, help="base-year HDF (buildings, parcels)")
    ap.add_argument("--vacancy-source", choices=["b25002", "acs_hh"], default="b25002",
                    help="BG vacancy rate: ACS B25002 vacant/HU, or 1 - ACS HH / BG units")
    ap.add_argument("--vacancy-shrink", type=float, default=100.0,
                    help="b25002 only: pseudo-units K pulling BG vacancy toward the tract rate, and tract toward county")
    args = ap.parse_args()

    est = load_estimate(args.estimate)
    city = pd.read_excel(CITY_TABLE)
    city = city[city.CITYID != DETROIT]  # Detroit total row; neighborhoods carry it
    acs = pd.read_csv(ACS_CONTROLS, dtype={"BLKGRPID": str})[["BLKGRPID", "HHBASE", "POPBASE"]]

    # ---- areas and household population (option B)
    areas = build_areas(est, city)
    areas["hhpop0"] = areas["pop"] - areas["gq"]
    if (areas.hhpop0 < 0).any():
        raise ValueError("areas with GQ > POP:\n%s" % areas[areas.hhpop0 < 0])
    gq_by_la = areas.groupby("large_area_id").gq.sum()
    rsqe_hhpop = load_rsqe_hhpop(gq_by_la)
    la_scale = rsqe_hhpop / areas.groupby("large_area_id").hhpop0.sum()
    areas["hhpop"] = 0
    for la, idx in areas.groupby("large_area_id").groups.items():
        areas.loc[idx, "hhpop"] = largest_remainder(
            areas.loc[idx, "hhpop0"] * la_scale[la], round(rsqe_hhpop[la]))
    if (areas.hhpop < areas.hh).any():
        raise ValueError("areas with HH pop < HH:\n%s" % areas[areas.hhpop < areas.hh])

    # ---- pieces and weights
    pieces = load_pieces(args.base_hdf)
    # keep the PopulationSim geography; buildings in other BG codes get no HH
    xw = pd.read_csv(GEO_XWALK, dtype={"BLKGRPID": str}).BLKGRPID
    outside = pieces[~pieces.BLKGRPID.isin(xw)]
    if len(outside):
        print("excluded %d BG codes outside the PopulationSim geography (%d units): %s"
              % (outside.BLKGRPID.nunique(), outside.units.sum(), sorted(outside.BLKGRPID.unique())))
    pieces = pieces[pieces.BLKGRPID.isin(xw)].reset_index(drop=True)
    pieces = pieces.merge(acs, on="BLKGRPID", how="left")
    if args.vacancy_source == "b25002":
        b25002 = pd.read_csv(ACS_B25002, dtype={"BLKGRPID": str})
        # BG vacancy is noisy (39% of BGs and 12% of tracts report 0 vacant
        # units): shrink each BG rate toward its tract rate, and each tract rate
        # toward its county rate, with K pseudo-units
        k = args.vacancy_shrink
        csum = b25002.groupby(b25002.BLKGRPID.str[:5])[["acs_vac", "acs_hu"]].transform("sum")
        county_rate = csum.acs_vac / csum.acs_hu
        tsum = b25002.groupby(b25002.BLKGRPID.str[:11])[["acs_vac", "acs_hu"]].transform("sum")
        tract_rate = (tsum.acs_vac + k * county_rate) / (tsum.acs_hu + k)
        b25002["vac_shrunk"] = (b25002.acs_vac + k * tract_rate) / (b25002.acs_hu + k).where(b25002.acs_hu + k > 0)
        pieces = pieces.merge(b25002[["BLKGRPID", "acs_hu", "acs_vac", "vac_shrunk"]], on="BLKGRPID", how="left")
        vac = pieces.vac_shrunk
    else:
        bg_units = pieces.groupby("BLKGRPID").units.transform("sum")
        vac = (1 - pieces.HHBASE / bg_units).clip(lower=0)
    # pieces with no ACS rate: the area's unit-weighted mean vacancy rate
    area_vac = (vac * pieces.units).groupby(pieces.area_id).sum() / \
        pieces.units.where(vac.notna()).groupby(pieces.area_id).sum()
    pieces["vac_rate"] = vac.fillna(pieces.area_id.map(area_vac)).fillna(0).clip(0, 1)
    pieces["weight"] = pieces.units * pieces.vac_rate

    unknown = set(areas.area_id[areas.hh > 0]) - set(pieces.area_id)
    if unknown:
        raise ValueError(f"estimate areas with no residential units: {sorted(unknown)}")
    orphan = pieces[~pieces.area_id.isin(areas.area_id)]
    if len(orphan):
        print("units in areas absent from the estimate (target 0 HH):",
              orphan.groupby("area_id").units.sum().to_dict())

    # ---- area -> piece allocation
    pieces["hh"] = 0
    pieces["hhpop"] = 0
    acs_size = (pieces.POPBASE / pieces.HHBASE.where(pieces.HHBASE > 0))
    for a in areas.itertuples():
        idx = pieces.index[pieces.area_id == a.area_id]
        units = pieces.loc[idx, "units"].values
        if a.hh > units.sum():
            raise ValueError(f"area {a.area_id}: {a.hh} HH > {units.sum()} units")
        vacant = allocate_capped(units.sum() - a.hh, pieces.loc[idx, "weight"].values, units)
        hh = units - largest_remainder(vacant, units.sum() - a.hh, cap=units)
        pieces.loc[idx, "hh"] = hh
        if a.hh == 0:
            continue
        # one person per household, then the persons above one per household
        # by HH x (ACS size - 1); keeps every piece at size >= 1
        size = acs_size.loc[idx].fillna(a.hhpop / a.hh).clip(upper=SIZE_CAP).values
        extra_w = hh * (size - 1)
        if extra_w.sum() <= 0:
            extra_w = hh.astype(float)
        extra = a.hhpop - hh.sum()
        pieces.loc[idx, "hhpop"] = hh + largest_remainder(extra * extra_w / extra_w.sum(), extra)

    # ---- block-group targets
    bg = pieces.groupby("BLKGRPID").agg(HHBASE=("hh", "sum"), POPBASE=("hhpop", "sum"),
                                         units=("units", "sum"), n_areas=("area_id", "nunique"))
    bg = bg.join(acs.set_index("BLKGRPID").add_prefix("acs_"), how="outer")
    bg = bg.reindex(xw).fillna(0).rename_axis("BLKGRPID")
    bg = bg.join(pieces.groupby("BLKGRPID").area_id.first().map(
        areas.set_index("area_id").large_area_id).rename("large_area_id"))
    bg[["HHBASE", "POPBASE", "units", "n_areas"]] = bg[["HHBASE", "POPBASE", "units", "n_areas"]].astype(int)
    size = bg.POPBASE / bg.HHBASE.where(bg.HHBASE > 0)
    bg["flag"] = ""
    bg.loc[(bg.HHBASE == 0) & (bg.POPBASE > 0), "flag"] = "pop_without_hh"
    bg.loc[(bg.HHBASE > 0) & (bg.POPBASE < bg.HHBASE), "flag"] = "pop_below_hh"
    bg.loc[size.notna() & (size > SIZE_FLAG_MAX) & (bg.flag == ""), "flag"] = "size_above_max"
    bg.loc[bg.HHBASE > bg.units, "flag"] = "hh_above_units"

    # ---- checks
    by_area = pieces.groupby("area_id")[["hh", "hhpop"]].sum()
    chk = areas.set_index("area_id")[["hh", "hhpop"]].sub(by_area, fill_value=0).abs().to_numpy().max()
    assert chk == 0, "area totals not preserved"
    assert (pieces.hh <= pieces.units).all(), "HH above units in a piece"
    assert bg.HHBASE.sum() == est.hh.sum(), "region HH does not match the estimate"

    # ---- outputs
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    bg.reset_index().to_csv(OUT_DIR / "bg_targets_2025.csv", index=False)
    pieces.to_csv(OUT_DIR / "bg_area_pieces_2025.csv", index=False)
    areas.to_csv(OUT_DIR / "area_targets_2025.csv", index=False)

    la = areas.groupby("large_area_id")[["hh", "pop", "gq", "hhpop0", "hhpop"]].sum()
    la["scale"] = la_scale.round(4)
    la["avg_size"] = (la.hhpop / la.hh).round(3)
    la.loc["region"] = la.sum()
    la.loc["region", ["scale", "avg_size"]] = [np.nan, round(la.loc["region", "hhpop"] / la.loc["region", "hh"], 3)]
    print("\nvacancy source:", args.vacancy_source)
    print(la.to_string())
    print("\nblock groups:", len(bg), "| with HH:", int((bg.HHBASE > 0).sum()),
          "| split across areas:", int((bg.n_areas > 1).sum()))
    print("flags:", bg.flag.value_counts().drop("", errors="ignore").to_dict())
    occupied = bg[bg.units > 0]
    print("region HH / units: %.3f; BGs with units: %d, at 100%% occupancy: %d"
          % (bg.HHBASE.sum() / bg.units.sum(), len(occupied), int((occupied.HHBASE == occupied.units).sum())))
    print("wrote", OUT_DIR)


if __name__ == "__main__":
    main()
