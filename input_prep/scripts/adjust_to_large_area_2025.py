# Adjust 2024 ACS PopulationSim control marginals to 2025 large-area household
# population targets.
#
# Method: one uniform scale ratio per large area,
#     ratio[LA] = HHPOP_2025[LA] / sum(current POPBASE over block groups in LA)
# is applied to every target column (block-group population + household controls
# and tract employment controls). A single per-area scalar preserves, exactly,
# every within-block-group identity (each category group already sums to its
# base POPBASE/HHBASE) and makes each large area's POPBASE total match the 2025
# household-population target. Results are integerized with a largest-remainder
# rule that keeps the per-BG category==base identities and the per-LA totals.
#
# Wayne County is split into Detroit (large_area_id 5) and OutWayne (3) via the
# block-group -> large_area crosswalk built from the land-use parcels table.
#
# Usage:
#     python input_prep/scripts/adjust_to_large_area_2025.py
#
# Inputs / outputs are the 2025_synthesis run folder (see paths below).

from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- paths
REPO_ROOT = Path(__file__).resolve().parents[2]
PROJECT_ROOT = REPO_ROOT.parent  # d_drive is a sibling of the SEMCOG_popsim repo
RUN_DIR = PROJECT_ROOT / "d_drive" / "popsim" / "runs" / "2025_synthesis"
DATA_DIR = RUN_DIR / "data"

BLKGRP_IN = DATA_DIR / "SEMCOG_2024_control_totals_blkgrp.csv"
TRACT_IN = DATA_DIR / "SEMCOG_2024_control_totals_tract.csv"
XWALK = DATA_DIR / "blockgroup_large_area_2025.csv"
LA_XLSX = (
    PROJECT_ROOT
    / "d_drive"
    / "forecast_inputs"
    / "group_quarters"
    / "data"
    / "LargeAreaControls_Pop_HHpop_GQ.xlsx"
)
TARGET_YEAR = 2025

BLKGRP_OUT = DATA_DIR / "SEMCOG_2025_control_totals_blkgrp.csv"
TRACT_OUT = DATA_DIR / "SEMCOG_2025_control_totals_tract.csv"

# ---------------------------------------------------------------- control groups
# each population category group sums (exactly) to POPBASE, each household group to HHBASE
POP_GROUPS = {
    "AGEP": ["AGEP1", "AGEP2", "AGEP3", "AGEP4"],
    "RACE": ["RACE1", "RACE2", "RACE3", "RACE4"],
    "SEX": ["SEX1", "SEX2"],
}
HH_GROUPS = {
    "HHAGE": ["HHAGE1", "HHAGE2", "HHAGE3", "HHAGE4"],
    "HHRACE": ["HHRACE1", "HHRACE2", "HHRACE3", "HHRACE4"],
    "HHHISP": ["HHHISP1", "HHHISP2"],
    "HHCHD": ["HHCHD1", "HHCHD2"],
    "HHINC": ["HHINC1", "HHINC2", "HHINC3", "HHINC4"],
    "HHCAR": ["HHCAR0", "HHCAR1", "HHCAR2"],
    "HHPERSONS": ["HHPERSONS%d" % i for i in range(1, 8)],
    "HHTENURE": ["HHTENURE1", "HHTENURE0"],
}
# tract employment controls: independent marginals (no shared base in the file)
TRACT_COLS = (
    ["HHWORKER0", "HHWORKER1", "HHWORKER2"]
    + ["INDUSTRY%d" % i for i in range(1, 15)]
    + ["EMPWORKER"]
)


def largest_remainder(values, target_total):
    """Round non-negative floats to integers summing exactly to target_total,
    preserving proportions via the largest-remainder (Hamilton) method."""
    values = np.asarray(values, dtype=float)
    target_total = int(round(target_total))
    floor = np.floor(values).astype(np.int64)
    diff = target_total - int(floor.sum())
    if diff == 0:
        return floor
    frac = values - np.floor(values)
    if diff > 0:  # hand out extra units to largest remainders
        idx = np.argsort(-frac, kind="stable")[:diff]
        floor[idx] += 1
    else:  # remove units from smallest remainders among cells with headroom
        order = np.argsort(frac, kind="stable")
        removed = 0
        for i in order:
            if removed >= -diff:
                break
            if floor[i] > 0:
                floor[i] -= 1
                removed += 1
    return floor


def integerize_base_by_la(scaled, la, target_by_la):
    """Integerize a scaled base series (indexed like the frame) within each large
    area so each LA sums to target_by_la[LA]."""
    out = pd.Series(0, index=scaled.index, dtype=np.int64)
    for la_id, pos in scaled.groupby(la).groups.items():
        out.loc[pos] = largest_remainder(scaled.loc[pos].values, target_by_la[la_id])
    return out


def allocate_group_to_base(df, cols, base_int):
    """Distribute the integer base across category cols per row, keeping original
    proportions; each row's cols sum exactly to base_int[row]."""
    orig = df[cols].values
    row_tot = orig.sum(axis=1)
    out = np.zeros_like(orig, dtype=np.int64)
    for r in range(orig.shape[0]):
        b = int(base_int.iloc[r])
        if b <= 0:
            continue
        if row_tot[r] > 0:
            shares = orig[r] * (b / row_tot[r])
        else:  # no original detail: put everything in the first category
            shares = np.zeros(len(cols))
            shares[0] = b
        out[r] = largest_remainder(shares, b)
    return pd.DataFrame(out, index=df.index, columns=cols)


def main():
    print("*** adjusting 2024 controls to %d large-area HH population ***\n" % TARGET_YEAR)

    xw = pd.read_csv(XWALK, dtype={"BLKGRPID": str})[["BLKGRPID", "large_area_id"]]
    bg = pd.read_csv(BLKGRP_IN, dtype={"BLKGRPID": str})
    tr = pd.read_csv(TRACT_IN, dtype={"TRACTID": str})

    # 2025 household-population target by large area
    hhpop = pd.read_excel(LA_XLSX, sheet_name="LargeAre_HHPOP").set_index("LargeArea")
    target = hhpop[TARGET_YEAR]  # index = large_area_id (3,5,93,...)

    # ---- block groups ----
    bg = bg.merge(xw, on="BLKGRPID", how="left")
    if bg["large_area_id"].isna().any():
        raise ValueError("block groups missing a large_area_id in crosswalk")

    cur_pop_by_la = bg.groupby("large_area_id")["POPBASE"].sum()
    la_ids = [3, 5, 93, 99, 115, 125, 147, 161]
    ratio = (target.loc[la_ids] / cur_pop_by_la.loc[la_ids]).rename("ratio")
    r_bg = bg["large_area_id"].map(ratio)

    # integer POPBASE (LA total == round(HHPOP_2025)) and HHBASE (LA total == round(scaled))
    pop_target_la = {la: round(target.loc[la]) for la in la_ids}
    hh_scaled = bg["HHBASE"] * r_bg
    hh_target_la = {la: round(hh_scaled[bg.large_area_id == la].sum()) for la in la_ids}

    popbase = integerize_base_by_la(bg["POPBASE"] * r_bg, bg["large_area_id"], pop_target_la)
    hhbase = integerize_base_by_la(hh_scaled, bg["large_area_id"], hh_target_la)

    out = bg[["BLKGRPID"]].copy()
    for grp in HH_GROUPS.values():
        for c, s in allocate_group_to_base(bg, grp, hhbase).items():
            out[c] = s
    out["HHBASE"] = hhbase.values
    out["POPBASE"] = popbase.values
    for grp in POP_GROUPS.values():
        for c, s in allocate_group_to_base(bg, grp, popbase).items():
            out[c] = s
    # restore original column order
    out = out[[c for c in bg.columns if c in out.columns]]

    # ---- tracts (employment) ----
    # tract -> large area by dominant block-group household population
    bg_la = bg[["BLKGRPID", "large_area_id", "HHBASE"]].copy()
    bg_la["TRACTID"] = bg_la["BLKGRPID"].str[:11]
    w = (
        bg_la.groupby(["TRACTID", "large_area_id"])["HHBASE"].sum().reset_index()
        .sort_values("HHBASE").drop_duplicates("TRACTID", keep="last")
        .set_index("TRACTID")["large_area_id"]
    )
    tr = tr.copy()
    tr["large_area_id"] = tr["TRACTID"].map(w)
    if tr["large_area_id"].isna().any():
        raise ValueError("tracts missing a large_area_id")
    r_tr = tr["large_area_id"].map(ratio)

    tout = tr[["TRACTID"]].copy()
    for c in TRACT_COLS:
        scaled = tr[c] * r_tr
        col = pd.Series(0, index=tr.index, dtype=np.int64)
        for la_id, pos in scaled.groupby(tr["large_area_id"]).groups.items():
            col.loc[pos] = largest_remainder(scaled.loc[pos].values, round(scaled.loc[pos].sum()))
        tout[c] = col
    tout = tout[[c for c in tr.columns if c in tout.columns]]

    out.to_csv(BLKGRP_OUT, index=False)
    tout.to_csv(TRACT_OUT, index=False)

    # ---- validation ----
    print("large-area scale ratios (HHPOP_2025 / current POPBASE):")
    print(ratio.round(4).to_string(), "\n")

    chk = out.merge(xw, on="BLKGRPID", how="left")
    post = chk.groupby("large_area_id")["POPBASE"].sum()
    v = pd.DataFrame(
        {"target_HHPOP_2025": target.loc[la_ids].round(), "post_POPBASE": post.loc[la_ids]}
    )
    v["diff"] = v.post_POPBASE - v.target_HHPOP_2025
    print("POPBASE reconciliation by large area:")
    print(v.to_string(), "\n")

    # per-BG identity check
    bad = 0
    for base, groups in (("POPBASE", POP_GROUPS), ("HHBASE", HH_GROUPS)):
        for name, cols in groups.items():
            d = (out[cols].sum(axis=1) - out[base]).abs()
            bad += int((d > 0).sum())
    print("per-BG category==base violations:", bad)
    print("all controls integer:",
          bool((out.select_dtypes("number") % 1 == 0).all().all()
               and (tout.select_dtypes("number") % 1 == 0).all().all()))
    print("\nwrote:\n  %s\n  %s" % (BLKGRP_OUT, TRACT_OUT))


if __name__ == "__main__":
    main()
