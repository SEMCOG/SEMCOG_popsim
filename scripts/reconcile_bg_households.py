# Make every block group's synthetic household count equal its HHBASE control.
#
# PopulationSim can miss HHBASE in a few BGs whose controls conflict (in the 2025
# MCD run: 7 BGs, +69 HH region-wide, 3 BGs above their residential units). The
# BG targets add up to the SEMCOG MCD estimates, so exact BG totals keep every
# MCD total exact for placement.
#
# Per BG that is off:
#   - over target: remove the households whose categories are most over their
#     controls (greedy, one at a time, residuals updated after each);
#   - under target: copy households of the same BG (the tract if the BG has fewer
#     than MIN_DONORS) whose categories are most under their controls.
# Persons follow their household. Only off-target BGs change.
#
# Usage:
#     python scripts/reconcile_bg_households.py --pass-dir .../pass2 \
#         --controls .../SEMCOG_2025_control_totals_blkgrp_hhsize_adj.csv

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

MIN_DONORS = 20

# household -> control category, as in configs/controls.csv
CATEGORIES = {
    "HHAGE": lambda h: np.select([h.AGEHOH <= 24, h.AGEHOH <= 44, h.AGEHOH <= 64], ["HHAGE1", "HHAGE2", "HHAGE3"], "HHAGE4"),
    "HHRACE": lambda h: "HHRACE" + h.HRACE.astype(int).astype(str),
    "HHHISP": lambda h: np.where(h.HHISP == 1, "HHHISP2", "HHHISP1"),
    "HHCHD": lambda h: np.where(h.R18 == 1, "HHCHD1", "HHCHD2"),
    "HHINC": lambda h: np.select([h.income <= 30000, h.income <= 60000, h.income <= 100000], ["HHINC1", "HHINC2", "HHINC3"], "HHINC4"),
    "HHCAR": lambda h: np.select([h.VEH == 0, h.VEH == 1], ["HHCAR0", "HHCAR1"], "HHCAR2"),
    "HHPERSONS": lambda h: "HHPERSONS" + h.NP.clip(upper=7).astype(int).astype(str),
    "HHTENURE": lambda h: np.where(h.TEN <= 2, "HHTENURE1", "HHTENURE0"),
}


def categorize(hh):
    return pd.DataFrame({g: f(hh) for g, f in CATEGORIES.items()}, index=hh.index)


def main():
    ap = argparse.ArgumentParser(description="Make each BG household count equal its HHBASE control")
    ap.add_argument("--pass-dir", required=True, type=Path, help="PopulationSim output folder (final pass)")
    ap.add_argument("--controls", required=True, type=Path, help="BG control file used by that pass")
    ap.add_argument("--out-dir", type=Path, help="default: <pass-dir>/../final")
    args = ap.parse_args()
    out_dir = args.out_dir or args.pass_dir.parent / "final"

    hh = pd.read_csv(args.pass_dir / "synthetic_households.csv", dtype={"BLKGRP": str, "TRACT": str})
    per = pd.read_csv(args.pass_dir / "synthetic_persons.csv", dtype={"BLKGRP": str, "TRACT": str})
    ctl = pd.read_csv(args.controls, dtype={"BLKGRPID": str}).set_index("BLKGRPID")

    count = hh.groupby("BLKGRP").size().reindex(ctl.index, fill_value=0)
    diff = count - ctl.HHBASE
    off = diff[diff != 0]
    print("BGs off target: %d | net %+d HH" % (len(off), off.sum()))

    cats = categorize(hh)
    cat_cols = sorted(set(cats.to_numpy().ravel()))
    remove_ids, add_rows, log = [], [], []
    next_id = int(hh.household_id.max()) + 1
    for bg, d in off.items():
        idx = hh.index[hh.BLKGRP == bg]
        # residual = synthesized - control, per category
        resid = cats.loc[idx].apply(pd.Series.value_counts).sum(axis=1).reindex(cat_cols, fill_value=0) \
            - ctl.loc[bg, cat_cols].astype(float)
        before = float(resid.abs().sum())
        if d > 0:
            pool = list(idx)
            for _ in range(int(d)):
                score = [resid[cats.loc[i]].sum() for i in pool]
                i = pool.pop(int(np.argmax(score)))
                remove_ids.append(hh.at[i, "household_id"])
                resid[cats.loc[i]] -= 1
        else:
            donors = idx if len(idx) >= MIN_DONORS else hh.index[hh.TRACT == bg[:11]]
            for _ in range(int(-d)):
                score = [resid[cats.loc[i]].sum() for i in donors]
                i = donors[int(np.argmin(score))]
                row = hh.loc[i].copy()
                row["BLKGRP"], row["household_id"] = bg, next_id
                add_rows.append((hh.at[i, "household_id"], row))
                resid[cats.loc[i]] += 1
                next_id += 1
        log.append({"BLKGRPID": bg, "control": int(ctl.loc[bg, "HHBASE"]), "synthesized": int(count[bg]),
                    "removed": int(max(d, 0)), "added": int(max(-d, 0)),
                    "abs_category_residual_before": before,
                    "abs_category_residual_after": float(resid.abs().sum())})

    new_hh = hh[~hh.household_id.isin(remove_ids)]
    new_per = per[~per.household_id.isin(remove_ids)]
    if add_rows:
        new_hh = pd.concat([new_hh, pd.DataFrame([r for _, r in add_rows])], ignore_index=True)
        copies = []
        for src_id, row in add_rows:
            p = per[per.household_id == src_id].copy()
            p["household_id"], p["BLKGRP"] = row["household_id"], row["BLKGRP"]
            copies.append(p)
        new_per = pd.concat([new_per] + copies, ignore_index=True)

    check = new_hh.groupby("BLKGRP").size().reindex(ctl.index, fill_value=0)
    assert (check == ctl.HHBASE).all(), "BG totals still off"
    assert set(new_per.household_id) == set(new_hh.household_id), "persons without household"

    out_dir.mkdir(parents=True, exist_ok=True)
    new_hh.to_csv(out_dir / "synthetic_households.csv", index=False)
    new_per.to_csv(out_dir / "synthetic_persons.csv", index=False)
    pd.DataFrame(log).to_csv(out_dir / "bg_reconciliation_log.csv", index=False)
    print(pd.DataFrame(log).to_string(index=False))
    print("households %d -> %d (control %d) | persons %d -> %d (control %d)"
          % (len(hh), len(new_hh), ctl.HHBASE.sum(), len(per), len(new_per), ctl.POPBASE.sum()))
    print("wrote", out_dir)


if __name__ == "__main__":
    main()
