# Build the 2025 PopulationSim run package (controls + configs) on the SEMCOG MCD
# targets from adjust_to_mcd_2025.py (step 2 of the October 2026 re-synthesis).
#
# Method:
#   - Totals: HHBASE / POPBASE per block group from the step-1 targets.
#   - Category controls: ACS 2020-2024 share x new base, integerized so every
#     group sums exactly to its base (one ratio per BG, the "hard update").
#     BG shares are shrunk toward the tract share, (count + K x tract share) /
#     (BG base + K), with K = sigma2/tau2 estimated per control group from the
#     data (ACS sampling noise vs real between-BG variation). A BG with no ACS
#     data gets the tract share; a tract with too little data, the PUMA share.
#   - Age of head: five bands; the 65+ band of B19037 is split into 65-74 and 75+
#     with ACS B25007 (input from fetch_acs_b25002_bg.py).
#   - Household size: NOT one ratio. The ACS size shares are tilted (shape-
#     preserving balancer) until implied persons = POPBASE; 7+ households count
#     at the PUMS mean size of 7+ households in the PUMA. BGs whose target mean
#     size is outside [1, MAX_TOP_BIN] or whose refit leaves a person gap are
#     flagged, not forced.
#   - Tract controls: worker groups = ACS tract shares x new tract HH; industry and
#     EMPWORKER = ACS tract rate per person x new tract HH population. Tracts on
#     PUMA rates derive EMPWORKER from the worker groups and industry from
#     EMPWORKER, so the three cannot contradict.
#   - Importance: totals hard; size high; workers lowered; income raised.
#
# Usage:
#     python input_prep/scripts/build_controls_mcd_2025.py

import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hh_size_balancer import rebalance_household_sizes_shape_preserving  # noqa: E402

# ---------------------------------------------------------------- paths
REPO_ROOT = Path(__file__).resolve().parents[2]
RUNS = REPO_ROOT.parent / "d_drive" / "popsim" / "runs"
RUN_DIR = RUNS / "2025_synthesis_mcd"
TARGETS = RUN_DIR / "data" / "targets" / "bg_targets_2025.csv"
ACS_DIR = RUNS / "2024_synthesis" / "data"          # unadjusted ACS 2020-2024 controls
ACS_B25007 = REPO_ROOT.parent / "d_drive" / "popsim" / "inputs" / "acs_2024_bg" / "acs2024_5yr_B25007_bg.csv"
TEMPLATE = RUNS / "2025_synthesis"                  # seed, crosswalk, configs

BLKGRP_OUT = RUN_DIR / "data" / "SEMCOG_2025_control_totals_blkgrp.csv"
TRACT_OUT = RUN_DIR / "data" / "SEMCOG_2025_control_totals_tract.csv"
DIAG_OUT = RUN_DIR / "data" / "targets" / "controls_build_diagnostics.csv"
COPY_DATA = ["SEMCOG_2024_seed_households.csv", "SEMCOG_2024_seed_persons.csv",
             "SEMCOG_2024_geo_cross_walk.csv"]

# ---------------------------------------------------------------- controls
HH_GROUPS = {
    "HHAGE": ["HHAGE1", "HHAGE2", "HHAGE3", "HHAGE4", "HHAGE5"],  # 4 = 65-74, 5 = 75+
    "HHRACE": ["HHRACE1", "HHRACE2", "HHRACE3", "HHRACE4"],
    "HHHISP": ["HHHISP1", "HHHISP2"],
    "HHCHD": ["HHCHD1", "HHCHD2"],
    "HHINC": ["HHINC1", "HHINC2", "HHINC3", "HHINC4"],
    "HHCAR": ["HHCAR0", "HHCAR1", "HHCAR2"],
    "HHTENURE": ["HHTENURE1", "HHTENURE0"],
}
SIZE_COLS = ["HHPERSONS%d" % i for i in range(1, 8)]
POP_GROUPS = {
    "AGEP": ["AGEP1", "AGEP2", "AGEP3", "AGEP4"],
    "RACE": ["RACE1", "RACE2", "RACE3", "RACE4"],
    "SEX": ["SEX1", "SEX2"],
}
WORKER_COLS = ["HHWORKER0", "HHWORKER1", "HHWORKER2"]
INDUSTRY_COLS = ["INDUSTRY%d" % i for i in range(1, 15)]

MIN_HH = 20      # tract with fewer ACS households than this: PUMA shares as the prior
MIN_POP = 50     # same, for person groups
MAX_K = 5000.0   # cap on the shrinkage strength (no real BG variation -> tract share)
MAX_TOP_BIN = 10.0  # the 7+ bin mean the pass-2 balancer allows (7-10)

# importance changes vs the July 2025 run (target name -> new importance). Age of
# head and children were lowered in run 1 to give the size refit room; the refit
# turned out mild, so they are back at the July values (5,000 / 10,000).
IMPORTANCE = {
    **{t: 1000 for t in ["hh_workers_0", "hh_workers_1", "hh_workers_2"]},
    **{t: 5000 for t in ["hh_inc_30", "hh_inc_30_60", "hh_inc_60_100", "hh_inc_100_plus"]},
}


# ---------------------------------------------------------------- helpers
def largest_remainder(values, total):
    """Round non-negative floats to integers summing to `total`."""
    values = np.asarray(values, dtype=float)
    out = np.floor(values).astype(np.int64)
    diff = int(round(total)) - int(out.sum())
    order = np.argsort(-(values - out), kind="stable")
    out[order[:diff]] += 1
    return out


def estimate_k(a, cols, base):
    """Shrinkage strength K = sigma2 / tau2 for one control group.

    For each BG, the squared gap between its share and the share of the OTHER BGs
    in its tract is tau2 (real between-BG variation) + sigma2 / n (sampling
    noise). A least-squares fit of the gaps on 1/n, pooled over the group's
    categories, separates the two. Returns (K, tau2, sigma2)."""
    b = a[(a[base] >= 10)].copy()
    tn = b.groupby("TRACT")[base].transform("sum")
    b = b[tn > b[base]]
    tn = b.groupby("TRACT")[base].transform("sum")
    n = b[base]
    d2, inv_n = [], []
    for c in cols:
        te = b.groupby("TRACT")[c].transform("sum")
        d2.append(((b[c] / n) - (te - b[c]) / (tn - n)) ** 2)
        inv_n.append(1 / n)
    y, x = np.concatenate(d2), np.concatenate(inv_n)
    (tau2, sigma2), *_ = np.linalg.lstsq(np.column_stack([np.ones_like(x), x]), y, rcond=None)
    k = sigma2 / tau2 if tau2 > 0 else MAX_K
    return float(np.clip(k, 0, MAX_K)), float(tau2), float(sigma2)


def share_table(acs, xw, cols, base, min_n):
    """Per-BG shares of `cols`, shrunk toward the tract share:
        share = (BG count + K x tract share) / (BG base + K)
    The tract share falls back to the PUMA share when the tract base < min_n.
    A BG with no ACS data gets the tract share."""
    a = acs.set_index("BLKGRPID")[cols + [base]].join(xw.set_index("BLKGRPID")[["TRACT", "PUMA"]])
    tract = a.groupby("TRACT")[cols + [base]].transform("sum")
    puma = a.groupby("PUMA")[cols + [base]].transform("sum")
    prior = tract[cols].where(tract[base] >= min_n, puma[cols])
    prior = prior.div(prior.sum(axis=1), axis=0)
    k, tau2, sigma2 = estimate_k(a, cols, base)
    shrunk = (a[cols] + prior.mul(k)).div(a[base] + k, axis=0)
    level = pd.Series("bg", index=a.index).where(a[base] > 0, "tract")
    level[(a[base] == 0) & (tract[base] < min_n)] = "puma"
    return shrunk.div(shrunk.sum(axis=1), axis=0), level, k


def allocate_group(shares, base):
    out = np.zeros(shares.shape, dtype=np.int64)
    for r, b in enumerate(base):
        if b > 0:
            out[r] = largest_remainder(shares[r] * b, b)
    return out


def workers_2plus_mean_by_puma():
    """Weighted PUMS mean number of workers in households with 2+ workers, by PUMA."""
    hh = pd.read_csv(TEMPLATE / "data" / COPY_DATA[0], usecols=["PUMA", "HWORKERS", "WGTP"])
    w2 = hh[hh.HWORKERS >= 2]
    m = (w2.HWORKERS * w2.WGTP).groupby(w2.PUMA).sum() / w2.WGTP.groupby(w2.PUMA).sum()
    return m.reindex(hh.PUMA.unique()).fillna(np.average(w2.HWORKERS, weights=w2.WGTP))


def split_head_age_65(acs, xw):
    """Split ACS age-of-head band HHAGE4 (65+, B19037) into HHAGE4 = 65-74 and
    HHAGE5 = 75+ with the B25007 householder ratio of the block group (the tract,
    then the PUMA, where the BG has no householders aged 65+).

    Without this split the 65-74 vs 75+ mix came from the PUMS seed and was biased
    toward 75+ in every large area (75+/65+ 0.43-0.50 vs ACS 0.39-0.43), which
    inflated the household forecast, because heads aged 75+ grow fastest."""
    b = pd.read_csv(ACS_B25007, dtype={"BLKGRPID": str}).set_index("BLKGRPID")
    b = b.join(xw.set_index("BLKGRPID")[["TRACT", "PUMA"]])
    young, old = b.acs_hoh_65_74, b.acs_hoh_75_plus
    ratio = young / (young + old)
    for level in ("TRACT", "PUMA"):
        y, o = young.groupby(b[level]).transform("sum"), old.groupby(b[level]).transform("sum")
        ratio = ratio.fillna(y / (y + o))
    ratio = ratio.fillna(young.sum() / (young.sum() + old.sum()))
    out = acs.copy()
    r = out.BLKGRPID.map(ratio).fillna(young.sum() / (young.sum() + old.sum()))
    total65 = out["HHAGE4"]
    out["HHAGE4"] = total65 * r
    out.insert(out.columns.get_loc("HHAGE4") + 1, "HHAGE5", total65 * (1 - r))
    print("head age 65+ split from B25007: region 75+/65+ = %.3f" % (out.HHAGE5.sum() / total65.sum()))
    return out


def split_head_age_control_rows(ctl):
    """Replace the control hh_age_65_plus by hh_age_65_74 (HHAGE4) and hh_age_75_plus (HHAGE5)."""
    i = ctl.index[ctl.target == "hh_age_65_plus"]
    assert len(i) == 1, "template controls changed: hh_age_65_plus not found"
    row = ctl.loc[i[0]]
    rows = [row.copy(), row.copy()]
    rows[0][["target", "control_field", "expression"]] = [
        "hh_age_65_74", "HHAGE4", "(households.AGEHOH > 64) & (households.AGEHOH <= 74)"]
    rows[1][["target", "control_field", "expression"]] = [
        "hh_age_75_plus", "HHAGE5", "(households.AGEHOH > 74)"]
    return pd.concat([ctl.loc[:i[0] - 1], pd.DataFrame(rows), ctl.loc[i[0] + 1:]], ignore_index=True)


def top_bin_mean_by_puma():
    """Weighted PUMS mean size of 7+ person households, by PUMA."""
    hh = pd.read_csv(TEMPLATE / "data" / COPY_DATA[0], usecols=["PUMA", "NP", "WGTP"])
    big = hh[hh.NP >= 7]
    m = (big.NP * big.WGTP).groupby(big.PUMA).sum() / big.WGTP.groupby(big.PUMA).sum()
    return m.reindex(hh.PUMA.unique()).fillna(np.average(big.NP, weights=big.WGTP))


# ---------------------------------------------------------------- main
def main():
    tgt = pd.read_csv(TARGETS, dtype={"BLKGRPID": str})[["BLKGRPID", "HHBASE", "POPBASE"]]
    acs = pd.read_csv(ACS_DIR / "SEMCOG_2024_control_totals_blkgrp.csv", dtype={"BLKGRPID": str})
    acs_tr = pd.read_csv(ACS_DIR / "SEMCOG_2024_control_totals_tract.csv", dtype={"TRACTID": str})
    xw = pd.read_csv(ACS_DIR / "SEMCOG_2024_geo_cross_walk.csv", dtype={"BLKGRPID": str, "TRACTID": str})
    xw = xw.rename(columns={"TRACTID": "TRACT"})
    xw["PUMA"] = xw.PUMA.astype(int)

    acs = split_head_age_65(acs, xw)
    acs = acs.set_index("BLKGRPID").reindex(tgt.BLKGRPID).fillna(0).reset_index()
    out = tgt.set_index("BLKGRPID")
    hh, pop = out.HHBASE.to_numpy(), out.POPBASE.to_numpy()

    # ---- household and person category groups
    diag = pd.DataFrame(index=out.index)
    k_by_group = {}
    for name, cols in HH_GROUPS.items():
        shares, level, k_by_group[name] = share_table(acs, xw, cols, "HHBASE", MIN_HH)
        out[cols] = allocate_group(shares.loc[out.index].to_numpy(), hh)
    diag["hh_share_level"] = level.loc[out.index]
    for name, cols in POP_GROUPS.items():
        shares, level, k_by_group[name] = share_table(acs, xw, cols, "POPBASE", MIN_POP)
        out[cols] = allocate_group(shares.loc[out.index].to_numpy(), pop)
    diag["pop_share_level"] = level.loc[out.index]

    # ---- household size: tilt to POPBASE, flag what cannot be met
    shares, _, k_by_group["HHPERSONS"] = share_table(acs, xw, SIZE_COLS, "HHBASE", MIN_HH)
    shares = shares.loc[out.index].to_numpy()
    top_mean = xw.set_index("BLKGRPID").PUMA.map(top_bin_mean_by_puma()).loc[out.index].to_numpy()
    size_counts = np.zeros(shares.shape, dtype=np.int64)
    gap = np.zeros(len(out))
    flag = np.array([""] * len(out), dtype=object)
    for r in range(len(out)):
        if hh[r] == 0:
            continue
        weights = np.array([1, 2, 3, 4, 5, 6, top_mean[r]])
        x0 = largest_remainder(shares[r] * hh[r], hh[r])
        target_mean = pop[r] / hh[r]
        if not 1 <= target_mean <= MAX_TOP_BIN:
            flag[r] = "infeasible_mean"
            size_counts[r] = x0
            continue
        x = rebalance_household_sizes_shape_preserving(x0, pop[r], weights)
        size_counts[r] = x
        gap[r] = float(np.dot(x, weights)) - pop[r]
        if abs(gap[r]) > max(1.0, top_mean[r] / 2):
            flag[r] = "person_gap"
    out[SIZE_COLS] = size_counts
    diag["acs_mean_size"] = (acs.set_index("BLKGRPID").POPBASE / acs.set_index("BLKGRPID").HHBASE
                             .where(lambda s: s > 0)).loc[out.index].round(3)
    diag["target_mean_size"] = (out.POPBASE / out.HHBASE.where(out.HHBASE > 0)).round(3)
    diag["top_bin_mean"] = top_mean.round(2)
    diag["size_person_gap"] = gap.round(2)
    diag["size_flag"] = flag

    # ---- tract controls
    tr_new = out.groupby(out.index.str[:11])[["HHBASE", "POPBASE"]].sum()
    tr_acs = acs.groupby(acs.BLKGRPID.str[:11])[["HHBASE", "POPBASE"]].sum()
    tr = acs_tr.set_index("TRACTID").reindex(tr_new.index).fillna(0)
    tr = tr.join(tr_acs.add_prefix("acs_")).join(xw.drop_duplicates("TRACT").set_index("TRACT").PUMA)
    puma = tr.groupby("PUMA")[WORKER_COLS + INDUSTRY_COLS + ["EMPWORKER", "acs_HHBASE", "acs_POPBASE"]].transform("sum")
    use_puma_hh = tr.acs_HHBASE < MIN_HH
    use_puma_pop = tr.acs_POPBASE < MIN_POP
    w_src = tr[WORKER_COLS].where(~use_puma_hh, puma[WORKER_COLS])
    w_shares = w_src.div(w_src.sum(axis=1), axis=0).to_numpy()
    rate_src = tr[INDUSTRY_COLS + ["EMPWORKER"]].where(~use_puma_pop, puma[INDUSTRY_COLS + ["EMPWORKER"]])
    rate_pop = tr.acs_POPBASE.where(~use_puma_pop, puma.acs_POPBASE)
    rates = rate_src.div(rate_pop, axis=0)

    tout = pd.DataFrame(index=tr.index)
    tout[WORKER_COLS] = allocate_group(w_shares, tr_new.HHBASE.to_numpy())
    ind = rates[INDUSTRY_COLS].mul(tr_new.POPBASE, axis=0)
    tout[INDUSTRY_COLS] = np.vstack([largest_remainder(v, round(v.sum())) for v in ind.to_numpy()])
    tout["EMPWORKER"] = (rates.EMPWORKER * tr_new.POPBASE).round().astype(int)
    # Tracts on PUMA rates: rounding the three groups separately can contradict
    # (e.g. a 0-worker household with EMPWORKER 1), which makes the integerizer
    # infeasible. Derive EMPWORKER from the worker households, then industry from
    # EMPWORKER, so one real household can meet all of them.
    m2 = tr.PUMA.map(workers_2plus_mean_by_puma())
    emp = (tout.HHWORKER1 + m2 * tout.HHWORKER2).round().astype(int)
    ind_shares = rate_src[INDUSTRY_COLS].div(rate_src[INDUSTRY_COLS].sum(axis=1), axis=0).fillna(0)
    for t in tr.index[use_puma_hh | use_puma_pop]:
        tout.loc[t, "EMPWORKER"] = emp[t]
        tout.loc[t, INDUSTRY_COLS] = largest_remainder(ind_shares.loc[t].to_numpy() * emp[t], emp[t])
    print("tracts: %d | worker shares from PUMA: %d | industry rates from PUMA: %d"
          % (len(tr), int(use_puma_hh[tr_new.HHBASE > 0].sum()), int(use_puma_pop[tr_new.POPBASE > 0].sum())))

    # ---- checks
    for cols in HH_GROUPS.values():
        assert (out[cols].sum(axis=1) == out.HHBASE).all()
    assert (out[SIZE_COLS].sum(axis=1) == out.HHBASE).all()
    for cols in POP_GROUPS.values():
        assert (out[cols].sum(axis=1) == out.POPBASE).all()
    assert (tout[WORKER_COLS].sum(axis=1) == tr_new.HHBASE).all()
    puma_tr = tr.index[use_puma_hh | use_puma_pop]
    assert (tout.loc[puma_tr, INDUSTRY_COLS].sum(axis=1) == tout.loc[puma_tr, "EMPWORKER"]).all()
    assert (out.select_dtypes("number") >= 0).all().all() and (tout >= 0).all().all()

    # ---- write the run package
    (RUN_DIR / "configs").mkdir(parents=True, exist_ok=True)
    (RUN_DIR / "output").mkdir(exist_ok=True)
    col_order = [c for c in acs.columns if c != "BLKGRPID"] + ["BLKGRPID"]
    out.reset_index()[col_order].to_csv(BLKGRP_OUT, index=False)
    tcols = [c for c in acs_tr.columns if c != "TRACTID"]
    tout[tcols].rename_axis("TRACTID").reset_index()[tcols + ["TRACTID"]].to_csv(TRACT_OUT, index=False)
    diag.to_csv(DIAG_OUT)
    for f in COPY_DATA:
        shutil.copy2(TEMPLATE / "data" / f, RUN_DIR / "data" / f)
    # OR-Tools integerizer: ~50-90x faster than CVXPY/GLPK with the same quality, and
    # GLPK has no working time limit. Needs the isolated popsim conda env (in the
    # base env, ortools + conda pyarrow segfault).
    settings = (TEMPLATE / "configs" / "settings.yaml").read_text()
    assert "USE_CVXPY: true\n" in settings, "template settings changed: check USE_CVXPY"
    settings = settings.replace("USE_CVXPY: true\n", "USE_CVXPY: false  # OR-Tools; run in the popsim conda env\n")
    (RUN_DIR / "configs" / "settings.yaml").write_text(settings)
    ctl = split_head_age_control_rows(pd.read_csv(TEMPLATE / "configs" / "controls.csv"))
    unknown = set(IMPORTANCE) - set(ctl.target)
    assert not unknown, unknown
    ctl["importance"] = ctl.target.map(IMPORTANCE).fillna(ctl.importance).astype(int)
    ctl.to_csv(RUN_DIR / "configs" / "controls.csv", index=False)

    # ---- report
    n = (out.HHBASE > 0).sum()
    print("BGs with HH: %d | HH share level: %s | pop share level: %s"
          % (n, diag.hh_share_level[out.HHBASE > 0].value_counts().to_dict(),
             diag.pop_share_level[out.POPBASE > 0].value_counts().to_dict()))
    print("shrinkage K by group:", {g: round(k) for g, k in k_by_group.items()})
    print("size refit flags:", diag.size_flag[diag.size_flag != ""].value_counts().to_dict(),
          "| max |person gap| %.2f" % np.abs(gap).max())
    d = diag[out.HHBASE > 0]
    print("mean size ACS -> target: median %.3f -> %.3f" % (d.acs_mean_size.median(), d.target_mean_size.median()))
    reg = out[SIZE_COLS].sum() / out.HHBASE.sum()
    acs_reg = acs[SIZE_COLS].sum() / acs.HHBASE.sum()
    print("region size shares ACS -> new:\n%s" % pd.DataFrame({"acs": acs_reg, "new": reg}).round(4).T.to_string())
    print("wrote", RUN_DIR)


if __name__ == "__main__":
    main()
