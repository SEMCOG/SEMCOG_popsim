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
#   - Industry controls count employed persons only (ESR 1-2; armed forces ESR 4-5),
#     as the ACS tables do; the seed copy codes the armed forces as industry 14.
#   - Worker controls count workers at work (ESR 1, 4), as B08014 / B08202 do; the
#     seed copy adds HWORKERS_AT_WORK for the household worker groups.
#   - Seed reweighting (SEED_REWEIGHT): household weights calibrated per PUMA to the employed and
#     unemployed rates by age of the recent survey years (the five-year PUMS includes the COVID labor
#     market); tract worker controls rebuilt from the reweighted seed, each with its own definition.
#   - Optional PUMA controls (USE_PUMA_PROFILE): person age inside the BG bands 25-64 and 65+, and
#     workers by age, at the PUMS seed shares of the PUMA, scaled to the BG / tract
#     totals of the PUMA. Without them the balancing moved workers to ages 20-34 and
#     persons from 65-69 to 70-74 (review 2026-10-13).
#   - Importance: totals hard; size high; workers lowered; income raised.
#
# Usage:
#     python input_prep/scripts/build_controls_mcd_2025.py

import os
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
# SYNTH_RUN_DIR: a pilot run folder (its data/targets/ needs bg_targets_2025.csv)
RUN_DIR = Path(os.environ.get("SYNTH_RUN_DIR") or RUNS / "2025_synthesis_mcd")
TARGETS = RUN_DIR / "data" / "targets" / "bg_targets_2025.csv"
ACS_DIR = RUNS / "2024_synthesis" / "data"          # unadjusted ACS 2020-2024 controls
ACS_B25007 = REPO_ROOT.parent / "d_drive" / "popsim" / "inputs" / "acs_2024_bg" / "acs2024_5yr_B25007_bg.csv"
TEMPLATE = RUNS / "2025_synthesis_07_08"            # July run: seed, crosswalk, configs

BLKGRP_OUT = RUN_DIR / "data" / "SEMCOG_2025_control_totals_blkgrp.csv"
TRACT_OUT = RUN_DIR / "data" / "SEMCOG_2025_control_totals_tract.csv"
PUMA_OUT = RUN_DIR / "data" / "SEMCOG_2025_control_totals_puma.csv"
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
# PUMA controls: BG age band -> {column: (min age, max age)}; workers by age
PUMA_AGE = {
    "AGEP3": {"PAGE25_34": (25, 34), "PAGE35_54": (35, 54), "PAGE55_64": (55, 64)},
    "AGEP4": {"PAGE65_69": (65, 69), "PAGE70_74": (70, 74), "PAGE75P": (75, 200)},
}
PUMA_WORKER = {"PWORK16_24": (16, 24), "PWORK25_54": (25, 54), "PWORK55_64": (55, 64), "PWORK65P": (65, 200)}
PUMA_IMPORTANCE = 1000   # moderate: below the age-of-head, income and size controls
USE_PUMA_PROFILE = False  # off in run 5: test the industry fix alone first
# Seed reweighting to recent-year employment: the five-year PUMS includes the COVID labor market of
# 2020-21. Household weights are calibrated per PUMA so that employed and unemployed persons by age
# band match the RECENT_YEARS rates, holding persons by age, households by size and by head age; the
# tract worker controls are then rebuilt from the reweighted seed (reweight_seed, rebuild_worker_controls).
SEED_REWEIGHT = True
RECENT_YEARS = ["2022", "2023", "2024"]
RW_AGE = [(0, 15), (16, 24), (25, 54), (55, 64), (65, 200)]
RW_EMP = [(16, 24), (25, 54), (55, 64), (65, 200)]
RW_HEAD = [(15, 24), (25, 44), (45, 64), (65, 74), (75, 200)]

# Age controls aged to 2025 (run 9, pilot A2 of 2026-10-21). The seed and the ACS age shares
# describe about 2022, and the region aged fast in 2022-25. Targets by large area x band from
# age_targets_2025.py (RSQE 2025 - GQ by age; needs data/targets/age_targets_2025.csv). The BG
# person bands are raked to them per large area; finer bands go to PUMA controls built from the
# BG pieces (PUMAs cross large areas), with 75+ split into 75-84 and 85+ (one 75+ band kept the
# seed's 2020-like split: 85+ +22%); head-age bands are scaled by target / current persons at the
# same ages. Pilots A (one 75+ band) and B (A2 + seed reweighted to the 2025 ages: same ages,
# 8,400 fewer employed) were not adopted.
AGE_2025 = True
AGE_TARGETS = RUN_DIR / "data" / "targets" / "age_targets_2025.csv"
AGE_BROAD = {"AGEP1": ["0_17"], "AGEP2": ["18_24"], "AGEP3": ["25_34", "35_54", "55_64"],
             "AGEP4": ["65_74", "75_84", "85_plus"]}
AGE_FINE = {"25_34": (25, 34), "35_54": (35, 54), "55_64": (55, 64), "65_74": (65, 74),
            "75_84": (75, 84), "85_plus": (85, 200)}
HEAD_BAND = {"HHAGE1": "15_24", "HHAGE2": "25_44", "HHAGE3": "45_64", "HHAGE4": "65_74", "HHAGE5": "75_plus"}

MIN_HH = 20      # tract with fewer ACS households than this: PUMA shares as the prior
MIN_POP = 50     # same, for person groups
MAX_K = 5000.0   # cap on the shrinkage strength (no real BG variation -> tract share)
MAX_TOP_BIN = 10.0  # the 7+ bin mean the pass-2 balancer allows (7-10)
MAX_IND_PER_WORKER = 1.3  # ACS industry / EMPWORKER above this: group-quarters tract (region p99 = 1.24)
MIN_HH_GQ = 100  # group-quarters tract with fewer ACS households: PUMA worker-household shares

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


def at_work_by_household():
    """Number of persons at work (ESR 1 or 4) per seed household id."""
    p = pd.read_csv(TEMPLATE / "data" / COPY_DATA[1], usecols=["hh_id", "ESR"])
    return p[p.ESR.isin([1, 4])].groupby("hh_id").size()


def workers_2plus_mean_by_puma():
    """Weighted PUMS mean number of workers at work in households with 2+ of them, by PUMA."""
    hh = pd.read_csv(TEMPLATE / "data" / COPY_DATA[0], usecols=["PUMA", "hh_id", "WGTP"])
    hh["w"] = hh.hh_id.map(at_work_by_household()).fillna(0)
    w2 = hh[hh.w >= 2]
    m = (w2.w * w2.WGTP).groupby(w2.PUMA).sum() / w2.WGTP.groupby(w2.PUMA).sum()
    return m.reindex(hh.PUMA.unique()).fillna(np.average(w2.w, weights=w2.WGTP))


def employed_per_at_work_by_puma():
    """Weighted PUMS ratio of employed (ESR 1, 2, 4, 5) to at-work (ESR 1, 4) persons, by PUMA:
    industry (employed) = EMPWORKER (at work) x this ratio in the fallback tracts."""
    p = pd.read_csv(TEMPLATE / "data" / COPY_DATA[1], usecols=["PUMA", "ESR", "PWGTP"])
    emp = (p.PWGTP * p.ESR.isin([1, 2, 4, 5])).groupby(p.PUMA).sum()
    atw = (p.PWGTP * p.ESR.isin([1, 4])).groupby(p.PUMA).sum()
    return emp / atw


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


def puma_profile_controls(out, tout, xw):
    """PUMA controls for person age inside the BG bands AGEP3 (25-64) and AGEP4 (65+),
    and for workers by age. Shares are the PUMS seed shares of the PUMA (household
    persons, PWGTP); totals are the BG / tract control sums of the PUMA, so the new
    controls agree with the existing ones."""
    p = pd.read_csv(TEMPLATE / "data" / COPY_DATA[1], usecols=["PUMA", "AGEP", "ESR", "PWGTP"])
    worker = p.ESR.isin([1, 2, 4, 5])
    bg_puma = xw.set_index("BLKGRPID").PUMA.loc[out.index]
    tr_puma = xw.drop_duplicates("TRACT").set_index("TRACT").PUMA.loc[tout.index]
    res = pd.DataFrame(index=pd.Index(sorted(bg_puma.unique()), name="PUMA"))

    def fill(masks, total):
        seed = pd.DataFrame({c: p.PWGTP.where(m, 0).groupby(p.PUMA).sum() for c, m in masks.items()})
        seed = seed.reindex(res.index).fillna(0)
        res[list(masks)] = allocate_group(seed.div(seed.sum(axis=1), axis=0).to_numpy(),
                                          total.reindex(res.index).to_numpy())

    for broad, bands in PUMA_AGE.items():
        fill({c: p.AGEP.between(lo, hi) for c, (lo, hi) in bands.items()}, out[broad].groupby(bg_puma).sum())
    fill({c: worker & p.AGEP.between(lo, hi) for c, (lo, hi) in PUMA_WORKER.items()},
         tout.EMPWORKER.groupby(tr_puma).sum())

    for broad, bands in PUMA_AGE.items():
        assert (res[list(bands)].sum(axis=1) == out[broad].groupby(bg_puma).sum()).all()
    assert (res[list(PUMA_WORKER)].sum(axis=1) == tout.EMPWORKER.groupby(tr_puma).sum()).all()
    assert (res.PWORK25_54 <= res.PAGE25_34 + res.PAGE35_54).all()
    assert (res.PWORK55_64 <= res.PAGE55_64).all()
    assert (res.PWORK65P <= res.PAGE65_69 + res.PAGE70_74 + res.PAGE75P).all()
    return res


def rake(prior, rows, cols, iters=500, tol=1e-9):
    """Scale a non-negative matrix to the row and column totals (cols rescaled to sum(rows))."""
    m = np.asarray(prior, dtype=float).copy()
    rows = np.asarray(rows, dtype=float)
    cols = np.asarray(cols, dtype=float) * rows.sum() / max(float(np.sum(cols)), 1e-12)
    for _ in range(iters):
        m *= np.divide(rows, m.sum(axis=1), out=np.zeros_like(rows), where=m.sum(axis=1) > 0)[:, None]
        m *= np.divide(cols, m.sum(axis=0), out=np.zeros_like(cols), where=m.sum(axis=0) > 0)[None, :]
        if np.abs(m.sum(axis=1) - rows).max() < tol * max(rows.max(), 1):
            break
    return m


def age_controls_2025(out, bg_la, xw):
    """Rake the BG person-age and head-age controls to the 2025 targets per large area, and
    build the PUMA controls for the finer bands from the BG pieces. Changes `out` in place;
    returns the PUMA control table."""
    t = pd.read_csv(AGE_TARGETS).set_index(["large_area_id", "band"])
    tgt, cur = t.target_hh.unstack(), t.current_hh.unstack()
    p = pd.read_csv(TEMPLATE / "data" / COPY_DATA[1], usecols=["PUMA", "AGEP", "PWGTP"])
    bg_puma = xw.set_index("BLKGRPID").PUMA.loc[out.index]
    agep, hhage = POP_GROUPS["AGEP"], HH_GROUPS["HHAGE"]
    fine = pd.DataFrame(0.0, index=out.index, columns=list(AGE_FINE))
    for la, idx in out.groupby(bg_la.loc[out.index]).groups.items():
        sub = out.loc[idx]
        # person bands: rake to the large area's targets, BG totals fixed
        cols = [tgt.loc[la, AGE_BROAD[c]].sum() for c in agep]
        m = rake(sub[agep], sub.POPBASE, cols)
        out.loc[idx, agep] = allocate_group(m / np.maximum(m.sum(axis=1, keepdims=True), 1e-12), sub.POPBASE.to_numpy())
        # head age: scale by target / current persons at the same ages (constant headship)
        f = np.array([tgt.at[la, HEAD_BAND[c]] / cur.at[la, HEAD_BAND[c]] for c in hhage])
        prior = sub[hhage].to_numpy() * f
        m = rake(prior, sub.HHBASE, prior.sum(axis=0))
        out.loc[idx, hhage] = allocate_group(m / np.maximum(m.sum(axis=1, keepdims=True), 1e-12), sub.HHBASE.to_numpy())
        # finer bands: PUMA seed shares inside the BG broad band, raked to the large area
        for broad in ("AGEP3", "AGEP4"):
            bands = AGE_BROAD[broad]
            seed = pd.DataFrame({b: p.PWGTP.where(p.AGEP.between(*AGE_FINE[b]), 0).groupby(p.PUMA).sum() for b in bands})
            shares = seed.div(seed.sum(axis=1), axis=0).reindex(bg_puma.loc[idx]).to_numpy()
            prior = shares * out.loc[idx, broad].to_numpy()[:, None]
            fine.loc[idx, bands] = rake(prior, out.loc[idx, broad], tgt.loc[la, bands].to_numpy())
    res = fine.groupby(bg_puma).sum()
    pout = pd.DataFrame(index=res.index.rename("PUMA"))
    for broad in ("AGEP3", "AGEP4"):
        bands = AGE_BROAD[broad]
        total = out[broad].groupby(bg_puma).sum().loc[res.index].to_numpy()
        pout[["PAGE" + b.upper() for b in bands]] = allocate_group(
            res[bands].div(res[bands].sum(axis=1), axis=0).fillna(0).to_numpy(), total)
        assert (pout[["PAGE" + b.upper() for b in bands]].sum(axis=1).to_numpy() == total).all()
    return pout


def age_puma_control_rows():
    """controls.csv rows for the 2025 PUMA age controls."""
    rows = []
    for b, (lo, hi) in AGE_FINE.items():
        age = "(persons.AGEP >= %d)" % lo + ("" if hi >= 200 else " & (persons.AGEP <= %d)" % hi)
        rows.append({"target": "persons_age_" + b, "geography": "PUMA", "seed_table": "persons",
                     "importance": PUMA_IMPORTANCE, "control_field": "PAGE" + b.upper(), "expression": age})
    return pd.DataFrame(rows)


def puma_control_rows():
    """controls.csv rows for the PUMA profile controls."""
    rows = []
    bands = {c: r for b in PUMA_AGE.values() for c, r in b.items()}
    for c, (lo, hi) in {**bands, **PUMA_WORKER}.items():
        age = "(persons.AGEP >= %d)" % lo + ("" if hi >= 200 else " & (persons.AGEP <= %d)" % hi)
        work = c.startswith("PWORK")
        rows.append({
            "target": ("persons_worker_" + c[5:] if work else "persons_age_" + c[4:]).lower(),
            "geography": "PUMA", "seed_table": "persons", "importance": PUMA_IMPORTANCE,
            "control_field": c,
            "expression": age + (" & persons.ESR.isin([1,2,4,5])" if work else "")})
    return pd.DataFrame(rows)


def fix_industry_rows(ctl):
    """Count only employed persons in the industry controls, as the ACS tables do.

    PUMS gives an industry code to everyone who worked in the last 5 years, so the
    template expression (persons.industry == k) also counted ~581,000 persons who are
    not employed. To meet the employed-only controls, the balancing removed 81% of
    them (recent retirees, teens between jobs) and moved the age profile.
    Industries 1-13 (C24030): civilian employed, ESR 1-2. Industry 14 (B08126,
    armed forces): ESR 4-5."""
    ind = ctl.target.str.match(r"Persons_industry_\d+$")
    assert ind.sum() == 14, "template controls changed: industry rows"
    ctl = ctl.copy()
    ctl.loc[ind, "expression"] = ctl.loc[ind, "expression"].str.strip() + " & persons.ESR.isin([1,2])"
    ctl.loc[ctl.target == "Persons_industry_14", "expression"] = "persons.ESR.isin([4,5])"
    return ctl


def fix_worker_rows(ctl):
    """Count workers AT WORK in the reference week (ESR 1, 4), as the ACS sources of
    EMPWORKER (B08014) and HHWORKER (B08202) do. The template counted all employed
    persons (ESR 1, 2, 4, 5); the controls were ~2% lower, so the balancing removed
    ~50,000 employed persons and every worker rate came out ~1 point low."""
    ctl = ctl.copy()
    w = ctl.target == "persons_worker"
    hw = ctl.target.str.match(r"hh_workers_\d$")
    assert w.sum() == 1 and hw.sum() == 3, "template controls changed: worker rows"
    ctl.loc[w, "expression"] = "persons.ESR.isin([1,4])"
    ctl.loc[hw, "expression"] = ctl.loc[hw, "expression"].str.replace(
        "households.HWORKERS ", "households.HWORKERS_AT_WORK ", regex=False)
    assert ctl.loc[hw, "expression"].str.contains("HWORKERS_AT_WORK").all()
    return ctl


def write_seed_households(dst):
    """Copy the seed households with HWORKERS_AT_WORK = persons with ESR 1 or 4 (used
    by the hh_workers controls). HWORKERS (all employed) is kept for the output."""
    h = pd.read_csv(TEMPLATE / "data" / COPY_DATA[0], low_memory=False)
    h["HWORKERS_AT_WORK"] = h.hh_id.map(at_work_by_household()).fillna(0).astype(int)
    assert (h.HWORKERS_AT_WORK <= h.HWORKERS).all()
    h.to_csv(dst, index=False)
    print("seed households: HWORKERS_AT_WORK added (differs from HWORKERS in %d records)"
          % (h.HWORKERS_AT_WORK != h.HWORKERS).sum())


def write_seed_persons(dst):
    """Copy the seed persons with industry = 14 for the armed forces (ESR 4-5); the
    seed codes them 13, because its NAICSP test missed the 92811P* codes."""
    p = pd.read_csv(TEMPLATE / "data" / COPY_DATA[1], low_memory=False)
    af = p.ESR.isin([4, 5])
    p.loc[af, "industry"] = 14
    p.to_csv(dst, index=False)
    print("seed persons: industry -> 14 for %d armed-forces records" % af.sum())


def _rw_incidence(h, p):
    """Household x constraint matrix (counts) for the seed reweighting."""
    cols = {}
    for lo, hi in RW_AGE:
        cols["persons_%d_%d" % (lo, hi)] = p.AGEP.between(lo, hi)
    for lo, hi in RW_EMP:
        cols["employed_%d_%d" % (lo, hi)] = p.AGEP.between(lo, hi) & p.ESR.isin([1, 2, 4, 5])
        cols["unemployed_%d_%d" % (lo, hi)] = p.AGEP.between(lo, hi) & (p.ESR == 3)
    A = pd.DataFrame(cols).groupby(p.hh_id).sum().reindex(h.hh_id).fillna(0)
    for k in range(1, 6):
        A["size_%d" % k] = (h.NP.clip(upper=5) == k).to_numpy()
    for lo, hi in RW_HEAD:
        A["head_%d_%d" % (lo, hi)] = h.AGEHOH.between(lo, hi).to_numpy()
    return A.astype(float)


def _calibrate(w0, A, T, tol=1e-10, max_sweeps=5000):
    """Minimum-entropy calibration w = w0 * exp(A @ lam) by cyclic dual coordinate ascent:
    for each constraint k a scalar Newton step on sum_h w_h exp(d a_hk) a_hk = T_k."""
    w = w0.copy()
    for sweep in range(max_sweeps):
        for k in range(A.shape[1]):
            a = A[:, k]
            nz = a > 0
            d = 0.0
            for _ in range(50):
                e = w[nz] * np.exp(d * a[nz])
                f, g = (e * a[nz]).sum() - T[k], (e * a[nz] ** 2).sum()
                step = f / g
                d -= step
                if abs(step) < 1e-14:
                    break
            w[nz] *= np.exp(d * a[nz])
        err = np.max(np.abs(A.T @ w - T) / np.maximum(T, 1))
        if err < tol:
            return w, sweep + 1, err
    return w, max_sweeps, err


def reweight_seed(data_dir, diag_out):
    """Calibrate the seed household weight WGTP per PUMA to the RECENT_YEARS employed and unemployed rates
    by age band (household weights, as PopulationSim counts persons). Rewrites the seed copies in data_dir;
    the pooled weight is kept as WGTP_pooled. Record attributes do not change."""
    h = pd.read_csv(data_dir / COPY_DATA[0], low_memory=False)
    p = pd.read_csv(data_dir / COPY_DATA[1], low_memory=False)
    year = p.groupby("hh_id").SERIALNO.first().astype(str).str[:4]
    h["survey_year"] = h.hh_id.map(year)
    assert h.survey_year.notna().all()
    p["hw"] = p.hh_id.map(h.set_index("hh_id").WGTP)
    p["yr"] = p.hh_id.map(year)
    mult = pd.Series(1.0, index=h.hh_id)
    rows = []
    for puma, hp in h.groupby("PUMA"):
        pp = p[p.hh_id.isin(hp.hh_id)]
        A = _rw_incidence(hp, pp)
        w0 = hp.WGTP.to_numpy(float)
        T = A.T.to_numpy() @ w0
        rec = pp.yr.isin(RECENT_YEARS)
        for lo, hi in RW_EMP:
            band = pp.AGEP.between(lo, hi)
            n = T[A.columns.get_loc("persons_%d_%d" % (lo, hi))]
            rate = (pp.hw * (band & pp.ESR.isin([1, 2, 4, 5]) & rec)).sum() / (pp.hw * (band & rec)).sum()
            T[A.columns.get_loc("employed_%d_%d" % (lo, hi))] = rate * n
            u = (pp.hw * (band & (pp.ESR == 3) & rec)).sum() / (pp.hw * (band & rec)).sum()
            T[A.columns.get_loc("unemployed_%d_%d" % (lo, hi))] = u * n
        w, sweeps, err = _calibrate(w0, A.to_numpy(), T)
        mult.loc[hp.hh_id] = w / w0
        rows.append({"PUMA": puma, "households": len(hp), "sweeps": sweeps, "max_rel_error": err,
                     "mult_min": (w / w0).min(), "mult_max": (w / w0).max()})
    diag = pd.DataFrame(rows)
    assert (diag.max_rel_error < 1e-8).all(), "seed reweighting did not converge"
    diag.to_csv(diag_out, index=False)
    h["WGTP_pooled"] = h.WGTP
    h["WGTP"] = h.WGTP * h.hh_id.map(mult)
    p["PWGTP"] = p.PWGTP * p.hh_id.map(mult)
    h.drop(columns=["survey_year"]).to_csv(data_dir / COPY_DATA[0], index=False)
    p.drop(columns=["hw", "yr"]).to_csv(data_dir / COPY_DATA[1], index=False)
    m = h.hh_id.map(mult)
    print("seed reweighted to %s employment: multiplier p1 %.3f, median %.3f, p99 %.3f"
          % ("-".join(RECENT_YEARS[::len(RECENT_YEARS) - 1]), m.quantile(.01), m.median(), m.quantile(.99)))


def _rw_puma_stats(h, p, wcol):
    """Per PUMA: at-work and civilian employed persons, at-work worker-group shares, m2, employed/at-work."""
    hw = p.hh_id.map(h.set_index("hh_id")[wcol])
    g = pd.DataFrame({
        "atw": hw * p.ESR.isin([1, 4]), "civ": hw * p.ESR.isin([1, 2]),
        "emp_p": p.PWGTP * p.ESR.isin([1, 2, 4, 5]), "atw_p": p.PWGTP * p.ESR.isin([1, 4])}).groupby(p.PUMA).sum()
    k = h.HWORKERS_AT_WORK.clip(upper=2)
    for i in range(3):
        g["s%d" % i] = (h[wcol] * (k == i)).groupby(h.PUMA).sum()
    g[["s0", "s1", "s2"]] = g[["s0", "s1", "s2"]].div(g[["s0", "s1", "s2"]].sum(axis=1), axis=0)
    w2 = h[h.HWORKERS_AT_WORK >= 2]
    g["m2"] = (w2.HWORKERS_AT_WORK * w2[wcol]).groupby(w2.PUMA).sum() / w2[wcol].groupby(w2.PUMA).sum()
    g["emp_per_atw"] = g.emp_p / g.atw_p
    return g


def rebuild_worker_controls(data_dir, fallback):
    """Rebuild the tract worker controls from the reweighted seed, each with its own definition:
    EMPWORKER (at work) x PUMA at-work ratio; HHWORKER0-2 tract shares x PUMA share ratio, renormalized;
    INDUSTRY1-13 (civilian employed) x PUMA ratio, INDUSTRY14 unchanged. Fallback tracts are rebuilt with
    the fallback formula and reweighted PUMA values. Ratios = reweighted / pooled seed."""
    workers, civ = WORKER_COLS, ["INDUSTRY%d" % i for i in range(1, 14)]
    h = pd.read_csv(data_dir / COPY_DATA[0], low_memory=False)
    p = pd.read_csv(data_dir / COPY_DATA[1], usecols=["hh_id", "PUMA", "ESR", "PWGTP"])
    p0 = p.copy()
    p0["PWGTP"] = p.PWGTP / p.hh_id.map(h.set_index("hh_id").eval("WGTP / WGTP_pooled"))
    old, new = _rw_puma_stats(h, p0, "WGTP_pooled"), _rw_puma_stats(h, p, "WGTP")
    xw = pd.read_csv(data_dir / COPY_DATA[2], dtype={"TRACTID": str})
    t = pd.read_csv(TRACT_OUT, dtype={"TRACTID": str}).set_index("TRACTID")
    puma = t.index.map(xw.drop_duplicates("TRACTID").set_index("TRACTID").PUMA)
    out = t.copy()
    hh = t[workers].sum(axis=1).to_numpy()
    share_ratio = (new[["s0", "s1", "s2"]] / old[["s0", "s1", "s2"]]).loc[puma].to_numpy()
    shares = t[workers].to_numpy() / np.maximum(hh, 1)[:, None] * share_ratio
    shares = shares / np.maximum(shares.sum(axis=1), 1e-12)[:, None]
    out[workers] = allocate_group(shares, hh)
    r_atw = (new.atw / old.atw).loc[puma].to_numpy()
    r_civ = (new.civ / old.civ).loc[puma].to_numpy()
    out["EMPWORKER"] = np.round(t.EMPWORKER.to_numpy() * r_atw).astype(int)
    civ_total = np.round(t[civ].sum(axis=1).to_numpy() * r_civ).astype(int)
    civ_shares = t[civ].div(t[civ].sum(axis=1).replace(0, 1), axis=0).to_numpy()
    out[civ] = np.vstack([largest_remainder(s * n, n) for s, n in zip(civ_shares, civ_total)])
    ind_shares = t[INDUSTRY_COLS].div(t[INDUSTRY_COLS].sum(axis=1).replace(0, 1), axis=0)
    for x in [x for x in t.index if x in fallback]:
        pu = puma[t.index.get_loc(x)]
        emp = int(round(out.at[x, "HHWORKER1"] + new.at[pu, "m2"] * out.at[x, "HHWORKER2"]))
        tot = int(round(emp * new.at[pu, "emp_per_atw"]))
        out.at[x, "EMPWORKER"] = emp
        out.loc[x, INDUSTRY_COLS] = largest_remainder(ind_shares.loc[x].to_numpy() * tot, tot)
    assert (out[workers].sum(axis=1) == t[workers].sum(axis=1)).all()
    assert (out.select_dtypes("number") >= 0).all().all()
    same = [c for c in t.columns if c not in workers + ["EMPWORKER"] + INDUSTRY_COLS]
    assert (out[same] == t[same]).all().all()
    out.reset_index()[t.reset_index().columns].to_csv(TRACT_OUT, index=False)
    print("worker controls rebuilt from the reweighted seed: EMPWORKER %d -> %d, industry %d -> %d"
          % (t.EMPWORKER.sum(), out.EMPWORKER.sum(), t[INDUSTRY_COLS].sum().sum(), out[INDUSTRY_COLS].sum().sum()))


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
    age_pout = None
    if AGE_2025:
        assert not USE_PUMA_PROFILE, "AGE_2025 and USE_PUMA_PROFILE both write PUMA controls"
        from age_targets_2025 import bg_large_area
        age_pout = age_controls_2025(out, bg_large_area(TARGETS).large_area_id, xw)

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
    gq_like = tr[INDUSTRY_COLS].sum(axis=1) > MAX_IND_PER_WORKER * tr.EMPWORKER
    # a group-quarters tract with few ACS households (U-M campus: 65 vs 697 target)
    # also takes its worker-household shares from the PUMA
    use_puma_hh = (tr.acs_HHBASE < MIN_HH) | (gq_like & (tr.acs_HHBASE < MIN_HH_GQ))
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
    # Fallback tracts: derive EMPWORKER from the worker households, then industry
    # from EMPWORKER, so one real household can meet all of them.
    #  - Tracts on PUMA rates: rounding the three groups separately can contradict
    #    (e.g. a 0-worker household with EMPWORKER 1) -> integerizer infeasible.
    #  - Group-quarters tracts (ACS industry > MAX_IND_PER_WORKER x EMPWORKER, e.g.
    #    campus tracts): ACS person tables count dorm residents but the rate is per
    #    household person, so industry came out up to 9x EMPWORKER. B08202 counts
    #    households only. Without this, the balancing added "employed, not at work"
    #    persons (Washtenaw 3.2x the seed; review 2026-10-13).
    m2 = tr.PUMA.map(workers_2plus_mean_by_puma())
    emp = (tout.HHWORKER1 + m2 * tout.HHWORKER2).round().astype(int)
    ind_total = (emp * tr.PUMA.map(employed_per_at_work_by_puma())).round().astype(int)
    ind_shares = rate_src[INDUSTRY_COLS].div(rate_src[INDUSTRY_COLS].sum(axis=1), axis=0).fillna(0)
    fallback = tr.index[use_puma_hh | use_puma_pop | gq_like]
    for t in fallback:
        tout.loc[t, "EMPWORKER"] = emp[t]
        tout.loc[t, INDUSTRY_COLS] = largest_remainder(ind_shares.loc[t].to_numpy() * ind_total[t], ind_total[t])
    print("tracts: %d | worker shares from PUMA: %d | industry rates from PUMA: %d | group-quarters tracts: %s"
          % (len(tr), int(use_puma_hh[tr_new.HHBASE > 0].sum()), int(use_puma_pop[tr_new.POPBASE > 0].sum()),
             list(tr.index[gq_like & ~(use_puma_hh | use_puma_pop)])))

    # ---- checks
    for cols in HH_GROUPS.values():
        assert (out[cols].sum(axis=1) == out.HHBASE).all()
    assert (out[SIZE_COLS].sum(axis=1) == out.HHBASE).all()
    for cols in POP_GROUPS.values():
        assert (out[cols].sum(axis=1) == out.POPBASE).all()
    assert (tout[WORKER_COLS].sum(axis=1) == tr_new.HHBASE).all()
    assert (tout.loc[fallback, INDUSTRY_COLS].sum(axis=1) == ind_total[fallback]).all()
    assert (out.select_dtypes("number") >= 0).all().all() and (tout >= 0).all().all()
    pout = puma_profile_controls(out, tout, xw) if USE_PUMA_PROFILE else age_pout

    # ---- write the run package
    (RUN_DIR / "configs").mkdir(parents=True, exist_ok=True)
    (RUN_DIR / "output").mkdir(exist_ok=True)
    col_order = [c for c in acs.columns if c != "BLKGRPID"] + ["BLKGRPID"]
    out.reset_index()[col_order].to_csv(BLKGRP_OUT, index=False)
    tcols = [c for c in acs_tr.columns if c != "TRACTID"]
    tout[tcols].rename_axis("TRACTID").reset_index()[tcols + ["TRACTID"]].to_csv(TRACT_OUT, index=False)
    if pout is not None:
        pout.reset_index().to_csv(PUMA_OUT, index=False)
    else:
        PUMA_OUT.unlink(missing_ok=True)
    diag.to_csv(DIAG_OUT)
    for f in COPY_DATA:
        if f == COPY_DATA[0]:
            write_seed_households(RUN_DIR / "data" / f)
        elif f == COPY_DATA[1]:
            write_seed_persons(RUN_DIR / "data" / f)
        else:
            shutil.copy2(TEMPLATE / "data" / f, RUN_DIR / "data" / f)
    if SEED_REWEIGHT:
        reweight_seed(RUN_DIR / "data", RUN_DIR / "data" / "targets" / "seed_reweight_by_puma.csv")
        rebuild_worker_controls(RUN_DIR / "data", set(fallback))
    # OR-Tools integerizer: ~50-90x faster than CVXPY/GLPK with the same quality, and
    # GLPK has no working time limit. Needs the isolated popsim conda env (in the
    # base env, ortools + conda pyarrow segfault).
    settings = (TEMPLATE / "configs" / "settings.yaml").read_text()
    assert "USE_CVXPY: true\n" in settings, "template settings changed: check USE_CVXPY"
    settings = settings.replace("USE_CVXPY: true\n", "USE_CVXPY: false  # OR-Tools; run in the popsim conda env\n")
    if pout is not None:
        tract_entry = "  tablename: TRACT_control_data\n"
        assert settings.count(tract_entry) == 1, "template settings changed: check TRACT_control_data"
        settings = settings.replace(tract_entry, tract_entry + "- filename: %s\n  tablename: PUMA_control_data\n"
                                    % PUMA_OUT.name)
    (RUN_DIR / "configs" / "settings.yaml").write_text(settings)
    ctl = split_head_age_control_rows(pd.read_csv(TEMPLATE / "configs" / "controls.csv"))
    ctl = fix_worker_rows(fix_industry_rows(ctl))
    unknown = set(IMPORTANCE) - set(ctl.target)
    assert not unknown, unknown
    ctl["importance"] = ctl.target.map(IMPORTANCE).fillna(ctl.importance).astype(int)
    if USE_PUMA_PROFILE:
        ctl = pd.concat([ctl, puma_control_rows()], ignore_index=True)
    elif AGE_2025:
        ctl = pd.concat([ctl, age_puma_control_rows()], ignore_index=True)
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
    if USE_PUMA_PROFILE:
        print("PUMA profile controls, region:\n%s" % pout.sum().to_string())
    print("wrote", RUN_DIR)


if __name__ == "__main__":
    main()
