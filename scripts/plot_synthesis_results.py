# Charts of PopulationSim results vs controls, in the formats used for past runs
# (input_prep/scripts/output_plot.py and the validation frequency plots):
#   error_plot_<GEO>.png       percent error per control: mean (dot, value printed),
#                              +-1 std (thick bar), min-max (thin bar)
#   diff_frequency_<GEO>.png   number of zones by (synthesized - control), per control
#   error_summary_<GEO>.csv    the statistics behind the error plot
#   control_vs_synth_<GEO>.png control vs synthesized count, one panel per control group
#   fit_summary.png            (with --acs-blkgrp) fit to own controls by run, fit to ACS at
#                              tract level, and within-tract BG variation kept vs raw ACS
# GEO = BLKGRP and TRACT. With --compare-dir, the error plots put a second run
# (e.g. the previous synthesis) beside this one.
#
# Usage:
#     python scripts/plot_synthesis_results.py --summary-dir .../pass2 --out-dir .../validation/charts \
#         --label "2025 MCD run" [--compare-dir .../old/pass2 --compare-label "July 2025"] \
#         [--extra-run .../run1/pass2 "run 1"] [--acs-blkgrp .../SEMCOG_2024_control_totals_blkgrp.csv]

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

GEOS = ["BLKGRP", "TRACT"]
GEO_NAME = {"BLKGRP": "block group", "TRACT": "tract"}


def pct_error_stats(summary):
    """Same rule as output_plot.pct_dif: (result - control) / control x 100;
    0 where the difference is 0; zones with control 0 and result > 0 dropped."""
    names = [c[:-8] for c in summary.columns if c.endswith("_control")]
    rows = {}
    for n in names:
        c, r = summary[n + "_control"].astype(float), summary[n + "_result"].astype(float)
        d = r - c
        pct = (d / c * 100).where(d != 0, 0.0).replace([np.inf, -np.inf], np.nan).dropna()
        rows[n] = {"mean": pct.mean(), "std": pct.std(), "min": pct.min(), "max": pct.max(),
                   "zones": len(pct), "control_total": c.sum(), "result_total": r.sum(),
                   "zones_exact": int((d == 0).sum()), "abs_diff_sum": float(d.abs().sum())}
    return pd.DataFrame(rows).T


def error_panel(ax, stats, title, xlim):
    """Whiskers are clipped to xlim (as in past plots, -100..400%); a clipped
    max is marked with an arrow and its value."""
    stats = stats.iloc[::-1]
    y = np.arange(len(stats))
    lo, hi = stats["min"].clip(lower=xlim[0]), stats["max"].clip(upper=xlim[1] * 0.8)
    ax.axvline(0, color="black", lw=1)
    ax.errorbar(stats["mean"], y, xerr=[stats["mean"] - lo, hi - stats["mean"]],
                fmt="none", ecolor="gray", lw=1, capsize=4)
    ax.errorbar(stats["mean"], y, xerr=stats["std"].clip(upper=xlim[1] * 0.8), fmt="ok", lw=4, ms=5)
    for yi, mx in zip(y, stats["max"]):
        if mx > xlim[1] * 0.8:
            ax.annotate("max %.0f" % mx, (xlim[1] * 0.8, yi), xytext=(4, 0), textcoords="offset points",
                        va="center", fontsize=7, color="gray", arrowprops=None)
    for yi, m in zip(y, stats["mean"]):
        ax.text(xlim[1] * 0.99, yi, "%.2f" % m, ha="right", va="center", fontsize=9)
    ax.set_yticks(y)
    ax.set_yticklabels(stats.index, fontsize=10)
    ax.set_xlim(*xlim)
    ax.set_title(title, fontsize=12)
    ax.set_xlabel("% error (synthesized vs control)")
    ax.grid(axis="x", alpha=0.3)


def error_plot(panels, geo, year, path):
    xlim = (-105, 400)
    n = max(len(s) for _, s in panels)
    fig, axs = plt.subplots(1, len(panels), figsize=(7 * len(panels), 0.32 * n + 1.5), sharey=True, squeeze=False)
    for ax, (label, stats) in zip(axs[0], panels):
        error_panel(ax, stats, label, xlim)
    fig.suptitle("%d synthesis: percent error by control, %s level\n"
                 "dot = mean (value at right), thick bar = +-1 std, thin bar = min-max" % (year, GEO_NAME[geo]),
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def frequency_plot(summary, geo, year, label, path, span=25):
    names = [c[:-8] for c in summary.columns if c.endswith("_control")]
    cols = 6
    rows = int(np.ceil(len(names) / cols))
    fig, axs = plt.subplots(rows, cols, figsize=(3.2 * cols, 2.4 * rows), squeeze=False)
    for ax, n in zip(axs.ravel(), names):
        d = (summary[n + "_result"] - summary[n + "_control"]).round().astype(int)
        vc = d.clip(-span, span).value_counts().sort_index()
        ax.vlines(vc.index, 0, vc.values, color="#d9662b", lw=1.5)
        ax.plot(vc.index, vc.values, "o", color="#d9662b", ms=3)
        ax.axvline(0, color="#3b6ea5", lw=1)
        ax.set_xlim(-span - 1, span + 1)
        exact = (d == 0).mean()
        ax.set_title("%s\nexact %.0f%%" % (n, 100 * exact), fontsize=9)
        ax.tick_params(labelsize=7)
    for ax in axs.ravel()[len(names):]:
        ax.axis("off")
    fig.supxlabel("difference (households or persons)")
    fig.supylabel("number of zones")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.suptitle("%d synthesis, %s: zones by synthesized - control, %s level (|diff| > %d shown at +-%d)"
                 % (year, label, GEO_NAME[geo], span, span), fontsize=13)
    fig.savefig(path, dpi=110, bbox_inches="tight")
    plt.close(fig)


# reference palette (dataviz skill), validated light-mode all-pairs for 3 slots
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
GROUPS = {"age of head": "hh_age_", "race of head": "hh_race_", "hispanic": "hh_hisp",
          "children": "hh_children", "income": "hh_inc_", "cars": "hh_cars_", "size": "hh_persons_",
          "tenure": "hh_tenure_", "workers": "hh_workers_", "person age": "persons_age_",
          "person race": "persons_race_", "person sex": "persons_sex_", "industry": "Persons_industry_"}
ACS_GROUPS = {"age of head": ("HHAGE", "hh_age_"), "race of head": ("HHRACE", "hh_race_"),
              "children": ("HHCHD", "hh_children"), "income": ("HHINC", "hh_inc_"),
              "cars": ("HHCAR", "hh_cars_"), "tenure": ("HHTENURE", "hh_tenure_")}


def _style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def group_columns(summary, prefix):
    return [c[:-8] for c in summary.columns if c.endswith("_control") and c.startswith(prefix)]


def misfit(summary, prefix):
    """Share of HH or persons in a wrong category: sum |result - control| / 2 / total."""
    cs = group_columns(summary, prefix)
    if not cs:
        return np.nan
    c = summary[[x + "_control" for x in cs]].to_numpy()
    r = summary[[x + "_result" for x in cs]].to_numpy()
    return 100 * np.abs(r - c).sum() / 2 / c.sum()


def scatter_plot(summary, geo, year, label, path):
    groups = [(g, p) for g, p in GROUPS.items() if group_columns(summary, p)]
    cols = 4
    rows = int(np.ceil(len(groups) / cols))
    fig, axs = plt.subplots(rows, cols, figsize=(3.4 * cols, 3.3 * rows), squeeze=False, facecolor="#fcfcfb")
    for ax, (g, p) in zip(axs.ravel(), groups):
        cs = group_columns(summary, p)
        c = summary[[x + "_control" for x in cs]].to_numpy().ravel()
        r = summary[[x + "_result" for x in cs]].to_numpy().ravel()
        top = np.quantile(np.maximum(c, r), 0.999) * 1.05
        ax.plot([0, top], [0, top], color=INK2, lw=1, ls="--", zorder=1)
        ax.scatter(c, r, s=4, color=SERIES[0], alpha=0.25, lw=0, zorder=2)
        ax.set_xlim(0, top)
        ax.set_ylim(0, top)
        ax.set_aspect("equal")
        rmse = np.sqrt(np.mean((r - c) ** 2))
        within = np.mean(np.abs(r - c) <= 2)
        ax.set_title(g, fontsize=10, color=INK, loc="left")
        ax.text(0.03, 0.97, "RMSE %.1f\nwithin +-2: %.0f%%\ncells: %s" % (rmse, 100 * within, format(len(c), ",")),
                transform=ax.transAxes, va="top", fontsize=8, color=INK2)
        _style(ax)
    for ax in axs.ravel()[len(groups):]:
        ax.axis("off")
    fig.supxlabel("control (households or persons per %s and category)" % GEO_NAME[geo], color=INK2, fontsize=10)
    fig.supylabel("synthesized", color=INK2, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.suptitle("%d synthesis, %s: control vs synthesized, %s level (dashed = 1:1; top 0.1%% of cells beyond the axes)"
                 % (year, label, GEO_NAME[geo]), fontsize=12, color=INK)
    fig.savefig(path, dpi=110, bbox_inches="tight", facecolor="#fcfcfb")
    plt.close(fig)


def _acs_tract_fit_and_spread(summary, acs):
    """Fit to raw ACS at tract level, and within-tract BG variance of shares
    (BGs with 100+ HH), for the groups in ACS_GROUPS."""
    fit, spread = {}, {}
    for g, (a_pref, s_pref) in ACS_GROUPS.items():
        acols = [c for c in acs.columns if c.startswith(a_pref)]
        res = summary[[c + "_result" for c in group_columns(summary, s_pref)]]
        if res.shape[1] == len(acols) + 1:  # age of head split 65-74 / 75+; raw ACS has 65+
            res = pd.concat([res.iloc[:, :-2], res.iloc[:, -2:].sum(axis=1)], axis=1)
        res.columns = acols
        tr_a = acs[acols].groupby(acs.index.str[:11]).sum()
        tr_r = res.groupby(res.index.str[:11]).sum()
        ok = tr_a.index[tr_a.sum(axis=1) >= 50].intersection(tr_r.index)
        sa = tr_a.loc[ok].div(tr_a.loc[ok].sum(axis=1), axis=0)
        n = tr_r.loc[ok].sum(axis=1)
        sr = tr_r.loc[ok].div(n, axis=0)
        fit[g] = 100 * ((sr - sa).abs().sum(axis=1) / 2 * n).sum() / n.sum()

        def within_var(df):
            tot = df.sum(axis=1)
            sh = df[tot >= 100].div(tot[tot >= 100], axis=0)
            return ((sh - sh.groupby(sh.index.str[:11]).transform("mean")) ** 2).sum(axis=1).mean()
        spread[g] = within_var(res) / within_var(acs[acols])
    return pd.Series(fit), pd.Series(spread)


def fit_summary_plot(runs, acs, year, path):
    """runs: list of (label, BG summary), oldest first; the last is this run."""
    fig, axs = plt.subplots(1, 3, figsize=(16, 5.6), facecolor="#fcfcfb",
                            gridspec_kw={"width_ratios": [1.25, 1, 1]})
    names = [g for g in GROUPS if g not in ("workers", "industry")]
    y = np.arange(len(names))[::-1]
    ax = axs[0]
    for k, (label, s) in enumerate(runs):
        v = [misfit(s, GROUPS[g]) for g in names]
        ax.plot(v, y + (k - (len(runs) - 1) / 2) * 0.18, "o", ms=7, color=SERIES[(len(runs) - 1 - k) % 3],
                label=label, zorder=3)
        if k == len(runs) - 1:
            for vi, yi in zip(v, y):
                ax.text(vi, yi + (k - (len(runs) - 1) / 2) * 0.18 - 0.35, "%.2f" % vi, ha="center",
                        fontsize=7, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=9, color=INK)
    ax.set_xlim(left=0)
    ax.set_xlabel("% of households or persons in a wrong category", color=INK2, fontsize=9)
    ax.set_title("(a) Fit to own controls, block-group level", loc="left", fontsize=11, color=INK)
    _style(ax)
    handles, labels = ax.get_legend_handles_labels()

    stats = [(label, *_acs_tract_fit_and_spread(s, acs)) for label, s in runs[-2:]]
    gy = np.arange(len(ACS_GROUPS))[::-1]
    for idx, (panel, title, xlabel) in enumerate([
            (1, "(b) Fit to ACS, tract level", "% in a wrong category vs raw ACS tract shares"),
            (2, "(c) Within-tract BG variation kept", "share of raw ACS within-tract variance")]):
        ax = axs[panel]
        for k, (label, fit, spread) in enumerate(stats):
            v = (fit if panel == 1 else spread * 100)[list(ACS_GROUPS)]
            off = (k - 0.5) * 0.22
            color = SERIES[(len(stats) - 1 - k) % 3]
            ax.plot(v, gy + off, "o", ms=7, color=color, label=label, zorder=3)
            for vi, yi in zip(v, gy):
                dy = 0.16 if k == len(stats) - 1 else -0.3  # this run above its dot, older run below
                ax.text(vi, yi + off + dy, ("%.2f" if panel == 1 else "%.0f%%") % vi, ha="center",
                        fontsize=7, color=INK2)
        if panel == 2:
            ax.axvline(100, color=INK2, lw=1, ls="--")
            ax.text(99, gy[0] + 0.45, "raw ACS = 100%\n(includes sampling noise)", ha="right", va="bottom",
                    fontsize=7, color=INK2)
            ax.set_xlim(0, 110)
            ax.set_ylim(gy[-1] - 0.6, gy[0] + 1.1)
        else:
            ax.set_xlim(left=0)
        ax.set_yticks(gy)
        ax.set_yticklabels(list(ACS_GROUPS), fontsize=9, color=INK)
        ax.set_xlabel(xlabel, color=INK2, fontsize=9)
        ax.set_title(title, loc="left", fontsize=11, color=INK)
        _style(ax)
    fig.suptitle("%d synthesis: fit and spatial detail by run" % year, fontsize=13, color=INK, x=0.01, ha="left")
    fig.legend(handles, labels, loc="upper right", ncol=len(runs), frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=120, bbox_inches="tight", facecolor="#fcfcfb")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Charts of PopulationSim results vs controls")
    ap.add_argument("--summary-dir", required=True, type=Path, help="folder with final_summary_<GEO>.csv")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--label", default="this run")
    ap.add_argument("--year", type=int, default=2025)
    ap.add_argument("--compare-dir", type=Path, help="second run to show beside this one")
    ap.add_argument("--compare-label", default="previous run")
    ap.add_argument("--extra-run", nargs=2, action="append", metavar=("DIR", "LABEL"), default=[],
                    help="more runs for fit_summary.png, oldest first, between --compare-dir and this run")
    ap.add_argument("--acs-blkgrp", type=Path, help="raw ACS BG control file; enables fit_summary.png")
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    for geo in GEOS:
        s = pd.read_csv(args.summary_dir / ("final_summary_%s.csv" % geo), index_col="id")
        stats = pct_error_stats(s)
        stats.to_csv(args.out_dir / ("error_summary_%s.csv" % geo))
        panels = []
        if args.compare_dir:
            old = pd.read_csv(args.compare_dir / ("final_summary_%s.csv" % geo), index_col="id")
            panels.append((args.compare_label, pct_error_stats(old).reindex(stats.index)))
        panels.append((args.label, stats))
        error_plot(panels, geo, args.year, args.out_dir / ("error_plot_%s.png" % geo))
        frequency_plot(s, geo, args.year, args.label, args.out_dir / ("diff_frequency_%s.png" % geo))
        scatter_plot(s, geo, args.year, args.label, args.out_dir / ("control_vs_synth_%s.png" % geo))
        print("%s: %d controls, %d zones -> %s" % (geo, len(stats), len(s), args.out_dir))

    if args.acs_blkgrp:
        read = lambda d: pd.read_csv(Path(d) / "final_summary_BLKGRP.csv", dtype={"id": str}).set_index("id")
        runs = ([(args.compare_label, read(args.compare_dir))] if args.compare_dir else []) + \
            [(lab, read(d)) for d, lab in args.extra_run] + [(args.label, read(args.summary_dir))]
        acs = pd.read_csv(args.acs_blkgrp, dtype={"BLKGRPID": str}).set_index("BLKGRPID")
        fit_summary_plot(runs[-3:], acs, args.year, args.out_dir / "fit_summary.png")
        print("fit_summary.png ->", args.out_dir)


if __name__ == "__main__":
    main()
