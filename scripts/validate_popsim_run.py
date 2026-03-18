from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize PopulationSim run quality and suggest next adjustments."
    )
    parser.add_argument(
        "--run-dir",
        required=True,
        help="PopulationSim run folder containing configs/, data/, and output/.",
    )
    parser.add_argument(
        "--top-controls",
        type=int,
        default=10,
        help="Number of worst controls to print.",
    )
    parser.add_argument(
        "--top-geos",
        type=int,
        default=5,
        help="Number of worst geographies to print per control.",
    )
    parser.add_argument(
        "--write-csv",
        action="store_true",
        help="Write validation CSVs into output/validation/.",
    )
    return parser.parse_args()


def find_summary_files(output_dir: Path) -> list[Path]:
    patterns = ["summary_*.csv", "final_summary_*.csv"]
    files = []
    for pattern in patterns:
        files.extend(sorted(output_dir.glob(pattern)))
    return files


def summary_geo_name(path: Path) -> str:
    stem = path.stem
    for prefix in ("summary_", "final_summary_"):
        if stem.startswith(prefix):
            return stem[len(prefix) :]
    return stem


def pct_error(result: pd.Series, control: pd.Series) -> pd.Series:
    control = control.astype(float)
    result = result.astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        pct = (result - control) / control.replace(0, np.nan) * 100.0
    return pct


def control_columns(df: pd.DataFrame) -> list[str]:
    return sorted(col for col in df.columns if col.endswith("_control"))


def id_column(df: pd.DataFrame) -> str | None:
    for candidate in ["id", "geography", "BLKGRP", "BLKGRPID", "TRACT", "TRACTID", "TAZ", "ZONE", "zone_id"]:
        if candidate in df.columns:
            return candidate
    return None


def summarize_file(path: Path) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    df = pd.read_csv(path)
    geography = summary_geo_name(path)
    identifier = id_column(df)
    metrics = []
    worst_geos: dict[str, pd.DataFrame] = {}

    for control_col in control_columns(df):
        base = control_col[: -len("_control")]
        result_col = f"{base}_result"
        if result_col not in df.columns:
            continue

        control = df[control_col].fillna(0)
        result = df[result_col].fillna(0)
        abs_diff = (result - control).abs()
        pct = pct_error(result, control)
        abs_pct = pct.abs()
        nonzero = control > 0

        metrics.append(
            {
                "summary_file": path.name,
                "geography": geography,
                "control": base,
                "n_nonzero": int(nonzero.sum()),
                "target_total": float(control.sum()),
                "result_total": float(result.sum()),
                "total_diff": float((result - control).sum()),
                "total_pct_diff": float(((result.sum() - control.sum()) / control.sum() * 100.0) if control.sum() else np.nan),
                "mean_abs_pct_diff": float(abs_pct[nonzero].mean() if nonzero.any() else np.nan),
                "rmse_pct_diff": float(np.sqrt(np.nanmean(np.square(pct[nonzero]))) if nonzero.any() else np.nan),
                "stdev_pct_diff": float(pct[nonzero].std() if nonzero.any() else np.nan),
                "max_abs_pct_diff": float(abs_pct[nonzero].max() if nonzero.any() else np.nan),
                "mean_abs_diff": float(abs_diff.mean()),
            }
        )

        if identifier is not None:
            worst = pd.DataFrame(
                {
                    "geography_id": df[identifier],
                    "control": control,
                    "result": result,
                    "abs_diff": abs_diff,
                    "pct_diff": pct,
                    "abs_pct_diff": abs_pct,
                }
            )
            worst = worst.loc[control > 0].sort_values("abs_pct_diff", ascending=False)
            worst_geos[base] = worst

    metrics_df = pd.DataFrame(metrics)
    return metrics_df, worst_geos


def load_controls(configs_dir: Path) -> pd.DataFrame | None:
    path = configs_dir / "controls.csv"
    if not path.exists():
        return None
    return pd.read_csv(path)


def suggest_adjustments(metrics: pd.DataFrame, controls_df: pd.DataFrame | None) -> list[str]:
    suggestions: list[str] = []
    if metrics.empty:
        return ["No summary metrics available yet. Run the synthesis to completion so summary_<geo>.csv files are written."]

    worst = metrics.sort_values(["rmse_pct_diff", "mean_abs_pct_diff"], ascending=False).head(8)
    seen = set()
    for _, row in worst.iterrows():
        control = str(row["control"])
        geography = str(row["geography"])
        key = (control, geography)
        if key in seen:
            continue
        seen.add(key)

        if control.startswith("hh_persons_") or control in {"persons_num", "num_hh"}:
            suggestions.append(
                f"{control} at {geography}: inspect household-size controls first; if large-household fit is poor, rebalance the hh size targets before changing other controls."
            )
        elif control.startswith("hh_"):
            suggestions.append(
                f"{control} at {geography}: if the miss is systematic across many zones, consider modestly increasing its importance in controls.csv or reviewing the target totals for that control."
            )
        elif control.startswith("persons_") or control.startswith("Persons_"):
            suggestions.append(
                f"{control} at {geography}: persistent misses usually point to person-control importance, target inconsistency, or a seed sample that lacks enough matching households."
            )
        else:
            suggestions.append(
                f"{control} at {geography}: review the control definition and whether the seed data has enough support in the affected zones."
            )

    if controls_df is not None and not controls_df.empty:
        high_importance = controls_df.sort_values("importance", ascending=False).head(5)
        top_names = ", ".join(high_importance["target"].astype(str).tolist())
        suggestions.append(
            f"Current highest-importance controls are: {top_names}. If one of the worst controls is missing from this set, that is a candidate for a small importance increase."
        )

    suggestions.append(
        "If only a handful of geographies are failing badly, check those zones for inconsistent control totals before changing global weights or importances."
    )
    return suggestions


def maybe_write_csv(output_dir: Path, metrics: pd.DataFrame, worst_geos: dict[str, pd.DataFrame]) -> list[Path]:
    validation_dir = output_dir / "validation"
    validation_dir.mkdir(parents=True, exist_ok=True)
    written = []

    metrics_path = validation_dir / "control_fit_summary.csv"
    metrics.to_csv(metrics_path, index=False)
    written.append(metrics_path)

    all_worst = []
    for control, df in worst_geos.items():
        if df.empty:
            continue
        temp = df.copy()
        temp.insert(0, "control_name", control)
        all_worst.append(temp)
    if all_worst:
        worst_path = validation_dir / "worst_geographies.csv"
        pd.concat(all_worst, ignore_index=True).to_csv(worst_path, index=False)
        written.append(worst_path)

    return written


def print_section(title: str) -> None:
    print(f"\n{title}")
    print("-" * len(title))


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir).resolve()
    configs_dir = run_dir / "configs"
    output_dir = run_dir / "output"

    summary_files = find_summary_files(output_dir)
    if not summary_files:
        print("No summary output files found yet.")
        print(f"Checked: {output_dir}")
        print("Expected files like summary_BLKGRP.csv or final_summary_BLKGRP.csv.")
        print("Finish the PopulationSim run first, then rerun this validator.")
        return

    metrics_frames = []
    worst_by_control: dict[str, pd.DataFrame] = {}
    for path in summary_files:
        metrics_df, worst_geos = summarize_file(path)
        if not metrics_df.empty:
            metrics_frames.append(metrics_df)
        for control, df in worst_geos.items():
            worst_by_control[f"{summary_geo_name(path)}::{control}"] = df

    metrics = pd.concat(metrics_frames, ignore_index=True) if metrics_frames else pd.DataFrame()
    controls_df = load_controls(configs_dir)

    print_section("Run")
    print(run_dir)
    print("summary files:")
    for path in summary_files:
        print(f"- {path.name}")

    if metrics.empty:
        print_section("Status")
        print("Summary files were found, but no *_control / *_result pairs could be read.")
        return

    worst_controls = metrics.sort_values(["rmse_pct_diff", "mean_abs_pct_diff"], ascending=False).head(args.top_controls)
    print_section("Worst Controls")
    cols = [
        "geography",
        "control",
        "n_nonzero",
        "total_pct_diff",
        "mean_abs_pct_diff",
        "rmse_pct_diff",
        "max_abs_pct_diff",
    ]
    print(worst_controls.loc[:, cols].to_string(index=False, float_format=lambda x: f"{x:,.2f}"))

    print_section("Worst Geographies")
    for _, row in worst_controls.head(min(5, len(worst_controls))).iterrows():
        key = f"{row['geography']}::{row['control']}"
        df = worst_by_control.get(key)
        if df is None or df.empty:
            continue
        print(f"{row['control']} ({row['geography']})")
        print(df.head(args.top_geos).to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
        print()

    print_section("Suggested Adjustments")
    for item in suggest_adjustments(metrics, controls_df):
        print(f"- {item}")

    if args.write_csv:
        written = maybe_write_csv(output_dir, metrics, worst_by_control)
        print_section("Wrote CSVs")
        for path in written:
            print(f"- {path}")


if __name__ == "__main__":
    main()
