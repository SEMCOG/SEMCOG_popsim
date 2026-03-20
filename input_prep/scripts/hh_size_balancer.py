# %% [markdown]
# Household size control balancer
#
# This script adjusts household-size control totals after a first PopulationSim run.
# It preserves total households and total persons while recalibrating the household
# size distribution, especially the open-ended 7+ bucket.
#
# It is designed to work with the run-based folder structure created by
# `input_prep/scripts/popsim_input_maker.py`.

# %%
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    import oyaml as yaml
except ModuleNotFoundError:
    import yaml


# %%
SIZE_COLS = [f"HHPERSONS{i}" for i in range(1, 8)]
SUMMARY_CONTROL_COLS = [f"hh_persons_{i}_control" for i in range(1, 8)]
SUMMARY_RESULT_COLS = [f"hh_persons_{i}_result" for i in range(1, 8)]
PERSONS_CONTROL_COL = "persons_num_control"
PERSONS_RESULT_COL = "persons_num_result"
DEFAULT_RANDOM_SEED = 42
DEFAULT_MIN_TOP_BIN_SIZE = 7.0
DEFAULT_MAX_TOP_BIN_SIZE = 10.0
DEFAULT_SUFFIX = "_hhsize_adj"
DEFAULT_METHOD = "shape_preserving"

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_PREP_DIR = SCRIPT_DIR.parent
REPO_ROOT = INPUT_PREP_DIR.parent
DEFAULT_RUN_ROOT = REPO_ROOT.parent / "d_drive" / "popsim" / "runs"


# %%
def resolve_config_path(path_str: str | None, config_dir: Path) -> Path | None:
    if not path_str:
        return None
    path = Path(path_str)
    if path.is_absolute():
        return path
    candidate = config_dir / path
    if candidate.exists():
        return candidate
    candidate = INPUT_PREP_DIR / path
    if candidate.exists():
        return candidate
    return REPO_ROOT / path


def find_single_file(base_dir: Path, pattern: str) -> Path:
    matches = sorted(base_dir.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No files matched pattern {pattern!r} in {base_dir}")
    if len(matches) > 1:
        raise FileExistsError(
            f"Expected one match for pattern {pattern!r} in {base_dir}, found {len(matches)}: {matches}"
        )
    return matches[0]


def inferred_top_bin_size_from_summary(summary_row: pd.Series) -> float:
    hh7 = float(summary_row[SUMMARY_RESULT_COLS[-1]])
    if hh7 <= 0:
        return np.nan
    lower_pop = sum(float(summary_row[col]) * size for size, col in zip(range(1, 7), SUMMARY_RESULT_COLS[:-1]))
    return (float(summary_row[PERSONS_RESULT_COL]) - lower_pop) / hh7


def implied_top_bin_size_from_controls(control_row: pd.Series) -> float:
    hh7 = float(control_row[SIZE_COLS[-1]])
    if hh7 <= 0:
        return np.nan
    lower_pop = sum(float(control_row[col]) * size for size, col in zip(range(1, 7), SIZE_COLS[:-1]))
    return (float(control_row["POPBASE"]) - lower_pop) / hh7


def bounded_top_bin_size(value: float, minimum: float, maximum: float) -> float:
    if pd.isna(value):
        return minimum
    return float(np.clip(value, minimum, maximum))


def rebalance_household_sizes(
    x0: np.ndarray,
    target_population: float,
    size_weights: np.ndarray,
    rng: np.random.Generator,
    max_iterations: int = 100000,
) -> np.ndarray:
    """
    Reallocate household counts across size bins while preserving total households.

    This preserves the original logic direction of the legacy script: if the current
    weighted household-size distribution understates the person target, it shifts
    households from lower bins to higher bins until the target is met.
    """
    x = np.rint(np.asarray(x0, dtype=float)).astype(int).copy()
    total_households = int(x.sum())
    if total_households == 0:
        return x

    weights = np.asarray(size_weights, dtype=float)
    current_population = float(np.dot(x, weights))
    iterations = 0

    while current_population < target_population and x.max() < total_households and iterations < max_iterations:
        eligible_i = np.flatnonzero(x[:-1] > 0)
        if eligible_i.size == 0:
            break
        p_i = x[eligible_i].astype(float)
        p_i /= p_i.sum()
        i = int(rng.choice(eligible_i, p=p_i))

        eligible_j = np.arange(i + 1, len(x))
        if eligible_j.size == 0:
            break
        p_j = x[eligible_j].astype(float) + 1.0
        p_j /= p_j.sum()
        j = int(rng.choice(eligible_j, p=p_j))

        x[i] -= 1
        x[j] += 1
        current_population = float(np.dot(x, weights))
        iterations += 1

    return x


def _shape_penalty(counts: np.ndarray, original_counts: np.ndarray) -> float:
    scale = np.maximum(original_counts, 1.0)
    return float(np.sum(((counts - original_counts) ** 2) / scale))


def _round_preserving_sum(values: np.ndarray, total: int) -> np.ndarray:
    floored = np.floor(values).astype(int)
    remainder = total - int(floored.sum())
    if remainder > 0:
        order = np.argsort(-(values - floored))
        floored[order[:remainder]] += 1
    elif remainder < 0:
        order = np.argsort(values - floored)
        for idx in order:
            if remainder == 0:
                break
            if floored[idx] > 0:
                floored[idx] -= 1
                remainder += 1
    return floored


def _tilted_distribution(
    original_counts: np.ndarray,
    lambda_value: float,
    total_households: int,
) -> np.ndarray:
    base = np.asarray(original_counts, dtype=float) + 1e-9
    steps = np.arange(len(base), dtype=float)
    log_weights = np.log(base) + lambda_value * steps
    log_weights -= np.max(log_weights)
    tilted = np.exp(log_weights)
    tilted *= total_households / tilted.sum()
    return tilted


def _refine_population_gap(
    counts: np.ndarray,
    original_counts: np.ndarray,
    target_population: float,
    size_weights: np.ndarray,
    max_iterations: int = 5000,
) -> np.ndarray:
    x = counts.copy()
    current_population = float(np.dot(x, size_weights))
    iterations = 0

    while iterations < max_iterations:
        gap = target_population - current_population
        if abs(gap) < 0.5:
            break

        best_move = None
        best_score = None
        best_population_change = 0.0
        current_penalty = _shape_penalty(x, original_counts)

        for i in range(len(x) - 1):
            direction = 1 if gap > 0 else -1
            src = i if direction > 0 else i + 1
            dst = i + 1 if direction > 0 else i
            if x[src] <= 0:
                continue

            candidate = x.copy()
            candidate[src] -= 1
            candidate[dst] += 1
            population_change = float(np.dot(candidate - x, size_weights))

            if gap > 0 and population_change <= 0:
                continue
            if gap < 0 and population_change >= 0:
                continue

            penalty_change = _shape_penalty(candidate, original_counts) - current_penalty
            score = penalty_change - 0.001 * min(abs(population_change), abs(gap))
            if best_score is None or score < best_score:
                best_score = score
                best_move = (src, dst)
                best_population_change = population_change

        if best_move is None:
            break

        src, dst = best_move
        x[src] -= 1
        x[dst] += 1
        current_population += best_population_change
        iterations += 1

    return x


def rebalance_household_sizes_shape_preserving(
    x0: np.ndarray,
    target_population: float,
    size_weights: np.ndarray,
    max_iterations: int = 5000,
) -> np.ndarray:
    """
    Reallocate household counts by smoothly tilting the original household-size
    curve, then making small adjacent-bin corrections to close any remaining
    person-gap after rounding.
    """
    x = np.rint(np.asarray(x0, dtype=float)).astype(int).copy()
    total_households = int(x.sum())
    if total_households == 0:
        return x

    weights = np.asarray(size_weights, dtype=float)
    base_population = float(np.dot(x, weights))
    if abs(base_population - target_population) < 0.5:
        return x

    def population_for_lambda(lambda_value: float) -> float:
        tilted = _tilted_distribution(x, lambda_value, total_households)
        return float(np.dot(tilted, weights))

    lower, upper = -4.0, 4.0
    pop_lower = population_for_lambda(lower)
    pop_upper = population_for_lambda(upper)
    while target_population < pop_lower:
        upper = lower
        pop_upper = pop_lower
        lower *= 2.0
        pop_lower = population_for_lambda(lower)
        if lower < -64:
            break
    while target_population > pop_upper:
        lower = upper
        pop_lower = pop_upper
        upper *= 2.0
        pop_upper = population_for_lambda(upper)
        if upper > 64:
            break

    for _ in range(60):
        mid = (lower + upper) / 2.0
        pop_mid = population_for_lambda(mid)
        if pop_mid < target_population:
            lower = mid
        else:
            upper = mid

    tilted = _tilted_distribution(x, (lower + upper) / 2.0, total_households)
    rounded = _round_preserving_sum(tilted, total_households)
    return _refine_population_gap(rounded, x.astype(float), target_population, weights, max_iterations=max_iterations)


def load_yaml_config(config_path: Path) -> dict[str, Any]:
    with open(config_path, "r") as stream:
        return yaml.load(stream, Loader=yaml.Loader)


def derive_run_dir(config_path: Path, conf: dict[str, Any]) -> Path:
    config_dir = config_path.parent
    run_name = conf["project"].get("run_name", f"{conf['project']['acs_year']}_synthesis")
    paths_conf = conf.get("paths", {})
    run_root = resolve_config_path(paths_conf.get("run_root", str(DEFAULT_RUN_ROOT)), config_dir)
    return run_root / run_name


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Rebalance household-size controls after a PopulationSim run.")
    parser.add_argument("--config", help="input_prep prepare.yaml file used to derive the run folder and optional defaults")
    parser.add_argument("--run-dir", help="PopulationSim run folder containing data/ and output/")
    parser.add_argument("--control-file", help="Explicit block-group control totals CSV to adjust")
    parser.add_argument("--summary-file", help="Explicit PopulationSim summary CSV to use as the first-pass result")
    parser.add_argument("--min-top-bin-size", type=float, default=None, help="Lower bound for the effective 7+ household size")
    parser.add_argument("--max-top-bin-size", type=float, default=None, help="Upper bound for the effective 7+ household size")
    parser.add_argument("--random-seed", type=int, default=None, help="Random seed used by the rebalancing heuristic")
    parser.add_argument("--suffix", default=None, help="Suffix for the adjusted control file")
    parser.add_argument("--replace-original", action="store_true", help="Replace the original control file after writing a backup copy")
    parser.add_argument("--diagnostics-file", help="Optional explicit diagnostics CSV path")
    parser.add_argument(
        "--method",
        choices=["shape_preserving", "legacy"],
        default=None,
        help="Balancing method to use (default: shape_preserving)",
    )
    return parser.parse_args()


def resolve_inputs(args: argparse.Namespace) -> tuple[Path, Path, Path, dict[str, Any]]:
    config_options: dict[str, Any] = {}
    run_dir: Path | None = Path(args.run_dir).resolve() if args.run_dir else None

    if args.config:
        config_path = Path(args.config).resolve()
        conf = load_yaml_config(config_path)
        config_options = conf.get("postprocess", {}).get("hh_size_balancer", {}) or {}
        if run_dir is None:
            run_dir = derive_run_dir(config_path, conf)

    if run_dir is None:
        raise ValueError("Either --run-dir or --config is required.")

    data_dir = run_dir / "data"
    output_dir = run_dir / "output"

    control_file = Path(args.control_file).resolve() if args.control_file else None
    if control_file is None:
        configured = config_options.get("control_file")
        if configured:
            control_file = Path(configured).resolve() if Path(configured).is_absolute() else (data_dir / configured)
        else:
            pattern = config_options.get("control_file_pattern", "*_control_totals_blkgrp.csv")
            control_file = find_single_file(data_dir, pattern)

    summary_file = Path(args.summary_file).resolve() if args.summary_file else None
    if summary_file is None:
        configured = config_options.get("summary_file")
        if configured:
            summary_file = Path(configured).resolve() if Path(configured).is_absolute() else (output_dir / configured)
        else:
            for pattern in [config_options.get("summary_file_pattern", "final_summary_BLKGRP.csv"), "summary_BLKGRP.csv"]:
                matches = sorted(output_dir.glob(pattern))
                if matches:
                    summary_file = matches[0]
                    break
            if summary_file is None:
                raise FileNotFoundError(f"Could not find summary file in {output_dir}")

    return run_dir, control_file, summary_file, config_options


def build_adjusted_controls(
    control_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    minimum_top_bin_size: float,
    maximum_top_bin_size: float,
    random_seed: int,
    method: str = DEFAULT_METHOD,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    control_df = control_df.copy()
    summary_df = summary_df.copy()

    if "BLKGRPID" not in control_df.columns:
        raise KeyError("Control file must contain BLKGRPID")
    if "id" not in summary_df.columns:
        raise KeyError("Summary file must contain id")

    control_df = control_df.set_index("BLKGRPID", drop=False)
    summary_df = summary_df.set_index("id", drop=False)

    missing_control_cols = [c for c in ["POPBASE", *SIZE_COLS] if c not in control_df.columns]
    missing_summary_cols = [c for c in [PERSONS_CONTROL_COL, PERSONS_RESULT_COL, *SUMMARY_RESULT_COLS] if c not in summary_df.columns]
    if missing_control_cols:
        raise KeyError(f"Control file is missing columns: {missing_control_cols}")
    if missing_summary_cols:
        raise KeyError(f"Summary file is missing columns: {missing_summary_cols}")

    common_ids = control_df.index.intersection(summary_df.index)
    if common_ids.empty:
        raise ValueError("Control file and summary file have no overlapping geography ids")

    adjusted = control_df.copy()
    diagnostics = []

    for idx in sorted(common_ids):
        c_row = control_df.loc[idx]
        s_row = summary_df.loc[idx]
        original_counts = c_row[SIZE_COLS].to_numpy(dtype=float)
        target_population = float(c_row["POPBASE"])
        synthesized_top_bin = inferred_top_bin_size_from_summary(s_row)
        bounded_top_bin = bounded_top_bin_size(synthesized_top_bin, minimum_top_bin_size, maximum_top_bin_size)
        weights = np.array([1, 2, 3, 4, 5, 6, bounded_top_bin], dtype=float)
        if method == "legacy":
            local_rng = np.random.default_rng(random_seed + int(idx) % 1000003)
            rebalanced_counts = rebalance_household_sizes(original_counts, target_population, weights, local_rng)
        elif method == "shape_preserving":
            rebalanced_counts = rebalance_household_sizes_shape_preserving(original_counts, target_population, weights)
        else:
            raise ValueError(f"Unsupported method: {method}")

        adjusted.loc[idx, SIZE_COLS] = rebalanced_counts
        diagnostics.append(
            {
                "BLKGRPID": idx,
                "method": method,
                "target_population": target_population,
                "original_households": float(np.rint(original_counts).sum()),
                "original_implied_top_bin_size": implied_top_bin_size_from_controls(c_row),
                "synthesized_top_bin_size": synthesized_top_bin,
                "bounded_top_bin_size": bounded_top_bin,
                "original_persons_from_bounded_weights": float(np.dot(np.rint(original_counts), weights)),
                "adjusted_persons_from_bounded_weights": float(np.dot(rebalanced_counts, weights)),
                "original_hhpersons7": float(original_counts[-1]),
                "adjusted_hhpersons7": float(rebalanced_counts[-1]),
                "original_hhpersons6": float(original_counts[-2]),
                "adjusted_hhpersons6": float(rebalanced_counts[-2]),
            }
        )

    adjusted[SIZE_COLS] = adjusted[SIZE_COLS].round().astype(int)
    diagnostics_df = pd.DataFrame(diagnostics)
    return adjusted.reset_index(drop=True), diagnostics_df


def main() -> None:
    args = parse_args()
    run_dir, control_file, summary_file, config_options = resolve_inputs(args)

    minimum_top_bin_size = (
        args.min_top_bin_size
        if args.min_top_bin_size is not None
        else float(config_options.get("min_top_bin_size", DEFAULT_MIN_TOP_BIN_SIZE))
    )
    maximum_top_bin_size = (
        args.max_top_bin_size
        if args.max_top_bin_size is not None
        else float(config_options.get("max_top_bin_size", DEFAULT_MAX_TOP_BIN_SIZE))
    )
    random_seed = (
        args.random_seed if args.random_seed is not None else int(config_options.get("random_seed", DEFAULT_RANDOM_SEED))
    )
    suffix = args.suffix or config_options.get("adjusted_suffix", DEFAULT_SUFFIX)
    replace_original = args.replace_original or bool(config_options.get("replace_original", False))
    method = args.method or config_options.get("method", DEFAULT_METHOD)

    control_df = pd.read_csv(control_file)
    summary_df = pd.read_csv(summary_file)
    adjusted_df, diagnostics_df = build_adjusted_controls(
        control_df,
        summary_df,
        minimum_top_bin_size=minimum_top_bin_size,
        maximum_top_bin_size=maximum_top_bin_size,
        random_seed=random_seed,
        method=method,
    )

    adjusted_path = control_file.with_name(f"{control_file.stem}{suffix}{control_file.suffix}")
    diagnostics_path = (
        Path(args.diagnostics_file).resolve()
        if args.diagnostics_file
        else run_dir / "output" / "validation" / f"{control_file.stem}{suffix}_diagnostics.csv"
    )
    diagnostics_path.parent.mkdir(parents=True, exist_ok=True)

    adjusted_df.to_csv(adjusted_path, index=False)
    diagnostics_df.to_csv(diagnostics_path, index=False)

    if replace_original:
        backup_path = control_file.with_name(f"{control_file.stem}_pre_balancer{control_file.suffix}")
        if not backup_path.exists():
            control_file.replace(backup_path)
        else:
            control_file.unlink()
        adjusted_df.to_csv(control_file, index=False)

    print(f"run_dir: {run_dir}")
    print(f"control_file: {control_file}")
    print(f"summary_file: {summary_file}")
    print(f"adjusted_file: {adjusted_path}")
    print(f"diagnostics_file: {diagnostics_path}")
    print(f"min_top_bin_size: {minimum_top_bin_size}")
    print(f"max_top_bin_size: {maximum_top_bin_size}")
    print(f"random_seed: {random_seed}")
    print(f"replace_original: {replace_original}")
    print(f"method: {method}")
    print()
    print("diagnostic summary:")
    print(diagnostics_df[["original_implied_top_bin_size", "synthesized_top_bin_size", "bounded_top_bin_size"]].describe().to_string())


if __name__ == "__main__":
    main()
