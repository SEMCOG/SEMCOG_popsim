#!/usr/bin/env python3
"""Run a two-pass PopulationSim workflow with household-size rebalancing.

The workflow keeps the base input package immutable:
1. Pass 1 uses a copied settings folder that points to the base BLKGRP controls.
2. The household-size balancer writes a separate *_hhsize_adj.csv control file.
3. Pass 2 uses a copied settings folder that points to that adjusted control file.

Outputs are written under output/<YYYY-MM-DD>_<HH>_two_pass/ so each run has
its own logs, pass outputs, copied configs, and validation artifacts.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

try:
    import oyaml as yaml
except ModuleNotFoundError:
    import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]
INPUT_PREP_DIR = REPO_ROOT / "input_prep"
BALANCER_SCRIPT = INPUT_PREP_DIR / "scripts" / "hh_size_balancer.py"
DEFAULT_RUN_NAME = "2024_synthesis"
DEFAULT_CONFIG = INPUT_PREP_DIR / "configs" / DEFAULT_RUN_NAME / "prepare.yaml"
DEFAULT_METHOD = "shape_preserving"
DEFAULT_RUN_ROOT = REPO_ROOT.parent / "d_drive" / "popsim" / "runs"


def resolve_config_path(path_str: str | None, config_dir: Path) -> Path | None:
    """Resolve a config path relative to the prepare config, input_prep, or repo root."""
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


def default_config_for_run_name(run_name: str) -> Path:
    """Return the default prepare.yaml path for a named run."""
    return INPUT_PREP_DIR / "configs" / run_name / "prepare.yaml"


def default_run_dir_for_run_name(run_name: str) -> Path:
    """Return the default generated PopulationSim package directory for a named run."""
    return DEFAULT_RUN_ROOT / run_name


def load_yaml_config(config_path: Path) -> dict[str, Any]:
    """Load a YAML config while preserving the repo's oyaml preference when available."""
    with open(config_path, "r") as stream:
        return yaml.load(stream, Loader=yaml.Loader)


def derive_run_dir(config_path: Path, conf: dict[str, Any]) -> Path:
    """Derive the generated run package path from prepare.yaml."""
    config_dir = config_path.parent
    run_name = conf["project"].get("run_name", f"{conf['project']['acs_year']}_synthesis")
    run_root = resolve_config_path(conf.get("paths", {}).get("run_root"), config_dir)
    if run_root is None:
        raise ValueError("Could not derive run_root from config")
    return run_root / run_name


def find_blockgroup_control_file(data_dir: Path) -> Path:
    """Find the base BLKGRP control file, excluding adjusted and backup copies."""
    candidates = []
    for path in sorted(data_dir.glob("*_control_totals_blkgrp.csv")):
        name = path.name
        if "_hhsize_" in name or name.endswith("_pre_balancer.csv"):
            continue
        candidates.append(path)
    if not candidates:
        raise FileNotFoundError(f"No base block-group control file found in {data_dir}")
    if len(candidates) > 1:
        raise FileExistsError(f"Expected one base block-group control file in {data_dir}, found {candidates}")
    return candidates[0]


def adjusted_path_for(control_file: Path, suffix: str = "_hhsize_adj") -> Path:
    """Return the path where the balancer should write adjusted BLKGRP controls."""
    return control_file.with_name(f"{control_file.stem}{suffix}{control_file.suffix}")


def copy_config_dir(source_config_dir: Path, destination_config_dir: Path) -> None:
    """Copy the run config directory so pass-specific settings can be edited safely."""
    destination_config_dir.mkdir(parents=True, exist_ok=True)
    shutil.copytree(source_config_dir, destination_config_dir, dirs_exist_ok=True)


def update_blkgrp_control_filename(settings_path: Path, control_filename: str) -> None:
    """Point BLKGRP_control_data in a copied settings.yaml to a chosen control file."""
    settings = load_yaml_config(settings_path)
    for table_item in settings.get("input_table_list", []):
        if table_item.get("tablename") == "BLKGRP_control_data":
            table_item["filename"] = control_filename
            break
    else:
        raise KeyError(f"Could not find BLKGRP_control_data in {settings_path}")

    with open(settings_path, "w") as stream:
        yaml.dump(settings, stream, default_flow_style=False, sort_keys=False)


def resolve_summary_file(pass1_output: Path) -> Path:
    """Find the BLKGRP summary produced by pass 1 for the household-size balancer."""
    for name in ["final_summary_BLKGRP.csv", "summary_BLKGRP.csv"]:
        candidate = pass1_output / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Could not find BLKGRP summary in {pass1_output}")


def write_status(path: Path, status: str, exit_code: int = 0) -> None:
    """Write a small status file for workflow and step monitoring."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"status={status}\nexit_code={exit_code}\n")


def run_stamp(dt: datetime | None = None) -> str:
    """Return the hourly run stamp used in output folder names."""
    if dt is None:
        dt = datetime.now().astimezone()
    return dt.strftime("%Y-%m-%d_%H")


def log(message: str, workflow_log: Path) -> None:
    """Write a workflow message to the terminal and workflow log."""
    line = f"[{datetime.now().astimezone().strftime('%Y-%m-%d %H:%M:%S')}] {message}"
    print(line)
    with open(workflow_log, "a") as stream:
        stream.write(line + "\n")


def run_command(cmd: list[str], log_path: Path, status_path: Path, workflow_log: Path, cwd: Path | None = None) -> int:
    """Run a subprocess, stream output to a log, and write a step status file."""
    log(f"Running: {' '.join(cmd)}", workflow_log)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w") as log_stream:
        process = subprocess.Popen(
            cmd,
            cwd=str(cwd) if cwd else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            log_stream.write(line)
        process.wait()
    if process.returncode == 0:
        write_status(status_path, "success", 0)
    else:
        write_status(status_path, "failed", process.returncode)
    return int(process.returncode)


def parse_args() -> argparse.Namespace:
    """Parse the minimal CLI for choosing the run package and balancer method."""
    parser = argparse.ArgumentParser(
        description="Run PopulationSim in a two-pass workflow with household-size rebalancing between passes."
    )
    target_group = parser.add_mutually_exclusive_group()
    target_group.add_argument(
        "--run-name",
        default=DEFAULT_RUN_NAME,
        help="Short run name under input_prep/configs/ and d_drive/popsim/runs/ (default: 2024_synthesis)",
    )
    target_group.add_argument(
        "--config",
        help="prepare.yaml used to derive the run folder and balancer defaults",
    )
    target_group.add_argument(
        "--run-dir",
        help="PopulationSim run folder; overrides config and run-name lookup",
    )
    parser.add_argument(
        "--method",
        default=DEFAULT_METHOD,
        choices=["shape_preserving", "legacy"],
        help="Household-size balancer method",
    )
    return parser.parse_args()


def resolve_run_context(args: argparse.Namespace) -> tuple[Path | None, dict[str, Any], Path]:
    """Resolve config metadata and the generated run directory from CLI arguments."""
    if args.run_dir:
        return None, {}, Path(args.run_dir).resolve()
    if args.config:
        config_path = Path(args.config).resolve()
        conf = load_yaml_config(config_path)
        return config_path, conf, derive_run_dir(config_path, conf)

    config_path = default_config_for_run_name(args.run_name)
    if config_path.exists():
        conf = load_yaml_config(config_path)
        return config_path, conf, derive_run_dir(config_path, conf)
    return None, {}, default_run_dir_for_run_name(args.run_name)


def main() -> int:
    """Coordinate pass 1, household-size balancing, and pass 2."""
    args = parse_args()
    config_path, config, run_dir = resolve_run_context(args)

    config_dir = run_dir / "configs"
    data_dir = run_dir / "data"
    output_root = run_dir / "output"
    control_file = find_blockgroup_control_file(data_dir)

    config_options = config.get("postprocess", {}).get("hh_size_balancer", {}) if config else {}
    adjusted_suffix = config_options.get("adjusted_suffix", "_hhsize_adj")
    adjusted_control_file = adjusted_path_for(control_file, adjusted_suffix)

    workflow_dir = output_root / f"{run_stamp()}_two_pass"
    pass1_output = workflow_dir / "pass1"
    pass2_output = workflow_dir / "pass2"
    logs_dir = workflow_dir / "logs"
    workflow_log = logs_dir / "workflow.log"
    overall_status = logs_dir / "workflow.status"
    validation_dir = workflow_dir / "validation"
    workflow_config_dir = workflow_dir / "configs"
    pass1_config_dir = workflow_config_dir / "pass1"
    pass2_config_dir = workflow_config_dir / "pass2"
    balancer_diagnostics = validation_dir / f"{control_file.stem}_{args.method}_diagnostics.csv"

    # Pass 1 and pass 2 get independent config copies. Only pass 2 is edited to
    # use the adjusted BLKGRP controls produced after pass 1.
    logs_dir.mkdir(parents=True, exist_ok=True)
    validation_dir.mkdir(parents=True, exist_ok=True)
    pass1_output.mkdir(parents=True, exist_ok=True)
    pass2_output.mkdir(parents=True, exist_ok=True)
    copy_config_dir(config_dir, pass1_config_dir)
    copy_config_dir(config_dir, pass2_config_dir)
    update_blkgrp_control_filename(pass2_config_dir / "settings.yaml", adjusted_control_file.name)
    workflow_log.write_text("")

    log(f"Two-pass workflow started for {run_dir}", workflow_log)
    log(f"Workflow dir: {workflow_dir}", workflow_log)
    log(f"Control file: {control_file}", workflow_log)
    log(f"Adjusted control file: {adjusted_control_file}", workflow_log)
    log(f"Pass 1 config dir: {pass1_config_dir}", workflow_log)
    log(f"Pass 2 config dir: {pass2_config_dir}", workflow_log)
    log(f"Balancer method: {args.method}", workflow_log)

    log("Base control file is left unchanged; pass 2 uses a separate adjusted control file.", workflow_log)

    pass1_status = run_command(
        ["populationsim", "-c", str(pass1_config_dir), "-d", str(data_dir), "-o", str(pass1_output)],
        logs_dir / "pass1.log",
        logs_dir / "pass1.status",
        workflow_log,
    )
    if pass1_status != 0:
        log(f"Pass 1 failed with exit code {pass1_status}", workflow_log)
        write_status(overall_status, "failed_pass1", pass1_status)
        return pass1_status

    summary_file = resolve_summary_file(pass1_output)
    balancer_cmd = [
        sys.executable,
        str(BALANCER_SCRIPT),
        "--run-dir",
        str(run_dir),
        "--summary-file",
        str(summary_file),
        "--method",
        args.method,
        "--suffix",
        adjusted_suffix,
        "--diagnostics-file",
        str(balancer_diagnostics),
    ]
    if config_path is not None:
        balancer_cmd.extend(["--config", str(config_path)])

    balancer_status = run_command(
        balancer_cmd,
        logs_dir / "balancer.log",
        logs_dir / "balancer.status",
        workflow_log,
    )
    if balancer_status != 0:
        log(f"Household-size balancer failed with exit code {balancer_status}", workflow_log)
        write_status(overall_status, "failed_balancer", balancer_status)
        return balancer_status

    if not adjusted_control_file.exists():
        raise FileNotFoundError(f"Expected adjusted control file was not created: {adjusted_control_file}")

    review_adjusted_control = validation_dir / adjusted_control_file.name
    shutil.copy2(adjusted_control_file, review_adjusted_control)
    log(f"Saved adjusted control review copy: {review_adjusted_control}", workflow_log)

    pass2_status = run_command(
        ["populationsim", "-c", str(pass2_config_dir), "-d", str(data_dir), "-o", str(pass2_output)],
        logs_dir / "pass2.log",
        logs_dir / "pass2.status",
        workflow_log,
    )
    if pass2_status != 0:
        log(f"Pass 2 failed with exit code {pass2_status}", workflow_log)
        write_status(overall_status, "failed_pass2", pass2_status)
        return pass2_status

    log("Two-pass workflow completed successfully.", workflow_log)
    log(f"Pass 1 output: {pass1_output}", workflow_log)
    log(f"Pass 2 output: {pass2_output}", workflow_log)
    log(f"Balancer diagnostics: {balancer_diagnostics}", workflow_log)
    log(f"Adjusted control review copy: {review_adjusted_control}", workflow_log)
    write_status(overall_status, "success", 0)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
