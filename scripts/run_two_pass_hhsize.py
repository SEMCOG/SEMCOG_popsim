#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from datetime import UTC, datetime
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
    return INPUT_PREP_DIR / "configs" / run_name / "prepare.yaml"


def default_run_dir_for_run_name(run_name: str) -> Path:
    return DEFAULT_RUN_ROOT / run_name


def load_yaml_config(config_path: Path) -> dict[str, Any]:
    with open(config_path, "r") as stream:
        return yaml.load(stream, Loader=yaml.Loader)


def derive_run_dir(config_path: Path, conf: dict[str, Any]) -> Path:
    config_dir = config_path.parent
    run_name = conf["project"].get("run_name", f"{conf['project']['acs_year']}_synthesis")
    run_root = resolve_config_path(conf.get("paths", {}).get("run_root"), config_dir)
    if run_root is None:
        raise ValueError("Could not derive run_root from config")
    return run_root / run_name


def find_blockgroup_control_file(data_dir: Path) -> Path:
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


def backup_path_for(control_file: Path) -> Path:
    return control_file.with_name(f"{control_file.stem}_pre_balancer{control_file.suffix}")


def resolve_summary_file(pass1_output: Path) -> Path:
    for name in ["final_summary_BLKGRP.csv", "summary_BLKGRP.csv"]:
        candidate = pass1_output / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Could not find BLKGRP summary in {pass1_output}")


def write_status(path: Path, status: str, exit_code: int = 0) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"status={status}\nexit_code={exit_code}\n")


def archive_existing_two_pass_outputs(workflow_dir: Path) -> Path | None:
    tracked = ["pass1", "pass2", "logs", "validation"]
    existing = [workflow_dir / name for name in tracked if (workflow_dir / name).exists()]
    if not existing:
        return None

    run_index = 1
    while (workflow_dir / f"run_{run_index}").exists():
        run_index += 1
    archive_dir = workflow_dir / f"run_{run_index}"
    archive_dir.mkdir(parents=True, exist_ok=False)

    for source in existing:
        shutil.move(str(source), str(archive_dir / source.name))

    return archive_dir


def log(message: str, workflow_log: Path) -> None:
    line = f"[{datetime.now(UTC).strftime('%Y-%m-%d %H:%M:%S')}] {message}"
    print(line)
    with open(workflow_log, "a") as stream:
        stream.write(line + "\n")


def run_command(cmd: list[str], log_path: Path, status_path: Path, workflow_log: Path, cwd: Path | None = None) -> int:
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
    parser.add_argument(
        "--skip-restore-original",
        action="store_true",
        help="Do not restore the base control file from *_pre_balancer before pass 1",
    )
    return parser.parse_args()


def resolve_run_context(args: argparse.Namespace) -> tuple[Path | None, dict[str, Any], Path]:
    if args.run_dir:
        return None, {}, Path(args.run_dir).resolve()
    if args.config:
        config_path = Path(args.config).resolve()
        return config_path, load_yaml_config(config_path), derive_run_dir(config_path, load_yaml_config(config_path))

    config_path = default_config_for_run_name(args.run_name)
    if config_path.exists():
        conf = load_yaml_config(config_path)
        return config_path, conf, derive_run_dir(config_path, conf)
    return None, {}, default_run_dir_for_run_name(args.run_name)


def main() -> int:
    args = parse_args()
    config_path, config, run_dir = resolve_run_context(args)

    config_dir = run_dir / "configs"
    data_dir = run_dir / "data"
    output_dir = run_dir / "output"
    control_file = find_blockgroup_control_file(data_dir)
    backup_file = backup_path_for(control_file)

    workflow_dir = output_dir / "two_pass"
    archived_run_dir = archive_existing_two_pass_outputs(workflow_dir)
    pass1_output = workflow_dir / "pass1"
    pass2_output = workflow_dir / "pass2"
    logs_dir = workflow_dir / "logs"
    workflow_log = logs_dir / "workflow.log"
    overall_status = logs_dir / "workflow.status"
    validation_dir = workflow_dir / "validation"
    balancer_diagnostics = validation_dir / f"{control_file.stem}_{args.method}_diagnostics.csv"

    logs_dir.mkdir(parents=True, exist_ok=True)
    validation_dir.mkdir(parents=True, exist_ok=True)
    pass1_output.mkdir(parents=True, exist_ok=True)
    pass2_output.mkdir(parents=True, exist_ok=True)
    workflow_log.write_text("")

    log(f"Two-pass workflow started for {run_dir}", workflow_log)
    if archived_run_dir is not None:
        log(f"Archived previous two-pass outputs to {archived_run_dir}", workflow_log)
    log(f"Control file: {control_file}", workflow_log)
    log(f"Balancer method: {args.method}", workflow_log)

    if backup_file.exists() and not args.skip_restore_original:
        shutil.copy2(backup_file, control_file)
        log(f"Restored original control file from backup: {backup_file}", workflow_log)
    elif not backup_file.exists():
        log("No pre-balancer backup found; pass 1 will use the current control file.", workflow_log)

    pass1_status = run_command(
        ["populationsim", "-c", str(config_dir), "-d", str(data_dir), "-o", str(pass1_output)],
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
        "--replace-original",
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

    review_adjusted_control = validation_dir / control_file.name
    shutil.copy2(control_file, review_adjusted_control)
    log(f"Saved adjusted control review copy: {review_adjusted_control}", workflow_log)

    if backup_file.exists():
        review_original_control = validation_dir / backup_file.name
        shutil.copy2(backup_file, review_original_control)
        log(f"Saved pre-balancer control review copy: {review_original_control}", workflow_log)

    pass2_status = run_command(
        ["populationsim", "-c", str(config_dir), "-d", str(data_dir), "-o", str(pass2_output)],
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
