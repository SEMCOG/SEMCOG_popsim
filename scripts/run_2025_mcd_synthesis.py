# Run the 2025 synthesis on the SEMCOG MCD household estimate, end to end.
#
# Steps (each is its own script; this runner calls them in order):
#   targets    input_prep/scripts/adjust_to_mcd_2025.py
#                MCD estimate -> BG HHBASE / POPBASE targets (data/targets/)
#   controls   input_prep/scripts/build_controls_mcd_2025.py
#                targets + ACS 2020-2024 -> BG / tract controls, configs, seed copies
#   synthesis  scripts/run_two_pass_hhsize.py
#                PopulationSim pass 1 -> size balancer -> pass 2 (output/<stamp>_two_pass/)
#   reconcile  scripts/reconcile_bg_households.py
#                exact HHBASE in every BG -> output/<stamp>_two_pass/final/
#
# Run it in the isolated `popsim` conda env (OR-Tools does not load in the base
# env). One-time input: ACS B25002 by BG (input_prep/scripts/fetch_acs_b25002_bg.py).
#
# Usage:
#     PYTHONNOUSERSITE=1 /opt/conda/envs/popsim/bin/python scripts/run_2025_mcd_synthesis.py \
#         --estimate .../July1_2025_Population_revised.xlsx --base-hdf .../main_100226.h5

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RUN_DIR = REPO.parent / "d_drive" / "popsim" / "runs" / "2025_synthesis_mcd"
ACS_B25002 = REPO.parent / "d_drive" / "popsim" / "inputs" / "acs_2024_bg" / "acs2024_5yr_B25002_bg.csv"
STEPS = ["targets", "controls", "synthesis", "reconcile"]


def run(cmd, env):
    print("\n>>> " + " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], cwd=REPO, env=env, check=True)


def check_env():
    """OR-Tools must load after pandas/pyarrow (true in the popsim env, not in base)."""
    probe = "import pandas; from ortools.linear_solver import pywraplp; pywraplp.Solver('t', pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)"
    if subprocess.run([sys.executable, "-c", probe], capture_output=True).returncode != 0:
        sys.exit("OR-Tools does not load in this Python (%s). Use /opt/conda/envs/popsim/bin/python." % sys.executable)


def latest_two_pass():
    runs = sorted((RUN_DIR / "output").glob("*_two_pass"), key=lambda p: p.stat().st_mtime)
    if not runs:
        sys.exit("no *_two_pass output under %s" % (RUN_DIR / "output"))
    return runs[-1]


def main():
    ap = argparse.ArgumentParser(description="2025 synthesis on the SEMCOG MCD estimate")
    ap.add_argument("--estimate", type=Path, help="SEMCOG MCD estimate workbook (step targets)")
    ap.add_argument("--base-hdf", type=Path, help="base-year HDF with buildings/parcels (step targets)")
    ap.add_argument("--from-step", choices=STEPS, default="targets", help="first step to run")
    ap.add_argument("--to-step", choices=STEPS, default="reconcile", help="last step to run")
    args = ap.parse_args()
    steps = STEPS[STEPS.index(args.from_step): STEPS.index(args.to_step) + 1]

    env = dict(os.environ, PYTHONNOUSERSITE="1")
    env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env.get("PATH", "")  # `populationsim` CLI
    py = sys.executable

    if "targets" in steps:
        if not (args.estimate and args.base_hdf):
            sys.exit("step targets needs --estimate and --base-hdf")
        if not ACS_B25002.exists():
            sys.exit("missing %s\nrun: CENSUS_API_KEY=... python input_prep/scripts/fetch_acs_b25002_bg.py" % ACS_B25002)
        run([py, "input_prep/scripts/adjust_to_mcd_2025.py", "--estimate", args.estimate, "--base-hdf", args.base_hdf], env)
    if "controls" in steps:
        run([py, "input_prep/scripts/build_controls_mcd_2025.py"], env)
    if "synthesis" in steps:
        check_env()
        run([py, "scripts/run_two_pass_hhsize.py", "--run-dir", RUN_DIR], env)
        # keep the exact inputs of this run next to its output
        out = latest_two_pass()
        snap = out / "inputs_snapshot"
        snap.mkdir(exist_ok=True)
        for f in list((RUN_DIR / "data").glob("SEMCOG_2025_control_totals_*.csv")) + \
                [RUN_DIR / "configs" / "controls.csv", RUN_DIR / "configs" / "settings.yaml",
                 RUN_DIR / "data" / "targets" / "controls_build_diagnostics.csv"]:
            shutil.copy2(f, snap / f.name)
    if "reconcile" in steps:
        out = latest_two_pass()
        run([py, "scripts/reconcile_bg_households.py", "--pass-dir", out / "pass2",
             "--controls", RUN_DIR / "data" / "SEMCOG_2025_control_totals_blkgrp_hhsize_adj.csv"], env)
        print("\nfinal synthetic population: %s" % (out / "final"))


if __name__ == "__main__":
    main()
