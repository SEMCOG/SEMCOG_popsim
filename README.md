# SEMCOG PopulationSim

A simple guide for preparing and running SEMCOG PopulationSim.

## Main Folders
- `input_prep/configs/<run_name>/`: prep configs such as `prepare.yaml`, `controls_pre.csv`, and optional adjustment configs
- `d_drive/popsim/runs/<run_name>/configs/`: generated run configs such as `settings.yaml` and `controls.csv`
- `d_drive/popsim/runs/<run_name>/data/`: generated synthesis inputs
- `d_drive/popsim/runs/<run_name>/output/`: one-pass and two-pass synthesis outputs
- `scripts/`: operational run helpers
- `d_drive/popsim/inputs/`: hand-placed inputs (PUMS, MCD estimate, Census key file) and fetched ACS tables

## Base-Year Synthesis on the SEMCOG MCD Estimate
Households are controlled to SEMCOG's MCD household estimate. Household population is the
RSQE total population minus group quarters, per large area. Attribute patterns come from the
ACS 5-year tables, and the households are copied from the ACS PUMS (the seed).

Run folder: `d_drive/popsim/runs/<run_name>/` (set as `RUN_DIR` at the top of the step scripts).
`d_drive/` is a sibling of this repository.

### Inputs
| Input | Where | Read by |
|---|---|---|
| ACS PUMS household and person files (Michigan) | `d_drive/popsim/inputs/pums/<year>/` | step 0 |
| Census API key (one line) | a file in `d_drive/popsim/inputs/` (never in the repository) | step 0, ACS fetch |
| SEMCOG MCD household estimate workbook | `d_drive/popsim/inputs/<estimate folder>/` | targets (`--estimate`) |
| Base-year HDF with buildings and parcels | `d_drive/forecast_inputs/base_year/` | targets (`--base-hdf`) |
| City table of population, group quarters, households (Detroit by neighborhood) | `d_drive/forecast_inputs/group_quarters/data/` | targets (`CITY_TABLE` in `adjust_to_mcd_2025.py`) |
| RSQE forecast workbook | `d_drive/forecast_inputs/forecast_controls/data/` | targets (`RSQE_WORKBOOK` in `adjust_to_mcd_2025.py`) |

### Steps
**Step 0. ACS run package (once per ACS vintage).** Builds the raw ACS block-group and tract
controls, the PUMS seed, the crosswalk and template configs in `d_drive/popsim/runs/<acs_run>/`:
```bash
python input_prep/scripts/popsim_input_maker.py "$(cat <key file>)" input_prep/configs/<acs_run>/prepare.yaml
```
The step scripts read this package through `ACS_DIR` / `ACS_CONTROLS` and `TEMPLATE` (set at the top of
`adjust_to_mcd_2025.py` and `build_controls_mcd_2025.py`). Keep the template folder.

**ACS fetch (once per ACS vintage).** Block-group occupancy (B25002) and age of householder (B25007):
```bash
CENSUS_API_KEY="$(cat <key file>)" python input_prep/scripts/fetch_acs_b25002_bg.py
```
The output goes to `d_drive/popsim/inputs/<acs folder>/`. The key is read from the environment only.

**Steps 1-4.** `scripts/run_2025_mcd_synthesis.py` runs them in order:

| Step | Script | What it does | Output |
|---|---|---|---|
| targets | `input_prep/scripts/adjust_to_mcd_2025.py` | Turns the MCD estimate into household and household-population targets per block group. Detroit is split into neighborhoods. Vacant units are spread over block group × area pieces by ACS vacancy rate, so households never exceed residential units. | `data/targets/` |
| controls | `input_prep/scripts/build_controls_mcd_2025.py` | ACS shares × targets → block-group and tract controls. Block-group shares are shrunk toward the tract share to reduce ACS noise; household size is refit to the population target. Writes the run's seed copy and configs. | `data/`, `configs/` |
| synthesis | `scripts/run_two_pass_hhsize.py` | PopulationSim pass 1, household-size balancer, pass 2 | `output/<stamp>_two_pass/pass1/`, `pass2/` |
| reconcile | `scripts/reconcile_bg_households.py` | Makes the household count exact in every block group | `output/<stamp>_two_pass/final/` |

Notes on the controls step:
- Each seed expression must count the same persons as its ACS table. Industry counts employed
  persons; workers count persons at work. The builder sets these expressions and adds the seed
  columns they need to the run's seed copy (the template seed is not changed).
- Tracts with large group quarters (for example, campus dorms) get worker and industry controls
  built from household worker counts, because ACS person tables include group-quarters residents.
- With `SEED_REWEIGHT` (top of `build_controls_mcd_2025.py`), the run's seed copy is reweighted so that
  employed and unemployed persons by age follow the recent survey years (`RECENT_YEARS`), and the tract worker
  controls are rebuilt from it. Diagnostics: `data/targets/seed_reweight_by_puma.csv`.

### Environment (one time)
OR-Tools, the fast integerizer, does not load in the base conda env (pip `ortools` and conda-forge
`pyarrow` bundle different abseil libraries). Use the isolated `popsim` env; the base env, which
urbansim uses, stays unchanged.
```bash
conda create -y -n popsim python=3.12 pip
PYTHONNOUSERSITE=1 /opt/conda/envs/popsim/bin/python -m pip install "populationsim==0.10.0" oyaml pyyaml openpyxl
```
Always set `PYTHONNOUSERSITE=1`, so packages in `~/.local` do not load into the env.

### Run
```bash
cd SEMCOG_popsim
PYTHONNOUSERSITE=1 /opt/conda/envs/popsim/bin/python scripts/run_2025_mcd_synthesis.py \
  --estimate <estimate workbook> --base-hdf <base-year HDF>
```
- `--from-step` / `--to-step` run part of the pipeline (`targets`, `controls`, `synthesis`,
  `reconcile`). A control change needs only `--from-step controls` (no `--estimate`, `--base-hdf`).
- A full run takes about 2 hours. In the background:
  `PYTHONNOUSERSITE=1 nohup /opt/conda/envs/popsim/bin/python scripts/run_2025_mcd_synthesis.py ... > <log> 2>&1 &`
- To stop a run, kill its PIDs. Do not use `pkill -f`, which can kill the calling shell.
- `<stamp>` in the output folder is the start date and hour in UTC; log timestamps are UTC too.

### Run-folder structure
```text
d_drive/popsim/runs/<run_name>/
├── configs/                  controls.csv, settings.yaml
├── data/
│   ├── targets/              block-group targets, block group × area pieces (input to placement),
│   │                         area targets, build diagnostics
│   ├── *_control_totals_*.csv   block-group and tract controls; size-adjusted block-group controls
│   └── *_seed_*.csv, *_geo_cross_walk.csv   the run's seed copy and crosswalk
└── output/<stamp>_two_pass/
    ├── configs/pass1/, pass2/   settings and controls of each pass
    ├── logs/                 pass1, balancer, pass2 logs and status files
    ├── pass1/, pass2/        PopulationSim outputs and summaries
    ├── final/                synthetic households and persons  <- the deliverable
    ├── inputs_snapshot/      control files and configs used by this run
    └── validation/           size-balancer review copy, charts/
```
Use `final/`, not `pass2/`.

### Charts
```bash
python scripts/plot_synthesis_results.py \
  --summary-dir <run>/output/<stamp>_two_pass/pass2 \
  --out-dir     <run>/output/<stamp>_two_pass/validation/charts \
  --final-dir   <run>/output/<stamp>_two_pass/final --seed-dir <run>/data \
  [--compare-dir <older run>/pass2 --compare-label "<label>"] \
  [--acs-blkgrp <acs_run>/data/<ACS block-group controls>]
```
Use `--final-dir`, so results are counted from the output population. PopulationSim's tract summary
reports the tract-level weights, which differ a little from the output.

### Integerizer
- Use `USE_CVXPY: false` (OR-Tools); the controls step writes it. With `true`, PopulationSim uses
  GLPK, which is much slower and has no working time limit.
- A few very small zones are INFEASIBLE in the integerizer, and PopulationSim falls back to smart
  rounding. This is expected; the totals stay exact.

### After the synthesis
1. Place the households into buildings (`final/` and the block group × area pieces in `data/targets/`).
2. Build the base-year HDF.
3. Re-run the household controls with the new HDF (`forecast_data_input`).
4. Realign the employed worker rates to the new base.

`input_prep/scripts/adjust_to_large_area_2025.py` is superseded: it scaled ACS controls by one ratio
per large area and kept the ACS household size. Do not use it for the base year.

## Other Workflows
The sections below describe the general run package (used as step 0 above), plain one-pass and
two-pass runs, and the refinement inputs.

## 1. Prepare a Run Package
Create the PopulationSim run package from Census and PUMS inputs.

```bash
python input_prep/scripts/popsim_input_maker.py <census_key> input_prep/configs/<run_name>/prepare.yaml
```

This produces a run package under:
```text
d_drive/popsim/runs/<run_name>/
├── configs/
│   ├── settings.yaml
│   └── controls.csv
├── data/
└── output/
```

## 2. Run One-Pass Synthesis
Run the standard one-pass synthesis.

Using the helper script (set `RUN_DIR` for another run package):
```bash
RUN_DIR=../d_drive/popsim/runs/<run_name> bash scripts/run_2024_synthesis.sh
```

Or directly with PopulationSim:
```bash
populationsim \
  -c ../d_drive/popsim/runs/<run_name>/configs \
  -d ../d_drive/popsim/runs/<run_name>/data \
  -o ../d_drive/popsim/runs/<run_name>/output/$(date +%Y-%m-%d_%H)_one_pass/run
```

One-pass output structure:
```text
output/
└── <YYYY-MM-DD>_<HH>_one_pass/
    ├── logs/
    │   ├── run.log
    │   └── run.status
    ├── run/
    │   ├── mem.csv
    │   ├── pipeline.h5
    │   ├── timing_log.csv
    │   ├── summary_*.csv
    │   ├── final_summary_*.csv
    │   ├── synthetic_households.csv
    │   └── synthetic_persons.csv
    └── validation/
```

Notes:
- rerunning within the same hour reuses the same one-pass folder
- overwriting within the same hour is expected

## 3. Run Two-Pass Synthesis
Run the two-pass workflow with household-size balancing between passes.

```bash
python scripts/run_two_pass_hhsize.py --run-name <run_name>
```

Useful options:
```bash
python scripts/run_two_pass_hhsize.py --run-name <run_name> --method legacy
python scripts/run_two_pass_hhsize.py --run-dir ../d_drive/popsim/runs/<run_name>
```

Two-pass output structure:
```text
output/
└── <YYYY-MM-DD>_<HH>_two_pass/
    ├── logs/
    │   ├── workflow.log
    │   ├── workflow.status
    │   ├── pass1.log
    │   ├── pass1.status
    │   ├── balancer.log
    │   ├── balancer.status
    │   ├── pass2.log
    │   └── pass2.status
    ├── configs/
    │   ├── pass1/
    │   │   └── settings.yaml
    │   └── pass2/
    │       └── settings.yaml
    ├── pass1/
    ├── pass2/
    └── validation/
```

Notes:
- `pass1/` is the first synthesis run
- `configs/pass1/settings.yaml` points to the base block-group control file
- the household-size balancer writes an adjusted `_hhsize_adj` control file without overwriting or restoring the base control file
- `configs/pass2/settings.yaml` points to the adjusted control file
- `pass2/` is the final synthesis output
- rerunning within the same hour reuses the same two-pass folder

## 4. Run in Background
If you want the run to continue while you keep using the terminal:

One-pass:
```bash
nohup bash scripts/run_2024_synthesis.sh >/dev/null 2>&1 &
```

Two-pass:
```bash
nohup python scripts/run_two_pass_hhsize.py --run-name <run_name> >/dev/null 2>&1 &
```

Check whether a PID is still running:
```bash
ps -fp <PID>
```

## 5. Validate a Run
Validate a standard run:

```bash
python scripts/validate_popsim_run.py --run-name <run_name> --write-csv
```

Validate a two-pass final run:
```bash
python scripts/validate_popsim_run.py --output-dir ../d_drive/popsim/runs/<run_name>/output/<YYYY-MM-DD>_<HH>_two_pass/pass2 --configs-dir ../d_drive/popsim/runs/<run_name>/configs --write-csv
```

Validation writes CSV summaries into a `validation/` folder under the target output folder.

## 6. Generate Refinement Inputs
Use the refinement input generator when you want to build a new synthesis package from existing UrbanSim outputs.

This workflow creates the package inputs only. It does not run the refinement synthesis itself.

```bash
python scripts/generate_refinement_inputs.py refinement/configs/template_refinement.yaml
```

The YAML config controls:
- project name
- year
- target geography
- sample geography
- input and target HDF paths
- output root
- optional county filter, weight column, and person-control toggle; omit county to keep all counties

Each run writes one target geography package under:
```text
d_drive/popsim/runs/<project_name>_<target_geo_lower>/
├── configs/
│   ├── controls.csv
│   ├── settings.yaml
│   └── refinement_input_config.yaml
├── data/
│   ├── <project_name>_geo_cross_walk.csv
│   ├── <project_name>_control_totals_<TARGET_GEO>.csv
│   ├── <project_name>_seed_households.csv
│   └── <project_name>_seed_persons.csv
└── output/
```

Notes:
- one target geography per run
- review and edit `data/<project_name>_control_totals_<TARGET_GEO>.csv`, then replace that file in-place before the later refinement/synthesis run
- `refinement/urbansim_refine_input.ipynb` remains reference material, but `scripts/generate_refinement_inputs.py` is the maintained entrypoint
