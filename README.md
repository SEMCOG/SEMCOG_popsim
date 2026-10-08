# SEMCOG PopulationSim

A simple guide for preparing and running SEMCOG PopulationSim.

## Main Folders
- `input_prep/configs/<run_name>/`: prep configs such as `prepare.yaml`, `controls_pre.csv`, and optional adjustment configs
- `d_drive/popsim/runs/<run_name>/configs/`: generated run configs such as `settings.yaml` and `controls.csv`
- `d_drive/popsim/runs/<run_name>/data/`: generated synthesis inputs
- `d_drive/popsim/runs/<run_name>/output/`: one-pass and two-pass synthesis outputs
- `scripts/`: operational run helpers

## Current Workflow: 2025 Base-Year Synthesis on the SEMCOG MCD Estimate
The 2025 base-year households are controlled to SEMCOG's July 1, 2025 MCD household
estimate. Household population is the RSQE 2025 total minus group quarters, per large
area. Attribute patterns come from ACS 2020-2024. Run folder:
`d_drive/popsim/runs/2025_synthesis_mcd/`.

### Pipeline
```text
fetch_acs_b25002_bg.py        (one time)  ACS B25002 occupancy by block group
        |
adjust_to_mcd_2025.py         step "targets"    MCD estimate -> BG HHBASE / POPBASE
        |                                       (data/targets/)
build_controls_mcd_2025.py    step "controls"   ACS shares x targets -> BG / tract controls,
        |                                       configs/, seed copies
run_two_pass_hhsize.py        step "synthesis"  PopulationSim pass 1 -> hh_size_balancer.py
        |                                       -> pass 2   (output/<stamp>_two_pass/)
reconcile_bg_households.py    step "reconcile"  exact HHBASE in every BG
                                                (output/<stamp>_two_pass/final/)
```
`scripts/run_2025_mcd_synthesis.py` runs the four steps in order.

| Step | Script | Main inputs | Main outputs |
|---|---|---|---|
| targets | `input_prep/scripts/adjust_to_mcd_2025.py` | MCD estimate workbook, base-year HDF (buildings, parcels), city GQ table, RSQE workbook, ACS B25002 | `data/targets/bg_targets_2025.csv`, `bg_area_pieces_2025.csv`, `area_targets_2025.csv` |
| controls | `input_prep/scripts/build_controls_mcd_2025.py` | targets, ACS 2020-2024 controls (`runs/2024_synthesis/data/`), PUMS seed | `data/SEMCOG_2025_control_totals_{blkgrp,tract}.csv`, `configs/`, `data/targets/controls_build_diagnostics.csv` |
| synthesis | `scripts/run_two_pass_hhsize.py` | run folder | `output/<stamp>_two_pass/pass1/`, `pass2/`, logs |
| reconcile | `scripts/reconcile_bg_households.py` | pass 2 output, size-adjusted BG controls | `output/<stamp>_two_pass/final/` (households, persons, log) |

Method notes:
- **targets:** Detroit is split into its 55 neighborhoods. Each area's vacant units are spread
  over its BG x area building pieces by the ACS vacancy rate (shrunk toward tract and county),
  so households never exceed residential units.
- **controls:** each BG category share is shrunk toward its tract share,
  `(count + K x tract share) / (BG base + K)`. K is estimated per control group from the data
  (ACS sampling noise vs real BG variation) and printed in the log. Household size is refit
  to POPBASE, not scaled.
- **reconcile:** it changes only BGs that miss HHBASE. It writes a new `final/` folder and does
  not change the PopulationSim output.

### Environment (one time)
OR-Tools, the fast integerizer, does not load in the base conda env, because pip `ortools`
and conda-forge `pyarrow` bundle different abseil libraries. Use the isolated `popsim` env.
The base env, which urbansim uses, stays unchanged.

```bash
conda create -y -n popsim python=3.12 pip
PYTHONNOUSERSITE=1 /opt/conda/envs/popsim/bin/python -m pip install "populationsim==0.10.0" oyaml pyyaml openpyxl
```

Always set `PYTHONNOUSERSITE=1`, so packages in `~/.local` do not load into the env.

### ACS input (one time)
```bash
CENSUS_API_KEY=<key> python input_prep/scripts/fetch_acs_b25002_bg.py
```
The key is read from the environment and is not written to disk.

### Run
```bash
cd SEMCOG_popsim
PYTHONNOUSERSITE=1 /opt/conda/envs/popsim/bin/python scripts/run_2025_mcd_synthesis.py \
  --estimate ../d_drive/popsim/inputs/2025_semcog_estimate/July1_2025_Population_revised.xlsx \
  --base-hdf ../d_drive/forecast_inputs/base_year/main_100226.h5
```

Useful options:
- `--from-step controls` or `--to-step controls`: run part of the pipeline. Steps:
  `targets`, `controls`, `synthesis`, `reconcile`.
- A full run takes about 2 hours. In the background:
  `PYTHONNOUSERSITE=1 nohup /opt/conda/envs/popsim/bin/python scripts/run_2025_mcd_synthesis.py ... > run.log 2>&1 &`

The final synthetic population is in `output/<stamp>_two_pass/final/`. The exact inputs of
each run are copied to `output/<stamp>_two_pass/inputs_snapshot/`.

### Integerizer notes
- Use `USE_CVXPY: false` (OR-Tools); `build_controls_mcd_2025.py` writes it. With
  `USE_CVXPY: true`, PopulationSim uses GLPK through CVXPY, which is about 50-90 times slower
  and has no working time limit (one solve cycled for 7 hours).
- A few one-household tracts are always INFEASIBLE in the integerizer, and PopulationSim falls
  back to smart rounding. This is expected and does not change the totals.

### Superseded
`input_prep/scripts/adjust_to_large_area_2025.py` (July 2026) scaled ACS controls by one
household-population ratio per large area. It kept the ACS household size, so it produced
about 57,000 fewer households than the MCD estimate. Do not use it for the 2025 base year.

## 1. Prepare a Run Package
Create the PopulationSim run package from Census and PUMS inputs.

```bash
python input_prep/scripts/popsim_input_maker.py <census_key> input_prep/configs/<run_name>/prepare.yaml
```

Example:
```bash
python input_prep/scripts/popsim_input_maker.py <census_key> input_prep/configs/2024_synthesis/prepare.yaml
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

Using the helper script:
```bash
bash scripts/run_2024_synthesis.sh
```

Or directly with PopulationSim:
```bash
populationsim \
  -c /home/da/RDF2055/d_drive/popsim/runs/2024_synthesis/configs \
  -d /home/da/RDF2055/d_drive/popsim/runs/2024_synthesis/data \
  -o /home/da/RDF2055/d_drive/popsim/runs/2024_synthesis/output/$(date +%Y-%m-%d_%H)_one_pass/run
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
python scripts/run_two_pass_hhsize.py --run-name 2024_synthesis
```

Useful options:
```bash
python scripts/run_two_pass_hhsize.py --run-name 2024_synthesis --method legacy
python scripts/run_two_pass_hhsize.py --run-dir /home/da/RDF2055/d_drive/popsim/runs/2024_synthesis
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
nohup python scripts/run_two_pass_hhsize.py --run-name 2024_synthesis >/dev/null 2>&1 &
```

Check whether a PID is still running:
```bash
ps -fp <PID>
```

## 5. Validate a Run
Validate a standard run:

```bash
python scripts/validate_popsim_run.py --run-name 2024_synthesis --write-csv
```

Validate a two-pass final run:
```bash
python scripts/validate_popsim_run.py --output-dir /home/da/RDF2055/d_drive/popsim/runs/2024_synthesis/output/<YYYY-MM-DD>_<HH>_two_pass/pass2 --configs-dir /home/da/RDF2055/d_drive/popsim/runs/2024_synthesis/configs --write-csv
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
- example run folders are `2024_refine_taz` and `2024_refine_mcd`
- review and edit `data/<project_name>_control_totals_<TARGET_GEO>.csv`, then replace that file in-place before the later refinement/synthesis run
- `refinement/urbansim_refine_input.ipynb` remains reference material, but `scripts/generate_refinement_inputs.py` is the maintained entrypoint
