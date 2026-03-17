# SEMCOG PopulationSim Package

## Overview
SEMCOG population synthesis package based on [RSG PopulationSim](https://github.com/ActivitySim/populationsim).

The repository is now organized around three layers:

1. `src/semcog_popsim/` contains reusable package code.
2. `scripts/` contains operational entrypoints.
3. compatibility wrappers remain in a few legacy locations where needed, but `scripts/` is the canonical entrypoint surface.

Supporting notes are available in:
- `docs/project_structure_review.md`
- `docs/target_structure_migration_plan.md`

---
## 1. Population Synthesis Preparation

### Preferred entrypoint
```bash
python scripts/prepare_inputs.py <key> <yml>
```

Preferred config location:
- `projects/<name>/prepare.yaml`

Examples:
```bash
python scripts/prepare_inputs.py <key> projects/2019/prepare.yaml
python scripts/prepare_inputs.py <key> projects/2022/prepare.yaml
python scripts/prepare_inputs.py <key> projects/baseline_a/prepare.yaml --output-dir data/interim/baseline_a
python scripts/prepare_inputs.py <key> projects/baseline_a/prepare.yaml --data-root /path/to/shared_data
```

Legacy compatibility entrypoint:
```bash
python input_prep/popsim_input_maker.py <key> <yml>
```

Legacy config locations under `input_prep/<year>/` still work during migration.

### Inputs
- `projects/<name>/prepare.yaml`: prep configuration
- `projects/<name>/settings.yaml`: project-specific PopulationSim settings template
- `projects/<name>/controls_pre.csv`: control-prep table
- `projects/geo/`: shared geographic crosswalk and equivalency tables
- PUMS input files referenced by the prep config

### Outputs
At the current migration stage, prep outputs land in `data/` relative to the chosen project folder by default, or can be redirected with `--output-dir` or `--data-root`.

Typical outputs:
- `[region]_[year]_geo_cross_walk.csv`
- `[region]_[year]_control_totals_[geo].csv`
- `[region]_[year]_seed_households.csv`
- `[region]_[year]_seed_persons.csv`
- `[region]_[year]_settings.yaml`

### Optional adjustments
Adjustment workflow example:
```bash
python input_prep/popsim_input_control_adj.py <key> projects/2019/control_adjustment.yaml
```

---
## 2. Run Population Synthesis

Preferred runner:
```bash
python scripts/run_popsim.py --run-config <name> -d <data_dir>
```

Compatibility runners:
```bash
python scripts/run_populationsim.py
python run_populationsim.py
```

The root-level runner files are compatibility wrappers only; use the `scripts/` versions for new work.

### Runtime config layout
The target runtime structure now exists in parallel with the legacy flat config directory.

Current migration state:
- legacy active config files still exist directly under `configs/`
- reusable config fragments are mirrored under `configs/base/`
- run-specific sets now exist for `configs/runs/2019/`, `configs/runs/2020/`, and `configs/runs/2022/`
- the runner can resolve `configs/runs/<name>/` plus `configs/base/` automatically with `--run-config <name>`
- `configs/runs/2022/` currently expects a prep-generated `controls.csv` to be staged before execution

The active runtime still primarily uses the legacy flat `configs/` path, so `configs/runs/<name>/` should be treated as the preferred emerging layout rather than the sole enforced path until the canonical run flow is finalized.

Preferred generic command shape:
```bash
python scripts/run_popsim.py --run-config 2020 -d data/2020_census_blkgrp
python scripts/run_popsim.py --run-config 2020 --data-root /path/to/shared_data
```

Override example:
```bash
python scripts/run_popsim.py \
  --run-config baseline_a \
  -d data/baseline_a \
  -o outputs/baseline_a_custom
```

Notes:
- `--run-config <name>` automatically uses `configs/runs/<name>/` and `configs/base/`
- the default data folder becomes `<data-root>/<name>` when `--run-config` is used without `-d`
- the default output folder becomes `outputs/<name>` when `--run-config` is used
- you can add a future `configs/runs/scenario_x/` and run it through the same entrypoint
- if you pass explicit `-c` values, those are used instead of the run-config defaults
- `--year <name>` is still accepted as a compatibility alias
- `configs/runs/2022/` will still need a staged `controls.csv` before execution

---
## 3. Forecast Refinement And Placement

Shared logic for refinement and placement now lives under:
- `src/semcog_popsim/forecast_refinement/`

Year-specific wrappers remain under `scripts/`:
- `scripts/run_pop_refinement_2015.py`
- `scripts/run_pop_refinement_2020.py`
- `scripts/run_placement_2015.py`
- `scripts/run_placement_2020.py`

These scripts still use environment-specific file paths and should currently be treated as project-local operational wrappers rather than portable command-line tools.

---
## 4. Current Migration State

Implemented so far:
- package-backed prep code
- package-backed PopulationSim run code
- package-backed refinement and placement orchestration
- packaged forecast-input helpers under `src/semcog_popsim/forecast_input/`
- legacy top-level `forecast_input/` kept only as a compatibility shim
- `projects/<year>/` structure for year-specific assets
- `configs/base/` and `configs/runs/` scaffolding
- `archive/` for low-risk legacy material

Still to finish:
- full activation of run-specific configs as the default runtime path
- validation packaging
- broader cleanup of remaining legacy notebook-derived material
- stricter separation of generated data into `data/raw`, `data/interim`, `data/processed`, and `outputs`
