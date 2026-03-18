# SEMCOG PopulationSim Package

## Overview
SEMCOG population synthesis package based on [RSG PopulationSim](https://github.com/ActivitySim/populationsim).

---
## 1. Population Synthesis Preparation

A major function of `SEMCOG_popsim` is to prepare PopulationSim run packages, including project settings, controls, demographic marginals, PUMS seed samples, and a geographic crosswalk table.

### Popsim Input Maker (`/input_prep/scripts/popsim_input_maker.py`)
##### Usage
```
  python input_prep/scripts/popsim_input_maker.py <census_key> input_prep/configs/<run_name>/prepare.yaml 

 - key: Census API key
 - yaml: input maker configuration (example: input_prep/configs/2024_synthesis/prepare.yaml)
```
##### Inputs
- `input_prep/configs/<run_name>/prepare.yaml`: input-maker configuration, including target year, run name, PUMS paths, and variable updates.
- `input_prep/configs/<run_name>/controls_pre.csv`: extended PopulationSim control specification. This adds an `acs_variables` field used to download and compile Census marginals into PopulationSim control fields.
- `d_drive/popsim/inputs/pums/`: shared PUMS source data for all runs.
- `input_prep/geo/`: shared geographic equivalency and tract-to-PUMA crosswalk files.

##### Outputs
All run-ready outputs are produced to `d_drive/popsim/runs/<run_name>/` with `configs/`, `data/`, and `output/` subfolders.
- `[region]_[year]_geo_cross_walk.csv`: crosswalk table for synthesis geographies.
- `[region]_[year]_control_totals_[geo].csv`: compiled control marginals by geography.
- `[region]_[year]_seed_households.csv`: seed households.
- `[region]_[year]_seed_persons.csv`: seed persons.
- `configs/settings.yaml`: generated PopulationSim runtime settings.
- `configs/controls.csv`: generated PopulationSim runtime control specification.

##### *(Optional)* Adjustments:
All margional controls could be scaled to a closer-to-reality totals. For example, adjusting 2019 5-year ACS BGs to 2019 1-year ACS County totals,so the results are closer to 2019 ground 'Truth'. A 2-step process is needed to accomplish this adjustment. Using county adjustment as example:
- Step 1. download county level control totals as adjustment targets
```
  python input_prep/scripts/popsim_input_control_adj.py <census_key> input_prep/configs/<run_name>/adjust.yaml 
 - key: Census API key
 - yml: adjustment configuration (example: input_prep/configs/2019_county_adjustment/adjust.yaml)
```
In additon, a new county-level control file is needed. Format is similar to `input_prep/configs/<run_name>/controls_pre.csv`.
- Step 2. adjust the control by county totals or category totals using `adjust_to_acs1_county.py`.

### Run Population Synthesis
Use the PopulationSim CLI against a generated run package under `d_drive/popsim/runs/<run_name>/`.

Validated run package structure:
- `configs/settings.yaml`
- `configs/controls.csv`
- `data/xxx_geo_cross_walk.csv`
- `data/xxx_control_totals_<geo>.csv`
- `data/xxx_seed_households.csv`
- `data/xxx_seed_persons.csv`
- `output/`

CLI example:
```bash
populationsim \
  -c /home/da/RDF2055/d_drive/popsim/runs/2024_synthesis/configs \
  -d /home/da/RDF2055/d_drive/popsim/runs/2024_synthesis/data \
  -o /home/da/RDF2055/d_drive/popsim/runs/2024_synthesis/output
```

Notes:
- the generated settings file uses `data_dir: data`, so the CLI should point `-d` to the run package data folder
- the SEMCOG helper script `run_populationsim.py` is legacy compatibility code; prefer the upstream `populationsim` CLI in Docker

##### *(Optional)* household size rebalance:
- To adjust household size and solve the over sized 7+ HHs issue, a rebalance process is needed.
- `hh_size_balancer.py` will need a household-size control file and the output summary file to create a new control file with new household size distribution.
- Rerun PopulationSim with new household sizes.
- Repeat this process as needed.

### Results and visualization
- `output/` has synthetic households, persons and one or more `summary_<geo>.csv` files.
- `output_stats_plots.ipynb` could be used to generate [error plot](https://raw.githubusercontent.com/SEMCOG/SEMCOG_popsim/master/validation/semcog_python/synpop_popsim_error_plot.png) and [histograms](https://github.com/SEMCOG/SEMCOG_popsim/blob/master/validation/semcog_python/popsim_oakland_BLKGRP__histograms.html).

---
## 2. Forecast Refinement
SEMCOG tested using PopulationSim as refinement tool for UrbanSim model.

`urbansim_refine_input("refinement/urbansim_refine_input.ipynb")` script uses model data to prepare inputs for the refinement process.

### Test Inputs
Similar inputs to population synthesis are expected
1. settings: with manual updates;
2. controls: `urbansim_refine_input` extracts and compiles information from annual household control totals from the forecast model;
3. geo_cross_walk: generated from model parcels and buildings;
4. control_totals: summarized from official SEMCOG forecast, or use reviewed indicators;
5. seed_households: use model output households, add weight and geographies;
6. seed_persons: use model output persons, add weight and geographies.
