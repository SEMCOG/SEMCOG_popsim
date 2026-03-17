# SEMCOG PopulationSim Package

## Overview
SEMCOG population synthesis package based on [RSG PopulationSim](https://github.com/ActivitySim/populationsim).

The repository is organized around two main stages:

1. `input_prep/` prepares ACS controls, seed households/persons, geographic crosswalks, and run settings.
2. `scripts/` contains operational entrypoints for running PopulationSim and related helper workflows.

A structure review and reorganization note is available in `docs/project_structure_review.md`.

---
## 1. Population Synthesis Preparation

A major function of SEMCOG_popsim is to prepare input configuration and dataset for PopulationSim, including project settings, controls, demographic marginals and samples for target geographies, and a geographic cross walk table.

### Popsim Input Maker (`input_prep/popsim_input_maker.py`)
##### Usage
```
python input_prep/popsim_input_maker.py <key> <yml>

- key: Census API key
- yml: input maker configuration (example: region_2019.yaml)
```

##### Inputs
- *[year]\_controls\_pre\_[year].csv*: an extended PopulationSim control file with an added `acs_variables` field used to download Census marginals and compile target control variables.
- *[year]\_region\_[year].yaml*: configuration for the input maker, including PUMS locations, geography files, and PUMS variable updates.
- *PUMS/*: household and person PUMS files for the region.
- *geo/*: geographic crosswalk and equivalency tables.

##### Outputs
All outputs are produced to `[year]/data`.
- *[region]\_[year]\_geo\_cross\_walk.csv*: crosswalk table for synthesis geographies.
- *[region]\_[year]\_control\_totals\_[geo].csv*: control marginals by geography.
- *[region]\_[year]\_seed\_households.csv*: seed households.
- *[region]\_[year]\_seed\_persons.csv*: seed persons.

##### Optional adjustments
All marginal controls could be scaled to closer-to-reality totals. For example, 2019 5-year ACS block group controls can be adjusted to 2019 1-year ACS county totals.

Step 1. Download county-level control totals as adjustment targets.
```
python input_prep/popsim_input_control_adj.py <key> <yml>

- key: Census API key
- yml: adjustment configuration (example: region_2019_control_adj.yaml)
```

A county-level control file is also needed. Its format is similar to *[year]\_controls\_pre\_[year].csv*.

Step 2. Adjust controls using `input_prep/adjust_to_acs1_county.py`.

### Run Population Synthesis
Both [PopulationSim](https://github.com/ActivitySim/populationsim) and [ActivitySim](https://github.com/ActivitySim/activitysim) are required for the synthesis process.

Primary runner:
```
python scripts/run_populationsim.py
```

The legacy top-level wrapper `python run_populationsim.py` is still available for compatibility.

Input structure:
- `configs/settings.yaml`: project settings
- `configs/controls.csv` or the control file defined in `settings.yaml`
- `data/xxx_geo_cross_walk.csv`
- `data/xxx_control_totals.csv`
- `data/xxx_seed_households.csv`
- `data/xxx_seed_persons.csv`

##### Optional household size rebalance
- `scripts/hh_size_balancer.py` uses a household size control file plus output summaries to create an adjusted control file.
- Rerun PopulationSim with the updated controls as needed.

### Results and visualization
- The output folder contains synthesized households, persons, and one or more `summary_<geo>.csv` files.
- `notebooks/output_stats_plots.ipynb` can be used to generate error plots and histograms.

---
## 2. Forecast Refinement
SEMCOG has also tested PopulationSim as a refinement tool for the UrbanSim model.

`refinement/urbansim_refine_input.ipynb` prepares inputs for the refinement process.

### Test inputs
Similar inputs to population synthesis are expected:
1. settings with manual updates
2. controls compiled from annual household control totals from the forecast model
3. geo crosswalks generated from model parcels and buildings
4. control totals summarized from official SEMCOG forecasts or reviewed indicators
5. seed households from model output households with weights and geographies
6. seed persons from model output persons with weights and geographies
