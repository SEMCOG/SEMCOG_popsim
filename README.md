# SEMCOG PopulationSim Package

## Overview
SEMCOG population synthesis package based on [RSG PopulationSim](https://github.com/ActivitySim/populationsim).

The repository is now organized around three layers:

1. `src/semcog_popsim/` contains reusable package code for input preparation and PopulationSim execution.
2. `scripts/` contains operational entrypoints.
3. compatibility wrappers remain in the legacy locations where needed.

Supporting notes are available in:
- `docs/project_structure_review.md`
- `docs/target_structure_migration_plan.md`

---
## 1. Population Synthesis Preparation

A major function of SEMCOG_popsim is to prepare input configuration and dataset for PopulationSim, including project settings, controls, demographic marginals and samples for target geographies, and a geographic cross walk table.

### Preferred entrypoint
```bash
python scripts/prepare_inputs.py <key> <yml>
```

Preferred config location:
- `projects/<year>/prepare.yaml`

Examples:
```bash
python scripts/prepare_inputs.py <key> projects/2019/prepare.yaml
python scripts/prepare_inputs.py <key> projects/2022/prepare.yaml
```

Legacy compatibility entrypoint:
```bash
python input_prep/popsim_input_maker.py <key> <yml>
```

Legacy config locations under `input_prep/<year>/` still work during migration.

### Inputs
- *controls_pre.csv*: an extended PopulationSim control file with an added `acs_variables` field used to download Census marginals and compile target control variables.
- *prepare.yaml*: project/year-specific prep configuration, including PUMS locations, geography files, and PUMS variable updates.
- *settings.yaml*: project/year-specific PopulationSim settings template.
- *PUMS/*: household and person PUMS files for the region.
- *projects/geo/*: shared geographic crosswalk and equivalency tables.

### Outputs
All outputs are still produced to `[year]/data` in the current migration stage.
- *[region]\_[year]\_geo\_cross\_walk.csv*: crosswalk table for synthesis geographies.
- *[region]\_[year]\_control\_totals\_[geo].csv*: control marginals by geography.
- *[region]\_[year]\_seed\_households.csv*: seed households.
- *[region]\_[year]\_seed\_persons.csv*: seed persons.
- *[region]\_[year]\_settings.yaml*: run-ready settings file.

### Optional adjustments
All marginal controls could be scaled to closer-to-reality totals. For example, 2019 5-year ACS block group controls can be adjusted to 2019 1-year ACS county totals.

Adjustment config example:
```bash
python input_prep/popsim_input_control_adj.py <key> projects/2019/control_adjustment.yaml
```

### Phase 1 and 2 migration note
The migration now includes:
- package-backed prep code in `src/semcog_popsim/input_prep/`
- a package-backed run entrypoint in `src/semcog_popsim/pipeline/`
- standardized project assets under `projects/<year>/`
- shared prep geography assets under `projects/geo/`

The original files under `input_prep/<year>/` are still present for backward compatibility.

---
## 2. Run Population Synthesis
Both [PopulationSim](https://github.com/ActivitySim/populationsim) and [ActivitySim](https://github.com/ActivitySim/activitysim) are required for the synthesis process.

Preferred runner:
```bash
python scripts/run_popsim.py
```

Compatibility runners:
```bash
python scripts/run_populationsim.py
python run_populationsim.py
```

Input structure:
- `configs/settings.yaml`: project settings
- `configs/controls.csv` or the control file defined in `settings.yaml`
- `data/xxx_geo_cross_walk.csv`
- `data/xxx_control_totals.csv`
- `data/xxx_seed_households.csv`
- `data/xxx_seed_persons.csv`

### Optional household size rebalance
- `scripts/hh_size_balancer.py` uses a household size control file plus output summaries to create an adjusted control file.
- Rerun PopulationSim with the updated controls as needed.

### Results and visualization
- The output folder contains synthesized households, persons, and one or more `summary_<geo>.csv` files.
- `notebooks/output_stats_plots.ipynb` can be used to generate error plots and histograms.

---
## 3. Forecast Refinement
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
