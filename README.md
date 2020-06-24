# SEMCOG PopulationSim Package

## Overview
SEMCOG population synthesis package based on [RSG PopulationSim](https://github.com/ActivitySim/populationsim).

---
### 1. Population Synthesis Preparation

A major function of SEMCOG_popsim is to prepare input configuration and dataset for PopulationSim, including project settings, controls, demographic marginals and samples for target geographies, and a geographic cross work table.

#### Popsim Input Maker ('/notebooks/popsim_input_maker.py')
usage:  python popsim_input_maker.py key yml 
 - key: Census API key
 - yml: yaml configuration file name (example: region.yml)

Data for input maker:
1. controls_pre.csv: an csv table extended from PopulationSim control file ("configs/controls.csv"). Control file is a customized CSV table containing synthesis variable definitions, geographies and sample query expression. "controls_pre.csv" table adds a new "acs_variables" field with ACS variables and pandas expressions. Input maker will use this field to download Census marginals and compile to "targe" variables. 
2. PUMS sample files: include both PUMS households and persons for the whole region. 

Outputs: 
1. geo_cross_walk.csv: cross walk table for all geographies to be used in synthesis (PUMS, Census Block Groups, Tract, TAZ, etc);
2. <geo>_control_totals: control marginals at single or multiple levels;
3. <region>_seed_households: seed households;
4. <region>_seed_persons: seed persons

#### Run Population Synthesis
Need both RSG[PopulationSim](https://github.com/ActivitySim/populationsim) and [ActivitySim](https://github.com/ActivitySim/activitysim) for the synthesis process.

Input structure:
 - /configs/settings.yaml(project settings)
 - /configs/controls.csv(popsim control)
 - /data/xxx_geo_cross_walk.csv
 - /data/xxx_control_totals.csv(could have multiple control totals by geography, such as TAZ_control_totals.csv)
 - /data/xxx_seed_households.csv
 - /data/xxx_seed_persons.csv

#### Results and visualization
 * output folder has synthetic HHs, persons and one or more summary_<geo>.csv files
 * output_stats_plots.ipynb could be used to generate [error plot]
(https://raw.githubusercontent.com/SEMCOG/SEMCOG_popsim/master/validation/semcog_python/synpop_popsim_error_plot.png) and [histograms]

---
### 2. Forecast Refinement
SEMCOG tested using PopulationSim as refinement tool for UrbanSim model. 

urbansim_refine_input("refinement/urbansim_refine_input.ipynb") script uses model data to prepare inputs for the refinement process. 

#### Test Inputs
Similar inputs to population synthesis are expected 
1. settings: with manual updates;
2. controls: urbansim_refine_input script extracts and compiles information from "annual households control totals" from forcast model;
3. geo_cross_walk: script generated using model parcels and buildings;
4. control_totals: summarized from official SEMCOG forecast, or use reviewed indicators;
5. seed_households: use model output households, add weight and geographies
6. see_persons: use model output persons, add weight and geographies
   
---
### 3. ABM test inputs
This package also contains the work to develop [test demographic inputs](https://github.com/SEMCOG/SEMCOG_popsim/tree/master/urbansim_to_abm) for SEMCOG ABM model. 
