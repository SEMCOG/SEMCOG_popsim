# SEMCOG PopulationSim Package

## Overview
SEMCOG population synthesis package based on [RSG PopulationSim](https://github.com/ActivitySim/populationsim).

---
## 1. Population Synthesis Preparation

A major function of SEMCOG_popsim is to prepare input configuration and dataset for PopulationSim, including project settings, controls, demographic marginals and samples for target geographies, and a geographic cross work table.

### Popsim Input Maker ('/input_prep/popsim_input_maker.py')
##### usage: 
```
  python popsim_input_maker.py key yml 

 - key: Census API key
 - yaml: input maker configuration (example: region_2019.yaml)
```
##### Inputs:
 - *[year]\_controls\_pre\_[year].csv*: &nbsp;&nbsp;&nbsp;&nbsp;an csv table extended from PopulationSim control file ("configs/controls.csv"). Control file is a customized CSV table containing synthesis variable definitions, geographies and sample query expression. "controls_pre_year.csv" table adds a new "acs_variables" field with ACS variables and pandas expressions. Input maker will use this field to download Census marginals and compile to "targe" variables. 
 - *[year]\_region\_[year]/yaml*: &nbsp;&nbsp;&nbsp;configuration file for input_maker, including necessary input, such as PUMS data, geo equivalency files location, PUMS variable updates, etc.
 - *PUMS/*: &nbsp;&nbsp;&nbsp;folder include both PUMS households and persons for the whole region. 
 - *geo/*: &nbsp;&nbsp;&nbsp;folder for geographic information tables.

##### Outputs:  
 All outputs are produced to '[year]/data' folder
 - *[region]\_[year]\_geo\_cross\_walk.csv*: &nbsp;&nbsp;&nbsp;cross walk table for all geographies to be used in synthesis (PUMS, Census Block Groups, Tract, TAZ, etc);
 - *[region]\_[year]\_control\_totals\_[geo].csv*: &nbsp;&nbsp;&nbsp;control marginals at single or multiple levels;
 - *[region]\_[year]\_seed\_households.csv*: &nbsp;&nbsp;&nbsp;seed households;
 - *[region]\_[year]\_seed\_persons.csv*: &nbsp;&nbsp;&nbsp;seed persons

##### *(Optional)* Adjustments: 
All margional controls could be scaled to a closer-to-reality totals. For example, adjusting 2019 5-year ACS BGs to 2019 1-year ACS County totals,so the results are closer to 2019 ground 'Truth'. A 2-step process is needed to accomplish this adjustment. Using county adjustment as example:
- Step 1. download county level control totals as adjustment targets
```
  python popsim_input_control_adj.py key yml 
 - key: Census API key
 - yml: adjustment configurations (example: region_2019_control_adj.yaml)
```
In additon, a new county-level control file is needed. Format is similar to *[year]\_controls\_pre\_[year].csv*
 - step 2. adjust the control by county totals or category totals using *adjust_to_acs1_county.py*.

### Run Population Synthesis
Need both RSG[PopulationSim](https://github.com/ActivitySim/populationsim) and [ActivitySim](https://github.com/ActivitySim/activitysim) for the synthesis process.

Input structure:
 - /configs/settings.yaml (project settings)
 - /configs/controls.csv or the control file name defined in settings.yaml (popsim control)
 - /data/xxx_geo_cross_walk.csv
 - /data/xxx_control_totals.csv(could have multiple control totals by geography, such as TAZ_control_totals.csv)
 - /data/xxx_seed_households.csv
 - /data/xxx_seed_persons.csv
##### *(Optional)* household size rebalance: 
 - To adjust household size and solve the over sized 7+ HHs issue, a rebalance process is needed.
 - hh_size_balancer.py will need hh size control file and the output summary file to create a new control file with new household size distribution.
 - Rerun Popsim with new household sizes
 - repeat this process as needed 

### Results and visualization
 * output folder has synthetic HHs, persons and one or more summary_<geo>.csv files
 * output_stats_plots.ipynb could be used to generate [error plot]
(https://raw.githubusercontent.com/SEMCOG/SEMCOG_popsim/master/validation/semcog_python/synpop_popsim_error_plot.png) and [histograms](https://github.com/SEMCOG/SEMCOG_popsim/blob/master/validation/semcog_python/popsim_oakland_BLKGRP__histograms.html). 

---
## 2. Forecast Refinement
SEMCOG tested using PopulationSim as refinement tool for UrbanSim model. 

urbansim_refine_input("refinement/urbansim_refine_input.ipynb") script uses model data to prepare inputs for the refinement process. 

### Test Inputs
Similar inputs to population synthesis are expected 
1. settings: with manual updates;
2. controls: urbansim_refine_input script extracts and compiles information from "annual households control totals" from forcast model;
3. geo_cross_walk: script generated using model parcels and buildings;
4. control_totals: summarized from official SEMCOG forecast, or use reviewed indicators;
5. seed_households: use model output households, add weight and geographies
6. see_persons: use model output persons, add weight and geographies
   
---
## 3. ABM test inputs
This package also contains the work to develop [test demographic inputs](https://github.com/SEMCOG/SEMCOG_popsim/tree/master/urbansim_to_abm) for SEMCOG ABM model. 
