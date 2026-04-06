# SEMCOG PopulationSim

A simple guide for preparing and running SEMCOG PopulationSim.

## Main Folders
- `input_prep/configs/<run_name>/`: prep configs such as `prepare.yaml`, `controls_pre.csv`, and optional adjustment configs
- `d_drive/popsim/runs/<run_name>/configs/`: generated run configs such as `settings.yaml` and `controls.csv`
- `d_drive/popsim/runs/<run_name>/data/`: generated synthesis inputs
- `d_drive/popsim/runs/<run_name>/output/`: one-pass and two-pass synthesis outputs
- `scripts/`: operational run helpers

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
