# SEMCOG PopulationSim Target Structure Migration Plan

## Goal

Migrate the repository from its current mixed layout into a more maintainable structure:

```text
SEMCOG_popsim/
  README.md
  pyproject.toml
  scripts/
    prepare_inputs.py
    run_popsim.py
  src/semcog_popsim/
    input_prep/
      maker.py
      control_adjustment.py
      utils.py
    pipeline/
      run.py
    forecast_refinement/
    validation/
  configs/
    base/
    runs/
  projects/
    2019/
      prepare.yaml
      settings.yaml
      controls_pre.csv
    2020/
    2022/
  data/
    raw/
    interim/
    processed/
  outputs/
  notebooks/
  docs/
  archive/
```

The intent is not only to move files, but to separate:

- reusable source code
- run entrypoints
- year/scenario-specific project assets
- generated data
- legacy or exploratory material

## Current-State Summary

The current repository already has a useful high-level workflow:

1. `input_prep/` prepares ACS controls, PUMS seed files, geographic crosswalks, and settings.
2. `run_populationsim.py` runs PopulationSim/ActivitySim.
3. `forecast_input/`, `refinement/`, `validation/`, and placement/refinement runners support additional downstream workflows.

However, the current structure mixes code, configs, generated outputs, notebooks, and older material in ways that create maintenance risk.

Examples of structural issues already identified:

- `input_prep/popsim_input_maker.py` hardcodes a working directory.
- `input_prep/popsim_input_maker.py` ignores its CLI yaml argument and loads a hardcoded config.
- `input_prep/` mixes code and year-specific project assets.
- Some runners and helpers use machine-specific absolute paths.
- Notebook-derived scripts duplicate logic that should eventually live in package modules.

## Recommended Delivery Strategy

Do this migration in phases. Avoid a single large cutover.

### Phase 1: Package Skeleton and Stable Entrypoints

Goal: create the future structure without breaking current usage.

Tasks:

1. Add `pyproject.toml`.
2. Create `src/semcog_popsim/` package skeleton.
3. Create these initial modules:
   - `src/semcog_popsim/input_prep/maker.py`
   - `src/semcog_popsim/input_prep/control_adjustment.py`
   - `src/semcog_popsim/input_prep/utils.py`
   - `src/semcog_popsim/pipeline/run.py`
4. Add new script entrypoints:
   - `scripts/prepare_inputs.py`
   - `scripts/run_popsim.py`
5. Keep compatibility wrappers for the current top-level runners.
6. Update README usage examples to prefer the new script paths.

Definition of done:

- old commands still work
- new script paths exist
- code can be imported from `src/semcog_popsim/`

Estimated effort:

- 1 to 2 days

Main risk:

- import and path assumptions from scripts written to run only from the repo root

### Phase 2: Input-Prep Refactor

Goal: move the prep logic from script-style code into reusable functions.

Tasks:

1. Refactor `input_prep/popsim_input_maker.py` into package code:
   - config loading
   - geography crosswalk creation
   - ACS marginal download and compilation
   - PUMS seed extraction
   - settings mutation and output writing
2. Refactor `input_prep/popsim_input_control_adj.py` into package code.
3. Move shared helper logic into `src/semcog_popsim/input_prep/utils.py`.
4. Replace wildcard imports such as `from input_utils import *` with explicit imports.
5. Remove hardcoded `os.chdir(...)` and resolve paths from the provided config or repository root.
6. Make the prep workflow callable via `main()` and importable functions.

Definition of done:

- `scripts/prepare_inputs.py` calls package code
- no required hardcoded working directory
- yaml/config path is honored from the command line

Estimated effort:

- 2 to 3 days

Main risks:

- hidden reliance on relative paths like `./geo/...`
- behavior drift while converting notebook-style script blocks into functions

### Phase 3: Project Asset Separation

Goal: separate reusable code from year/scenario-specific assets.

Tasks:

1. Create `projects/2017`, `projects/2019`, `projects/2020`, `projects/2022`.
2. Move year-specific files from `input_prep/<year>/` into `projects/<year>/`.
3. Standardize filenames where useful:
   - `prepare.yaml`
   - `settings.yaml`
   - `controls_pre.csv`
   - adjustment configs if needed
4. Update code so it treats a project folder as an explicit input, not an implicit relative location.
5. Document how to add a new year or scenario.

Definition of done:

- year/scenario assets no longer live inside the source-code folder
- prep scripts can target a project directory explicitly

Estimated effort:

- 1 to 2 days

Main risk:

- references embedded in yaml files or assumptions about old file names

### Phase 4: Runtime and Data Layout Cleanup

Goal: separate generated artifacts from source and configuration.

Tasks:

1. Create:
   - `data/raw/`
   - `data/interim/`
   - `data/processed/`
   - `outputs/`
2. Decide what belongs in each location.
   - raw: source downloads and immutable inputs
   - interim: temporary or partially processed prep outputs
   - processed: finalized synthesis-ready inputs
   - outputs: run results from PopulationSim and downstream scripts
3. Update prep and runtime code to write to the appropriate locations.
4. Update docs and examples.

Definition of done:

- generated files no longer accumulate beside source code or year configs
- data flow is visible from folder structure alone

Estimated effort:

- 1 to 2 days

Main risk:

- downstream scripts may assume current output locations

### Phase 5: Forecast Refinement and Validation Packaging

Goal: bring non-core workflows into the same structural model.

Tasks:

1. Move reusable code from `forecast_input/` into:
   - `src/semcog_popsim/forecast_refinement/`
2. Move validation helpers into:
   - `src/semcog_popsim/validation/`
3. Keep any highly project-specific scripts in `scripts/`.
4. Review placement and refinement runners for machine-specific paths and convert them to config-driven inputs.

Definition of done:

- reusable forecast/refinement code is importable
- script entrypoints are thin wrappers

Estimated effort:

- 2 to 4 days

Main risks:

- absolute paths in existing runners
- mixed assumptions about runtime environment and output destinations

### Phase 6: Archive and Legacy Cleanup

Goal: reduce noise without losing historical material.

Tasks:

1. Create `archive/`.
2. Move clearly inactive or exploratory material there, for example:
   - notebook-derived duplicate scripts
   - older preprocess assets if no longer part of the active workflow
   - `configs/old/`
   - older tests or one-off experiments such as `Census2010_DP_test/` if appropriate
3. Keep notebooks that remain useful in `notebooks/`.
4. Add notes about what is archival versus active.

Definition of done:

- active workflow is easier to discover
- historical materials remain available without cluttering the core repo

Estimated effort:

- 0.5 to 1 day

Main risk:

- moving something that is still quietly used

## Suggested Execution Order

Recommended order:

1. Phase 1: package skeleton and stable entrypoints
2. Phase 2: input-prep refactor
3. Phase 3: project asset separation
4. Phase 4: runtime and data layout cleanup
5. Phase 5: forecast refinement and validation packaging
6. Phase 6: archive and legacy cleanup

This order reduces risk because it stabilizes code imports and entrypoints before moving project data and outputs.

## What Can Be Deferred

The following items do not need to happen in the first migration pass:

- complete conversion of all notebooks into package code
- aggressive renaming of every file for stylistic consistency
- full validation-framework redesign
- immediate archival of every older folder

These can wait until the package structure and paths are stable.

## Effort Estimate

### Minimal near-target migration

Scope:

- add package skeleton
- add `pyproject.toml`
- add `scripts/prepare_inputs.py` and `scripts/run_popsim.py`
- refactor just enough code for imports and paths to work
- create `projects/` and move active year assets
- preserve compatibility wrappers

Estimated effort:

- 2 to 4 days

### Full clean target-structure migration

Scope:

- all phases above
- path cleanup
- forecast/refinement packaging
- archive cleanup
- documentation refresh
- basic smoke-test verification

Estimated effort:

- about 1 week for a strong first pass
- up to 1 to 2 weeks if we also want higher confidence testing and broader cleanup

## Critical Risks to Watch

1. Hardcoded working directories and relative paths.
2. Absolute machine-specific paths in placement/refinement runners.
3. Implicit references to `input_prep/<year>/` from yaml or scripts.
4. Duplicate notebook/script logic drifting apart.
5. Downstream expectations about where generated files appear.

## Success Criteria

The migration is successful when:

- the prep workflow runs from `scripts/prepare_inputs.py`
- the synthesis workflow runs from `scripts/run_popsim.py`
- reusable logic lives under `src/semcog_popsim/`
- year/scenario assets live under `projects/`
- generated artifacts live under `data/` and `outputs/`
- legacy materials are present but no longer mixed into the main operational path
- the README explains the workflow in terms of the new structure

## Suggested Immediate Next Step

Start with a narrow implementation pass that completes only Phase 1 and the highest-value parts of Phase 2:

1. add `pyproject.toml`
2. create `src/semcog_popsim/input_prep/` and `src/semcog_popsim/pipeline/`
3. move the current prep and run logic behind importable functions
4. add `scripts/prepare_inputs.py` and `scripts/run_popsim.py`
5. fix `popsim_input_maker.py` path handling and config loading during the move

That gets the repository onto the target trajectory with the best risk-to-value ratio.
