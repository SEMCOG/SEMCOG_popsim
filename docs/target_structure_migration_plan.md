# SEMCOG PopulationSim Target Structure Migration Plan

## Goal

Migrate the repository from its original mixed layout into a more maintainable structure:

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

The intent is to separate:
- reusable source code
- run entrypoints
- year/scenario-specific project assets
- generated data
- legacy or exploratory material

## Progress Snapshot

### Completed or substantially completed

- `pyproject.toml` added
- `src/semcog_popsim/` package skeleton added
- package-backed input-prep code added under `src/semcog_popsim/input_prep/`
- package-backed PopulationSim run entrypoint added under `src/semcog_popsim/pipeline/`
- packaged runner now owns its own standard CLI parsing instead of depending on ActivitySim's fragile CLI import path
- `scripts/prepare_inputs.py` and `scripts/run_popsim.py` added
- compatibility wrappers retained for legacy entrypoints
- year-specific project assets copied into `projects/<year>/`
- shared geography crosswalks copied into `projects/geo/`
- `configs/base/` plus run-specific sets for 2019, 2020, and 2022 created as parallel target structure
- `archive/` created and low-risk legacy content moved there
- duplicated `input_utils` logic replaced with a compatibility shim
- forecast refinement runner logic packaged under `src/semcog_popsim/forecast_refinement/`
- placement runner logic packaged under `src/semcog_popsim/forecast_refinement/`
- two legacy adjustment helpers converted from notebook-style code into valid standalone scripts

### Still incomplete

- the canonical runtime flow still depends on the legacy flat `configs/` layout
- `configs/runs/2019/`, `configs/runs/2020/`, and `configs/runs/2022/` now exist, but the runtime still leans on the legacy flat `configs/` path
- generated outputs still primarily land in legacy year/data locations rather than a fully enforced `data/raw|interim|processed` flow
- validation helpers have not yet been packaged into `src/semcog_popsim/validation/`
- several notebooks and legacy helper files still remain outside `archive/`
- path-heavy runner scripts still use machine-specific locations and have not yet been generalized into project-driven configs

## Current-State Summary

The repository now has a clearer high-level workflow than it started with:

1. `scripts/prepare_inputs.py` calls packaged input-prep code in `src/semcog_popsim/input_prep/`.
2. `scripts/run_popsim.py` calls packaged PopulationSim execution code in `src/semcog_popsim/pipeline/`.
3. project/year assets now also exist under `projects/<year>/`.
4. non-core refinement and placement runners have started moving toward package-backed shared logic.

The remaining work is less about initial structure creation and more about consolidating the active workflow around the new structure.

## Recommended Next Focus Areas

### 1. Canonical runtime configuration

Goal: move from “parallel target structure exists” to “target structure is the preferred operational path”.

Recommended tasks:
- document one canonical run recipe per supported year
- decide whether `configs/settings.yaml` remains a compatibility artifact or becomes generated from `configs/runs/<year>/`
- decide when `configs/runs/<year>/controls.csv` is required to be present before execution

Why this matters:
- the prep side now has a usable `projects/` structure, but the run side still primarily points to the legacy flat config layout

### 2. Validation packaging

Goal: move reusable validation logic into `src/semcog_popsim/validation/`.

Recommended tasks:
- identify reusable code in `validation/` and `input_prep/output_plot.py`
- package only the reusable parts
- leave notebook-heavy or one-off plotting workflows for later or archive

Why this matters:
- the main prep and run paths are increasingly package-backed, but validation remains mostly ad hoc

### 3. Archive review

Goal: continue reducing noise without moving anything still operational.

Recommended tasks:
- review duplicate notebook-derived scripts
- review older preprocess assets
- move only obviously inactive materials into `archive/`

Why this matters:
- the repo is cleaner, but there is still legacy material competing with the active path

## Phase-by-Phase Status

### Phase 1: Package Skeleton and Stable Entrypoints
Status: completed

Delivered:
- `pyproject.toml`
- `src/semcog_popsim/`
- `scripts/prepare_inputs.py`
- `scripts/run_popsim.py`
- compatibility wrappers

### Phase 2: Input-Prep Refactor
Status: mostly completed for the core prep path

Delivered:
- packaged input-prep logic in `maker.py`, `control_adjustment.py`, and `utils.py`
- removal of hardcoded working-directory dependency from the active prep flow
- CLI config path honored in the packaged prep path
- duplicated `input_utils` logic replaced with a shim

Remaining:
- fully retire remaining legacy helpers or keep them clearly transitional

### Phase 3: Project Asset Separation
Status: partially completed

Delivered:
- `projects/2017`, `projects/2019`, `projects/2020`, `projects/2022`
- standardized `prepare.yaml`, `settings.yaml`, and `controls_pre.csv` copies
- shared `projects/geo/`

Remaining:
- decide when the legacy `input_prep/<year>/` assets stop being the fallback path

### Phase 4: Runtime and Data Layout Cleanup
Status: scaffolding created, not yet fully enforced

Delivered:
- `data/raw/`, `data/interim/`, `data/processed/`
- `outputs/`
- `configs/base/`
- `configs/runs/2020/`

Remaining:
- make the active runtime prefer the new structure instead of just mirroring it
- finalize the 2022 control-file staging convention

### Phase 5: Forecast Refinement and Validation Packaging
Status: partially completed

Delivered:
- forecast refinement runner packaging
- placement runner packaging

Remaining:
- validation packaging
- further cleanup of machine-specific path assumptions

### Phase 6: Archive and Legacy Cleanup
Status: started

Delivered:
- `archive/`
- archived `Census2010_DP_test/`
- archived `configs/old/`

Remaining:
- review additional notebook-derived and experimental materials

## Suggested Execution Order From Here

Recommended next order:

1. define and document one canonical end-to-end run flow
2. decide how prep-generated controls are staged into `configs/runs/<year>/`
3. package validation helpers that are truly reusable
4. continue selective archive cleanup

## Success Criteria

The migration will feel operationally complete when:
- prep uses `scripts/prepare_inputs.py` with `projects/<year>/prepare.yaml`
- synthesis uses `scripts/run_popsim.py` with clearly defined `configs/runs/<year>/`
- reusable logic lives under `src/semcog_popsim/`
- year/scenario assets live under `projects/`
- generated artifacts have a clearer home under `data/` and `outputs/`
- legacy materials are still available but no longer visually dominate the active path

## Suggested Immediate Next Step

Document one canonical run recipe in the README for each supported year, then decide how prep-generated control files should be staged into `configs/runs/<year>/` before execution.

That is the most valuable next step because the structural scaffolding now exists and the remaining work is mostly about making the run path explicit and repeatable.
