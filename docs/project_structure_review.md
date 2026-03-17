# SEMCOG PopulationSim Structure Review

## Summary

The current repository already has the right high-level split:

1. `input_prep/` prepares ACS controls, geography crosswalks, seed households/persons, and derived PopulationSim settings.
2. `run_populationsim.py` runs the PopulationSim/ActivitySim synthesis once those prepared inputs are ready.

That separation of concerns is sound. The main structural issue is not the existence of two major scripts, but that the repository mixes source code, entrypoints, year-specific project assets, generated outputs, notebooks, and older experimental material in ways that make the workflow harder to maintain.

## Recommendation

The two major scripts should not both live as full implementations in the project root.

Recommended pattern:

- Keep input-prep logic under `input_prep/` or a future package module such as `src/semcog_popsim/input_prep/`.
- Keep run/entry scripts under a dedicated `scripts/` folder.
- Keep the repository root for high-level documentation, stable wrappers, and project-wide assets such as `configs/`.

For this repository, a practical near-term structure is:

```text
SEMCOG_popsim/
  README.md
  RunPopulationSim.bat
  docs/
  scripts/
    run_populationsim.py
    run_pop_refinement_2015.py
    run_pop_refinement_2020.py
    run_placement_2015.py
    run_placement_2020.py
    hh_size_balancer.py
  input_prep/
    popsim_input_maker.py
    popsim_input_control_adj.py
    input_utils.py
    geo/
    2017/
    2019/
    2020/
    2022/
  configs/
  notebooks/
  refinement/
  validation/
```

In a later cleanup, the year-specific inputs could move from `input_prep/<year>/` to something like `projects/<year>/`, and generated artifacts could be separated more clearly into `data/`, `outputs/`, or `interim/` locations.

## Why This Direction Helps

- It keeps the workflow understandable: `input_prep/` prepares, `scripts/` runs.
- It reduces clutter at the repository root.
- It makes room for future packaging without forcing a large refactor now.
- It preserves compatibility by allowing lightweight wrappers to remain at the old top-level paths.

## Issues Noticed During Review

These are more important than folder naming alone:

- `input_prep/popsim_input_maker.py` hardcodes a working directory with `os.chdir(...)`.
- `input_prep/popsim_input_maker.py` defines a CLI `yaml` argument but currently loads a hardcoded `2022/prepare_2022.yaml`.
- The repository contains duplicate or notebook-derived script logic, including `notebooks/popsim_input_maker.py`.
- Top-level helper and run scripts are mixed together in the root even though they are all operational entrypoints.
- `RunPopulationSim.bat` is environment-specific and references the root runner path directly.

## Improvements Recommended After This Reorganization

1. Refactor `popsim_input_maker.py` into functions with a `main()` entrypoint.
2. Remove hardcoded paths and always resolve files relative to the script or supplied config.
3. Consolidate duplicated notebook/script logic.
4. Separate year/scenario assets from generated outputs more explicitly.
5. Consider a future package layout such as `src/semcog_popsim/` if testability and reuse become priorities.

## Reorganization Applied

This reorganization focuses on low-risk structural cleanup:

- Added this review under `docs/`.
- Created `scripts/` for operational runners.
- Moved top-level runner/helper implementations into `scripts/`.
- Kept small root-level wrapper scripts so existing commands still work.
- Updated the README and batch launcher to point at the new script locations.
