# SEMCOG PopulationSim Structure Review

## Summary

The repository now has a substantially improved structure compared with its original state.

The core workflow is clearly separated into:

1. package-backed prep code under `src/semcog_popsim/input_prep/`
2. package-backed PopulationSim execution under `src/semcog_popsim/pipeline/`
3. operational entrypoints under `scripts/`
4. project/scenario assets under `projects/`
5. archived legacy material under `archive/`

The original design idea was sound from the start: one phase prepares inputs and one phase runs synthesis. The main work has been turning that idea into a clearer repository structure and reducing the amount of duplicated, notebook-derived, and path-coupled code.

## Recommendation

The two major scripts should not both live as full implementations in the project root.

That recommendation has now largely been implemented in practice:
- reusable logic lives under `src/semcog_popsim/`
- thin operational entrypoints live under `scripts/`
- legacy root-level wrappers remain only for compatibility and no longer carry primary implementation logic

A practical current structure is now closer to:

```text
SEMCOG_popsim/
  README.md
  pyproject.toml
  docs/
  scripts/
    prepare_inputs.py
    run_popsim.py
    run_pop_refinement_2015.py
    run_pop_refinement_2020.py
    run_placement_2015.py
    run_placement_2020.py
    hh_size_balancer.py
  src/semcog_popsim/
    input_prep/
    pipeline/
    forecast_refinement/
  projects/
    2017/
    2019/
    2020/
    2022/
    geo/
  configs/
    base/
    runs/
  data/
  outputs/
  archive/
  notebooks/
  refinement/
  validation/
```

## What Improved

- top-level operational runners were moved behind `scripts/`
- the remaining root-level runner files are explicitly compatibility-only
- input-prep code was packaged
- PopulationSim run logic was packaged
- project/scenario assets were copied into a clearer `projects/` layout
- duplicated utility logic was reduced
- legacy adjustment helpers were converted into valid scripts
- refinement and placement runners now share package-backed logic
- legacy `forecast_input/` now acts as a compatibility shim for the packaged forecast helpers
- low-risk legacy materials were moved into `archive/`

## Remaining Gaps

The reorganization is substantial, but not finished.

The biggest remaining gaps are:
- the canonical runtime config flow still leans on legacy `configs/` files
- the new `configs/runs/<name>/` directories still coexist with the legacy flat runtime config path
- validation utilities are not yet packaged under `src/semcog_popsim/validation/`
- some legacy notebook-derived material still remains visible in the active repo
- many non-core scripts still use machine-specific file paths

## Current Assessment

The repository structure is now sound enough that future cleanup can focus less on “where files live” and more on “which path is canonical”.

In other words, the project has moved from structural ambiguity to structural transition. That is a good place to be.

## Best Next Step

The highest-value next step is to make the run-side layout canonical, not just parallel.

Specifically:
1. document the generic `--run-config <name>` recipe
2. decide when `configs/runs/<name>/` becomes the only supported runtime layout
3. keep legacy flat config paths only as transitional compatibility until that cutover
