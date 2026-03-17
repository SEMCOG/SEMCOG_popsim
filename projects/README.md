# Projects

This directory is the first step toward separating year-specific assets from source code.

Current conventions:
- `projects/<year>/prepare.yaml` is the preferred prep config for that year.
- `projects/<year>/settings.yaml` is the standardized PopulationSim settings template for that year.
- `projects/<year>/controls_pre.csv` is the standardized control-prep table for that year.
- `projects/geo/` holds shared geography crosswalk files used by multiple project years.

The original files under `input_prep/<year>/` are still present for backward compatibility during migration.

Additional notes:
- project folder names do not need to be numeric years; the packaged prep flow now resolves files relative to the config folder itself.
- `python scripts/prepare_inputs.py <key> <yaml> --output-dir <path>` can be used when you want a fully explicit output location.
- `python scripts/prepare_inputs.py <key> <yaml> --data-root <path>` can be used when you want outputs grouped under a shared external data folder.
