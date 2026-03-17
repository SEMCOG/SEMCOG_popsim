import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from semcog_popsim.input_prep.utils import (  # noqa: F401
    CensusDownloader,
    combine_puma_data,
    group_pums_data,
    marginal_summary,
    preprocess_pums,
    pums_update,
    read_tract_puma_crosswalk,
)

# Backward-compatible name for legacy scripts.
Census_Downloader = CensusDownloader
