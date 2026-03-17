import runpy
from pathlib import Path


if __name__ == "__main__":
    runpy.run_path(Path(__file__).resolve().parent / "run_popsim.py", run_name="__main__")
