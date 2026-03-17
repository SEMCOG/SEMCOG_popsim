import runpy
from pathlib import Path


if __name__ == "__main__":
    runpy.run_path(Path(__file__).resolve().parent / "scripts" / "hh_size_balancer.py", run_name="__main__")
