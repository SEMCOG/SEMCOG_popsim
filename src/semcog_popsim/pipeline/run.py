import argparse
import logging
import os

from activitysim.core import config
from activitysim.core import inject
from activitysim.core import pipeline
from activitysim.core import tracing
from activitysim.core.config import setting
from activitysim.core.tracing import print_elapsed_time

logger = logging.getLogger("populationsim")
INJECTABLES = ["data_dir", "configs_dir", "output_dir", "settings_file_name"]


def add_run_args(parser, multiprocess=True):
    parser.add_argument("-w", "--working_dir", type=str, metavar="PATH", help="path to project directory")
    parser.add_argument("-c", "--config", type=str, action="append", metavar="PATH", help="path to config dir")
    parser.add_argument("-o", "--output", type=str, metavar="PATH", help="path to output dir")
    parser.add_argument("-d", "--data", type=str, action="append", metavar="PATH", help="path to data dir")
    parser.add_argument("-r", "--resume", type=str, metavar="STEPNAME", help="resume after step")
    parser.add_argument("-p", "--pipeline", type=str, metavar="FILE", help="pipeline file name")
    parser.add_argument("-s", "--settings_file", type=str, metavar="FILE", help="settings file name")
    parser.add_argument("-g", "--chunk_size", type=int, metavar="BYTES", help="chunk size")

    if multiprocess:
        parser.add_argument(
            "-m",
            "--multiprocess",
            default=False,
            const=-1,
            metavar="(N)",
            nargs="?",
            type=int,
            help="run multiprocess and optionally override the process count",
        )


def validate_injectable(name):
    try:
        dir_paths = inject.get_injectable(name)
    except RuntimeError:
        raise SystemExit(
            "Error: please specify either a --working_dir containing 'configs', 'data', and 'output' folders "
            "or all three of --config, --data, and --output"
        )

    dir_paths = [dir_paths] if isinstance(dir_paths, str) else dir_paths

    for dir_path in dir_paths:
        if not os.path.exists(dir_path):
            raise SystemExit(f"Could not find {name} '{os.path.abspath(dir_path)}'")

    return dir_paths


def handle_standard_args(args, multiprocess=True):
    def inject_arg(name, value, cache=False):
        assert name in INJECTABLES
        inject.add_injectable(name, value, cache=cache)

    if args.working_dir:
        os.chdir(args.working_dir)

    if args.settings_file:
        inject_arg("settings_file_name", args.settings_file, cache=True)

    if args.config:
        inject_arg("configs_dir", args.config)

    if args.data:
        inject_arg("data_dir", args.data)

    if args.output:
        inject_arg("output_dir", args.output)

    if multiprocess and args.multiprocess:
        config_paths = validate_injectable("configs_dir")

        if not os.path.exists("configs_mp"):
            logger.warning("could not find 'configs_mp'. skipping...")
        else:
            logger.info("adding 'configs_mp' to config_dir list...")
            config_paths.insert(0, "configs_mp")
            inject_arg("configs_dir", config_paths)

        config.override_setting("multiprocess", True)
        if args.multiprocess > 0:
            config.override_setting("num_processes", args.multiprocess)

    if args.chunk_size:
        config.override_setting("chunk_size", int(args.chunk_size))

    for injectable in ["configs_dir", "data_dir", "output_dir"]:
        validate_injectable(injectable)

    if args.pipeline:
        inject.add_injectable("pipeline_file_name", args.pipeline)

    if args.resume:
        config.override_setting("resume_after", args.resume)


def run_pipeline(argv=None):
    parser = argparse.ArgumentParser(description="Run SEMCOG PopulationSim")
    add_run_args(parser)
    args = parser.parse_args(argv)

    handle_standard_args(args)
    tracing.config_logger()
    t0 = print_elapsed_time()
    logger.info("GROUP_BY_INCIDENCE_SIGNATURE: %s", setting("GROUP_BY_INCIDENCE_SIGNATURE"))
    logger.info("INTEGERIZE_WITH_BACKSTOPPED_CONTROLS: %s", setting("INTEGERIZE_WITH_BACKSTOPPED_CONTROLS"))
    logger.info("SUB_BALANCE_WITH_FLOAT_SEED_WEIGHTS: %s", setting("SUB_BALANCE_WITH_FLOAT_SEED_WEIGHTS"))
    logger.info("meta_control_data: %s", setting("meta_control_data"))
    logger.info("control_file_name: %s", setting("control_file_name"))
    logger.info("USE_CVXPY: %s", setting("USE_CVXPY"))
    logger.info("USE_SIMUL_INTEGERIZER: %s", setting("USE_SIMUL_INTEGERIZER"))
    run_list_name = inject.get_injectable("run_list_name", "run_list")
    run_list = setting(run_list_name)
    assert "steps" in run_list, "Did not find steps in run_list"
    steps = run_list.get("steps")
    resume_after = run_list.get("resume_after", None)
    if resume_after:
        print("resume_after", resume_after)
    pipeline.run(models=steps, resume_after=resume_after)
    pipeline.close_pipeline()
    return ("all models", t0)


def main(argv=None):
    run_pipeline(argv)


if __name__ == "__main__":
    main()
