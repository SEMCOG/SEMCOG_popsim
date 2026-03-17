from activitysim.core import inject
from activitysim.core import pipeline
from activitysim.core import tracing
from activitysim.core.config import handle_standard_args
from activitysim.core.config import setting
from activitysim.core.tracing import print_elapsed_time
from populationsim import lp
from populationsim import multi_integerizer
import logging


def run_pipeline():
    handle_standard_args()
    tracing.config_logger()
    t0 = print_elapsed_time()
    logger = logging.getLogger("populationsim")
    logger.info("GROUP_BY_INCIDENCE_SIGNATURE: %s", setting("GROUP_BY_INCIDENCE_SIGNATURE"))
    logger.info("INTEGERIZE_WITH_BACKSTOPPED_CONTROLS: %s", setting("INTEGERIZE_WITH_BACKSTOPPED_CONTROLS"))
    logger.info("SUB_BALANCE_WITH_FLOAT_SEED_WEIGHTS: %s", setting("SUB_BALANCE_WITH_FLOAT_SEED_WEIGHTS"))
    logger.info("meta_control_data: %s", setting("meta_control_data"))
    logger.info("control_file_name: %s", setting("control_file_name"))
    logger.info("USE_CVXPY: %s", lp.use_cvxpy())
    logger.info("USE_SIMUL_INTEGERIZER: %s", multi_integerizer.use_simul_integerizer())
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


def main():
    run_pipeline()


if __name__ == "__main__":
    main()
