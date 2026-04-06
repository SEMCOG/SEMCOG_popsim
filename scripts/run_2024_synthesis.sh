#!/usr/bin/env bash
set -uo pipefail

# Run the standard one-pass PopulationSim synthesis for the 2024 package.
#
# This wrapper keeps the command predictable for daily use:
# - reads the base run configs from the run package configs folder
# - writes model outputs under output/<YYYY-MM-DD>_<HH>_one_pass/run
# - writes wrapper logs/status under output/<YYYY-MM-DD>_<HH>_one_pass/logs
# Environment variables can override paths when a one-off run needs a
# different package or output location.
RUN_DIR="${RUN_DIR:-/home/da/RDF2055/d_drive/popsim/runs/2024_synthesis}"
CONFIG_DIR="${CONFIG_DIR:-$RUN_DIR/configs}"
DATA_DIR="${DATA_DIR:-$RUN_DIR/data}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$RUN_DIR/output}"
RUN_STAMP="${RUN_STAMP:-$(date +%Y-%m-%d_%H)}"
RUN_FOLDER_NAME="${RUN_FOLDER_NAME:-${RUN_STAMP}_one_pass}"
OUTPUT_DIR="${OUTPUT_DIR:-$OUTPUT_ROOT/$RUN_FOLDER_NAME}"
RUN_OUTPUT_DIR="${RUN_OUTPUT_DIR:-$OUTPUT_DIR/run}"
LOG_DIR="${LOG_DIR:-$OUTPUT_DIR/logs}"
VALIDATION_DIR="${VALIDATION_DIR:-$OUTPUT_DIR/validation}"
LOG_FILE="${LOG_FILE:-$LOG_DIR/run.log}"
STATUS_FILE="${STATUS_FILE:-$LOG_DIR/run.status}"

mkdir -p "$RUN_OUTPUT_DIR" "$LOG_DIR" "$VALIDATION_DIR"

# Mirror each wrapper message to both the terminal and the run log.
log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG_FILE"
}

# Record a clear status file if the shell receives Ctrl-C or a terminate signal.
on_interrupt() {
  log "Run interrupted by signal."
  printf 'status=interrupted\n' > "$STATUS_FILE"
  exit 130
}

trap on_interrupt INT TERM

rm -f "$LOG_FILE" "$STATUS_FILE"
: > "$LOG_FILE"
log "Starting PopulationSim run"
log "Run directory : $RUN_DIR"
log "Config dir    : $CONFIG_DIR"
log "Data dir      : $DATA_DIR"
log "Output dir    : $OUTPUT_DIR"
log "Run output dir: $RUN_OUTPUT_DIR"
log "Log file      : $LOG_FILE"
log "Status file   : $STATUS_FILE"

populationsim \
  -c "$CONFIG_DIR" \
  -d "$DATA_DIR" \
  -o "$RUN_OUTPUT_DIR" \
  2>&1 | tee -a "$LOG_FILE"

POP_STATUS=${PIPESTATUS[0]}

if [ "$POP_STATUS" -eq 0 ]; then
  log "PopulationSim finished successfully."
  printf 'status=success\nexit_code=0\n' > "$STATUS_FILE"
else
  log "PopulationSim failed with exit code $POP_STATUS."
  printf 'status=failed\nexit_code=%s\n' "$POP_STATUS" > "$STATUS_FILE"
fi

exit "$POP_STATUS"
