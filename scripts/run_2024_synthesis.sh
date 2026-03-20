#!/usr/bin/env bash
set -uo pipefail

RUN_DIR="${RUN_DIR:-/home/da/RDF2055/d_drive/popsim/runs/2024_synthesis}"
CONFIG_DIR="${CONFIG_DIR:-$RUN_DIR/configs}"
DATA_DIR="${DATA_DIR:-$RUN_DIR/data}"
OUTPUT_DIR="${OUTPUT_DIR:-$RUN_DIR/output}"
LOG_DIR="${LOG_DIR:-$OUTPUT_DIR/logs}"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_FILE:-$LOG_DIR/popsim_${TIMESTAMP}.log}"
STATUS_FILE="${STATUS_FILE:-$LOG_DIR/popsim_${TIMESTAMP}.status}"

ARCHIVE_DIR="${ARCHIVE_DIR:-$OUTPUT_DIR/archive}"

archive_existing_output() {
  mkdir -p "$ARCHIVE_DIR"
  RUN_INDEX=1
  while [ -d "$ARCHIVE_DIR/run_$RUN_INDEX" ]; do
    RUN_INDEX=$((RUN_INDEX + 1))
  done
  RUN_ARCHIVE="$ARCHIVE_DIR/run_$RUN_INDEX"

  FOUND_EXISTING=0
  for ENTRY in "$OUTPUT_DIR"/*; do
    [ -e "$ENTRY" ] || continue
    BASENAME="$(basename "$ENTRY")"
    if [ "$BASENAME" = "archive" ] || [ "$BASENAME" = "two_pass" ]; then
      continue
    fi
    if [ "$FOUND_EXISTING" -eq 0 ]; then
      mkdir -p "$RUN_ARCHIVE"
      FOUND_EXISTING=1
    fi
    mv "$ENTRY" "$RUN_ARCHIVE/"
  done

  if [ "$FOUND_EXISTING" -eq 1 ]; then
    printf '%s' "$RUN_ARCHIVE"
  fi
}

mkdir -p "$OUTPUT_DIR"
ARCHIVED_RUN="$(archive_existing_output)"
mkdir -p "$LOG_DIR"

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG_FILE"
}

on_interrupt() {
  log "Run interrupted by signal."
  printf 'status=interrupted\n' > "$STATUS_FILE"
  exit 130
}

trap on_interrupt INT TERM

: > "$LOG_FILE"
log "Starting PopulationSim run"
if [ -n "$ARCHIVED_RUN" ]; then
  log "Archived previous one-pass output to $ARCHIVED_RUN"
fi
log "Run directory : $RUN_DIR"
log "Config dir    : $CONFIG_DIR"
log "Data dir      : $DATA_DIR"
log "Output dir    : $OUTPUT_DIR"
log "Log file      : $LOG_FILE"
log "Status file   : $STATUS_FILE"

populationsim \
  -c "$CONFIG_DIR" \
  -d "$DATA_DIR" \
  -o "$OUTPUT_DIR" \
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
