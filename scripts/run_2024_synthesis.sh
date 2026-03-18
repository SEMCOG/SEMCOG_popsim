#!/usr/bin/env bash
set -euo pipefail

RUN_DIR="${RUN_DIR:-/home/da/RDF2055/d_drive/popsim/runs/2024_synthesis}"
CONFIG_DIR="${CONFIG_DIR:-$RUN_DIR/configs}"
DATA_DIR="${DATA_DIR:-$RUN_DIR/data}"
OUTPUT_DIR="${OUTPUT_DIR:-$RUN_DIR/output}"
LOG_DIR="${LOG_DIR:-$OUTPUT_DIR/logs}"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_FILE:-$LOG_DIR/popsim_${TIMESTAMP}.log}"

mkdir -p "$LOG_DIR"

echo "Run directory : $RUN_DIR"
echo "Config dir    : $CONFIG_DIR"
echo "Data dir      : $DATA_DIR"
echo "Output dir    : $OUTPUT_DIR"
echo "Log file      : $LOG_FILE"

populationsim \
  -c "$CONFIG_DIR" \
  -d "$DATA_DIR" \
  -o "$OUTPUT_DIR" \
  2>&1 | tee "$LOG_FILE"
