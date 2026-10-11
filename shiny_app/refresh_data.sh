#!/usr/bin/env bash
# Copy the small summary tables from the pipeline output into the app.
set -e
ROOT="${MIGRATE_DATA_DIR:-$(cd "$(dirname "$0")/.." && pwd)/data}"
SRC="$ROOT/pm25_avoidance_lag3_olag1"; DST="$(dirname "$0")/data"
mkdir -p "$DST"; n=0
for f in "$SRC"/*.csv; do
  s=$(stat -f %z "$f" 2>/dev/null || stat -c %s "$f")
  if [ "$s" -lt 50000 ]; then cp "$f" "$DST/"; n=$((n+1)); fi
done
echo "copied $n tables to $DST"
