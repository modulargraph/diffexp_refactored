#!/usr/bin/env bash
# Opt-in correctness and timing gate for the equal-mass part of the original
# DiffExp Banana.nb example. It is intentionally not run on every change.

set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
cd "$repo_root"

max_seconds=${ORIGINAL_BANANA_DEADLINE_SECONDS:-180}
max_transport_seconds=${ORIGINAL_BANANA_MAX_TRANSPORT_SECONDS:-75}
data_directory=${ORIGINAL_BANANA_DATA:-"$repo_root/Examples/OriginalDiffExp/Data/Banana"}

for file in EqualMass/dt_0.m EqualMass/dt_1.m; do
  if [[ ! -f "$data_directory/$file" ]]; then
    echo "missing $data_directory/$file" >&2
    echo "run Scripts/fetch_original_banana_data.sh first" >&2
    exit 2
  fi
done

scratch=$(mktemp -d "${TMPDIR:-/tmp}/diffexp2-original-banana.XXXXXX")
trap 'rm -rf "$scratch"' EXIT

python3 Scripts/run_with_deadline.py \
  "$max_seconds" "$scratch/original-banana.log" -- \
  env \
    "ORIGINAL_BANANA_DATA=$data_directory" \
    "ORIGINAL_BANANA_WORKING_PRECISION=${ORIGINAL_BANANA_WORKING_PRECISION:-100}" \
    "ORIGINAL_BANANA_EXPANSION_ORDER=${ORIGINAL_BANANA_EXPANSION_ORDER:-50}" \
    wolframscript -file \
      Examples/OriginalDiffExp/BananaEqualMass.wl

if ! grep -q '^ORIGINAL_BANANA PASS$' "$scratch/original-banana.log"; then
  echo "original equal-mass Banana correctness marker is missing" >&2
  tail -80 "$scratch/original-banana.log" >&2
  exit 1
fi

if ! grep -q '^ORIGINAL_BANANA route=real-prescribed ' \
    "$scratch/original-banana.log"; then
  echo "original equal-mass Banana did not use the prescribed real singular route" >&2
  tail -80 "$scratch/original-banana.log" >&2
  exit 1
fi

transport_seconds=$(
  sed -n 's/.*transportSeconds=\([^ ]*\).*/\1/p' \
    "$scratch/original-banana.log" | tail -1
)
python3 - "$transport_seconds" "$max_transport_seconds" <<'PY'
import sys

elapsed = float(sys.argv[1])
ceiling = float(sys.argv[2])
if elapsed > ceiling:
    raise SystemExit(
        "original equal-mass Banana transport regression: "
        f"{elapsed:.3f}s > {ceiling:.3f}s"
    )
print(
    "original equal-mass Banana timing passed: "
    f"{elapsed:.3f}s <= {ceiling:.3f}s"
)
PY
