#!/usr/bin/env bash
# Opt-in correctness and timing gate for the unequal-mass part of the
# original DiffExp Banana.nb example. It is intentionally not run on every
# change.

set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
cd "$repo_root"

max_seconds=${ORIGINAL_BANANA_UNEQUAL_DEADLINE_SECONDS:-600}
max_total_seconds=${ORIGINAL_BANANA_UNEQUAL_MAX_TOTAL_SECONDS:-420}
data_directory=${ORIGINAL_BANANA_DATA:-"$repo_root/Examples/OriginalDiffExp/Data/Banana"}

required=(
  EqualMass/dt_0.m
  EqualMass/dt_1.m
  UnequalMass/dmm1_0.m
  UnequalMass/dmm1_1.m
  UnequalMass/dmm2_0.m
  UnequalMass/dmm2_1.m
  UnequalMass/dmm3_0.m
  UnequalMass/dmm3_1.m
)
for file in "${required[@]}"; do
  if [[ ! -f "$data_directory/$file" ]]; then
    echo "missing $data_directory/$file" >&2
    echo "run Scripts/fetch_original_banana_data.sh first" >&2
    exit 2
  fi
done

scratch=$(mktemp -d "${TMPDIR:-/tmp}/diffexp2-original-banana-unequal.XXXXXX")
trap 'rm -rf "$scratch"' EXIT

python3 Scripts/run_with_deadline.py \
  "$max_seconds" "$scratch/original-banana-unequal.log" -- \
  env \
    "ORIGINAL_BANANA_DATA=$data_directory" \
    "ORIGINAL_BANANA_UNEQUAL_WORKING_PRECISION=${ORIGINAL_BANANA_UNEQUAL_WORKING_PRECISION:-100}" \
    "ORIGINAL_BANANA_UNEQUAL_EXPANSION_ORDER=${ORIGINAL_BANANA_UNEQUAL_EXPANSION_ORDER:-50}" \
    wolframscript -file \
      Examples/OriginalDiffExp/BananaUnequalMass.wl

if ! grep -q '^ORIGINAL_BANANA_UNEQUAL PASS$' \
    "$scratch/original-banana-unequal.log"; then
  echo "original unequal-mass Banana correctness marker is missing" >&2
  tail -80 "$scratch/original-banana-unequal.log" >&2
  exit 1
fi

total_seconds=$(
  sed -n 's/.*totalSeconds=\([^ ]*\).*/\1/p' \
    "$scratch/original-banana-unequal.log" | tail -1
)
python3 - "$total_seconds" "$max_total_seconds" <<'PY'
import sys

elapsed = float(sys.argv[1])
ceiling = float(sys.argv[2])
if elapsed > ceiling:
    raise SystemExit(
        "original unequal-mass Banana timing regression: "
        f"{elapsed:.3f}s > {ceiling:.3f}s"
    )
print(
    "original unequal-mass Banana timing passed: "
    f"{elapsed:.3f}s <= {ceiling:.3f}s"
)
PY
