#!/usr/bin/env bash
# Opt-in correctness and timing gate for the original DiffExp/Henn 108-master
# nonplanar canonical pentagon. It is intentionally not run on every change.

set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
cd "$repo_root"

max_seconds=${HENN_NONPLANAR_DEADLINE_SECONDS:-120}
max_transport_seconds=${HENN_NONPLANAR_MAX_TRANSPORT_SECONDS:-30}
data_directory=${HENN_NONPLANAR_DATA:-"$repo_root/Examples/OriginalDiffExp/Data/HennNonplanar"}

for file in XB_Atilde.txt XB_Boundary_values_X0.txt XB_Boundary_values_X1.txt; do
  if [[ ! -f "$data_directory/$file" ]]; then
    echo "missing $data_directory/$file" >&2
    echo "run Scripts/fetch_henn_nonplanar_data.sh first" >&2
    exit 2
  fi
done

scratch=$(mktemp -d "${TMPDIR:-/tmp}/diffexp2-henn-timing.XXXXXX")
trap 'rm -rf "$scratch"' EXIT

python3 Scripts/run_with_deadline.py \
  "$max_seconds" "$scratch/henn.log" -- \
  env \
    "HENN_NONPLANAR_DATA=$data_directory" \
    "HENN_NONPLANAR_WORKING_PRECISION=${HENN_NONPLANAR_WORKING_PRECISION:-50}" \
    "HENN_NONPLANAR_EXPANSION_ORDER=${HENN_NONPLANAR_EXPANSION_ORDER:-25}" \
    wolframscript -file \
      Examples/OriginalDiffExp/HennNonplanarCanonical.wl

if ! grep -q '^HENN_NONPLANAR PASS$' "$scratch/henn.log"; then
  echo "Henn nonplanar correctness marker is missing" >&2
  tail -80 "$scratch/henn.log" >&2
  exit 1
fi

transport_seconds=$(
  sed -n 's/.*transportSeconds=\([^ ]*\).*/\1/p' "$scratch/henn.log" |
    tail -1
)
python3 - "$transport_seconds" "$max_transport_seconds" <<'PY'
import sys

elapsed = float(sys.argv[1])
ceiling = float(sys.argv[2])
if elapsed > ceiling:
    raise SystemExit(
        f"Henn nonplanar transport regression: {elapsed:.3f}s > {ceiling:.3f}s"
    )
print(f"Henn nonplanar timing passed: {elapsed:.3f}s <= {ceiling:.3f}s")
PY
