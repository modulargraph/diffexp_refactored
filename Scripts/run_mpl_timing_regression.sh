#!/usr/bin/env bash
# Opt-in correctness and timing gate for the original DiffExp
# MultiplePolylogarithms notebook, including its weight-20 stress case.

set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
cd "$repo_root"

max_seconds=${MPL_DEADLINE_SECONDS:-180}
max_weight20_seconds=${MPL_MAX_WEIGHT20_SECONDS:-120}

scratch=$(mktemp -d "${TMPDIR:-/tmp}/diffexp2-mpl-timing.XXXXXX")
trap 'rm -rf "$scratch"' EXIT

python3 Scripts/run_with_deadline.py \
  "$max_seconds" "$scratch/mpl.log" -- \
  env MPL_INCLUDE_WEIGHT20=1 \
    wolframscript -file \
      Examples/OriginalDiffExp/MultiplePolylogarithms.wl

if ! grep -q '^MPL PASS$' "$scratch/mpl.log"; then
  echo "multiple-polylogarithm correctness marker is missing" >&2
  tail -100 "$scratch/mpl.log" >&2
  exit 1
fi

weight20_seconds=$(
  sed -n \
    's/^MPL case=weight20.*wallSeconds=\([^ ]*\).*/\1/p' \
    "$scratch/mpl.log" |
    tail -1
)
python3 - "$weight20_seconds" "$max_weight20_seconds" <<'PY'
import sys

elapsed = float(sys.argv[1])
ceiling = float(sys.argv[2])
if elapsed > ceiling:
    raise SystemExit(
        f"MPL weight-20 regression: {elapsed:.3f}s > {ceiling:.3f}s"
    )
print(f"MPL weight-20 timing passed: {elapsed:.3f}s <= {ceiling:.3f}s")
PY
