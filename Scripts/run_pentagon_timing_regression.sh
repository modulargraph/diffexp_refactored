#!/usr/bin/env bash
# Opt-in end-to-end massless-pentagon performance gate.
#
# FIRE preparation is deliberately outside the timed region.  The measured
# run always starts with a fresh transport checkpoint directory and seeds the
# exact private matching halo discovered by the July 2026 recovery, so this
# gate measures the stable production transport path rather than discovery
# retries or stale checkpoint reuse.

set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
cd "$repo_root"

max_seconds=${PENTAGON_TIMING_MAX_SECONDS:-300}
working_precision=${PENTAGON_TIMING_WORKING_PRECISION:-300}
matching_digits=${PENTAGON_TIMING_MATCH_DIGITS:-8}
expansion_order=${PENTAGON_TIMING_EXPANSION_ORDER:-25}
boundary_extra_order=${PENTAGON_TIMING_BOUNDARY_EXTRA_ORDER:-16}
division_order=${PENTAGON_TIMING_DIVISION_ORDER:-3}
cpp_threads=${DE2_CPP_THREADS:-10}
prep_cache=${FT_PREP_CACHE_DIR:-${TMPDIR:-/tmp}/DiffExp2_FT_Prepared}
fire_path=${FT_FIRE_PATH:-$repo_root/Dependencies/fire/FIRE7/FIRE7}

case "$max_seconds" in
  ''|*[!0-9]*)
    echo "PENTAGON_TIMING_MAX_SECONDS must be a positive integer" >&2
    exit 2
    ;;
esac
if (( max_seconds < 1 )); then
  echo "PENTAGON_TIMING_MAX_SECONDS must be a positive integer" >&2
  exit 2
fi

if ! command -v wolframscript >/dev/null 2>&1; then
  echo "wolframscript is required for the pentagon timing regression" >&2
  exit 2
fi
if ! command -v python3 >/dev/null 2>&1; then
  echo "python3 is required to verify the pentagon timing result" >&2
  exit 2
fi
if [[ ! -e "$fire_path" ]]; then
  echo "FIRE7 was not found at $fire_path; set FT_FIRE_PATH" >&2
  exit 2
fi

library=${DE2_CPP_LIBRARY:-}
if [[ -z "$library" ]]; then
  for candidate in \
    "$repo_root/build/cpp/diffexp2_librarylink.dylib" \
    "$repo_root/build/cpp/diffexp2_librarylink.so" \
    "$repo_root/build-recovery/cpp/diffexp2_librarylink.dylib" \
    "$repo_root/build-recovery/cpp/diffexp2_librarylink.so"; do
    if [[ -f "$candidate" ]]; then
      library=$candidate
      break
    fi
  done
fi
if [[ -z "$library" || ! -f "$library" ]]; then
  echo "compiled DiffExp2 LibraryLink backend not found; set DE2_CPP_LIBRARY" >&2
  exit 2
fi

scratch=$(mktemp -d "${TMPDIR:-/tmp}/diffexp2-pentagon-timing.XXXXXX")
trap 'rm -rf "$scratch"' EXIT

common_environment=(
  env
  "DE2_CPP_LIBRARY=$library"
  "DE2_RECURRENCE_BACKEND=Cpp"
  "DE2_CPP_THREADS=$cpp_threads"
  "DE2_VALUE_TRANSPORT=1"
  "FT_CPP_BATCH_ENDPOINT_ARMS=1"
  "FT_EXAMPLES=pentagon"
  "FT_FIRE_PATH=$fire_path"
  "FT_PREP_CACHE_DIR=$prep_cache"
  "FT_WORKING_PRECISION=$working_precision"
  "FT_MATCH_DIGITS=$matching_digits"
  "FT_EXPANSION_ORDER=$expansion_order"
  "FT_EPS_ORDER=0"
  "FT_BOUNDARY_EXTRA_ORDER=$boundary_extra_order"
  "FT_LEVEL_EPS_HALOS=0"
  "FT_DIVISION_ORDER=$division_order"
  "FT_REBUILD_PREP=0"
  "FT_ALLOW_STALE_LADDER_CHECKPOINT=0"
  "FT_RESUME_LADDER_CHECKPOINT="
)

if [[ ${PENTAGON_TIMING_WARM_CACHE:-0} == 1 ]]; then
  echo "=== warming pentagon FIRE preparation cache (not timed)"
  "${common_environment[@]}" \
    "FT_STOP_AFTER_BOUNDARY_LEVEL=4" \
    "FT_LADDER_CHECKPOINT_DIR=$scratch/warm-checkpoints" \
    wolframscript -file Scripts/run_ft_stepwise2.m \
    >"$scratch/warm.log"
  if ! grep -Eq 'FTPREP CACHE (HIT|WRITE)' "$scratch/warm.log"; then
    echo "pentagon cache warmup did not produce or reuse a preparation cache" >&2
    tail -80 "$scratch/warm.log" >&2
    exit 1
  fi
fi

# The halo is private to the level-1 basis match.  Calling runExample directly
# keeps the public epsilon request at order zero while avoiding a deliberately
# expensive discovery replay.  A reservoir-regression test separately guards
# that widening this private window preserves the old physical prefix.
runner_code='Get[FileNameJoin[{Directory[], "Scripts", "run_ft_stepwise2.m"}]]; result = Global`runExample["pentagon", None, <|1 -> 2|>]; If[result === True, Quit[0], Print["FAILED pentagon: ", InputForm[result]]; Quit[1]]'

echo "=== pentagon timing regression"
echo "configuration: WP=$working_precision match=$matching_digits T=$expansion_order boundaryExtra=$boundary_extra_order division=$division_order threads=$cpp_threads ceiling=${max_seconds}s"
SECONDS=0
set +e
python3 Scripts/run_with_deadline.py \
  "$max_seconds" "$scratch/pentagon.log" -- \
  "${common_environment[@]}" \
  "FT_RUNNER_DEFINITIONS_ONLY=1" \
  "FT_STOP_AFTER_BOUNDARY_LEVEL=" \
  "FT_LADDER_CHECKPOINT_DIR=$scratch/timed-checkpoints" \
  wolframscript -code "$runner_code"
runner_status=$?
set -e
elapsed=$SECONDS

if (( runner_status != 0 )); then
  if (( runner_status == 124 )); then
    echo "pentagon timing hard deadline exceeded after ${max_seconds}s" >&2
  fi
  echo "pentagon timing runner failed with status $runner_status" >&2
  exit "$runner_status"
fi
if ! grep -q 'FTPREP CACHE HIT ' "$scratch/pentagon.log"; then
  echo "timed pentagon run was not a preparation-cache hit" >&2
  echo "rerun once with PENTAGON_TIMING_WARM_CACHE=1; FIRE preparation must not be charged to this gate" >&2
  exit 1
fi

python3 - "$scratch/pentagon.log" <<'PY'
from decimal import Decimal
import json
import pathlib
import sys

log_path = pathlib.Path(sys.argv[1])
lines = log_path.read_text(encoding="utf-8").splitlines()
final_lines = [line[len("FINAL "):] for line in lines if line.startswith("FINAL ")]
if len(final_lines) != 1:
    raise SystemExit(f"expected exactly one FINAL record, found {len(final_lines)}")

record = json.loads(final_lines[0], parse_float=Decimal, parse_int=int)
if record.get("Example") != "pentagon" or record.get("RawMinPower") != -2:
    raise SystemExit(f"unexpected pentagon FINAL record: {record}")
finite = record.get("Finite")
if not isinstance(finite, dict):
    raise SystemExit(f"pentagon finite value is not complex: {finite}")
real = finite.get("Re")
imag = abs(finite.get("Im", Decimal(0)))
if not isinstance(real, Decimal):
    real = Decimal(str(real))
if not isinstance(imag, Decimal):
    imag = Decimal(str(imag))

reference = Decimal("-0.025411779885306218280920810082412842534323133089462")
tolerance = Decimal("1e-12")
error = abs(real - reference)
if error > tolerance or imag > tolerance:
    raise SystemExit(
        f"pentagon result failed oracle: value={real} imag={imag} error={error}"
    )
if not any("FTLADDER NATIVE BATCH level=1" in line for line in lines):
    raise SystemExit("pentagon log has no completed native level-1 batch")
if any("FTLADDER NATIVE MATCH RETRY " in line for line in lines):
    raise SystemExit("pentagon timing gate unexpectedly entered matching-halo discovery")
print(f"pentagon oracle PASS: value={real} error={error}")
PY

if (( elapsed > max_seconds )); then
  echo "pentagon timing FAIL: ${elapsed}s exceeds ${max_seconds}s" >&2
  exit 1
fi

echo "pentagon timing PASS: ${elapsed}s <= ${max_seconds}s"
