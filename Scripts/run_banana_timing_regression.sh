#!/usr/bin/env bash
# Opt-in end-to-end banana performance gate.
#
# FIRE preparation is deliberately outside the timed region: its cache and
# external-process behavior are independent of the DiffExp2 transport
# regression this test guards. Set BANANA_TIMING_WARM_CACHE=1 once on a new
# machine to populate FT_PREP_CACHE_DIR before the measured run.

set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
cd "$repo_root"

max_seconds=${BANANA_TIMING_MAX_SECONDS:-300}
working_precision=${BANANA_TIMING_WORKING_PRECISION:-500}
matching_digits=${BANANA_TIMING_MATCH_DIGITS:-20}
expansion_order=${BANANA_TIMING_EXPANSION_ORDER:-50}
division_order=${BANANA_TIMING_DIVISION_ORDER:-3}
cpp_threads=${DE2_CPP_THREADS:-10}
prep_cache=${FT_PREP_CACHE_DIR:-${TMPDIR:-/tmp}/DiffExp2_FT_Prepared}
fire_path=${FT_FIRE_PATH:-$repo_root/Dependencies/fire/FIRE7/FIRE7}

case "$max_seconds" in
  ''|*[!0-9]*)
    echo "BANANA_TIMING_MAX_SECONDS must be a positive integer" >&2
    exit 2
    ;;
esac
if (( max_seconds < 1 )); then
  echo "BANANA_TIMING_MAX_SECONDS must be a positive integer" >&2
  exit 2
fi

if ! command -v wolframscript >/dev/null 2>&1; then
  echo "wolframscript is required for the banana timing regression" >&2
  exit 2
fi
if ! command -v python3 >/dev/null 2>&1; then
  echo "python3 is required to verify the banana timing result" >&2
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

scratch=$(mktemp -d "${TMPDIR:-/tmp}/diffexp2-banana-timing.XXXXXX")
trap 'rm -rf "$scratch"' EXIT

common_environment=(
  env
  "DE2_CPP_LIBRARY=$library"
  "DE2_RECURRENCE_BACKEND=Cpp"
  "DE2_CPP_THREADS=$cpp_threads"
  "FT_EXAMPLES=banana"
  "FT_FIRE_PATH=$fire_path"
  "FT_PREP_CACHE_DIR=$prep_cache"
  "FT_WORKING_PRECISION=$working_precision"
  "FT_MATCH_DIGITS=$matching_digits"
  "FT_EXPANSION_ORDER=$expansion_order"
  "FT_EPS_ORDER=0"
  "FT_BOUNDARY_EXTRA_ORDER=4"
  "FT_LEVEL_EPS_HALOS=0"
  "FT_DIVISION_ORDER=$division_order"
  "FT_REBUILD_PREP=0"
  "FT_ALLOW_STALE_LADDER_CHECKPOINT=0"
  "FT_RESUME_LADDER_CHECKPOINT="
)

if [[ ${BANANA_TIMING_WARM_CACHE:-0} == 1 ]]; then
  echo "=== warming FIRE preparation cache (not timed)"
  "${common_environment[@]}" \
    "FT_STOP_AFTER_BOUNDARY_LEVEL=2" \
    "FT_LADDER_CHECKPOINT_DIR=$scratch/warm-checkpoints" \
    wolframscript -file Scripts/run_ft_stepwise2.m \
    >"$scratch/warm.log"
  if ! grep -Eq 'FTPREP CACHE (HIT|WRITE)' "$scratch/warm.log"; then
    echo "banana cache warmup did not produce or reuse a preparation cache" >&2
    tail -80 "$scratch/warm.log" >&2
    exit 1
  fi
fi

if [[ ${BANANA_TIMING_SKIP_SINGULAR_GATE:-0} != 1 ]]; then
  echo "=== banana singular endpoint timing regression (diagnostic, not charged to full-run timer)"
  "${common_environment[@]}" \
    "BANANA_SINGULAR_MAX_SECONDS=${BANANA_SINGULAR_MAX_SECONDS:-60}" \
    wolframscript -file Tests/benchmark_banana_singular.m
fi

echo "=== banana timing regression"
echo "configuration: WP=$working_precision match=$matching_digits T=$expansion_order division=$division_order threads=$cpp_threads ceiling=${max_seconds}s"
SECONDS=0
set +e
"${common_environment[@]}" \
  "FT_STOP_AFTER_BOUNDARY_LEVEL=" \
  "FT_LADDER_CHECKPOINT_DIR=$scratch/timed-checkpoints" \
  wolframscript -file Scripts/run_ft_stepwise2.m \
  2>&1 | tee "$scratch/banana.log"
runner_status=${PIPESTATUS[0]}
set -e
elapsed=$SECONDS

if (( runner_status != 0 )); then
  echo "banana timing runner failed with status $runner_status" >&2
  exit "$runner_status"
fi
if ! grep -q 'FTPREP CACHE HIT ' "$scratch/banana.log"; then
  echo "timed banana run was not a preparation-cache hit" >&2
  echo "rerun once with BANANA_TIMING_WARM_CACHE=1; FIRE preparation must not be charged to this gate" >&2
  exit 1
fi

python3 - "$scratch/banana.log" <<'PY'
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
if record.get("Example") != "banana" or record.get("RawMinPower") != 0:
    raise SystemExit(f"unexpected banana FINAL record: {record}")
finite = record.get("Finite")
if isinstance(finite, dict):
    real = finite.get("Re")
    imag = abs(finite.get("Im", Decimal(0)))
else:
    real = finite
    imag = Decimal(0)
if not isinstance(real, Decimal):
    real = Decimal(str(real))
if not isinstance(imag, Decimal):
    imag = Decimal(str(imag))

reference = Decimal("8.2681045358689687315430153454799888687")
tolerance = Decimal("1e-20")
error = abs(real - reference)
if error > tolerance or imag > tolerance:
    raise SystemExit(
        f"banana result failed 20-digit oracle: value={real} imag={imag} error={error}"
    )
if not any("FTLADDER NATIVE BATCH level=1" in line for line in lines):
    raise SystemExit("banana log has no completed native level-1 batch")
print(f"banana oracle PASS: value={real} error={error}")
PY

if (( elapsed > max_seconds )); then
  echo "banana timing FAIL: ${elapsed}s exceeds ${max_seconds}s" >&2
  exit 1
fi

echo "banana timing PASS: ${elapsed}s <= ${max_seconds}s"
