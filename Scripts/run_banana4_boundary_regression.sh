#!/usr/bin/env bash
# Opt-in production-chain regression through the difficult banana4 level-2
# boundary.  This is intentionally excluded from ordinary per-change tests:
# it exercises Wolfram manifest construction, LibraryLink JSON forwarding,
# exact Rational-shadow ownership/checkpointing, the singular physical ODE,
# and the native terminal composed-adjoint certificate at expansion order 50.
#
# FIRE preparation is outside the measured transport contract.  Populate the
# ordinary preparation cache once with BANANA4_BOUNDARY_WARM_CACHE=1.

set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
cd "$repo_root"

max_seconds=${BANANA4_BOUNDARY_MAX_SECONDS:-900}
working_precision=${BANANA4_BOUNDARY_WORKING_PRECISION:-300}
matching_digits=${BANANA4_BOUNDARY_MATCH_DIGITS:-80}
expansion_order=${BANANA4_BOUNDARY_EXPANSION_ORDER:-50}
adjoint_order=${BANANA4_BOUNDARY_ADJOINT_ORDER:-100}
target_level=${BANANA4_BOUNDARY_TARGET_LEVEL:-2}
cpp_threads=${DE2_CPP_THREADS:-10}
prep_cache=${FT_PREP_CACHE_DIR:-${TMPDIR:-/tmp}/DiffExp2_FT_Prepared}
fire_path=${FT_FIRE_PATH:-$repo_root/Dependencies/fire/FIRE7/FIRE7}

case "$max_seconds" in
  ''|*[!0-9]*)
    echo "BANANA4_BOUNDARY_MAX_SECONDS must be a positive integer" >&2
    exit 2
    ;;
esac
if (( max_seconds < 1 )); then
  echo "BANANA4_BOUNDARY_MAX_SECONDS must be a positive integer" >&2
  exit 2
fi
case "$adjoint_order" in
  ''|*[!0-9]*)
    echo "BANANA4_BOUNDARY_ADJOINT_ORDER must be an integer at least as large as the expansion order" >&2
    exit 2
    ;;
esac
if (( adjoint_order < expansion_order )); then
  echo "BANANA4_BOUNDARY_ADJOINT_ORDER must be an integer at least as large as the expansion order" >&2
  exit 2
fi
if [[ "$target_level" != 2 && "$target_level" != 3 ]]; then
  echo "BANANA4_BOUNDARY_TARGET_LEVEL must be 2 or 3" >&2
  exit 2
fi

if ! command -v wolframscript >/dev/null 2>&1; then
  echo "wolframscript is required for the banana4 boundary regression" >&2
  exit 2
fi
if ! command -v python3 >/dev/null 2>&1; then
  echo "python3 is required to verify the banana4 boundary regression" >&2
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

scratch=$(mktemp -d "${TMPDIR:-/tmp}/diffexp2-banana4-boundary.XXXXXX")
trap 'rm -rf "$scratch"' EXIT

common_environment=(
  env
  "DE2_CPP_LIBRARY=$library"
  "DE2_RECURRENCE_BACKEND=Cpp"
  "DE2_CPP_THREADS=$cpp_threads"
  "DE2_VALUE_TRANSPORT=1"
  "DE2_DIAGNOSTIC_TERMINAL_COMPOSED_ADJOINT=require"
  "DE2_DIAGNOSTIC_TERMINAL_COMPOSED_TAYLOR_ORDER=$adjoint_order"
  "FT_CPP_BATCH_ENDPOINT_ARMS=1"
  "FT_EXAMPLES=banana4"
  "FT_FIRE_PATH=$fire_path"
  "FT_PREP_CACHE_DIR=$prep_cache"
  "FT_WORKING_PRECISION=$working_precision"
  "FT_MATCH_DIGITS=$matching_digits"
  "FT_EXPANSION_ORDER=$expansion_order"
  "FT_EPS_ORDER=4"
  "FT_BOUNDARY_EXTRA_ORDER=4"
  "FT_LEVEL_EPS_HALOS=20,14,8,0"
  "FT_DIVISION_ORDER=3"
  "FT_REBUILD_PREP=0"
  "FT_ALLOW_STALE_LADDER_CHECKPOINT=0"
  "FT_RESUME_LADDER_CHECKPOINT="
  "FT_DISABLE_MATCHING_HALO_PROFILE=1"
)

if [[ ${BANANA4_BOUNDARY_WARM_CACHE:-0} == 1 ]]; then
  echo "=== warming banana4 FIRE preparation cache (not timed)"
  "${common_environment[@]}" \
    "FT_STOP_AFTER_BOUNDARY_LEVEL=4" \
    "FT_LADDER_CHECKPOINT_DIR=$scratch/warm-checkpoints" \
    wolframscript -file Scripts/run_ft_stepwise2.m \
    >"$scratch/warm.log"
  if ! grep -Eq 'FTPREP CACHE (HIT|WRITE)' "$scratch/warm.log"; then
    echo "banana4 cache warmup did not produce or reuse a preparation cache" >&2
    tail -80 "$scratch/warm.log" >&2
    exit 1
  fi
fi

# These are matching-only reservoirs.  The public epsilon request remains 4.
# Calling runExample directly prevents discovery retries from changing the
# experiment being guarded.
runner_code='Get[FileNameJoin[{Directory[], "Scripts", "run_ft_stepwise2.m"}]]; result = Global`runExample["banana4", None, <|1 -> 20, 2 -> 14, 3 -> 8, 4 -> 0|>]; If[result === True, Quit[0], Print["FAILED banana4 boundary: ", InputForm[result]]; Quit[1]]'

echo "=== banana4 level-$target_level boundary regression"
echo "configuration: WP=$working_precision match=$matching_digits T=$expansion_order adjointT=$adjoint_order threads=$cpp_threads ceiling=${max_seconds}s"
SECONDS=0
set +e
python3 Scripts/run_with_deadline.py \
  "$max_seconds" "$scratch/banana4-boundary.log" -- \
  "${common_environment[@]}" \
  "FT_RUNNER_DEFINITIONS_ONLY=1" \
  "FT_STOP_AFTER_BOUNDARY_LEVEL=$target_level" \
  "FT_LADDER_CHECKPOINT_DIR=$scratch/checkpoints" \
  wolframscript -code "$runner_code"
runner_status=$?
set -e
elapsed=$SECONDS

if (( runner_status != 0 )); then
  if (( runner_status == 124 )); then
    echo "banana4 boundary hard deadline exceeded after ${max_seconds}s" >&2
  fi
  echo "banana4 boundary runner failed with status $runner_status" >&2
  exit "$runner_status"
fi
if ! grep -q 'FTPREP CACHE HIT ' "$scratch/banana4-boundary.log"; then
  echo "banana4 boundary run was not a preparation-cache hit" >&2
  echo "rerun once with BANANA4_BOUNDARY_WARM_CACHE=1; FIRE preparation is not part of this gate" >&2
  exit 1
fi

python3 - "$scratch/banana4-boundary.log" "$target_level" <<'PY'
import pathlib
import sys

lines = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8").splitlines()
target_level = int(sys.argv[2])

def require(fragment, description):
    if not any(fragment in line for line in lines):
        raise SystemExit(f"missing {description}: {fragment}")

segment_begin = 0
for source_level in range(4, target_level, -1):
    batch_fragment = f"FTLADDER NATIVE BATCH level={source_level}"
    batch_index = next(
        (index for index in range(segment_begin, len(lines))
         if batch_fragment in lines[index]),
        None,
    )
    if batch_index is None:
        raise SystemExit(f"missing completed native level-{source_level} batch")
    if not any(
        "status=certified-composed-tail" in line
        for line in lines[segment_begin:batch_index]
    ):
        raise SystemExit(
            f"missing center-ending composed-adjoint certificate in level {source_level}"
        )
    segment_begin = batch_index + 1

require(
    f"STOPPED_AFTER_BOUNDARY_LEVEL {target_level}",
    "intentional boundary stop",
)

for forbidden in (
    "FTLADDER NATIVE BATCH FAIL",
    "FTLADDER NATIVE MATCH RETRY",
    "status=unsupported",
    "FAILED banana4 boundary",
):
    if any(forbidden in line for line in lines):
        raise SystemExit(f"forbidden failure marker in banana4 boundary log: {forbidden}")

# A noncenter-ending terminal tile and a row whose forcing reaches t^0 are
# currently handled by the established production contraction.  The latter
# needs a Laurent/log adjoint plus an endpoint pairing before it belongs to
# the center-anchored one-contraction theorem.  Both cases must be classified
# explicitly; no other composed-adjoint failure may be silently skipped.
noncenter = [line for line in lines if "status=not-applicable" in line]
if noncenter and not all(
    (
        "detail=tile-does-not-end-at-chart-center" in line
        or (
            "detail=row-requires-laurent-log-center-adjoint" in line
            and "forcing_power=" in line
        )
    )
    for line in noncenter
):
    raise SystemExit("unclassified non-applicable terminal composed-adjoint tile")

print("banana4 boundary production-chain PASS")
PY

if (( elapsed > max_seconds )); then
  echo "banana4 boundary timing FAIL: ${elapsed}s exceeds ${max_seconds}s" >&2
  exit 1
fi

echo "banana4 boundary timing PASS: ${elapsed}s <= ${max_seconds}s"
