#!/usr/bin/env bash
# Opt-in end-to-end planar-double-box performance gate.
#
# FIRE preparation is deliberately outside the timed region. The measured
# run always starts with a fresh transport checkpoint directory and seeds the
# exact private matching halos learned by the July 2026 recovery, so this gate
# measures the stable production transport path rather than halo discovery or
# stale checkpoint reuse.

set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
cd "$repo_root"

max_seconds=${DOUBLE_BOX_TIMING_MAX_SECONDS:-1500}
working_precision=${DOUBLE_BOX_TIMING_WORKING_PRECISION:-150}
matching_digits=${DOUBLE_BOX_TIMING_MATCH_DIGITS:-8}
expansion_order=${DOUBLE_BOX_TIMING_EXPANSION_ORDER:-25}
boundary_extra_order=${DOUBLE_BOX_TIMING_BOUNDARY_EXTRA_ORDER:-16}
division_order=${DOUBLE_BOX_TIMING_DIVISION_ORDER:-3}
cpp_threads=${DE2_CPP_THREADS:-10}
prep_cache=${FT_PREP_CACHE_DIR:-${TMPDIR:-/tmp}/DiffExp2_FT_Prepared}
fire_path=${FT_FIRE_PATH:-$repo_root/Dependencies/fire/FIRE7/FIRE7}

case "$max_seconds" in
  ''|*[!0-9]*)
    echo "DOUBLE_BOX_TIMING_MAX_SECONDS must be a positive integer" >&2
    exit 2
    ;;
esac
if (( max_seconds < 1 )); then
  echo "DOUBLE_BOX_TIMING_MAX_SECONDS must be a positive integer" >&2
  exit 2
fi

if ! command -v wolframscript >/dev/null 2>&1; then
  echo "wolframscript is required for the double-box timing regression" >&2
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

scratch=$(mktemp -d "${TMPDIR:-/tmp}/diffexp2-double-box-timing.XXXXXX")
trap 'rm -rf "$scratch"' EXIT

common_environment=(
  env
  "DE2_CPP_LIBRARY=$library"
  "DE2_RECURRENCE_BACKEND=Cpp"
  "DE2_CPP_THREADS=$cpp_threads"
  "DE2_VALUE_TRANSPORT=1"
  "FT_CPP_BATCH_ENDPOINT_ARMS=1"
  "FT_EXAMPLES=double_box_planar"
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
  "FT_DISABLE_MATCHING_HALO_PROFILE=1"
)

if [[ ${DOUBLE_BOX_TIMING_WARM_CACHE:-0} == 1 ]]; then
  echo "=== warming double-box FIRE preparation cache (not timed)"
  "${common_environment[@]}" \
    "FT_STOP_AFTER_BOUNDARY_LEVEL=4" \
    "FT_LADDER_CHECKPOINT_DIR=$scratch/warm-checkpoints" \
    wolframscript -file Scripts/run_ft_stepwise2.m \
    >"$scratch/warm.log"
  if ! grep -Eq 'FTPREP CACHE (HIT|WRITE)' "$scratch/warm.log"; then
    echo "double-box cache warmup did not produce or reuse a preparation cache" >&2
    tail -80 "$scratch/warm.log" >&2
    exit 1
  fi
fi

# These halos belong to the private level-2 and level-3 basis matches. Calling
# runExample directly keeps the public epsilon request at order zero and
# avoids charging an already-learned matching-halo discovery retry to the
# regression timer.
runner_code='Get[FileNameJoin[{Directory[], "Scripts", "run_ft_stepwise2.m"}]]; result = Global`runExample["double_box_planar", None, <|1 -> 2, 2 -> 5, 3 -> 1|>]; If[result === True, Quit[0], Print["FAILED double_box_planar: ", InputForm[result]]; Quit[1]]'

echo "=== double-box timing regression"
echo "configuration: WP=$working_precision match=$matching_digits T=$expansion_order boundaryExtra=$boundary_extra_order division=$division_order threads=$cpp_threads ceiling=${max_seconds}s"
SECONDS=0
set +e
python3 Scripts/run_with_deadline.py \
  "$max_seconds" "$scratch/double-box.log" -- \
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
    echo "double-box timing hard deadline exceeded after ${max_seconds}s" >&2
  fi
  echo "double-box timing runner failed with status $runner_status" >&2
  exit "$runner_status"
fi
if ! grep -q 'FTPREP CACHE HIT ' "$scratch/double-box.log"; then
  echo "timed double-box run was not a preparation-cache hit" >&2
  echo "rerun once with DOUBLE_BOX_TIMING_WARM_CACHE=1; FIRE preparation must not be charged to this gate" >&2
  exit 1
fi
if grep -q 'FTLADDER NATIVE MATCH RETRY ' "$scratch/double-box.log"; then
  echo "double-box timing gate unexpectedly entered matching-halo discovery" >&2
  exit 1
fi

echo "=== double-box numerical oracle (8 digits)"
wolframscript -file Scripts/verify_double_box_planar.m \
  "$scratch/double-box.log" 8

if (( elapsed > max_seconds )); then
  echo "double-box timing FAIL: ${elapsed}s exceeds ${max_seconds}s" >&2
  exit 1
fi

echo "double-box timing PASS: ${elapsed}s <= ${max_seconds}s"
