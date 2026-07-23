#!/usr/bin/env bash
# Opt-in production-chain regression through a selected banana4 boundary, or
# through final publication with an independent Bessel-value check.  This is
# intentionally excluded from ordinary per-change tests:
# it exercises Wolfram manifest construction, LibraryLink JSON forwarding,
# exact Rational-shadow ownership/checkpointing, the singular physical ODE,
# and authoritative native terminal composed-adjoint selection at expansion
# order 50.  The latter must adapt its private adjoint Taylor order until its
# rigorous enclosure meets the retained matching tolerance; containment alone
# is not sufficient for a value consumed by the next ladder level.
#
# FIRE preparation is outside the measured transport contract.  Populate the
# ordinary preparation cache once with BANANA4_BOUNDARY_WARM_CACHE=1.

set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
cd "$repo_root"

target_level=${BANANA4_BOUNDARY_TARGET_LEVEL:-2}
if [[ "$target_level" == final || "$target_level" == 1 ]]; then
  max_seconds=${BANANA4_BOUNDARY_MAX_SECONDS:-2400}
  matching_digits=${BANANA4_BOUNDARY_MATCH_DIGITS:-15}
else
  max_seconds=${BANANA4_BOUNDARY_MAX_SECONDS:-420}
  matching_digits=${BANANA4_BOUNDARY_MATCH_DIGITS:-25}
fi
working_precision=${BANANA4_BOUNDARY_WORKING_PRECISION:-500}
expansion_order=${BANANA4_BOUNDARY_EXPANSION_ORDER:-50}
adjoint_order=${BANANA4_BOUNDARY_ADJOINT_ORDER:-100}
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
if [[ "$target_level" != 1 && "$target_level" != 2 &&
      "$target_level" != 3 &&
      "$target_level" != final ]]; then
  echo "BANANA4_BOUNDARY_TARGET_LEVEL must be 1, 2, 3, or final" >&2
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
preserved_log=${BANANA4_BOUNDARY_LOG_PATH:-}
keep_on_failure=${BANANA4_BOUNDARY_KEEP_ON_FAILURE:-0}
case "$keep_on_failure" in
  0|1) ;;
  *)
    echo "BANANA4_BOUNDARY_KEEP_ON_FAILURE must be 0 or 1" >&2
    exit 2
    ;;
esac
cleanup() {
  status=$?
  if [[ -n "$preserved_log" &&
        -f "$scratch/banana4-boundary.log" ]]; then
    cp "$scratch/banana4-boundary.log" "$preserved_log"
  fi
  if [[ "$keep_on_failure" == 1 && "$status" != 0 ]]; then
    echo "preserved failed banana4 artifacts: $scratch" >&2
    return
  fi
  rm -rf "$scratch"
}
trap cleanup EXIT

common_environment=(
  env
  "DE2_CPP_LIBRARY=$library"
  "DE2_RECURRENCE_BACKEND=Cpp"
  "DE2_CPP_THREADS=$cpp_threads"
  "DE2_VALUE_TRANSPORT=1"
  "DE2_NATIVE_STAGE_TIMING=1"
  "DE2_DIAGNOSTIC_TERMINAL_COMPOSED_ADJOINT=authoritative"
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
  # Public output halos are zero.  The larger values below are private
  # matching reservoirs and must not be counted a second time in the public
  # epsilon request.
  "FT_LEVEL_EPS_HALOS=0,0,0,0"
  "FT_DIVISION_ORDER=3"
  "FT_REBUILD_PREP=0"
  "FT_ALLOW_STALE_LADDER_CHECKPOINT=0"
  "FT_RESUME_LADDER_CHECKPOINT="
  "FT_DISABLE_MATCHING_HALO_PROFILE=1"
)

if [[ ${BANANA4_BOUNDARY_WARM_CACHE:-0} == 1 ]]; then
  echo "=== warming banana4 FIRE preparation cache (not timed)"
  "${common_environment[@]}" \
    "FT_PREPARATION_ONLY=1" \
    wolframscript -file Scripts/run_ft_stepwise2.m \
    >"$scratch/warm.log"
  if ! grep -Eq 'FTPREP CACHE (HIT|WRITE)' "$scratch/warm.log" ||
      ! grep -q 'FTPREP ONLY COMPLETE ' "$scratch/warm.log"; then
    echo "banana4 cache warmup did not produce or reuse a preparation cache" >&2
    tail -80 "$scratch/warm.log" >&2
    exit 1
  fi
fi

# These are matching-only reservoirs.  The public epsilon request remains 4.
# Calling runExample directly prevents discovery retries from changing the
# experiment being guarded.
runner_code='Get[FileNameJoin[{Directory[], "Scripts", "run_ft_stepwise2.m"}]]; result = Global`runExample["banana4", None, <|1 -> 20, 2 -> 14, 3 -> 8, 4 -> 0|>]; If[result === True, Quit[0], Print["FAILED banana4 boundary: ", InputForm[result]]; Quit[1]]'

if [[ "$target_level" == final ]]; then
  stop_after_boundary=
  gate_label="final-value"
else
  stop_after_boundary=$target_level
  gate_label="level-$target_level boundary"
fi

echo "=== banana4 $gate_label regression"
echo "configuration: WP=$working_precision match=$matching_digits T=$expansion_order adjointT=$adjoint_order threads=$cpp_threads ceiling=${max_seconds}s"
SECONDS=0
set +e
python3 Scripts/run_with_deadline.py \
  "$max_seconds" "$scratch/banana4-boundary.log" -- \
  "${common_environment[@]}" \
  "FT_RUNNER_DEFINITIONS_ONLY=1" \
  "FT_STOP_AFTER_BOUNDARY_LEVEL=$stop_after_boundary" \
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
from decimal import Decimal
import json
import pathlib
import sys

lines = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8").splitlines()
target = sys.argv[2]
target_level = 0 if target == "final" else int(target)

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
    if not any(
        "status=authoritative-certified-value" in line
        for line in lines[segment_begin:batch_index]
    ):
        raise SystemExit(
            f"missing accuracy-qualified authoritative composed value in level {source_level}"
        )
    segment_begin = batch_index + 1

if target == "final":
    final_lines = [
        line[len("FINAL "):] for line in lines if line.startswith("FINAL ")
    ]
    if len(final_lines) != 1:
        raise SystemExit(
            f"expected exactly one banana4 FINAL record, found {len(final_lines)}"
        )
    record = json.loads(
        final_lines[0], parse_float=Decimal, parse_int=int
    )
    if record.get("Example") != "banana4" or record.get("RawMinPower") != 0:
        raise SystemExit(f"unexpected banana4 FINAL record: {record}")
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
    reference = Decimal(
        "39.655526834297652529992823046933581156446060218710"
    )
    tolerance = Decimal("1e-8")
    error = abs(real - reference)
    if error > tolerance or imag > tolerance:
        raise SystemExit(
            "banana4 result failed eight-digit Bessel oracle: "
            f"value={real} imag={imag} error={error}"
        )
    if not any("FTLADDER NATIVE BATCH level=1" in line for line in lines):
        raise SystemExit("banana4 log has no completed native level-1 batch")
    print(f"banana4 Bessel oracle PASS: value={real} error={error}")
else:
    require(
        f"STOPPED_AFTER_BOUNDARY_LEVEL {target_level}",
        "intentional boundary stop",
    )

for forbidden in (
    "FTLADDER NATIVE BATCH FAIL",
    "FTLADDER NATIVE MATCH RETRY level=",
    "FTLADDER NATIVE MATCH TAYLOR RETRY STALLED",
    "FTLADDER NATIVE MATCH TAYLOR RETRY EXHAUSTED",
    "status=authoritative-fallback-insufficient-accuracy",
    "status=unsupported",
    "FAILED banana4 boundary",
):
    if any(forbidden in line for line in lines):
        raise SystemExit(f"forbidden failure marker in banana4 boundary log: {forbidden}")

# Crossing level 3 is the regression that exposed locally successful Acb
# columns losing their exact right-frame correlation only at the downstream
# singular match.  The producer certificate must therefore select and retain
# the Rational shadow before any speculative Acb column solve.
if target_level <= 2:
    require(
        "capability=producer-certified-proactive-rational-shadow-v1",
        "proactive exact Rational-shadow singular basis",
    )

# Physical value transport is deliberately ineligible at a singular chart,
# where a basis crossing remains necessary.  An ordinary basis solve can lose
# the private epsilon reservoir needed by a later regular value hop and
# previously turned one tail miss into a long fallback cascade.  It is
# harmless at the last ordinary chart directly before a singular receiver:
# that receiver necessarily starts a fresh basis crossing, with no intervening
# value hop consuming the old frame.  Reject every other regular fallback.
import re

basis_start = re.compile(
    r"DE2 NATIVE STAGE stream-basis-start arm=(\S+) index=(\d+)"
)
basis_done = re.compile(
    r"DE2 NATIVE STAGE stream-basis-done arm=(\S+) index=(\d+).* regular=(True|False)"
)
pending_basis = {}
basis_completions = []
for line in lines:
    start_match = basis_start.search(line)
    if start_match:
        key = start_match.groups()
        pending_basis[key] = pending_basis.get(key, 0) + 1
        continue
    done_match = basis_done.search(line)
    if not done_match:
        continue
    arm, index, regular = done_match.groups()
    key = (arm, index)
    if pending_basis.get(key, 0) < 1:
        raise SystemExit(
            f"basis completion has no matching start: arm={arm} index={index}"
        )
    pending_basis[key] -= 1
    basis_completions.append((arm, index, regular))
if any(count for count in pending_basis.values()):
    raise SystemExit("banana4 log ended with an unfinished basis fallback")

completion_set = set(basis_completions)
for arm, index, regular in basis_completions:
    if regular == "True" and (
        arm, str(int(index) + 1), "False"
    ) not in completion_set:
        raise SystemExit(
            "ordinary banana4 chart used basis fallback: "
            f"arm={arm} index={index}"
        )

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
        or (
            "detail=q-epsilon0-is-not-center-unit" in line
            and "q_t_valuation=" in line
        )
    )
    for line in noncenter
):
    raise SystemExit("unclassified non-applicable terminal composed-adjoint tile")

print(f"banana4 {target} production-chain PASS")
PY

if (( elapsed > max_seconds )); then
  echo "banana4 boundary timing FAIL: ${elapsed}s exceeds ${max_seconds}s" >&2
  exit 1
fi

echo "banana4 boundary timing PASS: ${elapsed}s <= ${max_seconds}s"
