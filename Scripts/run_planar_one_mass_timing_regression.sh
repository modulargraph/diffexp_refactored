#!/usr/bin/env bash
# Opt-in correctness and timing gate for all four canonical planar one-mass
# pentagon systems from the original DiffExp repository.

set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
cd "$repo_root"

max_seconds=${PLANAR_ONE_MASS_DEADLINE_SECONDS:-900}
data_directory=${PLANAR_ONE_MASS_DATA:-"$repo_root/Examples/OriginalDiffExp/Data/PlanarOneMass"}

required=(
  alphabet.m
  1loop/diffEq-1loop.m
  1loop/numIntegrals-1loop.m
  zmz/diffEq-zmz.m
  zmz/numIntegrals-zmz.m
  mzz/diffEq-mzz.m
  mzz/numIntegrals-mzz.m
  zzz/diffEq-zzz.m
  zzz/numIntegrals-zzz.m
)
for file in "${required[@]}"; do
  if [[ ! -f "$data_directory/$file" ]]; then
    echo "missing $data_directory/$file" >&2
    echo "run Scripts/fetch_planar_one_mass_data.sh first" >&2
    exit 2
  fi
done

scratch=$(mktemp -d "${TMPDIR:-/tmp}/diffexp2-planar-one-mass-timing.XXXXXX")
trap 'rm -rf "$scratch"' EXIT

python3 Scripts/run_with_deadline.py \
  "$max_seconds" "$scratch/planar-one-mass.log" -- \
  env \
    "PLANAR_ONE_MASS_DATA=$data_directory" \
    PLANAR_ONE_MASS_FAMILY=all \
    "PLANAR_ONE_MASS_WORKING_PRECISION=${PLANAR_ONE_MASS_WORKING_PRECISION:-50}" \
    "PLANAR_ONE_MASS_EXPANSION_ORDER=${PLANAR_ONE_MASS_EXPANSION_ORDER:-25}" \
    wolframscript -file \
      Examples/OriginalDiffExp/PlanarOneMassCanonical.wl

if ! grep -q '^PLANAR_ONE_MASS PASS$' "$scratch/planar-one-mass.log"; then
  echo "planar one-mass correctness marker is missing" >&2
  tail -100 "$scratch/planar-one-mass.log" >&2
  exit 1
fi

python3 - "$scratch/planar-one-mass.log" <<'PY'
import os
import re
import sys

ceilings = {
    "1loop": float(os.environ.get("PLANAR_ONE_MASS_MAX_1LOOP_SECONDS", "60")),
    "zmz": float(os.environ.get("PLANAR_ONE_MASS_MAX_ZMZ_SECONDS", "300")),
    "mzz": float(os.environ.get("PLANAR_ONE_MASS_MAX_MZZ_SECONDS", "300")),
    "zzz": float(os.environ.get("PLANAR_ONE_MASS_MAX_ZZZ_SECONDS", "350")),
}
pattern = re.compile(
    r"^PLANAR_ONE_MASS family=(\S+).* comparableSeconds=(\S+)"
)
measured = {}
with open(sys.argv[1], encoding="utf-8") as stream:
    for line in stream:
        match = pattern.match(line)
        if match:
            measured[match.group(1)] = float(match.group(2))

missing = sorted(set(ceilings) - set(measured))
if missing:
    raise SystemExit(f"missing timing records for: {', '.join(missing)}")

for family, ceiling in ceilings.items():
    elapsed = measured[family]
    if elapsed > ceiling:
        raise SystemExit(
            f"planar one-mass {family} regression: "
            f"{elapsed:.3f}s > {ceiling:.3f}s"
        )
    print(
        f"planar one-mass {family} timing passed: "
        f"{elapsed:.3f}s <= {ceiling:.3f}s"
    )
PY
