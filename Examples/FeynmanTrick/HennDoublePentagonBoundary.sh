#!/usr/bin/env sh
# Reconstruct Henn et al.'s double-pentagon X0 canonical boundary with FT.
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
example="$repo_root/Examples/FeynmanTrick/HennDoublePentagonBoundary.wl"
runner="$repo_root/Examples/FeynmanTrick/run_henn_boundary_manifest.py"

if [ -n "${HENN_FT_WOLFRAMSCRIPT:-}" ]; then
  wolfram_script=$HENN_FT_WOLFRAMSCRIPT
elif [ -x /Applications/Wolfram.app/Contents/MacOS/wolframscript ]; then
  wolfram_script=/Applications/Wolfram.app/Contents/MacOS/wolframscript
elif [ -x "/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/wolframscript" ]; then
  wolfram_script="/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/wolframscript"
else
  wolfram_script=$(command -v wolframscript || true)
fi

if [ -z "$wolfram_script" ]; then
  echo "wolframscript is required" >&2
  exit 2
fi

for argument in "$@"; do
  if [ "$argument" = "--plan" ]; then
    exec "$wolfram_script" -script "$example" "$@"
  fi
done

scratch=$(mktemp -d "${TMPDIR:-/tmp}/diffexp2-henn-ft.XXXXXX")
cleanup() {
  status=$?
  trap - 0
  if [ "$status" -eq 0 ]; then
    rm -rf "$scratch"
  else
    echo "Henn FT run artifacts preserved at $scratch" >&2
  fi
  exit "$status"
}
trap cleanup 0
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
manifest="$scratch/manifest.json"

HENN_FT_WOLFRAMSCRIPT="$wolfram_script" \
  "$wolfram_script" -script "$example" \
    --write-manifest "$manifest" "$@"
if [ -n "${HENN_FT_DEADLINE_SECONDS:-}" ]; then
  python3 "$repo_root/Scripts/run_with_deadline.py" \
    "$HENN_FT_DEADLINE_SECONDS" "$scratch/runner.log" -- \
    python3 "$runner" "$manifest"
else
  python3 "$runner" "$manifest"
fi
"$wolfram_script" -script "$example" --compare-manifest "$manifest"
