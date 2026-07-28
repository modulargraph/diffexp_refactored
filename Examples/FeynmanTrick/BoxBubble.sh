#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)

# Massless two-loop three-point box-bubble subfamily, s=-1, t=-1/3.
exec wolframscript -script \
  "$repo_root/Examples/FeynmanTrick/RunExample.m" "$@" box_bubble
