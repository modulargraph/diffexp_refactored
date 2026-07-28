#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)

# Fully massive equal-mass five-propagator kite, p^2=-1, D=2-2 eps.
exec wolframscript -script \
  "$repo_root/Examples/FeynmanTrick/RunExample.m" "$@" kite
