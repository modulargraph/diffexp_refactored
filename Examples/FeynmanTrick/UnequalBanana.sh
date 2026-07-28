#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)

# Squared masses {2,3/2,4/3,1}, p^2=-1, D=2-2 eps.
exec wolframscript -script \
  "$repo_root/Examples/FeynmanTrick/RunExample.m" "$@" banana_unequal
