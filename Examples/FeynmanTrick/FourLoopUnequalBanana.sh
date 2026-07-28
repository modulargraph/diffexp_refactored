#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)

# Experimental/high-cost example.  Squared masses
# {2,3/2,4/3,5/4,1}, p^2=-1, D=2-2 eps.  This example does not claim a
# completed source-controlled four-loop ladder result.
exec wolframscript -script \
  "$repo_root/Examples/FeynmanTrick/RunExample.m" "$@" banana4_unequal
