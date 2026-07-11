#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cache_root=${DIFFEXP2_CACHE_DIR:-"$HOME/.cache/diffexp2"}
threads=${DE2_CPP_THREADS:-4}

mkdir -p "$cache_root/fire" "$cache_root/ladder/banana4-unequal"
cd "$repo_root"

# Experimental/high-cost example.  Squared masses
# {2,3/2,4/3,5/4,1}, p^2=-1, D=2-2 eps.  This prototype does not claim a
# completed source-controlled four-loop ladder result.
exec env \
  DE2_RECURRENCE_BACKEND=Cpp \
  DE2_CPP_THREADS="$threads" \
  FT_EXAMPLES=banana4_unequal \
  FT_WORKING_PRECISION=500 \
  FT_EXPANSION_ORDER=70 \
  FT_EPS_ORDER=0 \
  FT_BOUNDARY_EXTRA_ORDER=10 \
  FT_LEVEL_EPS_HALOS=0,4,7,7 \
  FT_DIVISION_ORDER=3 \
  FT_RADIUS_OF_CONVERGENCE=1 \
  FT_PREP_CACHE_DIR="$cache_root/fire" \
  FT_LADDER_CHECKPOINT_DIR="$cache_root/ladder/banana4-unequal" \
  wolframscript -file Scripts/run_ft_stepwise2.m
