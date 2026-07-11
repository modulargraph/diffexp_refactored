#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cache_root=${DIFFEXP2_CACHE_DIR:-"$HOME/.cache/diffexp2"}
threads=${DE2_CPP_THREADS:-4}

mkdir -p "$cache_root/fire" "$cache_root/ladder/bubble"
cd "$repo_root"

exec env \
  DE2_RECURRENCE_BACKEND=Cpp \
  DE2_CPP_THREADS="$threads" \
  DE2_VALUE_TRANSPORT=1 \
  FT_CPP_BATCH_ENDPOINT_ARMS=1 \
  FT_EXAMPLES=bubble \
  FT_WORKING_PRECISION=300 \
  FT_EXPANSION_ORDER=40 \
  FT_EPS_ORDER=0 \
  FT_BOUNDARY_EXTRA_ORDER=10 \
  FT_LEVEL_EPS_HALOS=0 \
  FT_DIVISION_ORDER=3 \
  FT_RADIUS_OF_CONVERGENCE=1 \
  FT_PREP_CACHE_DIR="$cache_root/fire" \
  FT_LADDER_CHECKPOINT_DIR="$cache_root/ladder/bubble" \
  wolframscript -file Scripts/run_ft_stepwise2.m
