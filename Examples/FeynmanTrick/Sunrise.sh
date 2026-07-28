#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)

exec wolframscript -script \
  "$repo_root/Examples/FeynmanTrick/RunExample.m" "$@" sunrise
