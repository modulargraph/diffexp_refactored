#!/bin/sh
# Minimal example: compute one literature family with DiffExp 2.2.
#
#   sh compute_one.sh families/hvp_chpt_e5e6_d2.json
#
# The executable defaults to build/diffexp of this checkout (override: DIFFEXP=...).
#
# The family file is all the input DiffExp needs: propagators, kinematics,
# dimension and the requested integrals.  DiffExp derives the IBP reductions
# (built-in IBP Solver, no FIRE or Mathematica), the Feynman-parameter
# differential equations, the boundary values and the transport.  Exact
# reductions are cached, so a second run only repeats the numerics.
set -eu

family=${1:?usage: compute_one.sh FAMILY.json}
here=$(cd "$(dirname "$0")" && pwd)
diffexp=${DIFFEXP:-$here/../../build/diffexp}
threads=${THREADS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}
cache=${DIFFEXP_CACHE:-$HOME/.cache/diffexp-literature-2026}

OMP_NUM_THREADS=$threads exec "$diffexp" ft "$family" \
  --ibp-provider ibp-solver \
  --ibp-threads "$threads" \
  --cache "$cache" \
  --json
# The JSON response holds "coefficients": one row per requested integral,
# epsilon orders epsilon_low..epsilon_high, each an Arb ball
# "[midpoint +/- radius]"; "timings" splits preparation and numerics.
# Faster for finite integrals at ~12 digits:  --epsilon-strategy sample --samples 8 --leaf-digits 12
