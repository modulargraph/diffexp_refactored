# CHESS DP: matched baseline and standalone spectral experiment

The production adaptive Taylor transport did **not** beat CHESS on this case.
A separate C++ implementation of CHESS's ordinary sequential-epsilon
Chebyshev collocation method was about 3.7 times faster for the numerical stage.
Including input preparation, it was slower at the 20-digit target and 1.36 times
faster at the 40-digit target. It is an experiment, not the default solver.

## Provenance and matched problem

- CHESS method: Yuanche Liu and Yang Zhang. Public source and DP data:
  <https://github.com/Alice-Shimada/CHESS/tree/4a701fc1332f29f6237d14427336e60615b966e3>.
- The released physical-region DP path runs from parameter `t=0` to `t=1`.
  Both solvers use the same released `Atilde`, scalar letter derivatives,
  boundary at `pr1`, and independent AMFlow reference at `pr2`.
- All 267 components and seven coefficients are compared: 1,869 values.
  The released powers `ep^-4` through `ep^2` are shifted uniformly to indices
  zero through six, leaving the differential equation epsilon-linear.
- CHESS uses `SequentialEpsilon`, one kernel, no endpoint regularization,
  `Precision=60`, and `WorkingPrecisionA=80`. The prototype uses 200 binary
  bits for the nodal recurrence and 266 bits for matrix samples (rounding
  60/80 decimal digits upward). It computes the same Lobatto collocation
  equations, including the boundary row, at 16 or 32 nodes plus the left node.
- Runs are sequential with one numerical thread. Results below are medians
  of three alternating CHESS/prototype repetitions on the same machine.
  Every individual solve has a 600-second process limit.

The converter retains the exact scalar one-forms and rational sparse weights.
It rewrites rational radicands into eight independent polynomial square roots.
Their denominators are positive constants times polynomial squares on `[0,1]`;
exact root counts establish no radicand zero on this interval. Expressions
agree exactly at seven sample points, including both endpoints. The prototype
also continues roots between ordered nodes using the native branch tracker.

## Measurements, 2026-09-06

Times are seconds. “Preparation + numerical” excludes process startup,
JSON serialization, data download, compilation, and the one-off Mathematica
source conversion. The conversion took **1.532 seconds** and is recorded
separately; including it would remove the prototype's first-use total-time win.
CHESS setup includes package/data loading and its cached sparse decomposition.
The prototype preparation includes JSON parsing, exact expression compilation,
and the scalar collocation inverse; numerical time includes matrix sampling.

| Target | Nodes | Method | Numerical | Preparation + numerical | Max absolute error against AMFlow |
|---|---:|---|---:|---:|---:|
| 20 digits | 16 | CHESS | 1.551666 | 1.801432 | 3.22e-21 |
| 20 digits | 16 | C++ spectral experiment | 0.425032 | 2.080441 | 3.22e-21 |
| 40 digits | 32 | CHESS | 3.145165 | 3.395047 | 6.74e-42 |
| 40 digits | 32 | C++ spectral experiment | 0.854430 | 2.493518 | 6.74e-42 |

The maximum normalized error, `abs(error)/max(1,abs(reference))`, is the same
as the maximum absolute error in these runs. The summaries also preserve the
native accuracy-check normalization `abs(error)/(1+abs(reference))`.
The matching errors confirm the same retained collocation result; this is a
change of implementation, not an accuracy advantage or a new AMFlow run.

Earlier production trials are preserved in the summary. An ordinary-entry
N40/goal20 trial took 57.07 seconds (24.61 preparation, 32.46 numerical), with
3.86e-57 maximum error. Supplied-form N16/goal20 and N32/goal40 trials took
64.50 and 42.40 seconds respectively, with 7.55e-35 and 3.79e-54 errors.
Both form trials reached the finite-lag preparation time budget and used the
series fallback. These are exploratory single trials, not repetition medians;
the default ordinary-entry interface remains unchanged by the explicit-form addition.

## What the experiment does and does not establish

`spectral.cpp` independently constructs the standard barycentric Lobatto
differentiation matrix, replaces its left row with the boundary equation,
and inverts that small scalar matrix. At each epsilon order it applies the
sampled sparse connection to the preceding order, then applies the scalar
inverse to all components. It does not copy CHESS implementation code, solve a
large combined component/node system, or call Mathematica during the solve.

A polynomial check (`y'=epsilon*y`) recovered all seven endpoint coefficients
`1/k!` with maximum error 7.47e-60 at eight nodes. The released AMFlow endpoint
then validates every coefficient of the realistic DP case.

This fixed-node experiment has no adaptive spectral truncation estimator.
Reported radii enclose arithmetic for the finite collocation calculation only,
not the omitted spectral tail. It rejects asymptotic starts and endpoint basis
prefactors. Its fixed 200/266-bit settings are explicit; request precision and
accuracy-goal fields do not control it. None of these timings justify replacing
the production path/accuracy machinery with this experiment.

The main remaining first-use cost is exact input compilation (~1.64 seconds),
not the FLINT ball operations in the numerical stage (~0.43/0.85 seconds).
These measurements support reusing prepared scalar forms and optimizing the
expression/preparation layer before considering replacement of FLINT.

## Reproduce

From a configured DiffExp 2.1 checkout with FLINT and Boost JSON installed:

```sh
python3 tools/performance/literature/chess/run_benchmark.py /tmp/chess-reproduction
```

The runner downloads only the six pinned public source/data files, converts
them with Mathematica, builds the standalone C++ source, alternates three
repetitions of each method at both node counts, and validates all outputs.
It writes `summary.json` in the requested work directory. Add `--include-native`
to include the two production-form trials using `build/diffexp`; those are
kept separate from the collocation comparisons. `--kernel`, `--prefix`,
`--cxx`, `--repo`, and `--repeats` override the local defaults. The Homebrew
library prefix and macOS Wolfram kernel path are defaults, not requirements.

`results/2026-09-06-summary.json` records all measured repetitions, source
hashes, conversion checks, and prior production trials. Compressed first-repeat
CHESS values and complete prototype responses retain all 1,869 cells without
checking in duplicated 17-MB converted requests. Upstream source/data are fetched
at the pinned commit rather than redistributed here.

## Adaptive production follow-up

The fixed-resolution prototype above is now complemented by the native adaptive
transport backend. See [production results and reproduction](../../../../docs/adaptive-spectral.md).
Its automatic resolution checks add work and should not be conflated with the
fixed-node numerical timings above.
