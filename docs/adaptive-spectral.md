# Adaptive spectral transport

DiffExp 2.1 now integrates ordinary epsilon-linear differential equations with
an adaptive Chebyshev solver as well as local finite-lag/Taylor recurrences.
The implementation follows the sequential-epsilon collocation method of
[Yuanche Liu and Yang Zhang, CHESS](https://arxiv.org/abs/2606.26691).
It is an independent C++ implementation, not a new mathematical method.

## Measured improvement

Apple M4 laptop, FLINT 3.4, sequential processes. Native totals include exact
preparation, numerical transport and adaptive refinement. They exclude input
conversion, executable compilation and reference validation.

| System | Requested digits | Local transport | Adaptive auto | Speedup |
| --- | ---: | ---: | ---: | ---: |
| CHESS double pentagon, 267 masters × 7 epsilon coefficients | 20 | 51.49 s | 1.66 s | 31.1× |
| Same complete system | 40 | 49.04 s | 3.22 s | 15.2× |
| Canko–Pozzoli RL1, 27 masters × 7 coefficients | 40 | 0.338 s | 0.290 s | 1.16× |
| Canko–Pozzoli RL2, 25 masters × 7 coefficients | 40 | 4.668 s | 0.511 s | 9.13× |

Double-pentagon adaptive timings are medians of three fresh processes. The
matched local controls are single observations using the same request, 200-bit
precision and order-40 local expansion. Both controls produce about 56-digit
agreement even at the lower requested goal; adaptive transport avoids that
unnecessary work. Process-wall medians are 1.96 s and 3.52 s. RL timings are
single observations at 384 bits, order 80 and epsilon order 6; their local
baselines are the earlier [retained benchmark](literature-benchmarks.md).
The small RL1 improvement should not be overinterpreted.

All 1,869 double-pentagon coefficients pass comparison with the published
AMFlow reference distributed by CHESS. Maximum absolute differences are
1.37e-31 and 1.18e-56 at the two requested goals. RL1/RL2 comparisons with the
original DiffExp transport pass 40 digits, with maximum normalized differences
6.83e-50 and 2.47e-49. The solver chooses its resolution without these references.

The ordinary Mathematica matrix-file workflow also improves: `TransportTo`
takes **0.654 s for RL1** and **1.203 s for RL2**, compared with the previous
2.881 s and 3.619 s wrapper observations (about 4.4× and 3.0×). Both full
coefficient comparisons pass 40 digits. These are single observations; complete
fresh-kernel process times are 3.22 s and 3.51 s. A direct-kernel launch initially
stalled in Wolfram activation; the accepted runs use the usual `wolframscript`
launcher and include no concurrent builds or other benchmarks.

This supersedes the fixed-resolution prototype as a usable solver option. The
[earlier CHESS comparison](../tools/performance/literature/chess/README.md)
measures a different fixed-resolution workload. The adaptive production totals
are close to CHESS's retained fixed-resolution totals (1.80 s and 3.40 s),
while doing additional refinement checks. This is not evidence of a universal
speed advantage over CHESS or AMFlow, nor an end-to-end integral reconstruction
comparison. One-time input conversion remains excluded.

## Selection and interface

Native JSON accepts `"recurrence":"auto"` (default), `"spectral"`, `"taylor"`
or `"series"`. Mathematica exposes the same choice:

```wolfram
UpdateConfiguration[Recurrence -> "auto", AccuracyGoal -> 40];
result = TransportTo[boundary, target];
DiffExpLastTimings[]
```

Automatic selection considers systems with at least 16 components, a positive
accuracy goal, an ordinary boundary and no native Taylor-segment export. A
nearby-singularity estimate screens out unsuitable global intervals. Spectral
acceptance requires an epsilon-linear shared connection; ordinary derivative
matrices, dlog matrices and explicitly supplied scalar forms are supported.
There is no dispatch by family name. `"spectral"` also allows smaller systems
but reports an error if unsupported or unresolved; `"auto"` falls back.
`"taylor"` retains local finite-lag acceleration; `"series"` forces convolution.

The engine uses resolutions 8, 12, 16, 24, 32, 48, 64, 96 and 128, with at
least three levels checked. Nodes and shared scalar samples are cached across
levels; sparse contributions sharing a source and scalar form are grouped.
Repeated exact expressions are parsed once during preparation. The speedups
include that preparation improvement and preservation of shared scalar forms,
as well as the numerical solver change; they are not a pure inner-loop ratio.
FLINT remains the arithmetic and exact-algebra backend. Independent
polynomial square roots are continued along the path, including sheet changes,
and endpoint basis prefactors are respected.

The scalar collocation inverse is shared over masters and epsilon orders.
With m nodes, d masters, k epsilon orders and S sparse contributions, the
straightforward numerical work is O(m^3 + k m S + k d m^2), plus scalar sampling
and exact preparation. This is not the O(N) finite-lag recurrence: it wins here
by reaching the requested accuracy with few nodes and reusing work. It avoids
a dense solve whose dimension multiplies the master count by the node count.

## Accuracy and limits

Geometric refinement differences estimate the omitted spectral tail. For tiny
nondecaying components, two absolute differences with a safety factor may fit
the guarded tolerance instead; diagnostics report how many components used
this fallback. Incoming boundary boxes and arithmetic errors are also checked.
These are **estimated errors, not rigorous bounds on the infinite tail**.
Small differences alone can miss unresolved structure, so near-singularity
screening and reference/analytic regression tests remain important.
`omitted_tails_certified` is always false for this method.

Automatic attempts have an eight-second checked budget; explicit attempts have
thirty seconds. Node, arithmetic-work and storage budgets also apply. Time is
checked between operations, not by forcibly interrupting an individual FLINT
operation. Rejected automatic attempts preserve the ordinary solver fallback.
Native diagnostics include selection, resolutions, estimated geometry and
rejection reasons. On fallback requests the preparation interval also includes
rejected speculative work. The CLI reports preparation, numerical and total timings;
Mathematica additionally reports wrapper/process overhead.

This does not replace singular-start Frobenius matching, general epsilon-zero
coupled transport, nested-root support or the separate Feynman-trick recursion
backend. High-precision/local problems can still favor finite-lag recurrence.
The frozen heavy Henn and banana4 reconstructions were not restarted.

## Reproduction

Use the pinned CHESS data converter described in the
[benchmark guide](../tools/performance/literature/chess/README.md), then:

```sh
python3 tools/performance/literature/chess/run_production.py \
  --request /path/to/dp-shared-request.json --executable "$PWD/build/diffexp" \
  --output /tmp/diffexp-adaptive --repeats 3 --control
python3 tools/performance/literature/chess/validate_production.py \
  --reference /path/to/reference.wl --results /tmp/diffexp-adaptive
```

Validation requires Python `mpmath`; it reads the prepared numeric reference
without launching Mathematica. RL reproduction uses the existing
`tools/performance/literature/run_cp.py` workflow; native auto selection now
includes the adaptive method.

[Raw outputs, timings and validation](../tools/performance/literature/chess/results/adaptive-2026-09-06/summary.json)
retain every measured DP repetition and the RL endpoints. Forced local controls
precede the final small-component acceptance refinement; their algorithm is
unchanged. Executable hashes distinguish these runs.

The [regression record](validation/adaptive-spectral-tests.json) contains 80
passing tests, including seven original reference checks and the complete
Mathematica compatibility suite. Four targeted checks pass again after the
final kinematic-name fix. The new spectral tests include polynomial
cancellation, analytic epsilon coefficients, boundary uncertainty, complex
contours, root winding, narrow-feature rejection and solver fallback.
