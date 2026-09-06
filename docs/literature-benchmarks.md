# Literature benchmarks: Canko–Pozzoli, CHESS, and epsilon sampling

These are sequential laptop observations, not a universal solver ranking.
Every solver process is bounded to 600 seconds. Boundary generation and IBP
reconstruction are excluded from transport comparisons: the published boundary
values are inputs to both programs. Target accuracy is checked on every retained
component, using `abs(value-reference)/max(1,abs(reference))`.

## Canko–Pozzoli three-loop vector-boson integrals

Sources: [paper](https://arxiv.org/abs/2412.06972),
[ancillary archive](https://doi.org/10.5281/zenodo.14284044).
The complete RL1 and RL2 systems have 27 and 25 masters, respectively. All seven
coefficients epsilon^0 through epsilon^6 are transported from
`(x,y,z)=(3/2,1/5,1/2)` to `(5/3,1/11,1/7)`. Native protocol requests rename the
kinematic x to u; the Mathematica compatibility interface now performs its own
internal renaming when the user selects a different line parameter.

The authors' differential equations and supplied numerical boundary values are
used unchanged. Their runner is adapted only for local paths, output capture,
explicit working precision and time limits. Original DiffExp uses expansion
order 50, 115 decimal working digits, no Mobius transformation or Pade, and
40-digit AccuracyGoal. Native dlog requests use order 80 and 384 bits. The
Mathematica wrapper uses the same ordinary derivative matrices, order 50 and
115 working digits as the original. These are matched-output-accuracy tests,
not identical internal arithmetic or chart divisions.

Cross-method agreement is a validation check using the same supplied boundary;
it is not an independent derivation of the boundary constants. The native
omitted-tail estimates are not rigorous full-tail enclosures.


| Family | Original DiffExp transport | Native dlog request, preparation + transport | Mathematica wrapper transport | Max native/reference scaled difference |
|---|---:|---:|---:|---:|
| RL1 | 10.609s | 0.338s | 2.881s | 1.04e-54 |
| RL2 | 13.325s | 4.668s | 3.619s | 3.53e-49 |

All coefficients pass the 40-digit threshold. These are individual accepted
observations, not statistical medians. Native and wrapper columns use different
input representations and orders; the wrapper preserves the authors’ ordinary
derivative workflow. Mathematica kernel startup is excluded from the transport
columns; the accepted complete wrapper processes took 5.19s and 5.83s.

## Controlled epsilon experiment

This experiment isolates direct Laurent/Taylor coefficients versus Cauchy
sampling using the same C++ transport implementation and the actual RL1 matrix.
It does **not** run AMFlow, its optimized numerical kernel, its sampling-node
selection, or the full Feynman-trick recursion.

The boundary is truncated to its published coefficients through epsilon^2.
This defines an analytic boundary polynomial. Its evolved coefficients through
that order equal those from the full physical boundary, but its higher
coefficients are not the unknown full physical boundary. This distinction is
necessary to avoid pretending that finite published Laurent data determine an
all-orders physical boundary suitable for AMFlow.

Direct expansion and every sample use 384 bits, order 64, and 40-digit requested
accuracy. Sampling uses 12 and 16 complex nodes on a circle of radius 1e-4,
then a discrete Fourier/Cauchy reconstruction. Sampling-node coordinates are
rounded to exact rational numbers at 110 decimal places; input uncertainty
is included separately. Both sample counts are checked against the direct
coefficient solution at a common 40-digit threshold.

At this small epsilon depth, direct propagation wins. This is not a theorem:
sampling is naturally parallel, can use optimized fixed-epsilon solvers, and
may behave differently for deep Laurent windows, epsilon poles, large coupled
blocks, cancellation, or different boundary-generation costs. A full FT study
must also measure the intermediate epsilon depths required by recursive poles.


| Method | Numerical seconds | Sum of process wall times | Maximum scaled difference from direct series |
|---|---:|---:|---:|
| Direct epsilon coefficients, series | 0.380 | 0.394 | reference |
| Direct epsilon coefficients, auto | 0.448 | 0.507 | not separately analyzed |
| 12 Cauchy samples | 12.443 | 12.945 | 2.10e-47 |
| 16 Cauchy samples | 17.290 | 18.119 | 2.10e-47 |

The roughly 33x numerical ratio for 12 samples applies only to this controlled
transport experiment. Cauchy reconstruction and request-generation time are
excluded from its process-time sum. No conclusion is drawn about AMFlow’s own
optimized sampling solver from this ratio.

## CHESS comparison and spectral experiment

The production Taylor transport loses on the released physical-region
double-pentagon benchmark. A separate implementation of the ordinary
sequential-epsilon Chebyshev collocation method from
[Liu and Zhang’s CHESS paper](https://arxiv.org/abs/2606.26691) gives a numerical
speedup of about 3.7x. At 40 digits its median preparation-plus-numerical time is
2.494s versus CHESS’s 3.395s; at 20 digits it is slower overall (2.080s versus
1.801s). The one-time 1.532s input conversion removes the first-call total win.

These are implementation improvements on the published method, not a claim
of a new integration algorithm. The prototype is not a default or production
backend: spectral truncation is checked against external references in these
experiments, with no built-in adaptive tail estimator. See the
[complete reproducible comparison](../tools/performance/literature/chess/README.md)
for three alternating repetitions, all 1,869 checked coefficients, raw results,
source hashes, analytic smoke test and limitations.

## Reproduction

Download the authors' `ancillary.tar.xz`, verify MD5
`3e4c31002fcd2dd2e66e0aa123b46dcb`, and extract it. Download the original DiffExp
package from `https://gitlab.com/hiddingm/diffexp/-/raw/master/DiffExp.m`;
the recorded validation metadata pins the exact file hash used here.

```
python3 tools/performance/literature/run_cp.py \
  --data /path/to/ancillary --output /tmp/cp-benchmark \
  --kernel /path/to/WolframKernel --package /path/to/DiffExp2.1 \
  --original /path/to/original/DiffExp.m
DIFFEXP_HOME=/path/to/DiffExp2.1 python3 /tmp/cp-benchmark/epsilon_sampling.py
```

The Python analyzers require mpmath. The runner preserves logs and fails if a
solver fails or an observation would overwrite an existing output directory.
The initial archive download and installation/build times are outside the
solver timing columns. Numerical runs must be sequential to compare timings.

## Arithmetic backend implications

All native experiments use FLINT. The improvement from the spectral experiment
therefore comes from changing the algorithm and data organization, not from
replacing the arithmetic library. Exact preparation currently uses FLINT
multivariate rational functions (`fmpz_mpoly_q`) with explicit reductions for
algebraic roots; numerical work uses Arb/Acb. Specialized representations for
univariate rational functions and algebraic extensions, expression sharing,
and avoiding repeated normalization are plausible optimization targets. These
observations do not establish that FLINT is faster than every alternative CAS;
no matched CAS microbenchmark was performed here.

## Integrated fixes and checks

The published examples exposed two fixed input issues: Mathematica absolute
accuracy annotations on uncertain zeros, and kinematic `x` colliding with the
native path parameter. The wrapper now renames reserved symbols internally.
Explicit scalar forms preserve sparse epsilon-linear structure without forcing
large per-entry radical sums. Six transport, frontend, CLI and Mathematica
regressions pass in 25.02s. Benchmark outputs are cross-method validated; they
are not rigorous full-tail certificates.
