# Bounded chart performance experiments

The [high-precision zzz algebraic recurrence experiment](algebraic-zzz-experiment.md) removes the roots with an exact diagonal basis change and measures a finite-lag O(N) recurrence: 1.80–1.86× faster at order 230 on the first chart. The [production follow-up](finite-lag-default.md) integrates it as a guarded default.

The grouped `acb_dot` experiment accelerates this retained Taylor-chart workload, but changes interval radii. A follow-up shadow implementation passes five targeted transport tests using the existing precision/radius guard and scalar fallback; fallback-heavy cases can still slow down. Neither raw Acb storage nor in-place Horner evaluation alone shows a compelling benefit. `acb_dot_precise` is substantially slower here. The guarded candidate was subsequently integrated into the production rational chart, with the scalar fallback and existing accuracy controls retained. The five targeted production tests passed again after integration.

## Scope and reproducibility

This synthetic workload models four observable rows over ten masters, plus one constant-source component: dimension 41, 124 sparse matrix entries, ε orders 0–5, Taylor order 80, and 384-bit complex ball arithmetic. Each block has three rational connections per row with different ε shifts and one rational source. The chart center is `1/8 + I/16` and the step is `1/32 + I/64`. Sparse and dense initial boundaries exercise different amounts of structural zero skipping. The rational adjoint chart is a reference/cross-check path; this experiment does not measure the default polynomial recurrence or an entire integral family.

Machine: Apple M4, arm64 macOS; Apple clang 21.0.0; FLINT 3.4.0. Flags: `-std=c++20 -O3 -DNDEBUG`. One compiler/run worker was used. Each benchmark process has a 60-second timeout and the runner has a finite overall budget. The frozen large reconstruction examples were not run.

From the repository root, compile and run:

```sh
/usr/bin/c++ -std=c++20 -O3 -DNDEBUG -Iinclude -I/opt/homebrew/include tools/performance/benchmark_charts.cpp -L/opt/homebrew/lib -lflint -lmpfr -lgmp -lboost_json -o /tmp/diffexp-performance-charts
python3 tools/performance/run_charts.py /tmp/diffexp-performance-charts
```

Homebrew include/library locations are specific to this machine. The runner records all observations in `docs/validation/performance-charts.json`. Three fresh processes are used for each method/input combination, rotating method order. Each process measures its first chart, then the mean of five further charts. Summary values are medians of those three observations. “First” means no prior chart invocation in that process, after exact input preparation; it does **not** mean flushed operating-system or hardware caches. Input preparation, compilation, and verification are excluded from chart timings. Warm timings still include coefficient expansion, allocation, recurrence, evaluation, and cleanup on every call.

## Measured results

Times below are seconds per chart. CPU speedups use warm medians; first-call and warm wall times are shown separately.

| Input | Method | First wall | Warm wall | Warm CPU | CPU speedup |
|---|---|---:|---:|---:|---:|
| Sparse | baseline | 0.0794 | 0.0773 | 0.0772 | 1.00× |
| Sparse | inplace_horner | 0.0859 | 0.0786 | 0.0780 | 0.99× |
| Sparse | raw_scalar | 0.0794 | 0.0775 | 0.0771 | 1.00× |
| Sparse | dot | 0.0407 | 0.0401 | 0.0400 | 1.93× |
| Sparse | precise_dot | 0.7248 | 0.7483 | 0.7416 | 0.10× |
| Dense | baseline | 0.3741 | 0.3759 | 0.3743 | 1.00× |
| Dense | inplace_horner | 0.3669 | 0.3670 | 0.3662 | 1.02× |
| Dense | raw_scalar | 0.3705 | 0.3669 | 0.3657 | 1.02× |
| Dense | dot | 0.1498 | 0.1484 | 0.1482 | 2.53× |
| Dense | precise_dot | 2.6581 | 2.4902 | 2.4673 | 0.15× |

All 30 processes passed their output checks. Ordinary dot retained minimum relative accuracies of 378 bits (sparse) and 382 bits (dense), matching baseline minima. Its largest radius ratios were 1.088 (sparse) and 6.204 (dense); only 231/246 and 172/246 balls, respectively, were contained in the baseline balls. Thus matching minimum accuracy does not mean identical conditioning behavior.

The later trial ran more slowly for several methods, so these medians are a local workload observation, not a tightly controlled hardware performance guarantee. In particular, percent-level layout/Horner differences are not persuasive. The complete observations preserve that variation.

## Mathematical checks

All methods use the same exact system, boundary balls, precision, Taylor order, ε window, center, and step. The recurrence is

\[
 y_{n+1,r,k}=\frac{1}{n+1}\sum_e\sum_{(r,c)\in P_e}\sum_{m=0}^{n} a_{e,m}\,y_{n-m,c,k-\epsilon_e}.
\]

The dot variants group only the innermost lag sum. FLINT operates on genuine contiguous Acb arrays with documented signed strides, without reinterpreting C++ wrapper storage. No coefficient, retained order, or nonzero summand is discarded. The in-place Horner variant keeps the original wrapper storage and operation order and changes only temporary allocation during polynomial evaluation.

Every run checks all 246 output balls for finiteness and overlap with the existing recurrence. Raw scalar and in-place Horner must additionally match all 246 output balls bit for bit. Radius ratios, containment counts, and minimum relative accuracy are recorded for each method. Grouped arithmetic computes the same retained polynomial with different ball rounding; overlap alone is not an independent omitted-tail certificate. No full-solution accuracy or end-to-end speed claim follows from these measurements.

## Guarded transport follow-up

`tools/performance/guarded-dot.patch` records the reviewed change affecting only the rational chart kernel. It retains the old scalar chart, adds the grouped-dot implementation, and dispatches as follows: evaluate the grouped result; call the existing `needs_rational_cross_check(input, output)` predicate; if it fires, return the scalar reference result. All existing polynomial/rational intersections, centered actions, and conditioning subdivisions remain in place. Precision, Taylor order, the half-precision reserve, the factor-256 growth threshold, and the accuracy floor are unchanged. There is no new numerical acceptance rule and no requirement that accepted ball radii match bitwise.

The initial experiment applied the patch to a generated shadow include tree. Five existing tests were compiled and run against both baseline and shadow headers: `test_conditioned_adjoint.cpp`, `test_centered_adjoint.cpp`, `test_polynomial_transport.cpp`, `test_adjoint_transport.cpp`, and `test_adjoint_scaling.cpp`. All ten runs passed, including the independent exact adversarial retained-polynomial check, broad-input reserve check, nonfinite rejection, complex and Laurent solutions, centered cancellation, and batching. The polynomial-only test is a regression check; it does not itself demonstrate the new rational kernel. These are targeted tests, not a full-suite or cached-family validation.

A direct shadow-chart probe checks that the guarded chart equals the dot result bitwise when the unchanged predicate accepts it, and equals the old scalar result bitwise when the predicate rejects it. Sparse and dense ordinary inputs are accepted, despite the previously recorded radius differences. A dense boundary with an added real radius of `2^-190` in each observable coefficient is rejected; the exact constant source stays unchanged. In all three repetitions, the fallback result matches all 246 scalar output balls bitwise. This broad-input case demonstrates that the guard declines the dot candidate; it does not claim the fallback itself restores already-lost input precision. The surrounding transport retains its existing final acceptance/subdivision logic.

Reproduce the shadow tests and guarded probe with:

```sh
python3 tools/performance/run_guarded_tests.py
/usr/bin/c++ -std=c++20 -O3 -DNDEBUG -Itools/performance/shadow -Iinclude -I/opt/homebrew/include tools/performance/benchmark_guarded_charts.cpp -L/opt/homebrew/lib -lflint -lmpfr -lgmp -lboost_json -o /tmp/diffexp-performance-guarded
python3 tools/performance/run_guard_probe.py
```

Each compile and test process has a 60-second timeout; compilation and execution are sequential. The runner reconstructs a scalar baseline and a guarded shadow from the production header, then regenerates the comparison patch. `docs/validation/performance-guarded-dot.json` records the test names, output, compile/run times, and all guarded chart observations. The probe uses the same machine/compiler/flags and 384-bit/order-80 controls as above. Its first call is process-first after input preparation; warm samples average three further calls and are summarized over three fresh processes.

## Guarded timings

All times are wall seconds for single targeted-test runs, including process startup; these are validation timings, not stable microbenchmark estimates.

| Test | Baseline | Guarded shadow |
|---|---:|---:|
| conditioned_adjoint | 0.338 | 0.292 |
| centered_adjoint | 0.362 | 0.609 |
| polynomial_transport | 0.279 | 0.267 |
| adjoint_transport | 0.269 | 0.277 |
| adjoint_scaling | 0.585 | 0.563 |

Guarded chart warm CPU medians (seconds per chart):

| Input | Baseline | Guarded dot | Speedup |
|---|---:|---:|---:|
| Sparse | 0.07771 | 0.04009 | 1.94× |
| Dense | 0.36943 | 0.14882 | 2.48× |
| Dense, broad input | 0.36370 | 0.50311 | 0.72× |

The broad-input slowdown is expected: the candidate is evaluated and rejected before the unchanged scalar chart runs. No lower precision or looser reserve was used to avoid that cost.

## Production changes

The guarded grouped-dot chart is now integrated with its scalar fallback and
outer acceptance controls. It is not an unconditional speedup: a fallback runs
both kernels, and the centered-adjoint test was slower in the shadow run.
The raw-layout and Horner-only candidates were not adopted.

The generic matrix-transport engine also shares canonical letter expansions.
For the original one-loop planar transport at AccuracyGoal 15, the development
comparison fell from 80.15 seconds to 11.90 seconds. The runs used different
amounts of excess numerical work, so this is a same-requested-goal comparison,
not a claim of identical achieved digits.

Chained noncanonical transport needs additional care with uncertainty. Applying
an input interval repeatedly inside a Taylor recurrence can greatly widen it,
and grouped ball arithmetic can then lower midpoint accuracy. For significant
input uncertainty, small noncanonical systems with feedback in the epsilon-zero
dependency graph compute the retained local homogeneous map alongside the midpoint. The
incoming uncertainty is applied once through that map. Last-term estimates and
saved coefficients still include the uncertain contribution; precision and
acceptance thresholds are unchanged. A finite workspace cap retains the
original enclosing recurrence for maps that would exceed the budget.

The diagnostic equal-banana step from 2 to 32, with its original carried input
errors, took 8.81 seconds at order 70 and 499 bits. Its radius decreased from
398 to 4.11e-8, and its midpoint agreed with a separate zero-input-error run
within 1.9e-39. This diagnostic isolates interval overestimation; the zero-error
run is not an independent physical reference. The full wrapper route and
published-reference comparisons are reported separately. Canonical transport
keeps its original fast path; a regression verifies its retained balls are
unchanged.


## Shared canonical source convolutions

For an epsilon-linear canonical system, several matrix rows often need the
same convolution of one letter with one source component. The native transport
engine computes that convolution once per Taylor and epsilon order and reuses
its ball value. Each row retains its original rational weight and accumulation
order. A regression with multiple rows and signed, nonunit weights compares
all retained values, errors and saved coefficients bit for bit against the
unshared recurrence.

The original `zzz` matrix has 4,241 entries but 1,426 distinct letter/source
pairs; Henn has 18,079 entries but 1,873 pairs. These counts describe avoided
repeated work, not overall runtime speedups. Observed wrapper times were:

| Original workflow | Before sharing source convolutions | Final run |
| --- | ---: | ---: |
| `zmz` | 402.29 s | 165.18 s |
| `mzz` | 262.75 s | 120.56 s |
| `zzz` | 550.86 s | 192.95 s |
| Henn supplied-boundary transport | 151.12 s | 59.75 s |
| Planar one-loop | 11.90 s | 21.03 s |

These runs used the same requested accuracy and orders, but machine load was
not isolated. The one-loop slowdown is retained in the report. See the
[interface results](validation/interface-2.1.1.json) for final preparation,
numerical and wrapper timings and the independent reference comparisons.

## Selecting uncertainty propagation and plotting reuse

Strictly triangular iterated-integral systems use the direct enclosing
recurrence: computing a dense homogeneous map there adds work without solving
the feedback problem seen in the banana system. This selection preserves
input errors and the final accuracy criterion. The original weight-20 MPL
workflow fell from 34.31 s with the overbroad map selection to 1.65 s with the
triangular selector, and agrees with its independent value within 1.85e-57.

For Mathematica plots that vary one physical coordinate affinely, saved
functions reuse the nearest computed physical sample. All components and
epsilon coefficients share that sample. The banana point at 20.1 took 0.393 s
from a nearby sample, compared with 12.32 s from the original boundary, and
agreed within 4.35e-38. The bounded original `ReImPlot` check took 1.86 s. General
curves retain replay from the original boundary; a failed nearby continuation
falls back once. These measurements do not imply that arbitrary dense plots
are inexpensive.


## Whole-layer polynomial experiment

A further bounded experiment computes canonical epsilon layers in sequence,
using `acb_poly_mullow` for the complete letter/source products and polynomial
integration for the Taylor coefficients. This is possible because a canonical
epsilon layer depends only on the preceding layer. It is an experimental
benchmark, not a production solver option.

| First-chart input | Production shared convolution | Polynomial layers |
| --- | ---: | ---: |
| Synthetic 3-component system, order 80, 384 bits | 0.00228 s | 0.00172 s |
| Original `zzz`, 86 components, order 230, 931 bits | 23.49 s | 33.26 s |

All 15 synthetic and 430 original output balls were finite and overlapped the
production retained-chart values. The maximum radius ratios were 1.022 and
1.137 respectively; neither calculation tests omitted tails. The candidate
was **not adopted**, because the relevant large chart was slower. These are
single observations under concurrent load, with exact preparation excluded
and candidate coefficient extraction included. The large negative relative
accuracy sentinel in the raw report comes from balls containing zero and is
not a meaningful finite precision estimate.

The benchmark source is `tools/performance/benchmark_canonical_polynomials.cpp`.
Compile it with the same include/link flags used above, then pass a canonical
transport request JSON file and an exact rational first-chart step. The original
measurement used the exported high-precision request and step `1/100000`.
Each of the recorded processes had a 120-second limit. The
[raw observations](validation/performance-canonical-polynomials.json) preserve
the result and the decision.
