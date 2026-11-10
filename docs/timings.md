# Timings

## Current code

Wall times of the current code on a 4-core x86-64 Linux container (GCC 13,
FLINT 3.4.0, Release build), single observations. Feynman-trick runs use the
built-in IBP solver, a cold cache and `--no-numerical-cache`; "before" is the
same machine at the start of the transport and recursion rework (commit
`d68d2dd`).

| Feynman-trick example | Before | Now | Preparation / numerical now |
| --- | ---: | ---: | --- |
| `sunrise` | 1.46 s | 0.29 s | 0.05 / 0.24 s |
| `box` | 0.52 s | 0.16 s | 0.05 / 0.11 s |
| `pentagon` | 3.71 s | 1.24 s | 0.16 / 1.07 s |
| `box_bubble` | 3.02 s | 0.73 s | 0.65 / 0.08 s |
| `pentagon_massive` | 25.1 s | 3.20 s | 0.58 / 2.61 s |
| `box_triangle` | 17.2 s | 3.20 s | 1.88 / 1.31 s |
| `banana` | 30.1 s | 4.24 s | 0.48 / 3.74 s |
| `banana_unequal` | 50.7 s | 6.25 s | 0.63 / 5.59 s |
| `kite` | failed | 9.7 s | 0.9 / 8.8 s |
| `double_box_planar` | not run | 29.4 s | 5.3 / 24.1 s |
| `banana4` (four loops, equal masses) | not run | 85 s | 30.5 / 53.9 s |
| `banana4_unequal` (four loops, four masses) | not run | 940 s | 33.5 / 906 s |

`kite`, `double_box_planar` and `banana4` were measured with the threaded
chart expansion and probe propagation of uncertainty (four threads); the rows
above them, and `banana4_unequal`, before.
The `banana4_unequal` run shared the machine with two other jobs. Its finite
part and that of `banana4` reproduce the stored acceptance references
(28.647167500233570163... and 39.655526834297652529...).

`henn_double_pentagon_x0` (two loops, eight propagators, experimental): the
exact preparation of all seven levels takes 299 s cold with the settings in
its family file (level 1: 69 masters, 153 s; before, it did not get past
level 1). The epsilon-expanded numerical recursion completed levels 7 to 3 in
24 minutes on four threads; levels 2 and 1 run for hours. By epsilon sampling
the family completes: one sample takes 227 s on four threads, and twelve
samples (four at a time) gave its Laurent coefficients through epsilon^0 in
48 minutes. [What limits the recursion](recursion-cost-experiments.md) has
the method, the measurements and the coefficients.

Nonzero coefficients of the nine examples that ran before agree between the
two versions to 1e-33 relative or better. The JSON report splits the
numerical time into `endpoint_series_seconds`, `endpoint_operator_seconds`,
`basis_selection_seconds`, `leaf_seconds` and `ordinary_seconds`; for `kite`
the ordinary transport is now about three quarters of it.

| Generic transport request | Before | Now | Check |
| --- | ---: | ---: | --- |
| `tools/performance/planar-one-loop-request.json.gz` (13 canonical integrals, epsilon 0-4, order 40, AccuracyGoal 15) | 24.3 s | 3.9 s | versions agree to 7e-26 |
| Original equal-mass banana system, t=-1 to 20 on a complex path, order 50, 384 bits, AccuracyGoal 30 | 13.4 s | 1.8 s | published reference to 6e-28 |
| Same, order 120, 768 bits, AccuracyGoal 60 | 82.5 s | 4.0 s | published reference to 6e-28 |

The remaining sections record earlier measurements on other hardware and
earlier code; they are kept for reference and are not comparable with the
table above.

## Earlier measurements

These are observed wall times on an Apple M4, macOS ARM64, with a Release
build and FLINT 3.4.0. Other compilation and numerical checks sometimes ran
concurrently. These are development measurements, not isolated hardware
benchmarks. Wrapper times below exclude Wolfram kernel startup.

| Calculation | Observed time | Scope |
| --- | ---: | --- |
| New one-loop configured family | 0.0105 s | Masses squared 2, momentum squared -3; independent elementary formula |
| New two-loop configured family, empty cache | 3.15 s | Three distinct masses, changed momentum; finite part |
| Same two-loop family, reused cache | 2.67 s | Exact preparation reused; numerical evaluation included |
| Original MPL examples | 0.053–1.65 s each | Includes the weight-20 example; reference differences below 1.4e-37 |
| Original planar one-loop wrapper | 21.03 s | PH1 to PH6, epsilon 0–4, AccuracyGoal 15 |
| Original planar `zmz` wrapper | 2m 45s | PH1 to PH6, epsilon 0–4, AccuracyGoal 15 |
| Original planar `mzz` wrapper | 2m 1s | PH1 to PH6, epsilon 0–4, AccuracyGoal 15 |
| Original planar `zzz` wrapper | 3m 13s | PH1 to PH6, epsilon 0–4, AccuracyGoal 15 |
| Original 128-digit `zzz` demonstration | 41m 27s | PH1 to PH2, epsilon 0–4, 145 charts; endpoint basis conversion 0.71 s, reference replay 3.46 s separately |
| Original Henn wrapper | 59.75 s | Transport of 108 supplied boundary components, AccuracyGoal 30 |
| Original banana partial initialization | 27.69 s + 2.16 s | Mathematica gamma/epsilon preparation, then native matching and transport through the wrapper |
| Original equal-banana transport to 32 | 6.37 s | From the derived boundary at -1, with saved metadata |
| Full equal/unequal banana route workflow | 9m 34s | Both routes and reference comparisons; includes several separate transports |
| Banana plot sample at 20.1 from nearby sample | 0.393 s | Versus 12.32 s from the original boundary; difference 4.35e-38 |
| Original banana `ReImPlot` check | 1.86 s | Bounded plotting range and sample count in the regression script |
| Equal four-loop banana FT | 20m 14s | Previously accepted full finite-part calculation and reference work |
| Unequal four-loop banana FT | 46m 24s | Previously accepted full finite-part calculation and reference work |

The five planar/Henn workflows took 9m 43s together, including five kernel
startups and reference validation. The unequal-banana route workflow is
particularly uneven: changing the masses first took 44.8 s and then transporting
the momentum to 50 took 435.3 s; the reverse order took 8.56 s and 49.3 s.
The two endpoints agree within 4.53e-38. The published unequal-banana reference
agrees within 4.80e-28 at its supplied precision.

The separate 128-digit `zzz` PH1-to-PH2 demonstration completed at Taylor order
230 and WorkingPrecision 280. Its maximum difference from the published values
is 5.70e-134. The computation used continuous root sheets; the published
integrand prefactors then converted the endpoint to the principal basis in
0.71 s using the same helper as native transport. The complete request was
verified unchanged except for those supplied prefactors, so no ODE rerun was
needed. The 3.46 s Wolfram replay includes startup and reference validation;
its short wrapper timing does not include the 41-minute native computation.
The [interface report](validation/interface-2.1.1.json) preserves the initial
basis mismatch and the earlier incomplete trials, including a 30-minute timeout.

The four-loop timings are preserved acceptance measurements. Those expensive
calculations were not repeated for this interface release. The full Henn FT
reconstruction is separate from the passing Henn boundary transport above;
see "Current code" for how far it runs now.

The native FT JSON reports `preparation_seconds`, `numerical_seconds` and
`total_seconds`. Mathematica's `DiffExpLastTimings[]` additionally reports total
wrapper time and individual native calls. Boundary preparation such as a
Mathematica gamma-function epsilon expansion is measured separately.

The [performance experiments](performance-experiments.md) measured a guarded
rational-chart speedup of 1.94–2.48 times on their tested workload. Cases that
need the scalar fallback can be slower. Sharing identical canonical source
convolutions reduced the observed `zmz`, `mzz`, `zzz` and Henn wrapper times by
roughly 2.2–2.9 times in development comparisons. The one-loop measurement was
slower in the loaded final run (21.03 s versus an earlier 11.90 s), so these
figures should not be read as universal end-to-end speedup guarantees.

See [validation](validation.md) for accuracy, compatibility and omitted-tail
limitations, [interface results](validation/interface-2.1.1.json) for detailed
workflow timings and errors, and [generic-family results](validation/generic-families.json)
for the new-family configurations and separate cold/reused-cache times.

## Comparison with AMFlow 2.0

The [small-integral comparison](amflow-comparison.md) measures bubble, sunrise
and massless-box FT evaluations at a common independently checked 20-digit
threshold. It reports cold and reused-preparation timings, the isolated LiteRed
compatibility adaptation needed for sunrise, and an optimization experiment
that was not adopted after it slowed the cached sunrise matrix.

## Further cold FT laptop examples

The [four-example run](ft-laptop-examples.md) measured 7.85 s for a massive
pentagon through epsilon squared, 18.88 s for the box-triangle, 23.86 s for an
unequal three-loop banana, and 7m 45s for the planar double-box. All passed
independent reference checks. These runs use fresh exact caches and disable
persistent numerical checkpoints; the report records controls, tolerances and
separate preparation/reference times.
