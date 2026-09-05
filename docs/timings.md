# Timings

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

The separate original 128-digit `zzz` PH1-to-PH2 demonstration remains under
validation. Its 30-minute trial reached chart 100, at path parameter 0.787,
before the timeout. A longer bounded run is recorded separately in the
[interface report](validation/interface-2.1.1.json); no completed timing or
128-digit numerical acceptance is inferred from that partial run.

The four-loop timings are preserved acceptance measurements. Those expensive
calculations were not repeated for this interface release. Full Henn FT
reconstruction remains frozen and is separate from the passing Henn boundary
transport above.

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
