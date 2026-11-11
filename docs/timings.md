# Timings

Wall times on a 4-core x86-64 Linux container (GCC 13, FLINT 3.4.0, Release
build), single observations.

## Feynman-trick recursion

Complete runs from a family file: exact preparation with the built-in IBP
solver on a cold cache, then the numerical recursion, on four threads with
default options and `--no-numerical-cache`.

| Example | Total | Preparation / numerical |
| --- | ---: | --- |
| `sunrise` | 0.16 s | 0.05 / 0.10 s |
| `box` | 0.11 s | 0.06 / 0.05 s |
| `pentagon` | 0.70 s | 0.17 / 0.52 s |
| `box_bubble` | 0.49 s | 0.46 / 0.03 s |
| `pentagon_massive` | 1.43 s | 0.56 / 0.84 s |
| `box_triangle` | 2.12 s | 1.43 / 0.67 s |
| `banana` (three loops) | 1.76 s | 0.46 / 1.27 s |
| `banana_unequal` | 2.54 s | 0.54 / 1.97 s |
| `kite` | 4.6 s | 1.0 / 3.5 s |
| `double_box_planar` | 14.2 s | 3.7 / 10.3 s |
| `banana4` (four loops, equal masses) | 47.6 s | 19.4 / 27.9 s |
| `banana4_unequal` (four loops, four masses) | 87.1 s | 27.4 / 59.5 s |

The finite parts of the four-loop bananas reproduce the stored references
(39.655526834297652529... and 28.647167500233570163...).

Numerical phase on one thread: kite 4.1 s, planar double box 12.3 s,
four-loop banana 44.8 s; one epsilon sample of the four-loop banana 7.6 s.

`henn_double_pentagon_x0` (two loops, eight propagators, experimental):
preparation 158 s on four threads (298 s on one); one epsilon sample about
150 s on four threads; twelve samples give its Laurent coefficients through
epsilon^0. See [performance.md](performance.md).

`--epsilon-strategy sample` (preparation excluded, four cores): `box` with 8
samples 0.9 s, giving epsilon^-2..0 to 1e-12, 3e-12 and 5e-7; `banana4` with
14 samples 49.5 s (133 s of CPU), finite part to 4e-16. The double pentagon
needs twelve samples of about 150 s each on four threads, or of one thread
each side by side on a larger machine.

`higgs_jet_family_f` (non-planar, seven propagators, four massive; Euclidean
point, experimental): preparation 88 s, one epsilon sample 107 s on four
threads.

Preparation of the built-in IBP solver, four threads against one: planar
double box 3.4 s / 5.2 s, four-loop banana 15.3 s / 33.7 s.

## Generic transport

| Request | Time | Check |
| --- | ---: | --- |
| `tools/performance/planar-one-loop-request.json.gz` (13 canonical integrals, epsilon 0-4, order 40, AccuracyGoal 15) | 3.9 s | |
| Equal-mass banana system, t=-1 to 20 on a complex path, order 50, 384 bits, AccuracyGoal 30 | 1.8 s | published reference to 6e-28 |
| Same, order 120, 768 bits, AccuracyGoal 60 | 4.0 s | published reference to 6e-28 |

Comparisons with other programs are in
[amflow-comparison.md](amflow-comparison.md) and
[literature-benchmarks.md](literature-benchmarks.md); they were measured on
other hardware and are not comparable with the tables above.
