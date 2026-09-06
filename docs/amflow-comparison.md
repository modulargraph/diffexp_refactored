# Laptop comparison with AMFlow 2.0

DiffExp 2.1 is faster on the accepted small-integral comparisons below. These are
end-to-end measurements at a common **20 decimal digit validation threshold**,
not a universal ranking of differential-equation solvers. The bubble uses
DiffExp's specialized certified leaf; the other cases use generic FT recursion.

AMFlow 2.0 has a Feynman-trick recursion mode based on the Hidding–Usovitsch
approach and a C++ differential-equation backend. We used the official upstream,
`RecursionMode -> "FT"`, `DESolver -> "CPP"`, and FIRE7 with LiteRed.
See the [official README](https://gitlab.com/multiloop-pku/amflow/-/blob/326bbb420615220602ea697afa3355065f483f91/README.md)
and [AMFlow 2.0 paper](https://arxiv.org/abs/2607.08477).

## Timings

| Case | DiffExp cold | AMFlow FT cold | DiffExp warm | AMFlow FT warm |
| --- | ---: | ---: | ---: | ---: |
| Massive bubble | 0.0381 | 61.86 | 0.0157 | 24.34 |
| Equal-mass sunrise* | 1.245 | 176.11 | 0.791 | 48.41 |
| Massless box | 0.754 | 228.31 | 0.238 | 30.95 |

\* AMFlow sunrise uses the isolated LiteRed compatibility adaptation described
below. Stock AMFlow/LiteRed returned an incorrect zero. The cold speed ratios
are about 1,626×, 141× and 303× respectively; these are small-example end-to-end
ratios, strongly affected by orchestration and preparation costs.

Times are wall-clock seconds including program/kernel startup and shutdown,
initial reduction, equation preparation, numerical computation and result
collection. “Cold” means an unused package preparation-cache directory; OS and
hardware caches were not flushed. “Warm” starts a new process using the prepared
equations, and still recomputes numerical results. AMFlow also repeats its initial
reduction on these warm calls. There is one accepted cold and one warm observation
per case, not a statistical benchmark sample. Runs were sequential, with no other
benchmark or compilation running concurrently. Every integral invocation had a
600-second process-group timeout.

The initial setup/build attempts are excluded. A later 15.06 s warm sunrise attempt failed during licensed Mathematica child startup before numerical evaluation; it is preserved separately, and the accepted warm retry took 48.41 s. An uncorrected AMFlow/LiteRed
sunrise evaluation returned zero in 144.02 s and **fails accuracy validation**;
it is not a successful performance result. See the compatibility investigation
below. The bubble and box use the unchanged upstream algorithms and dependency.

## Inputs and accuracy

All integrals have unit physical propagator powers. The same family JSON drives
both programs. AMFlow's complete two-loop scalar-product basis additionally has
two zero-power irreducible numerators, which do not change the integral.

| Case | Kinematics | Dimension | Compared epsilon coefficients |
| --- | --- | --- | --- |
| Massive bubble | Both squared masses 1; p² = −1 | 2 − 2ε | 0, 1, 2 |
| Equal-mass sunrise | Three squared masses 1; p² = −1 | 2 − 2ε | 0; spurious poles also checked |
| Massless on-shell box | s = −1; t = −1/3 | 4 − 2ε | −2, −1, 0 |

Normalize AMFlow's Minkowski denominators to DiffExp's positive Euclidean scalar
convention by multiplying by `(-1)^(sum of propagator powers)`. Neither result
includes an `exp(L EulerGamma epsilon)` prefactor. The benchmark uses AMFlow's
standard `SolveIntegrals` interface; its order is relative to the conservative
leading power `−2L`, so the requested orders are 4, 4 and 2 respectively.

Independent references use Feynman-parameter quadrature for the bubble,
coordinate-space Bessel quadrature for sunrise, and the analytic massless-box
Laurent coefficients. The reference calculation takes 1.65 s, excluded from all
solver timings. Acceptance requires every compared coefficient to satisfy
`abs(computed-reference)/max(1,abs(reference)) <= 1e-20`; imaginary native midpoints
are checked too. All twelve accepted cold/warm observations pass. The largest scaled discrepancy is 8.17e-34 for AMFlow and 3.43e-46 for the printed DiffExp midpoints; these are measured discrepancies, not certified accuracy estimates.

The general DiffExp FT results do not certify omitted endpoint/transport tails.
Their printed ball radii cover retained arithmetic, not the full truncation error.
The independent comparisons validate these numerical observations; they do not
turn them into rigorous full-integral enclosures. The specialized bubble does
provide certified coefficient enclosures.

## What accounts for the gap here?

DiffExp computes the epsilon coefficients directly with series recurrences and
prepares its reductions in C++. AMFlow's standard interface samples epsilon and
reconstructs its Laurent series. Here it uses 12/14/7 epsilon samples for
bubble/sunrise/box, at 126/180/130 working decimal digits and Taylor orders
252/360/260. DiffExp uses the published family controls: 384 bits, endpoint order
32, ordinary order 80 and leaf goal 28 digits. These are different amounts of
excess work, compared at the same independently checked output threshold.

AMFlow's repeated Mathematica launches and its candidate-pair reductions are
substantial on small examples. The first box merge screens six pairs in 67 s;
the chosen pair is also DiffExp's configured default. Importing that exhaustive
search would add cost here. A future bounded search needs to recover its own
preparation cost on asymmetric families before becoming a default.

The numerical phases also differ, although their instrumentation is not identical:
DiffExp reports approximately 0.017/1.020/0.339 s on the cold calls, while AMFlow's
logs round its FT numerical phases up to 8/33/15 s. AMFlow's phase includes
numerical child-process work and epsilon sampling; it is not a pure C++ kernel
timing. The raw report retains both kinds of measurements.

This does not isolate performance of the two C++ arithmetic kernels at identical
orders, precision, master bases, requested intermediates or truncation policy.
The larger examples remain frozen; these results say nothing about their
relative AMFlow performance. AMFlow's default auxiliary-mass mode and alternative
IBP reducers were not benchmarked.

## LiteRed compatibility investigation

The bundled LiteRed 1.83 initially marks every sector with fewer denominators than
loops as scaleless, before its actual Feynman-parametric or IBP scaleless test:

```wl
st=(Plus@@#<nloops)&/@sectors;
```

This is invalid after FT merging. A single merged quadratic can involve two
independent loop momenta. For the recorded sunrise intermediate, the determinant
of its loop quadratic form is `x1 (4-3 x1)/4`; at `x1=1/2` it is `5/16`, and
completing the square leaves constant `−11/10`. It is a massive, full-rank Gaussian
integral, not a scaleless sector.

An isolated copy changes that preliminary filter to classify only the empty
sector as automatically zero. The unchanged actual scaleless test then retains
one nonzero merged sector and continues to reject the ordinary one-propagator,
two-loop control. The installed FIRE/LiteRed tree is untouched.

Exposing that sector also exposes optional internal-symmetry rules with free
continuous transformation parameters and square roots that FIRE's rational
coefficient reader cannot parse (`strange symbol in internal symmetries
coefficient`). The isolated adaptation skips these optional internal symmetries
for sectors with fewer denominators than loops. Full IBP identities and actual
scaleless checks remain. This is a dependency compatibility adaptation, not a
change to AMFlow's FT algorithm or numerical solver. The corrected cold and warm evaluations agree with the independent sunrise reference. Their timings are labeled separately from the stock result.

The small copy-and-patch tool refuses to change an installed source or overwrite
an existing destination and checks the expected source patterns. Its applicability
to other LiteRed versions is deliberately not assumed. No upstream issue or
message was sent automatically.

## Recurrence experiment

A small prototype keeps individual rational products `f=(p/q)y`, avoiding a large
common denominator for forcing terms. This is a limited probe of the rational
forcing idea used by AMFlow, not a port of its whole SCC solver.

| Retained chart | Production recurrence | Product prototype | Outcome |
| --- | ---: | ---: | --- |
| Manufactured, 10 components, N=128, 256 bits | 0.00964 s | 0.00705 s | 1.37× faster |
| Cached sunrise matrix, 3 components, N=80, 384 bits | 0.000556 s | 0.001864 s | 3.35× slower |

Warm times average five calls. Exact rational coefficients, carried input
uncertainty, baseline ball overlap and an independent dense recurrence were
checked. The cached-matrix case uses an exact coordinate translation and identical
cached test balls; it is not an original physical-chart replay. Additional edge
buffers and divisions outweigh the smaller convolution count in this case.
**The candidate was not adopted; production recurrence code is unchanged.**
The [source and bounded runner](../tools/performance/ft_product_experiment/README.md)
and [observations](validation/amflow-product-experiment.json) preserve the result.

## Reproduce

Machine: Apple M4, macOS 26.2 arm64; Wolfram 15.0.0; FLINT 3.4.0. Versions:

- DiffExp checkout `8ea296a` (2.1.1 code plus README update).
- AMFlow `326bbb420615220602ea697afa3355065f483f91`.
- FIRE7 `b038d5de256ff881c32f6e7345de39a2edcab836` (7.1).
- MPSolve `282b77c25667257b20a0d4d680e82da61ffb49a1` (3.2.3).

Build AMFlow's `desolver` and WSTP `link` following its official instructions,
and configure its FIRE installation path. This run used the upstream release
flags (`-O3 -ffast-math -march=native`), with OpenMP disabled. MPSolve was built
in a local prefix; the Apple yacc-generated C source needed
`-Wno-error=implicit-function-declaration` during compilation. No mathematical
MPSolve source was changed. A small executable wrapper passes `-j 1` to MPSolve;
AMFlow and FIRE also use one worker. Set `OMP_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1`, and `VECLIB_MAXIMUM_THREADS=1` (the runner sets these).

```sh
export AMFLOW_HOME=/path/to/amflow
export WOLFRAM_KERNEL=/path/to/WolframKernel
export FIRE7_BINARY=/path/to/FIRE7/bin/FIRE7
export AMFLOW_BENCH_ROOT=/tmp/diffexp-amflow-comparison

python3 tools/performance/amflow_compare.py diffexp box --group final
python3 tools/performance/amflow_compare.py diffexp box --group final --label warm
python3 tools/performance/amflow_compare.py amflow box --group final
python3 tools/performance/amflow_compare.py amflow box --group final --label warm
```

Use `bubble` and `sunrise` to select the other two bounded fixtures. The runner
refuses to overwrite recorded observations or call a reused preparation cache
“cold”. It never installs or builds dependencies. Do not run the commands
concurrently when measuring timings.

For the isolated sunrise dependency adaptation:

```sh
python3 tools/performance/apply_litered_ft_fix.py \
  /path/to/FIRE7/extra/LiteRed/Setup "$AMFLOW_BENCH_ROOT/litered-ft-setup"
export AMFLOW_LITERED_SETUP="$AMFLOW_BENCH_ROOT/litered-ft-setup"
python3 tools/performance/amflow_compare.py amflow sunrise --group corrected_v2
python3 tools/performance/amflow_compare.py amflow sunrise --group corrected_v2 --label warm
"$WOLFRAM_KERNEL" -noprompt -script tools/performance/amflow_litered_diagnostic.wl
unset AMFLOW_LITERED_SETUP
```

Copy `tools/performance/amflow_references.wl` into the benchmark root and run it
with the same Wolfram kernel. Then run:

```sh
python3 tools/performance/analyze_amflow.py "$AMFLOW_BENCH_ROOT" /tmp/comparison.json
```

The analyzer requires the cold/warm observations for all three cases, plus the
preserved stock sunrise failure and the sector diagnostic. Keep the generated
`ft-compatibility-patch.json` in `$AMFLOW_BENCH_ROOT/litered-ft-setup` (use that
location when creating the isolated copy). The checked report is
[validation/amflow-comparison.json](validation/amflow-comparison.json).
