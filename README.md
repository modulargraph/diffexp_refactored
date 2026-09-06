# DiffExp 2.1

DiffExp 2.1 computes dimensionally regulated Feynman integrals with series
recurrences and recursive Feynman-parameter integration. The mathematical
backend is a standalone C++20 package built on FLINT/Arb.

**Mathematica users can use the familiar DiffExp1 workflow:** load a
configuration, prepare boundary conditions, transport the integrals and plot
the results. The Mathematica interface calls the C++ backend through
`RunProcess`; you can work with Wolfram expressions and the original function
names directly.

[Mathematica quick start](#mathematica-diffexp1-style-workflow) ·
[Build and install](#build) · [Feynman-trick recursion](#feynman-trick-recursion) ·
[Examples and timings](#examples-and-validation) ·
[AMFlow 2.0 comparison](docs/amflow-comparison.md) ·
[Built-in IBP Solver](docs/ibp-solver.md)

## Mathematica: DiffExp1-style workflow

[Build the executable once](#build), then load `DiffExp.m` in Mathematica.
The wrapper finds `build/diffexp` automatically. For an installed executable,
set `$DiffExpExecutable = "/your/prefix/bin/diffexp"`.

This self-contained example solves `dy/dt = y/(1-t)` with `y(0) = 1`.
For your own system, point `MatrixDirectory` at your existing matrix files
and supply your boundary values:

```wl
Get["/path/to/DiffExp2.1/DiffExp.m"];

matrixDirectory = CreateDirectory[];
Put[{{1/(1-t)}}, FileNameJoin[{matrixDirectory, "dt_0.m"}]];

LoadConfiguration[{
  MatrixDirectory -> matrixDirectory,
  EpsilonOrder -> 0, ExpansionOrder -> 40,
  WorkingPrecision -> 80, AccuracyGoal -> 20
}];

boundary = PrepareBoundaryConditions[{1}, {t -> 0}];
result = TransportTo[boundary, {t -> 1/2}];
timings = DiffExpLastTimings[];
result[[2]]  (* {{2}}: integral values, ordered by component and epsilon power *)
```

`TransportTo` returns `{point, values, estimatedErrors}`. Its result can be
used as the boundary for another transport. Save a transport to obtain
numerical functions for evaluation and plotting:

```wl
saved = TransportTo[boundary, {t -> x}, 1/2, True];
functions = ToPiecewise[saved];
value = functions[[1, 1]][1/4];  (* 4/3 *)
Plot[functions[[1, 1]][x], {x, 0, 1/2}]
```

The familiar `UpdateConfiguration`, `CurrentConfiguration` and
`IntegrateSystem` calls are also available. Supported inputs include ordinary
`dVARIABLE_ORDER.m` and canonical `d_1.m` matrices, `SparseArray`, epsilon
expansions, partial asymptotic boundaries and causal prescriptions.

**All nine original numerical example groups have been checked through this
Mathematica interface.** Some API differences remain: symbolic general
solutions are not provided, and algebraic bases defined with principal
endpoint roots can require explicit `BasisPrefactors`. See the
[Mathematica guide](docs/mathematica.md) for the complete compatibility scope,
the [runnable example](examples/Mathematica.wl), and
[measured timings](docs/timings.md).

**Adaptive spectral transport is available through the same Mathematica interface.**
`Recurrence -> "auto"` selects it for suitable larger epsilon-linear systems,
with local series as a fallback. Use `UpdateConfiguration[Recurrence -> "spectral"]`
to request it explicitly with a positive `AccuracyGoal`; `"taylor"` selects
local finite-lag/Taylor transport. See [adaptive transport](docs/adaptive-spectral.md)
for accuracy checks, timings and applicability. In the complete RL1/RL2
40-digit examples, Mathematica `TransportTo` now takes 0.65 s / 1.20 s on the
test laptop, versus 2.88 s / 3.62 s previously (single observations).

**Feynman-trick recursion now also has automatic spectral acceleration.**
`DiffExpFeynmanTrick` uses the same generic family configuration and selects the
new backend from the equation structure. You can compare it with
`"--ft-transport", "taylor"` in the wrapper argument list. See
[FT algorithms, Mathematica usage and matched timings](docs/ft-spectral.md).
For the massive-pentagon FT example, matched runs show **2.9–3.6× faster
ordinary transport and about 1.5× faster complete calculations**. Other tested
families retain local recurrence where the spectral approach does not help.

**The built-in IBP Solver also works through Mathematica.**
Use `DiffExpFeynmanTrick[family, {"--ibp-provider", "ibp-solver",
"--cache", "/path/to/cache"}]` with a JSON filename or family Association.
It prepares coefficients depending on both the merge parameter and dimension,
without launching FIRE. See [provider selection and validation](docs/ibp-solver.md).

The installed wrapper is `share/diffexp/Mathematica/DiffExp.wl`. Communication
uses ordinary processes and JSON, with no LibraryLink or WSTP setup.
The [direct native operations](docs/mathematica.md#direct-native-operations),
including `DiffExpSeries` and `DiffExpFeynmanTrick`, are available from the
same package.

## Build

Requires a C++20 compiler, CMake 3.24+, FLINT 3.4+, Boost.JSON 1.80+ and pkg-config.
The tested platform is macOS; the recursive FIRE process provider uses POSIX
facilities. Mathematica is optional. A pinned core of the separate
[IBP Solver](https://github.com/modulargraph/ibp-solver) package is included for
cold finite-field reductions and reconstruction. FIRE7 and FIRE7p remain
selectable providers for larger or unsupported reductions.

On macOS, dependencies can be installed with Homebrew:

```sh
brew install cmake pkgconf flint boost
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
cmake --build build -j2
build/diffexp --version
build/diffexp series examples/logarithm.json
build/diffexp bubble
```

`bubble` uses the specialized certified recurrence. General recursion uses `ft`.
Large Feynman-trick examples are opt-in and are never started by the default
test suite. To build and run tests, use the commands in
[validation](docs/validation.md#reproduce-the-release-checks).

Install the executable, C++ headers and Mathematica wrapper with:

```sh
cmake --install build --prefix /your/prefix
```

CMake consumers use `find_package(DiffExp 2.1 CONFIG REQUIRED)` and link
`DiffExp::core`. Public headers are under `diffexp/`, with namespace `diffexp`.
See [the consumer example](tests/consumer) and [architecture](docs/architecture.md).

## Feynman-trick recursion

Editable [family configurations](docs/feynman-families.md) define propagators,
kinematics and requested powers for both CLI and Mathematica. The name is a
label; changing the geometry changes the integral without registering a family.

```sh
build/diffexp ft examples/feynman/sunrise.json --fire /path/to/FIRE7 --cache /path/to/cache --json
build/diffexp ft examples/feynman/banana4.json --fire /path/to/FIRE7 --cache /path/to/cache --json
build/diffexp ft examples/feynman/banana4_unequal.json --fire /path/to/FIRE7 --cache /path/to/cache --json
```

The same command from Mathematica:

```wl
family = DiffExpFamilyTemplate["sunrise"];
(* Edit the propagators, masses, scalar products and requested powers. *)
DiffExpFeynmanTrick[family, {
  "--fire", "/path/to/FIRE7", "--cache", "/path/to/cache",
  "--epsilon-order", "0"
}]
```

`--json` returns one JSON object on stdout, with progress on stderr. Coefficients
are ordered by component, then epsilon power, beginning at `epsilon_low`.
Each complex coefficient has `real` and `imaginary` Arb interval strings.
`omitted_tails_certified` reports whether omitted series tails are enclosed.
General recursive results currently return `false`: their arithmetic intervals
alone are not certified bounds for the full integral.

Use `prepare family.json` to prepare the exact recursion without numerical evaluation.
Use `--fire-prime /path/to/FIRE7p` for native modular reconstruction. Completed
exact reductions and endpoint columns are stored under `--cache` and verified
on reuse. Persistent numerical continuation also retains completed arms.

`--method adjoint|factored|auto|values` selects the transport method; the default
is `adjoint`. `--ft-transport auto|spectral|taylor` selects ordinary transport
inside the adjoint/factored calculation, with `auto` as the default. Numerical defaults are endpoint order 32, ordinary order 80,
384 working bits and 28 leaf digits. Preparation has finite budgets, adjustable
with `--fire-seconds`, `--level-seconds` and `--total-seconds`. These bound exact
preparation, not the duration of the full numerical calculation.
Run `build/diffexp --help` for the remaining options.

Compatible algebraic differential equations now automatically use a finite-lag
O(N) recurrence, with the existing series solver retained for unsuitable charts.
This also applies through Mathematica. The [production comparison](docs/finite-lag-default.md)
records a complete smaller-family test and high-precision chart timings.

## Examples and validation

Four [cold FT laptop examples](docs/ft-laptop-examples.md) pass independent
checks: massive pentagon (7.85 s), box-triangle (18.88 s), unequal three-loop
banana (23.86 s), and planar double-box (7m 45s). They use generic family JSON
and `--no-numerical-cache` to avoid large checkpoint files.

The numerical checks include both four-loop banana examples and all nine
original DiffExp example groups. The full 108-component Henn Feynman-trick
boundary reconstruction is **experimental and frozen**; it is not a completed
example. Transport from the published Henn boundary and reconstruction of its
first canonical component have passed separately.

See [timings](docs/timings.md) for measured runtimes and
[validation and limitations](docs/validation.md) for the tested families,
precision scope and Banana4 reference comparisons. The [example guide](docs/examples.md)
describes the original-data downloads and the opt-in FT acceptance runner.

The default build uses the new native core. The extracted prepared-input
compatibility runtime is available with `-DDIFFEXP_BUILD_KERNEL_RUNTIME=ON`
for inherited regression tests; the Mathematica wrapper does not use it.

DiffExp 2.1 replaces the previous application at the repository root. Earlier
versions remain available in Git history.

## License and citation

GPL-3.0-or-later. See [LICENSE](LICENSE), [CITATION.cff](CITATION.cff), and the
[method references](docs/citation.md).

Recent [literature benchmarks](docs/literature-benchmarks.md) cover the complete
Canko–Pozzoli RL1/RL2 systems through epsilon^6, Mathematica-wrapper timings, a
controlled epsilon-sampling study, and the original fixed-resolution
[CHESS-style C++ spectral comparison](tools/performance/literature/chess/README.md).
