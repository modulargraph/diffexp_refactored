# DiffExp 2.2

DiffExp computes dimensionally regulated Feynman integrals numerically, to
tens of digits, in two ways:

- **Feynman-trick recursion.** From a list of propagators and kinematics it
  derives the integrals itself: IBP reduction with a built-in solver,
  differential equations in Feynman parameters, boundary values and
  transport. No boundary conditions or external programs are needed.
- **Differential-equation transport.** Given a system of differential
  equations and boundary values it transports the integrals along a path,
  with the DiffExp 1 workflow in Mathematica.

The backend is a standalone C++20 package on FLINT/Arb with a command-line
tool and a Mathematica wrapper.

## How fast

Complete Feynman-trick calculations from a propagator list (cold cache,
IBP reduction included, 28 requested digits) on a 4-core Linux machine:

| Integral family | Time |
| --- | ---: |
| One-loop box, pentagon | 0.1 s, 0.7 s |
| Two-loop sunrise, box-bubble | 0.2 s, 0.5 s |
| One-loop pentagon with masses | 1.4 s |
| Two-loop box-triangle | 2.1 s |
| Three-loop banana, equal / unequal masses | 1.8 s / 2.5 s |
| Two-loop kite | 4.6 s |
| Two-loop planar double box | 14 s |
| Four-loop banana, equal / four masses | 48 s / 87 s |

Transport of a supplied system: 13 one-loop canonical integrals through
epsilon^4 in 3.9 s; the equal-mass banana system at 60 digits in 4.0 s.
Details and conditions are in [timings](docs/timings.md); comparisons with
[AMFlow 2.0](docs/amflow-comparison.md) and with
[published benchmarks](docs/literature-benchmarks.md) are separate.

## Build

Requires a C++20 compiler, CMake 3.24+, FLINT 3.4+, Boost.JSON 1.80+ and
pkg-config. [Singular](https://www.singular.uni-kl.de) is optional: with it on the
PATH the recursion chooses its merges from the graph. OpenMP is used when available (`-DDIFFEXP_OPENMP=OFF` disables
it); results do not depend on the thread count. Mathematica is optional.

```sh
brew install cmake pkgconf flint boost        # macOS; use your package manager elsewhere
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
cmake --build build -j2
build/diffexp --version
build/diffexp ft examples/feynman/sunrise.json --json
```

The executables are compiled for the CPU of the build machine
(`-DDIFFEXP_NATIVE=OFF` gives a portable build). Most of the run time is spent
inside GMP and FLINT, and distribution packages of GMP are generic builds:
`tools/build-tuned-deps.sh PREFIX` builds both for this CPU and prints the
CMake flags to use them. Homebrew's libraries are already built per architecture.

Install the executable, C++ headers and Mathematica wrapper with
`cmake --install build --prefix /your/prefix`. CMake consumers use
`find_package(DiffExp 2.2 CONFIG REQUIRED)` and link `DiffExp::core`; see
[the consumer example](tests/consumer) and [architecture](docs/architecture.md).
To build and run the tests, see
[validation](docs/validation.md#reproduce-the-release-checks).

## Feynman-trick recursion

A [family file](docs/feynman-families.md) lists propagators, kinematics and
the requested powers. Edit one of the [examples](docs/examples.md) to compute
a different integral; nothing has to be registered.

```sh
build/diffexp ft examples/feynman/banana4.json --cache /path/to/cache --json
```

From Mathematica:

```wl
family = DiffExpFamilyTemplate["sunrise"];
(* Edit the propagators, masses, scalar products and requested powers. *)
DiffExpFeynmanTrick[family, {"--cache", "/path/to/cache", "--epsilon-order", "0"}]
```

Kinematics may be Euclidean or physical. Propagators are `1/(m^2 - q^2 - i0)`
and the contour that this prescription implies is derived for every level
from its Symanzik polynomial; above a threshold the result is complex
([the Feynman i0](docs/feynman-families.md#the-feynman-i0): what is derived
and what is refused).

`--json` prints one JSON object on stdout, with progress on stderr.
Coefficients are ordered by component, then epsilon power, beginning at
`epsilon_low`; each has `real` and `imaginary` Arb interval strings.

- **Accuracy.** The reported radii are estimates of the numerical error, not
  certified bounds for the integral (`omitted_tails_certified` is `false`
  for general recursions). `--uncertainty map` propagates uncertainty
  rigorously through the transports at several times the cost.
- **IBP reduction.** The built-in [IBP solver](docs/ibp-solver.md) is the
  default and runs on several threads. `--fire PATH` (FIRE7) and
  `--fire-prime PATH` (FIRE7p) select external reducers for larger systems.
  `--fire-seconds`, `--level-seconds` and `--total-seconds` bound the
  preparation.
- **Cache.** `--cache` stores the exact reductions (verified on reuse), the
  chosen master bases and transport checkpoints, so an interrupted run
  resumes. `--no-numerical-cache` keeps the exact part only. `prepare`
  instead of `ft` stops after the exact preparation.
- **Defaults.** 384 working bits, 28 leaf digits, Taylor order 80.
  `build/diffexp --help` lists every option;
  [performance](docs/performance.md) explains the ones that matter for speed.
- **Large families.** `--epsilon-strategy sample` evaluates the integral at
  shifted dimensions in parallel and fits the Laurent coefficients
  ([epsilon sampling](docs/performance.md#epsilon-sampling)). It is how the
  two-loop double pentagon (`henn_double_pentagon_x0`, experimental) is
  computed; that result has the expected leading poles but no independent
  reference.

## Mathematica

[Build the executable once](#build), then load `DiffExp.m` in Mathematica.
The wrapper finds `build/diffexp` automatically. For an installed executable,
set `$DiffExpExecutable = "/your/prefix/bin/diffexp"`.

This self-contained example solves `dy/dt = y/(1-t)` with `y(0) = 1`.
For your own system, point `MatrixDirectory` at your existing matrix files
and supply your boundary values:

```wl
Get["/path/to/DiffExp2.2/DiffExp.m"];

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

`TransportTo` returns `{point, values, estimatedErrors}` and its result can
be the boundary of another transport. `UpdateConfiguration`,
`CurrentConfiguration` and `IntegrateSystem` are available, as are
`dVARIABLE_ORDER.m` and canonical `d_1.m` matrices, `SparseArray`, epsilon
expansions, partial asymptotic boundaries and causal prescriptions.

- All nine original DiffExp example groups run through this interface.
  Symbolic general solutions are not provided, and algebraic bases defined
  with principal endpoint roots can require explicit `BasisPrefactors`.
- `Recurrence -> "auto"` chooses between local series (a finite-lag O(N)
  recurrence where the system allows it) and
  [adaptive spectral transport](docs/adaptive-spectral.md) for suitable
  larger epsilon-linear systems; `"taylor"` and `"spectral"` force one.
- The wrapper talks to `build/diffexp` through `RunProcess` and JSON, with
  no LibraryLink or WSTP setup. Set `$DiffExpExecutable` for an installed
  executable.

See the [Mathematica guide](docs/mathematica.md) and the
[runnable example](examples/Mathematica.wl). The same requests can be sent
from the command line: [transport protocol](docs/transport-protocol.md).

## Documentation

- [Timings](docs/timings.md) and [performance](docs/performance.md)
- [Family files](docs/feynman-families.md) and [examples](docs/examples.md)
- [Mathematica interface](docs/mathematica.md)
- [IBP solver](docs/ibp-solver.md)
- [Validation and limitations](docs/validation.md)
- [Architecture](docs/architecture.md),
  [transport backend](docs/backend-review-implementation.md),
  [epsilon scheduling and basis search](docs/epsilon-demand-design.md),
  [moving poles](docs/moving-poles.md)
- [Changelog](CHANGELOG.md)

## License and citation

GPL-3.0-or-later. See [LICENSE](LICENSE), [CITATION.cff](CITATION.cff) and the
[method references](docs/citation.md).
