# DiffExp 2.1

DiffExp 2.1 computes dimensionally regulated Feynman integrals with series
recurrences and recursive Feynman-parameter integration. The mathematical
backend is a standalone C++20 package built on FLINT/Arb. Mathematica users can
call the same executable through a small `RunProcess` wrapper.

This is the C++ rewrite developed under the working name **DiffExp3**. It replaces
the 2.0 application at the repository root. The former implementation remains
in Git history; its Mathematica API is not source-compatible with this version.

The numerical checks include both four-loop banana examples and all nine
original DiffExp example groups. The full 108-component Henn Feynman-trick
boundary reconstruction is **experimental and frozen**; it is not a completed
example. Transport from the published Henn boundary and reconstruction of its
first canonical component have passed separately.

## Build

Requires a C++20 compiler, CMake 3.24+, FLINT 3.4+, Boost.JSON 1.80+ and pkg-config.
The tested platform is macOS; the recursive FIRE process provider uses POSIX
facilities. Mathematica is optional. FIRE7 or FIRE7p is needed for larger cold
IBP reductions, but not for building or running the basic examples.

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

## Mathematica

Build once, then load the wrapper from the repository:

```wl
Get["/path/to/DiffExp2.1/DiffExp.m"];
DiffExpBackendInfo[]

result = DiffExpSeries[<|
  "matrix" -> {{"eps/(1-x)"}},
  "boundary" -> {{"1", "0", "0"}},
  "center" -> "0", "taylor_order" -> 8,
  "epsilon_low" -> 0, "epsilon_high" -> 2
|>];
result["coefficients"][[4, 1, 2]]  (* "1/3": x^3 eps^1 *)
```

The wrapper finds `build/diffexp` automatically. For another build or installed
executable, set `$DiffExpExecutable = "/your/prefix/bin/diffexp"`.
The installed wrapper is `share/diffexp/Mathematica/DiffExp.wl`.

Every call launches one process, passes ordinary command arguments and, when
needed, JSON on stdin. There are no LibraryLink libraries, persistent native
handles, or Mathematica-side reductions. Exact coefficients and arbitrary
precision intervals remain strings; returned text is never evaluated as code.
See [the wrapper guide](docs/mathematica.md) and [runnable example](examples/Mathematica.wl).

## Feynman-trick recursion

```sh
build/diffexp ft sunrise --fire /path/to/FIRE7 --cache /path/to/cache --json
build/diffexp ft banana4 --fire /path/to/FIRE7 --cache /path/to/cache --json
build/diffexp ft banana4_unequal --fire /path/to/FIRE7 --cache /path/to/cache --json
```

The same command from Mathematica:

```wl
DiffExpFeynmanTrick["sunrise", {
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

Use `prepare FAMILY` to prepare the exact recursion without numerical evaluation.
Use `--fire-prime /path/to/FIRE7p` for native modular reconstruction. Completed
exact reductions and endpoint columns are stored under `--cache` and verified
on reuse. Persistent numerical continuation also retains completed arms.

`--method adjoint|factored|auto|values` selects the transport method; the default
is `adjoint`. Numerical defaults are endpoint order 32, ordinary order 80,
384 working bits and 28 leaf digits. Preparation has finite budgets, adjustable
with `--fire-seconds`, `--level-seconds` and `--total-seconds`. These bound exact
preparation, not the duration of the full numerical calculation.
Run `build/diffexp --help` for the remaining options.

## Examples and validation

See [validation and limitations](docs/validation.md) for the tested families,
precision scope and Banana4 reference comparisons. The [example guide](docs/examples.md)
describes the original-data downloads and the opt-in FT acceptance runner.

The default build uses the new native core. The extracted prepared-input
compatibility runtime is available with `-DDIFFEXP_BUILD_KERNEL_RUNTIME=ON`
for inherited regression tests; the Mathematica wrapper does not use it.

Some serialized cache and diagnostic schema identifiers retain the old
`DiffExp3` working name to preserve saved artifacts. They are format versions,
not a dependency on another package.

## License and citation

GPL-3.0-or-later. See [LICENSE](LICENSE), [CITATION.cff](CITATION.cff), and the
[method references](docs/citation.md).
