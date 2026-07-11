# C++ Recurrence Backend

DiffExp2 has an optional compiled backend for the finite-width recurrence
inside `SolveChart` and `SolveParticular`. It is an accelerator for the
recurrence kernel, not a separate implementation of the whole package.
Chart construction, singularity analysis, analytic-continuation choices,
transport, integration, and result validation remain in Wolfram Language.

The backend is experimental and opt-in. The default is
`"RecurrenceBackend" -> "Wolfram"`.

## Build

The build requires:

- a C++20 compiler, CMake 3.24 or newer, and optionally Ninja;
- FLINT (with its GMP/MPFR dependencies) and Boost.JSON;
- a Wolfram installation containing the LibraryLink C headers.

On Apple Silicon with Homebrew and the standard Wolfram application path, a
release build is:

```sh
brew install cmake ninja flint boost

cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_OSX_ARCHITECTURES=arm64 \
  -DWOLFRAM_INSTALL_DIR=/Applications/Wolfram.app/Contents
cmake --build build -j4
ctest --test-dir build --output-on-failure
```

Set `WOLFRAM_INSTALL_DIR` to the installation's `Contents` directory (or its
platform equivalent). The expected header is
`SystemFiles/IncludeFiles/C/WolframLibrary.h` below that directory. Make sure
the compiled library and Wolfram kernel have the same architecture.

The normal output is:

```text
build/cpp/diffexp2_librarylink.dylib
```

with `.so` or `.dll` on other platforms. DiffExp2 looks there automatically.
To keep the build elsewhere, set `DE2_CPP_LIBRARY` to the full shared-library
path before starting `wolframscript` or the kernel.

Run the Wolfram-side parity suite after building:

```sh
DE2_REQUIRE_CPP=1 DE2_CPP_THREADS=4 \
  wolframscript -file Tests/test_cpp_backend.m
```

`DE2_REQUIRE_CPP=1` makes an unavailable library a test failure instead of a
skip; it is a test-runner switch, not a solver option.

## Enable the backend

For direct DiffExp2 use, select it in the normal configuration:

```mathematica
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];
DiffExp2`Config`LoadConfiguration[{
  "WorkingPrecision" -> 500,
  "ExpansionOrder" -> 50,
  "RecurrenceBackend" -> "Cpp"
}];
```

The Feynman-trick runner exposes the same setting through the environment:

```sh
DE2_RECURRENCE_BACKEND=Cpp \
DE2_CPP_THREADS=4 \
FT_EXAMPLES=banana \
FT_WORKING_PRECISION=500 \
FT_EXPANSION_ORDER=50 \
wolframscript -file Scripts/run_ft_stepwise2.m
```

For the chart micro-benchmark, use `BENCH_BACKEND=Cpp`. The selected backend
is included in benchmark rows, solve-cache keys, and Feynman-trick checkpoint
metadata. A checkpoint made with one backend is not silently resumed with the
other.

`DE2_CPP_THREADS` controls the batch of independent homogeneous basis columns.
It defaults to 4, is clamped to the number of tasks, and does not consume
additional Wolfram licenses. Use 1 or 2 when memory is tighter. Inhomogeneous
recurrences are not independent-column batches, and exact symbolic-regulator
batches currently run single-threaded.

## Architecture and correctness contract

The LibraryLink call is deliberately coarse grained: Wolfram prepares and
serializes one complete framed recurrence rather than crossing the library
boundary for every coefficient.

Wolfram remains responsible for all structural mathematical decisions:

- denominator clearing and the prepared finite-lag tensors;
- indicial data, Jordan blocks, resonance classification, and the explicit
  T/P/R solve schedule;
- epsilon work-window selection and adaptive lower-frame retries;
- chart/branch prescriptions and the final ODE residual certificate.

The C++ kernel applies the polynomial and grouped-rational lags, propagates
explicit coefficient-validity data, solves the serialized T/P/R schedule,
handles materialized inhomogeneous sources, and can assemble the spectral
`V` transform before returning compact sector slabs. Independent homogeneous
columns are batched inside one LibraryLink call.

This split is important for analytic regularization. Resonance and structural
zero tests are never inferred from floating-point or ball midpoints. Laurent
frames retain separate validity information, and a shift that would discard
known lower-frame content is an error. Numeric coefficients use complex Arb
balls (`acb`). Inexact Wolfram inputs cross as reliable midpoint digits plus
their tracked uncertainty, and returned midpoint accuracy is capped by a
certified Arb radius bound. Additional analytic regulators are supported exactly when all
coefficients are rational functions over `Q` in the declared regulator
symbols; FLINT rational-function arithmetic is used for that path. This
preserves regulator dependence through recurrence and assembly instead of
sampling it numerically.

Selecting `"Cpp"` is strict. An unsupported coefficient field, malformed
request, lower-frame underflow, load failure, or compiled-kernel failure is
reported as a DiffExp2 error. It never silently falls back to the Wolfram
recurrence. Choose `"Wolfram"` explicitly if that is the desired backend.

## Current scope and limitations

- This is a recurrence-kernel port, not yet a C++ port of segmentation,
  matching, transport, regularized integration, FIRE, or FeynmanTrick.
- Symbolic regulator coefficients must be exact rational functions over
  `Q`. Algebraic functions, transcendental dependence, symbolic complex
  coefficients, and mixed inexact symbolic expressions are rejected.
- Extra analytic-regulator symbols are currently collected from Wolfram's
  `Global` context. The canonical dimensional regulator epsilon is handled
  by the framed Laurent representation and is not one of those coefficient
  field variables.
- The independent symbolic-regulator ODE residual is a five-point rational
  diagnostic that skips coefficient poles; it is not an exact identity proof
  in the regulator field. The recurrence and returned coefficients themselves
  remain exact rational functions and are never produced by those samples.
- Symbolic rational-function batches are intentionally single-threaded until
  the FLINT multivariate allocator/context path is certified for concurrent
  use.
- Small systems and low expansion orders can be dominated by Wolfram
  preparation, validation, JSON serialization, and LibraryLink startup, so
  they need not show a large speedup.
- Each parallel task owns recurrence state. Reducing `DE2_CPP_THREADS` lowers
  peak memory for large systems or wide epsilon frames.

## Measured performance

These development measurements are from the committed `banana_L1` fixture
(`d = 7`) on the same Apple Silicon machine and are intended as regression
anchors, not portable timing promises:

| Case | Wolfram recurrence | C++ recurrence | Speedup |
|---|---:|---:|---:|
| WP 100, expansion 50, epsilon order 5 | 81.268 s | 4.380 s | 18.6x |
| WP 500, expansion 50, epsilon order 11 | about 246.7 s | 11.492 s | about 21.5x |

At expansion order 20, the measured pair was 3.766 s versus 1.527 s because
fixed bridge and validation work is a larger fraction. That parity run had
identical sector/window structure and a maximum numerical coefficient
difference of approximately `2.7e-118`.

Use `Scripts/bench_chart.m` for reproducible local comparisons rather than
extrapolating these figures to a different chart or machine.

## Troubleshooting

Backend availability and build information can be inspected from a loaded
kernel:

```mathematica
DiffExp2`CppBackend`BackendAvailableQ[]
DiffExp2`CppBackend`BackendInformation[]
```

Common problems:

- **Library not found:** build into `build/cpp`, or set
  `DE2_CPP_LIBRARY=/absolute/path/to/diffexp2_librarylink.dylib` before the
  kernel starts.
- **Library fails to load:** check CPU architecture and linked FLINT/Boost
  libraries. On macOS, `file` and `otool -L` are useful diagnostics.
- **CMake cannot find FLINT or Boost:** add the package-manager prefixes to
  `PKG_CONFIG_PATH` (FLINT) and `CMAKE_PREFIX_PATH` (Boost), then configure a
  fresh build directory.
- **Rebuilt library is not picked up:** evaluate
  `DiffExp2`CppBackend`ResetBackend[]` after all calls have finished, or
  restart the kernel, before loading the rebuilt binary.
- **E5 unsupported scalar:** keep regulator dependence exact and rational,
  remove undeclared/algebraic symbolic coefficients, or explicitly use the
  Wolfram backend. There is intentionally no automatic fallback.
- **Little speedup:** confirm a `Release` build, use a realistic expansion
  order, and set `DE2_CPP_THREADS` above 1 for a multi-column homogeneous
  solve. `DEBUG_CPP_RECURRENCE=1` prints request/bridge/decode and compiled
  kernel timings.
- **High memory use:** lower `DE2_CPP_THREADS`; this trades parallel speed for
  a smaller peak working set.
