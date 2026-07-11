# C++ Recurrence Backend

DiffExp2 has a compiled backend for the finite-width recurrence
inside `SolveChart` and `SolveParticular`. It is an accelerator for the
recurrence kernel, not a separate implementation of the whole package.
Chart construction, singularity analysis, analytic-continuation choices,
transport, integration, and result validation remain in Wolfram Language.

The release default is `"RecurrenceBackend" -> "Cpp"`. The Wolfram
recurrence remains available as an explicit diagnostic/reference selection;
the selected backend never silently falls back to the other one.
The low-level `DiffExp2/DiffExp2.m` implementation loader retains the
Wolfram-reference default for focused module tests; users load the root
`DiffExp2.m` umbrella, which injects the release default.

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
DE2_REQUIRE_CPP=1 DE2_CPP_THREADS=2 \
  wolframscript -file Tests/test_cpp_arm_batch.m
```

`DE2_REQUIRE_CPP=1` makes an unavailable library a test failure instead of a
skip; it is a test-runner switch, not a solver option.

## Select the backend

The root release loader selects C++ by default:

```mathematica
Get[FileNameJoin[{repoRoot, "DiffExp2.m"}]];
DiffExp2`LoadConfiguration[
  "WorkingPrecision" -> 500,
  "ExpansionOrder" -> 50
];
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

### Lower/upper endpoint arms in one Wolfram kernel

`run_ft_stepwise2.m` can also batch the boundary-independent homogeneous bases
for corresponding lower/upper endpoint charts. Enable the experimental path
with `FT_CPP_BATCH_ENDPOINT_ARMS=1` together with
`DE2_RECURRENCE_BACKEND=Cpp`.
There is still only one Wolfram evaluator: planning, matching, transport state,
endpoint limits, integration, and checkpoint writes remain sequential. Before
the marches, Wolfram prepares one chart from each arm and sends their recurrence
requests through the same native task pool. The responses are replayed through
the ordinary `SolveHomogeneous` assembly and ODE certificates and then retained
in its bounded memo cache. This fills otherwise idle workers for scalar or
small systems and reduces LibraryLink round trips; a chart that already has at
least `DE2_CPP_THREADS` independent columns generally cannot gain much from
cross-arm batching.

The runner submits pair-sized waves, so `DE2_CPP_THREADS` remains the total
native worker budget rather than multiplying per arm. It only enables this
schedule when the level dimension is smaller than that budget, and skips a
single unpaired tail chart or the identical shared anchor: those cases expose
no idle native worker and the request-collection pass would be pure overhead.
No subkernel or second Wolfram license is opened. Native
prewarming is pure cache state and never records a completed arm: the lower
transport is still written atomically before the upper march, and a resumed
checkpoint computes only its missing arm. Boundary-dependent value transport
(`DE2_VALUE_TRANSPORT=1`) intentionally uses the sequential route.

Exact analytic regulators retain the same contract. If a collected request
uses the symbolic rational-function field, the native dispatcher keeps that
batch single-threaded; it does not specialize the regulator or fall back to a
numeric/Wolfram solve.

The scheduling heuristic is measurement-driven. With two native workers, a
scalar two-endpoint fixture at WP 500, expansion order 50, and epsilon order 5
took 0.526 s for sequential chart bases and 0.309 s with paired prewarming
(1.70x). Conversely, forcing the same mechanism on the `banana_L1` fixture
(`d = 7`) at WP 100, expansion order 20, epsilon order 5, and four workers took
44.281 s sequentially versus 59.915 s batched (0.74x): each individual chart
already filled the pool, while request capture/replay added work. Both runs had
identical sector structure and zero displayed coefficient difference at 114
and 519 decimal places, respectively. `Scripts/bench_cpp_arm_batch.m`
reproduces these comparisons.

### Unequal-mass three-loop banana comparison

At the Euclidean point

```text
m_i^2 = {2, 3/2, 4/3, 1},  p^2 = -1,  d = 2 - 2 eps,
```

the legacy 15-master route and the Feynman-trick route agree on the scalar
`I[1,1,1,1]`.  Legacy component 11 uses

```text
J11 = eps (1 + 3 eps) (1 + 4 eps) Exp[3 EulerGamma eps] I_FT.
```

The independent two-dimensional Bessel oracle is

```text
5.83402729266214946740741989567969814964058746213209...
```

Measured warm-solve results on the development Mac were:

| Route | Numerical settings | Timed scope | Seconds | Finite scalar | Oracle agreement |
| --- | --- | --- | ---: | --- | ---: |
| Legacy 15-master C++ | WP250, EO50, eps through 4, CP34, division 3 | exact slice audit, matrix reconstruction, planning, transport | 46.770 | `5.8340272926621494708226551820` | 18.2 digits |
| Feynman trick C++ | WP500, EO70, finite eps order, extra order 10, halos `0,4,7`, division 3 | analytic deepest boundary and full transport/integration ladder from prepared cache | 390.286 | `5.8340272926621494674057011919652946916` | 20.8 digits |

These are deliberately scoped warm timings. The legacy row excludes FIRE
matrix generation and generation of its frozen 38-digit equal-mass seed. The
FT row excludes FIRE/IBP preparation (`FTPREP CACHE HIT`). The legacy run also
delivers four positive epsilon orders while the FT timing targets the finite
coefficient, so the rows are not identical workloads; they establish route
agreement and the practical compiled-transport timing rather than a claim of
perfectly normalized end-to-end benchmarking. The direct benchmark is
reproduced by `Scripts/bench_unequal_banana_cpp.m`.

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
  global symbol context. The canonical dimensional regulator epsilon is handled
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
  ``DiffExp2`CppBackend`ResetBackend[]`` after all calls have finished, or
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
