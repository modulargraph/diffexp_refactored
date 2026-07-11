# Installation

DiffExp 2 consists of Wolfram Language modules and a C++ recurrence library
loaded through LibraryLink.  The release profile expects the compiled backend.
Direct transport does not need FIRE; the Feynman-trick workflow does.

## Requirements

- a Wolfram Language installation with `wolframscript`, a licensed kernel,
  and the LibraryLink C headers;
- CMake 3.24 or newer and a C++20 compiler;
- FLINT, GMP/MPFR, and Boost.JSON 1.80 or newer;
- Ninja is optional but used in the commands below;
- FIRE6 only when preparing new Feynman-trick systems.

The minimum supported Wolfram product version has not yet been pinned in the
build metadata.  That version pin is a release blocker; do not infer one from
a successful local build.

## Clone and build the C++ backend

On Apple Silicon with Homebrew and Mathematica installed at the standard
location:

```sh
brew install cmake ninja flint boost

cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_OSX_ARCHITECTURES=arm64 \
  -DWOLFRAM_INSTALL_DIR=/Applications/Wolfram.app/Contents
cmake --build build -j4
ctest --test-dir build --output-on-failure
```

`WOLFRAM_INSTALL_DIR` must contain
`SystemFiles/IncludeFiles/C/WolframLibrary.h`.  The kernel and compiled library
must use the same CPU architecture.

The expected output is one of:

```text
build/cpp/diffexp2_librarylink.dylib
build/cpp/diffexp2_librarylink.so
build/cpp/diffexp2_librarylink.dll
```

DiffExp 2 searches that location automatically.  For an out-of-tree build,
set `DE2_CPP_LIBRARY` to the absolute shared-library path before starting the
kernel.

## Verify the installation

From the repository root:

```sh
DE2_REQUIRE_CPP=1 DE2_CPP_THREADS=4 \
  wolframscript -file Tests/test_cpp_backend.m

wolframscript -file Examples/Direct/MinimalTransport.wl
```

Or inspect the backend in a kernel:

```mathematica
Get["/absolute/path/to/DiffExp2/DiffExp2/DiffExp2.m"];
DiffExp2`CppBackend`BackendAvailableQ[]
DiffExp2`CppBackend`BackendInformation[]
```

Selecting C++ is strict.  Library loading failures, unsupported coefficient
fields, and recurrence failures are returned as DiffExp 2 errors; they do not
trigger a hidden Wolfram recurrence.

## Wolfram-only development mode

The implementation still contains a reference recurrence:

```mathematica
DiffExp2`Config`LoadConfiguration[{
  "RecurrenceBackend" -> "Wolfram"
}];
```

This is useful for parity and debugging.  It is not the intended release
default.  In the prototype snapshot, however, the configuration schema still
defaults to `"Wolfram"`; always select `"Cpp"` explicitly until that mismatch
is fixed.

## FIRE for Feynman-trick preparation

The current package default expects FIRE6 below:

```text
Dependencies/fire/FIRE6
```

That directory is ignored by Git.  Alternatively set the package option after
loading `FeynmanTrick/FeynmanTrick.m`:

```mathematica
FeynmanTrick`SetFTOption["FIREPath", "/absolute/path/to/FIRE6"];
```

The command-line driver uses a watchdog for each FIRE invocation.  Override it
with `FT_FIRE_TIMEOUT_SECONDS`; a timeout or nonzero FIRE exit is a failed
preparation, not a partial result.

## Cache locations

The Feynman-trick driver has two independent persistent stores:

- `FT_PREP_CACHE_DIR` contains completed FIRE preparation and reduction data;
- `FT_LADDER_CHECKPOINT_DIR` contains DiffExp 2 endpoint transports and level
  boundaries.

Keep both outside the Git checkout.  A typical layout is:

```sh
export FT_PREP_CACHE_DIR="$HOME/.cache/diffexp2/fire"
export FT_LADDER_CHECKPOINT_DIR="$HOME/.cache/diffexp2/ladder"
```

Use `FT_REBUILD_PREP=1` only when the prepared topology must be regenerated.
See [Feynman Trick](FeynmanTrick.md) for checkpoint resumption.

## Troubleshooting

- If the backend is not found, confirm its path and set `DE2_CPP_LIBRARY`
  before the kernel starts.
- On macOS, use `file` and `otool -L` to check architecture and shared-library
  dependencies.
- If CMake cannot find FLINT or Boost, set `PKG_CONFIG_PATH` and
  `CMAKE_PREFIX_PATH`, then configure a fresh build directory.
- After rebuilding an already loaded library, restart the kernel or call
  ``DiffExp2`CppBackend`ResetBackend[]`` when no calls are active.
- Lower `DE2_CPP_THREADS` if a high-dimensional recurrence exhausts memory.

The backend architecture, supported coefficient fields, and performance
measurements are documented in [C++ Recurrence Backend](CppBackend.md).
