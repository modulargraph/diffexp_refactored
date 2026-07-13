# Installation

DiffExp 2 consists of Wolfram Language modules and a C++ recurrence library
loaded through LibraryLink. The release loader selects the compiled backend.
Direct transport does not need FIRE; preparing a new Feynman-trick family
does.

## Requirements

- Wolfram Language 15.0 with `wolframscript`, a licensed kernel, and the
  LibraryLink C headers;
- CMake 3.24 or newer and a C++20 compiler;
- FLINT, GMP/MPFR, and Boost.JSON 1.80 or newer;
- Ninja is optional but used in the commands below;
- FIRE 7.1 only for Feynman-trick preparation. The default modular backend
  also requires an MPI launcher and a FIRE build with prime, multiprime, and
  reconstruction executables.

The first release is validated against Wolfram Language 15.0. Earlier versions
may work but are outside the initial compatibility promise until added to the
release test matrix.

## Build the C++ backend

On Apple Silicon with Homebrew and Wolfram installed at the standard
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
`SystemFiles/IncludeFiles/C/WolframLibrary.h`. The kernel and compiled library
must use the same CPU architecture.

The expected library is one of:

```text
build/cpp/diffexp2_librarylink.dylib
build/cpp/diffexp2_librarylink.so
build/cpp/diffexp2_librarylink.dll
```

DiffExp 2 searches that build location automatically. For an out-of-tree
build, set `DE2_CPP_LIBRARY` to the absolute shared-library path before
starting the kernel.

## Verify the installation

From the repository root:

```sh
DE2_REQUIRE_CPP=1 DE2_CPP_THREADS=4 \
  wolframscript -file Tests/test_cpp_backend.m
```

Or inspect the backend through the stable umbrella:

```mathematica
repo = "/absolute/path/to/diffexp2";
Get[FileNameJoin[{repo, "DiffExp2.m"}]];

DiffExp2`BackendInformation[]
```

The root loader installs a validated configuration with
`"RecurrenceBackend" -> "Cpp"`. Library-loading failures, unsupported
coefficient fields, and recurrence failures are reported as DiffExp 2
errors; they never trigger a hidden Wolfram recurrence.

## Optional paclet discovery

`PacletInfo.wl` declares both root contexts. To register a built source tree
for the current Wolfram session and load by context name:

```mathematica
repo = "/absolute/path/to/diffexp2";
PacletDirectoryLoad[repo];
Needs["DiffExp2`"];
Needs["FeynmanTrick`"];
```

The explicit root-file `Get` forms used elsewhere remain supported and make
the checkout location visible in reproducible scripts.

## Wolfram reference backend

The reference recurrence remains available for parity testing and debugging:

```mathematica
DiffExp2`LoadConfiguration[
  "RecurrenceBackend" -> "Wolfram"
];
```

Calling `LoadConfiguration` resets all other keys to their validated defaults.
Use `UpdateConfiguration` when only the backend should change in an existing
configuration.

## FIRE for Feynman-trick preparation

The subprocess facade looks for FIRE 7.1 at:

```text
Dependencies/fire/FIRE7/FIRE7
```

`Dependencies/` is ignored by Git. A standard source build can be placed at
that location by cloning the official FIRE repository into
`Dependencies/fire/FIRE7`, entering its `FIRE7` subdirectory, and following
FIRE's `./configure`, `make dep`, and `make` instructions. The modular path
expects these nonempty executables below the configured FIRE root:

```text
bin/FIRE7
bin/FIRE7p
bin/FIRE7mp
bin/FIRE7_MPI
bin/reconstruct
```

Install FIRE there, make that path point to the desired installation, or set
the facade option explicitly:

```mathematica
plan = FeynmanTrick`PipelinePlan[
  "bubble", "FIREPath" -> "/absolute/path/to/fire/FIRE7"
];
```

The CLI equivalent is `FT_FIRE_PATH`. Modular FIRE7 is the default and uses
finite-field Zippel sampling plus exact rational reconstruction. Set
`"FIREBackend" -> "Classical"` only when an exact Classical comparison is
wanted. The clean facade subprocess does not inherit in-memory package options
from its calling kernel.

FIRE 7.1 modular reconstruction has a fixed 16-variable exponent array. The
facade therefore rejects modular families with more than 16 FIRE variables,
including the dimension variable, instead of risking memory corruption.
Classical FIRE7 execution does not use that reconstruction path.

The facade executes the ladder in a clean `wolframscript` subprocess. Its
current inherited-environment launcher requires `/usr/bin/env`, so facade
execution is supported on macOS and Linux, not native Windows. This does not
limit the direct DiffExp 2 solver. The command-line runner may be invoked
separately with environment syntax appropriate to the host shell.

Each FIRE invocation has a watchdog. Set the facade option
`"FIRETimeoutSeconds"` or the CLI variable `FT_FIRE_TIMEOUT_SECONDS` to
override it. A timeout or nonzero FIRE exit is a failed preparation, not a
partial result.

## Cache locations

Feynman-trick runs use two independent stores:

- prepared FIRE reductions;
- DiffExp 2 endpoint transports and level-boundary checkpoints.

With the facade, set them explicitly when persistent locations are desired:

```mathematica
run = FeynmanTrick`RunIntegrationPipeline[
  "banana_unequal",
  "PreparedCacheDirectory" ->
    FileNameJoin[{$HomeDirectory, ".cache", "diffexp2", "fire"}],
  "CheckpointDirectory" ->
    FileNameJoin[{$HomeDirectory, ".cache", "diffexp2", "ladder"}]
];
```

The corresponding command-line variables are:

```sh
export FT_PREP_CACHE_DIR="$HOME/.cache/diffexp2/fire"
export FT_LADDER_CHECKPOINT_DIR="$HOME/.cache/diffexp2/ladder"
```

Use the facade option `"RebuildPreparation" -> True` or
`FT_REBUILD_PREP=1` only when the prepared topology must be regenerated. See
[Feynman Trick](FeynmanTrick.md) for checkpoint resumption.

## Troubleshooting

- If the backend is not found, confirm its path and set `DE2_CPP_LIBRARY`
  before the kernel starts.
- On macOS, use `file` and `otool -L` to check architecture and shared-library
  dependencies.
- If CMake cannot find FLINT or Boost, set `PKG_CONFIG_PATH` and
  `CMAKE_PREFIX_PATH`, then configure a fresh build directory.
- Restart the kernel after rebuilding a library that it already loaded.
- Lower `DE2_CPP_THREADS` or the facade option `"CppThreads"` if a
  high-dimensional recurrence exhausts memory.

The backend architecture, supported coefficient fields, and performance
measurements are documented in [C++ Backend](CppBackend.md).
