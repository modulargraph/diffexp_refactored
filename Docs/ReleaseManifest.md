# Release Tree Manifest

This branch is the user-facing DiffExp 2 source release. The complete
development history, legacy implementation, benchmark campaigns, design
notes, and large oracle collection remain preserved on `codex/develop` and at
the annotated tag `dev-snapshot-2026-07-11-pre-release`.

`Scripts/check_release_tree.py` enforces the release contents as an exact
allowlist. It rejects both missing files and accidental additions.

## Product surface

| Path | Purpose |
| --- | --- |
| `DiffExp2.m` | Stable root loader; selects the compiled recurrence backend by default |
| `DiffExp2/` | Exact recurrence, transport, matching, endpoint, integration, and public API modules |
| `FeynmanTrick.m` | Stable root loader for Feynman-trick preparation and execution |
| selected `FeynmanTrick/*.m` | Topology algebra, FIRE interface, iteration, exact level reductions, boundary data, and pipeline facade |
| `cpp/`, `CMakeLists.txt` | C++20/FLINT recurrence backend and native test |
| `Examples/` | Runnable direct and Feynman-trick examples |
| selected `Docs/*.md` | Installation, tutorials, API, continuation, results, migration, citation, and backend documentation |
| `Scripts/run_ft_stepwise2.m` | Transitional runtime used by the public Feynman-trick subprocess facade |
| `Scripts/FTExamples.m` | Built-in example/topology registry used by that runtime |
| `Scripts/run_release_tests.sh` | Sequential release verification |
| selected `Tests/*.m` and fixtures | Focused exact, compiled-backend, facade, cache, checkpoint, and oracle coverage |

The public Feynman-trick runner uses `FeynmanTrick/LevelReduction.m`. That
module keeps the exact FIRE reduction, analytic-regulator, and epsilon-budget
logic independent of DiffExp 1.

## Deliberately archived on development

The release does not track:

- the `DiffExp/` tree or the old root `DiffExp.m` loader;
- `FeynmanTrick/DiffExpIntegration.m`, the legacy DiffExp 1 bridge;
- the original monolithic source under `Reference/`;
- design specifications, agent reviews, performance campaign notes, and
  historical status reports;
- benchmark/profiling scripts, external-oracle generation environments, and
  large diagnostic logs;
- matrix fixtures and tests used only by the archived DiffExp 1/reference
  implementation.

No archived material was discarded or history-rewritten. The unsafe
persistent FIRE-storage experiment is retained separately at
`research/fire-storage-no-go-2026-07-11`; it is not part of the release.

## Verification

From a built source tree, run:

```sh
python3 Scripts/check_release_tree.py
Scripts/run_release_tests.sh
```

The release runner executes the exact Wolfram test allowlist, requires the C++
backend for compiled parity tests, runs both direct examples, and checks every
shipped shell wrapper. Set `DE2_RUN_FIRE_TESTS=1` only when a FIRE 6
installation is configured and the optional real-FIRE integration tests are
desired.

The FIRE reduction cache is invocation-local and content-addressed. DiffExp 2
does not reuse FIRE persistent databases across invocations. Prepared
topology snapshots and transport checkpoints validate their source,
configuration, runtime, dimension-variable, and exact request fingerprints
before reuse.

## External dependencies and scope

FIRE is not vendored. It is needed to prepare a new Feynman-trick family, but
not to transport an existing exact differential system. The public facade
currently starts one clean `wolframscript` child, which can occupy a second
Wolfram license seat; native C++ workers do not consume Wolfram seats.

`PacletInfo.wl` registers the built source tree and its two root contexts.
The compiled library is discovered below `build/cpp` or through the explicit
`DE2_CPP_LIBRARY` environment variable; this release does not claim a
precompiled relocatable binary paclet.

The exact input matrix must be rational in the line variable and epsilon.
Non-integer powers with a line-variable-dependent base are rejected. Exact
algebraic constants, algebraic singular locations, and algebraic local
exponents are not rejected merely for being algebraic. Indicial exponents must
be affine in epsilon, `a + b eps`.
