# DiffExp Refactored

DiffExp is a Mathematica package for integrating Feynman integrals in terms of series expansions, using differential equations.

This repository is an experimental, AI-assisted refactoring and upgrade project based on the original DiffExp code. The original reference repository remains at [gitlab.com/hiddingm/diffexp](https://gitlab.com/hiddingm/diffexp).

## Version

This refactored repository is version 2.0.

## What Is New

- Modular package layout under `DiffExp/` instead of a single generated package file.
- A stricter recurrence-based transport path enabled by `UseRationalRecurrence -> True`.
- Recursive finite-width coefficient solvers intended to avoid long truncated products.
- Support for ordinary points, regular singular points, resonant/logarithmic sectors, fractional powers, singular transport, infinity transport, and inhomogeneous multi-sector source terms.
- Focused regression tests comparing recurrence transport against the standard solver on equal-mass and unequal-mass banana examples.

The recurrence methods are still an active experimental upgrade. The older Wronskian/Frobenius strategies remain available, but recurrence failures should be reported directly rather than silently hidden by fallback behavior when `UseRationalRecurrence -> True`.

DiffExp2 also has an opt-in C++/FLINT recurrence backend. It keeps indicial,
resonance, epsilon-window, and analytic-continuation decisions in Wolfram
Language while accelerating the finite-width coefficient recurrence through
LibraryLink. See [Docs/CppBackend.md](Docs/CppBackend.md) for build instructions,
the analytic-regularization contract, tests, and current limitations.

## Loading

From Mathematica or Wolfram Engine:

```mathematica
Get["/path/to/diffexp_refactored/DiffExp.m"]
```

The package prints its version on load.

## DiffExp2 Feynman-Trick runner

Run the stepwise DiffExp2 ladder from the repository root.  Environment
assignments precede `wolframscript`; this banana invocation uses the default
numerical settings:

```sh
FT_EXAMPLES=banana \
FT_WORKING_PRECISION=500 \
FT_EXPANSION_ORDER=50 \
FT_EPS_ORDER=0 \
FT_BOUNDARY_EXTRA_ORDER=4 \
FT_LEVEL_EPS_HALOS=0 \
FT_DIVISION_ORDER=3 \
FT_RADIUS_OF_CONVERGENCE=1 \
wolframscript -file Scripts/run_ft_stepwise2.m
```

After building the optional compiled backend, add
`DE2_RECURRENCE_BACKEND=Cpp DE2_CPP_THREADS=4` to this command. The Wolfram
recurrence remains the default, and a requested C++ backend fails loudly
rather than silently falling back.

`FT_EXAMPLES` is a comma-separated list (default `bubble`); available names
are `bubble`, `sunrise`, `banana`, `box`, `pentagon`, `box_bubble`,
`box_triangle`, and `double_box_planar`.  The other variables set the working
precision, local expansion order, requested epsilon order, extra epsilon
orders retained at level boundaries, optional comma-separated per-level
epsilon halos (levels listed from 1 upward), the classic coupled
predivision/matching divisor, and the affine chart-coordinate radius,
respectively. Adjacent regular charts meet at
`+1/DivisionOrder` and `-1/DivisionOrder`; the runner overrides a differing
legacy `FT_STEP_DIVISION_ORDER` value. A halo computes extra internal
coefficients without raising the downstream requested order; it compensates
only certified finite-window losses such as analytic log/Laurent basis
width. `DE2_VALUE_TRANSPORT=1` opts
into the experimental regular-chart
value-transport path.  It is a prototype, is off by default, and retains its
safety fallbacks to basis transport.

The runner also supports persistent caches and restart checkpoints:

- `FT_PREP_CACHE_DIR` selects the FIRE/preparation cache directory (default:
  the system temporary directory under `DiffExp2_FT_Prepared`).
  `FT_REBUILD_PREP=1` ignores a cached preparation and rebuilds it.  This
  cache stores the completed Feynman-trick/FIRE preparation and reduction
  cache; it is separate from transport-ladder checkpoints and intentionally
  is not invalidated by runner or DiffExp2 edits.
- `FT_LADDER_CHECKPOINT_DIR=/path/to/checkpoints` writes a transport
  checkpoint after each level's expensive two-way endpoint transport and a
  boundary checkpoint after assembling the next level's boundary values.
  Resuming a transport checkpoint reuses that transport and replays boundary
  assembly; resuming a boundary checkpoint starts transport at the named
  lower level.
- `FT_RESUME_LADDER_CHECKPOINT=/path/to/file.mx` resumes either checkpoint
  kind.  When no checkpoint directory is supplied, its directory is also
  used for subsequent checkpoints.  Ladder checkpoints record their example,
  numerical settings, prepared-data key, level metadata, and a fingerprint of
  the runner plus DiffExp2 sources.  Stale or unversioned checkpoints are
  rejected by default; `FT_ALLOW_STALE_LADDER_CHECKPOINT=1` accepts one with a
  warning and marks its descendants tainted.

## Tests

Focused recurrence and singularity tests can be run from the repository root:

```sh
wolframscript -file Tests/test_package_loading.m
wolframscript -file Tests/test_recurrence_no_fallback.m
wolframscript -file Tests/test_local_series_edge_cases.m
wolframscript -file Tests/test_resonant_2f1.m
wolframscript -file Tests/test_singular_recurrence.m
wolframscript -file Tests/test_unequal_mass_banana_parity_speed.m
```

Additional implementation notes are in `Docs/`.
