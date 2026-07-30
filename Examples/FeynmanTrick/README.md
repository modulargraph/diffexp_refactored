# Feynman-trick examples

These examples use the public `FeynmanTrick` pipeline facade. The shared
`RunExample.m` launcher constructs a validated `PipelinePlan` and executes it
with `RunIntegrationPipeline`; the shell files are only convenient names for
those public calls. They require FIRE for a cold preparation and the compiled
recurrence backend. Run the commands below from the repository root with
`sh`. The scripts may also be invoked by absolute path from another
directory.

| Example | Command | Status in the documentation snapshot |
| --- | --- | --- |
| bubble | ``sh Examples/FeynmanTrick/Bubble.sh`` | small recommended first run |
| sunrise | ``sh Examples/FeynmanTrick/Sunrise.sh`` | completed parity fixture |
| unequal three-loop banana | ``sh Examples/FeynmanTrick/UnequalBanana.sh`` | completed C++/oracle comparison |
| box-bubble | ``sh Examples/FeynmanTrick/BoxBubble.sh`` | completed Euclidean ladder |
| massive kite | ``sh Examples/FeynmanTrick/Kite.sh`` | completed Euclidean ladder |
| unequal four-loop banana | ``sh Examples/FeynmanTrick/FourLoopUnequalBanana.sh`` | experimental; no source-controlled completed ladder result in this snapshot |
| Henn double-pentagon boundary | ``sh Examples/FeynmanTrick/HennDoublePentagonBoundary.sh`` | opt-in FT reconstruction of a published canonical boundary component (`--all` for 108) |

The numerical profiles are centralized in `RunExample.m`. Every profile:

- selects the C++ recurrence backend;
- accepts an existing `DE2_CPP_THREADS` value or uses four workers;
- enables regular-chart value transport and native value-hop execution;
- requests paired C++ endpoint-arm prewarming;
- keeps FIRE preparation under
  `${DIFFEXP2_CACHE_DIR:-$HOME/.cache/diffexp2}/fire`;
- keeps resumable ladder checkpoints in an example-specific subdirectory;
- records all precision, expansion, epsilon-lookahead, division, and radius
  settings in the pipeline plan.

Inspect one exact plan without running it:

```sh
sh Examples/FeynmanTrick/Kite.sh --plan
```

Validate every documented profile:

```sh
wolframscript -script Examples/FeynmanTrick/RunExample.m --check
```

The generic launcher is also usable directly:

```sh
wolframscript -script Examples/FeynmanTrick/RunExample.m bubble
```

Set a different cache root without editing a script:

```sh
DIFFEXP2_CACHE_DIR=/fast/local/cache \
DE2_CPP_THREADS=8 \
sh Examples/FeynmanTrick/UnequalBanana.sh
```

Set `FT_FIRE_PATH=/absolute/path/to/fire/FIRE7` when FIRE is not installed at
`Dependencies/fire/FIRE7/FIRE7`; the wrappers preserve that environment
variable.

To force a new FIRE preparation, add ``FT_REBUILD_PREP=1`` to the environment.
That is normally unnecessary: solver-source changes do not invalidate the
separate prepared-data cache.

The exact propagators and Euclidean points are centralized in
`Scripts/FTExamples.m`, which the pipeline implementation loads. Their conventions and all
controls are documented in [Feynman Trick](../../Docs/FeynmanTrick.md).

## Henn double-pentagon boundary

`HennDoublePentagonBoundary.sh` is an additional, deliberately heavier public
family example. It independently reconstructs the canonical boundary at
`X0={3,-1,1,1,-1}` from the scalar integrals in the ancillary dlog basis of
Chicherina et al., [arXiv:1812.11160](https://arxiv.org/abs/1812.11160), and
compares all five published coefficients through epsilon order 4. The
maintained default is canonical component 1; `--all` requests the complete
108-component boundary.

Fetch the hash-pinned ancillary files once, inspect the exact universal-family
plan, and run the comparison:

```sh
Scripts/fetch_henn_nonplanar_data.sh
sh Examples/FeynmanTrick/HennDoublePentagonBoundary.sh --plan
sh Examples/FeynmanTrick/HennDoublePentagonBoundary.sh
sh Examples/FeynmanTrick/HennDoublePentagonBoundary.sh --all
```

The complete 108-component mode contains 257 distinct scalar integrals, but it
does not run 257 recursions. One family with physical denominators
`D1,...,D8` discovers and transports all 108 L0 FIRE masters through one
seven-level ladder. One exact L0 FIRE reduction maps the 257 scalar
observables onto that basis using resumable sector scheduling, after which the
published canonical combinations are formed.
`D9,D10,D11` are the genuine irreducible-numerator slots; negative powers of
physical denominators are handled by the L0 reduction. Subsector masters with
an absent merged line use the recursion's typed lower-limit, upper-limit, or
direct boundary case instead of an invalid beta integral with zero exponent.

Any individual canonical component can be selected explicitly:

```sh
sh Examples/FeynmanTrick/HennDoublePentagonBoundary.sh --component 1
```

The launcher avoids nested Wolfram kernels: one invocation writes the exact
runner manifest, a small Python driver executes its plans sequentially, and a
final invocation performs the paper comparison. It prefers Mathematica's
`/Applications/Wolfram.app` installation, falls back to Wolfram Engine, and
honors an explicit `HENN_FT_WOLFRAMSCRIPT` path. Numerical controls are
`HENN_FT_WORKING_PRECISION`, `HENN_FT_MATCH_DIGITS`,
`HENN_FT_MATCHING_CERTIFICATION_DIGITS`,
`HENN_FT_EXPANSION_ORDER`, `HENN_FT_BOUNDARY_EXTRA_ORDER`,
`HENN_FT_MASTER_EPSILON_ORDER`, `HENN_FT_DIVISION_ORDER`,
`HENN_FT_FIRE_TIMEOUT_SECONDS`, `HENN_FT_L0_FIRE_BACKEND`,
`HENN_FT_L0_FIRE_THREADS`, `HENN_FT_L0_REDUCTION_BATCH_SIZE`,
`HENN_FT_CHECKPOINT_DIR`,
`HENN_FT_DEADLINE_SECONDS`, and `DE2_CPP_THREADS`.
The physical continuation is fixed by the example rather than exposed as a
performance knob: the deepest Symanzik value is evaluated on the lower rim,
and the level contours use `{+1,-1,...}` from level 1 upward. The default uses
working precision 80, Taylor order 25, an eight-digit publication contract,
and fifth-radius chart overlaps for the closely spaced post-singular geometry.
It uses distinct exact anchors
`{1/5,3/10,2/5,1/2,3/5,7/10,4/5}` because a repeated common anchor would put
the second-level system directly on the diagonal singularity `x2=x1`.
The runner propagates guarded accuracy through the full ladder and caches the
resolved master basis, exact FIRE work, request-specific epsilon width, and
producer-accuracy requirements discovered by the rigorous matching checks.
