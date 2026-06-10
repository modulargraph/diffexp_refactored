# Feynman Trick Banana Status

Date: 2026-06-09 (updated 2026-06-10: level 1 matches pySecDec; level 0 in
progress, four further local-solver bugs found and fixed)

This note records the Feynman trick banana debugging work against pySecDec.
The level-1 mismatches described in earlier revisions of this note are
RESOLVED: all banana level-1 reference masters match pySecDec through the
finite term. The root causes and fixes are documented below, followed by
the state of the level-0 (full banana) verification.

## Current Reference Setup

The pySecDec reference setup exists outside the repo in:

- `/tmp/pysecdec-venv`
- `/tmp/pysecdec-ft-family-bubble/results.json`
- `/tmp/pysecdec-ft-family-sunrise/results.json`
- `/tmp/pysecdec-ft-family-banana-deep/results.json`
- `/tmp/pysecdec-ft-family-banana-level1-pos-factored/results.json`

The comparison scripts use the Feynman trick parameter fixed at `11/23`.

The pySecDec convention used for the exported family specs is:

```text
prefactor = gamma(v - L + L*eps) / product_i gamma(v_i)
U power   = v - (L + 1) + (L + 1)*eps
F power   = -v + L - L*eps
```

## Verification Status (all match pySecDec)

Reproduce with:

```sh
env WolframKernel='/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/WolframKernel' \
FT_EXAMPLES=banana \
FT_STOP_AFTER_BOUNDARY_LEVEL=1 \
FT_WORKING_PRECISION=500 \
FT_EXPANSION_ORDER=50 \
FT_DIVISION_ORDER=4 \
FT_BOUNDARY_EXTRA_ORDER=8 \
wolframscript -file temp/feynmantrick_stepwise_pysecdec_check.m
```

Level 3 and level 2 (unchanged, matched before):

- Level 3 `{1,0,0,0}`: `4.54705/eps - 2.63571`
- Level 2 `{1,0,0,2}`: `-0.0527956607/eps^2 - 1.8287423769/eps + 3.096706646`
- Level 2 `{1,0,0,1}`: `-0.5435916861/eps^2 - 0.7399338561/eps + 3.819725487`
- Level 2 `{1,0,0,0}`: `5.5968978699/eps - 4.173153149`

Level 1 (all previously failing entries now match):

| Master      | DiffExp (this repo)                                              | pySecDec                                       |
|-------------|------------------------------------------------------------------|------------------------------------------------|
| `{2,0,1,1}` | `5.4025802965`                                                   | `5.4025802965`                                 |
| `{1,0,1,0}` | `-2.0037878788/eps^2 + 2.2439088015/eps - 13.5042252469`         | same                                           |
| `{1,0,1,1}` | `1/3/eps^3 - 0.1144868285/eps^2 + 1.3904392399/eps - 5.8738356912` | `0.3333333333/eps^3 - 0.114487/eps^2 + 1.39044/eps - 5.87384` |
| `{1,0,1,2}` | `0.5/eps^2 - 0.1717302427/eps - 0.9025092631`                    | `0.5/eps^2 - 0.171730/eps - 0.902509`          |
| `{1,0,2,1}` | `0.5/eps^2 - 0.1717302427/eps - 0.9025092631`                    | same as `{1,0,1,2}`                            |

All imaginary parts are at or below `1e-16` (numerical noise).

## Root Causes Found

The two remaining level-1 failures had distinct mechanisms, both inside the
multi-sector endpoint machinery flagged by earlier rounds of this note.
Multi-sector means: several local solutions `x^(a + b_i eps)` share the same
integer power `a` at a singular point (for banana segment 1 near
`xx2 = 0.048`: effective `b` values 2 and 1 over `a = -1`), which
`DecomposeSingularity` collapses into a single averaged exponent, leaving
residual `Logx` towers that `IntegrateAnalyticRegularizedBySubtractionLaurent`
must resum back into sectors.

### Root cause 1: insufficient epsilon lookahead (hit `{1,0,1,1}`, `{1,0,1,2}`)

`RequiredTransportEpsilonOrder` added only `integrationPoleAllowance = 1`
extra epsilon order per level. The residual-sector resummation needs more
input orders than a plain endpoint pole:

- each endpoint pole deepens the Laurent window by one order,
- solving the sector weights at epsilon offset `q` consumes offsets
  `q .. q + sectorCount - 1`,
- the top offsets serve as the truncation boundary.

Worse, `IntegrateCombinedMasters` silently dropped IBP-coefficient times
transport-order products whose shifted index `n + k_j - k` exceeded the
available transport orders, so the top output orders of the combined series
were *incomplete but present* - and the old fitter consumed them as data.
This is why `FT_BOUNDARY_EXTRA_ORDER` alone could not help: the per-level
required order capped the transport before the extra orders arrived.

The eps^-1 mismatch of `{1,0,1,1}` (`2.31409 + 1.99613 I` instead of
`1.39044`) came from sector coefficients solved out of those partially
combined top orders.

### Root cause 2: the `s0 = 0` bail-out (hit `{1,0,2,1}`)

For `{1,0,2,1}` the two sector weights cancel at their leading epsilon
order: the first nonzero tower coefficient is pure `Logx`, so the leading
moment `s0` vanishes. `fitResidualEndpointSectors` returned `$Failed`
whenever `s0 == 0` and fell back to explicit-log integration against the
(here complex, meaningless) averaged exponent from `DecomposeSingularity`.
The Hankel system is perfectly solvable in this case (`det = -s1^2 != 0`),
and the recovered sectors reproduce `0.5/eps^2` exactly.

## Fixes Implemented

### `DiffExp/RegularizedIntegration.m`

1. `fitResidualEndpointSectors` rewritten:
   - Works on absolute epsilon offsets (no first-nonzero shift, which
     misaligned the moment diagonal when leading weights cancel).
   - Scans the moment-ladder reference order `q0` so cancelling leading
     weight slices (`s0 = 0`) still yield Prony root candidates; the
     two-sector Hankel solve no longer bails on `s0 = 0`.
   - Candidate fits (one- and two-sector) are validated against the FULL
     reconstructed `Logx` content of every available offset, with a
     relative tolerance; the fit is only trusted on the leading run of
     validated offsets, retreating one extra order when validation fails
     mid-tower (polluted upstream data can hide in the by-construction
     rows for one offset).
   - Offsets beyond the trusted run have their non-log content salvaged
     against the plain branch exponent (the old truncation-boundary
     convention) and their `Logx` content dropped and counted.
2. Warnings now distinguish "dropped data only affects orders beyond the
   requested window" from "dropped data affects the reported Laurent window
   starting at eps^k; results at and above k are NOT trustworthy".
   On the stale (insufficient lookahead) dumps the flagged order is exactly
   where pySecDec disagreement started.
3. The `b = 0` unit-regulator substitution is restricted to coefficient
   lists without `Logx` towers, and towers now always take the subtraction
   path: a vanishing averaged exponent with surviving towers (perfectly
   cancelling weights) otherwise had every recovered sector shifted by one.
4. A hard warning fires when the sector fit fails and the explicit-log
   fallback runs on data that still carries `Logx` towers (the fallback is
   wrong for genuinely multi-sector data).

### `FeynmanTrick/FeynmanTrick.m`, `FeynmanTrick/DiffExpIntegration.m`

5. New config option `IntegrationPoleAllowance` (default 4 = 1 pole +
   2 sector-solve + 1 truncation guard), used by
   `RequiredTransportEpsilonOrder`; overridable with the
   `FT_INTEGRATION_POLE_ALLOWANCE` environment variable.
6. `IntegrateCombinedMasters` trims the combined epsilon window to the
   provably complete orders (every active master must still have transport
   data at the shifted index) instead of emitting partially combined top
   orders, and warns when the trim cuts into the requested window.
7. The level driver (`RunIntegrationPipeline` and the stepwise script)
   warns when the available boundary orders cap the transport below
   `RequiredTransportEpsilonOrder`.

### Tests (`Tests/test_regularized_integration_edge_cases.m`)

New regression tests (all 25 pass, as do the other suites):

- two sectors with equal `a`, different `b`, constant and
  epsilon-dependent weights, recovered exactly;
- cancelling leading weights (`s0 = 0`, the `{1,0,2,1}` mechanism);
- the same two-sector recovery on a negative local interval, validated
  against direct per-sector integration (branch-convention independent);
- zero residual root with epsilon-dependent weight (the historical `0^0`
  validation bug);
- a truncated top epsilon order must not corrupt lower Laurent orders
  (the stale-dump mechanism).

## Code Cleanup (2026-06-10)

Dead ends from earlier debugging rounds were removed once the live path
was verified: 18 unused private functions in
`DiffExp/RegularizedIntegration.m` (the abandoned `*LaurentPrepared`
integration path, the `reorient*`/`prepareBranchResolved*` branch-handling
experiments, the `stripReconstructedEndpointLogs` approach superseded by
the residual-sector fit, `ApplyRegularizationStepOld`), 6 unused helpers
and 3 dead debug blocks in `FeynmanTrick/DiffExpIntegration.m`, and the
duplicated non-finite-expression helpers now live once in
`DiffExp`Utilities` (`NonFiniteExpressionQ`, `FiniteAbsMax`,
`ReplaceSparseArrays`).

The canonical end-to-end check for all examples is
`temp/feynmantrick_stepwise_pysecdec_check.m`
(`FT_EXAMPLES=bubble,sunrise,banana`); the box example is covered by
`Tests/test_feynmantrick_pipeline.m`.

## Diagnostic Tools

- `temp/eval_dump_generic.m`: evaluate any Laurent dump
  (`DUMP_FILE=... [PRINT_SEGMENTS=0]`), prints totals and per-segment
  results.
- `temp/debug_banana_l1_sector_fit.m`: inspect the residual sector fit for
  a dump segment (`DUMP_FILE=... SEG_INDEX=... LOCAL_POWER=...`), prints
  moments, Hankel data, roots, weights, and per-q coefficient rows.
- Dumps are produced by setting `DIFFEXP_DUMP_LAURENT_DIR=<dir>` during a
  run; fresh known-good dumps from the final verification run are in
  `/tmp/diffexp_banana_l1_dumps_v2final/` (mapping: 0003 -> {2,0,1,1},
  0004 -> {1,-1,1,1}, 0005 -> {1,0,1,0}, 0006 -> {1,0,1,1},
  0007 -> {1,0,1,2}, 0008 -> {1,0,2,1}; in-run values carry an extra
  Gamma-prefactor factor relative to raw dump replays).

## Banana Level 0 (full banana): Status

Reference for the finite part (direct numerical Feynman-parameter
integration, see `/tmp/banana_direct_feynman_param_c0.m` history):

```text
banana {1,1,1,1} finite = 8.26810451329511583109184...   (no eps poles)
```

The level-0 run (`FT_EXAMPLES=banana`, no stop level) initially produced
spurious `eps^-2`/`eps^-1` poles and O(10^-3) imaginary contamination.
Diagnosing the level-1 -> 0 transport with per-segment ODE-residual checks
(`temp/check_transport_ode_residual.m`) exposed FOUR local-solver bugs at
the singular endpoints `xx1 = 0, 1` of the level-1 system in the
Feynman-trick variable (the level-1 system is the first place where a
nilpotent DOUBLE pole and 5-fold resonant residues appear; bubble/sunrise
and the deeper banana levels never exercise these paths):

1. `TrimFuchsianLattice`/`FuchsianizeLocal` ran the lattice saturation on
   NUMERIC matrices; exact pole-cancellation tests never fire on floats,
   so the trim pass ratcheted lattice columns by x once per pass up to
   x^50-scale entries while the precision collapsed to ~15 digits.  Fixed
   by rationalizing the local theta-matrix before the reduction and
   guarding against degenerate lattices (`BuildFuchsianizedRecurrenceData`).
2. Exact `JordanDecomposition` on (rationalized) noisy residue matrices
   sees noise-split "distinct" eigenvalues with catastrophically
   ill-scaled eigenvectors (the float noise splits degenerate eigenvalues
   by ~sqrt(noise)).  `ComputeResonanceStructure` now detects degenerate
   snapped spectra and builds the Jordan data with toleranced
   generalized-eigenvector chains from SVD null spaces
   (`NumericJordanData`); this also recovers the genuine size-3 Jordan
   block (log^2 solutions) that the old path missed entirely.
3. `PrepareSingularRecurrence` accepted degenerate snapped spectra (its
   rank test cannot see near-parallel eigenvectors and its resonance test
   only rejects positive-integer differences), so blocks could be solved
   by plain diagonalization with a 1e111-amplified basis.  Degenerate
   spectra are now rejected and routed to the general solver.
4. The resonant fundamental-matrix construction resolved kernel freedom of
   singular recurrence steps only for n >= 1; the n = 0 Jordan-chain
   descent used minimum-norm pseudo-inverse solutions whose missing kernel
   components made the next chain equation inconsistent (silently
   projected away).  The n = 0 descent now resolves each step's freedom
   from the next step's solvability, like the n >= 1 loop.

Supporting hardening from the same investigation:

- `DecomposeSingularity` only keeps a cleanly rational extracted `b`;
  a complex/non-rational average (collapsed multi-sector data) now keeps
  `b = 0` with the towers explicit.
- `CheckParticularResidual` evaluates both theta branches instead of
  silently treating symbolic residuals as zero.
- Gated debug tooling: `DiffExp`State`$LogStrategyDispatch` records the
  per-block strategy choices; `DiffExp`State`$DebugFuchsianizedCheck`
  residual-checks every solver output piece against the block ODE.

Current level-0 state after these fixes: the spurious `eps^-2` pole and
all imaginary contamination are gone and level 1 still matches pySecDec,
but the value still disagrees:

```text
current:   -0.3976/eps + 12.9242
reference:  0/eps      +  8.2681
```

All 33 transported segments now satisfy the level-1 ODE to ~1e-400,
match at the anchor (the combined needed-integrand at xx1 = 11/23
reproduces pySecDec's 5.4025802965 with vanishing poles), and are
continuous across segment boundaries.  Yet the combination's pointwise
`eps^-1` order does not vanish away from the anchor, although the IBP
identity says it must.  Since sunrise runs the identical level-0 machinery
to an exact final value, the remaining suspect is banana-specific input:
the level-1 differential matrix in the trick variable (`dxx1`, which has a
nilpotent double pole in the numerator-master row and is the only level
never validated against pySecDec) and/or the IBP tables involving the
numerator master `{1,-1,1,1}`.

## Remaining Notes / Next Steps

1. Validate the banana level-1 `dxx1` matrix / numerator-master IBP layer
   directly: pySecDec the level-1 family at a second parameter point
   (e.g. z0 = 1/4 or 3/4) and compare the transported combination
   pointwise (`temp/check_dump_anchor.m` with `ANCHOR=...`); or verify
   the matrix against an independent derivation.
2. `ShiftRawBoundariesToFinite` still zero-pads shorter boundary lists up
   to the longest master's window; with the per-level warnings and the
   combine-window trim this is detected downstream, but a per-master
   trusted-window tag threaded through transport would be cleaner.
3. The residual-sector fitter handles one and two pure sectors. Three or
   more sectors at one integer power will fail validation and fall back
   loudly.
4. After banana level 0, rerun the full pySecDec comparison matrix:
   bubble 2d/4d, sunrise 2d/4d, and banana
   (`temp/feynmantrick_stepwise_pysecdec_check.m` with
   `FT_EXAMPLES=bubble,sunrise,banana` and no stop level).
   Bubble (0.860817881928) and sunrise (2.2367927002126465) currently
   match through their final values with all fixes in place.
