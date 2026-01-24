# Feynman Trick End-to-End Pipeline - Implementation Status

## What Has Been Done

### Files Created/Modified

1. **`DiffExp/RegularizedIntegration.m`** (MODIFIED)
   - Added `DefiniteIntegralWithPrefactor` function declaration and implementation
   - This function integrates `∫ x^alpha * (upper-x)^beta * r(x) * f(x) dx` where f(x) is a DiffExp piecewise series
   - Handles singular power-law prefactors at boundaries by modifying the decomposition's `a` value
   - Series-expands smooth rational prefactors at each segment's expansion point
   - Added helper `IntegrateSegmentWithPrefactor`

2. **`FeynmanTrick/BoundaryConditions.m`** (NEW)
   - `ComputeSymanzikPolynomials` - uses FIRE6's UF function (multiloop-ready)
   - `RescaledFeynmanParameters` - eq. (2.11) of paper
   - `EvaluateTadpoleBoundary` - Gamma(v-Ld/2)/Gamma(v) * U^(...)/F^(...) expanded in eps
   - `EvaluateTadpoleNumerical` - fixed eps value version
   - `DeepestLevelBoundary` - high-level orchestrator

3. **`FeynmanTrick/DiffExpIntegration.m`** (NEW)
   - `TransportLevel` - loads DiffExp, configures, transports from fixedVal to cover [0,1]
   - `IntegrateLevelMaster` - integrates one master with Feynman trick prefactors
   - `ComputeLevelBoundary` - **HAS A BUG** (see below)
   - `RunIntegrationPipeline` - full bottom-up orchestrator

4. **`FeynmanTrick/FeynmanTrick.m`** (MODIFIED)
   - Added `Get` calls for BoundaryConditions.m and DiffExpIntegration.m
   - Added `WorkingPrecision` config key

5. **`Tests/test_feynmantrick_endtoend.m`** (NEW)
   - Tests BoundaryConditions (rescaled params, Symanzik, tadpole formula)
   - Tests 3-level iteration setup (box with combination {1,2},{1,3},{1,4})
   - Tests matrix computation at each level
   - Tests deepest level boundary
   - Tests DiffExp transport

## What Still Needs To Be Fixed

### Bug in `ComputeLevelBoundary` (DiffExpIntegration.m)

The function incorrectly calls `IdentifyNeededIntegrals[ftData, level, masterVec]` but that function only takes 2 arguments: `IdentifyNeededIntegrals[ftData, level]`.

More importantly, `IdentifyNeededIntegrals` goes the WRONG direction for our use case. It takes masters at level+1 and splits them to get integrals at level. We need the opposite: given a master at level, what integral at level+1 is needed.

**The correct mapping is:**

Given master M = {m_1, ..., m_n} at level l, and combination positions {i,j} for level l+1:

- **If m_i > 0 AND m_j > 0**: Integration case
  - Needed at l+1: M' where M'_i = m_i + m_j, M'_j = 0, rest unchanged
  - Feynman trick: v1 = m_i, v2 = m_j
  - Result: Γ(v1+v2)/(Γ(v1)Γ(v2)) * ∫₀¹ x^(v1-1)*(1-x)^(v2-1) * I^(l+1)_M' dx

- **If m_i > 0 AND m_j == 0**: Limit at x=1
  - Needed at l+1: M' = M (since m_j already 0)
  - Result: lim_{x→1} I^(l+1)_M

- **If m_i == 0 AND m_j > 0**: Limit at x=0
  - Needed at l+1: M' where M'_i = m_j, M'_j = 0, rest unchanged
  - Result: lim_{x→0} I^(l+1)_M'

- **If m_i == 0 AND m_j == 0**: Direct equality
  - I^(l)_M = I^(l+1)_M (no integration)

The needed integral M' at level+1 might NOT be a master. It needs IBP reduction via `ReduceIntegrals` to express in terms of masters at level+1. The IBP coefficients are rational functions of xx.

### Steps to Complete

1. **Fix `ComputeLevelBoundary`** - Rewrite to use the correct recursion mapping above directly, instead of calling IdentifyNeededIntegrals

2. **Fix IBP reduction interface** - The `ReduceIntegrals` return format from FIREInterface needs to be properly parsed. Currently returns a list of `{integral, coefficient}` pairs but the exact format needs verification.

3. **Test the pipeline** - Run `test_feynmantrick_endtoend.m` and debug

4. **Handle eps prefactors consistently** - When boundary has poles in eps, the prefactors from `DeepestLevelBoundary` need to be tracked through the transport and integration steps

5. **(Optional) Sunrise test** - If box works, try a 2-loop sunrise integral

## Architecture

```
Level 3 (tadpole, 1 master)
  ↓ [Tadpole formula: Gamma(v-Ld/2)/Gamma(v) * U^(...)/F^(...)]
  Boundary at xx=11/23
  ↓ [DiffExp TransportLevel: transport from 11/23 to cover [0,1]]
  Piecewise series in [0,1]
  ↓ [IntegrateLevelMaster: ∫ x^(v1-1)*(1-x)^(v2-1)*r(x)*f(x) dx]
Level 2 (bubble, ~2 masters)
  ↓ [Same: transport + integrate]
Level 1 (triangle, ~4 masters)
  ↓ [Same: transport + integrate]
Level 0 (box, original masters)
  = Final result!
```

## Key Technical Details

- DiffExp uses global state → reload config for each level
- Matrix files: `dxx_0.m`, `dxx_1.m` (eps-expanded, rational in xx)
- DiffExp auto-detects `xx` as kinematic variable from filenames
- LineParameter (default Global`x) differs from xx → no conflict
- `UseMobius -> False` required for DefiniteIntegral to work
- `SaveExpansions -> True` required to get segment data
- Fixed parameter value: 11/23 (avoids singularities)
- Euclidean point: s=-1, t=-1/3 (both negative → no thresholds)

## Test Commands

```bash
export WOLFRAMSCRIPT_KERNELPATH="/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/WolframKernel"
cd /Users/mhidding/Desktop/diffexp_refactored/Tests

# Quick test (existing, should still pass)
wolframscript -file test_feynmantrick_algebra.m

# End-to-end test (new)
wolframscript -file test_feynmantrick_endtoend.m
```

## Previous Working Tests (should still pass)

- `test_feynmantrick_algebra.m`: 21/21 PASS
- `test_feynmantrick_fire.m`: 15/15 PASS
- `test_feynmantrick_iteration.m`: 17/17 PASS
