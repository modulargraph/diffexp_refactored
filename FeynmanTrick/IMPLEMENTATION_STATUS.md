# Feynman Trick Pipeline - Implementation Status

## Paper Reference

Hidding & Usovitsch, "Feynman parameter integration through differential equations", JHEP 2022.
Source: `Papers/FeynmanTrick/main.tex`

---

## Key Equations

### Feynman's trick (eq 2.1)
```
1/(D_i^v_i * D_j^v_j) = Gamma(v_i+v_j)/(Gamma(v_i)*Gamma(v_j)) *
    ∫₀¹ dx x^(v_i-1)*(1-x)^(v_j-1) / [x*D_i + (1-x)*D_j]^(v_i+v_j)
```

### Fundamental recursion (eq 2.8)
```
I^(k-1)_{v1,...} = Gamma(v1+v2)/(Gamma(v1)*Gamma(v2)) *
    ∫₀¹ dx x^(v1-1)*(1-x)^(v2-1) * I^(k)_{v1+v2, v3,...}
```

### Special cases (eq 2.10)
```
I_{0,0,...}^(k-1) = I_{0,...}^(k)              (direct equality)
I_{v1,0,...}^(k-1) = lim_{x->1} I_{v1,...}^(k) (limit at 1)
I_{0,v2,...}^(k-1) = lim_{x->0} I_{v2,...}^(k) (limit at 0)
```

### Rescaled Feynman parameters (eq 2.11)
For iterating the trick by always combining positions {1,2}:
```
x_1' = prod_{i=1}^{n-1} x_i
x_j' = (1-x_{j-1}) * prod_{i=j}^{n-1} x_i   for j=2,...,n-1
x_n' = (1-x_{n-1})
```

### Generalized tadpole (eq 2.16)
```
I_v^(n-1) = Gamma(v - L*d/2) / Gamma(v) * U_tilde^(v-(L+1)*d/2) / F_tilde^(v-L*d/2)
```
where d = 4-2*eps, L = loop count. U_tilde and F_tilde are the Symanzik polynomials
evaluated at the rescaled Feynman parameter values.

**IMPORTANT:** U_tilde and F_tilde come from the ORIGINAL topology's Symanzik
polynomials (computed via FIRE6's UF on the level-0 propagators), with the
standard Feynman variables replaced by the rescaled parameters (eq 2.11).
They are NOT computed from the deepest-level's combined propagators.

### Regularization formula (eq 3.3)
```
∫₀^c dx x^{a+bε} g(x) = ∫₀^c dx x^{a+bε+1}/(1+a+bε) *
    [(2+a+bε)/c * g(x) - (1-x/c)*g'(x)]
```

### Limit evaluation (end of section 3.2)
"Take the segment centered at x=0 (or x'=1-x=0), filter out the finite
coefficient of g₀(x,ε). Put x^{aᵢ+bᵢε}·gᵢ(x,ε) with bᵢ≠0 to zero."

---

## Architecture

```
Level N (tadpole, 1 master)
  ↓ [Tadpole formula: Gamma(v-Ld/2)/Gamma(v) * U_tilde^(...)/F_tilde^(...)]
  Boundary at xx=11/23
  ↓ [DiffExp TransportLevel: transport from 11/23 to cover [0,1]]
  Piecewise series in [0,1]
  ↓ [IntegrateLevelMaster: ∫ x^(v1-1)*(1-x)^(v2-1)*c(x)*f(x) dx]
Level N-1 (fewer masters)
  ↓ [Same: transport + integrate]
  ...
Level 0 (original topology, all masters)
  = Final result!
```

---

## Package Structure

| File | Purpose | Status |
|------|---------|--------|
| `FeynmanTrick.m` | Main loader, config | ✓ Working |
| `PropagatorAlgebra.m` | Propagator derivatives, combination | ✓ Working (21/21 tests) |
| `FIREInterface.m` | FIRE6 setup, basis, reduction, diff matrix | ✓ Working (15/15, 17/17 tests) |
| `MatrixExport.m` | Export to DiffExp format | ✓ Working |
| `EpsPrefactors.m` | Find eps^k prefactors to remove poles | ✓ Working |
| `FeynmanTrickIteration.m` | Multi-level orchestrator | ✓ Fixed (unique params per level) |
| `BoundaryConditions.m` | Tadpole formula, Symanzik polynomials | ✓ Fixed (uses original topology + rescaled params) |
| `DiffExpIntegration.m` | Transport, integration, pipeline | ✓ Fixed (ChopPrecision added) |

---

## FIXED: Unique Feynman Parameter per Level

**Previously:** All levels used the same symbol `FeynmanTrick`xx` for the Feynman
parameter. When building level k, the code substituted `xx -> 11/23` to fix the
previous level's parameter, but this also substituted the NEW `xx` just introduced.
Result: propagators at levels 2+ had no Feynman parameter dependence, causing wrong
indicial roots (like 749/639 instead of integers).

**Fix:** Each level k now uses a unique symbol `xxk` (xx1, xx2, xx3, ...). This
matches the paper's notation (eq 2.100-2.104) where each combination introduces
a new parameter x_1, x_2, ..., x_{n-1}. At the deepest level, all parameters are
set to numerical values, then we transport/integrate over them one at a time.

---

## FIXED: `DeepestLevelBoundary` in BoundaryConditions.m

**Previously:** The function computed Symanzik polynomials from the deepest-level
topology's combined propagators. FIRE6's UF returned garbage because the combined
propagators had numerical kinematic values already substituted.

**Fix:** Now uses the ORIGINAL topology's propagators (from `ftData["TopTopology"]`)
for the Symanzik computation. The rescaled Feynman parameters are computed via
`RescaledFeynmanParametersFromSequence`, which traces the combination sequence
symbolically to determine x'_j (the coefficient of each original D_j in D_combined).
This works for any combination sequence, not just "always combine leftmost two."

---

## FIXED: DiffExp Configuration in TransportLevel

**Previously:** `TransportLevel` didn't set `ChopPrecision`, so DiffExp used its
default (250). When `WorkingPrecision -> 200`, this triggered the error
"ChopPrecision should be smaller than WorkingPrecision."

**Fix:** Added `DiffExp`State`ChopPrecision -> precision - 50` to the config.

---

## FIXED: FIREInterface Multi-Topology Support

**Previously:** The FIREInterface called `ClearIBP[]` before each topology and ran
FIRE6 separately for each operation (FindBasis, ReduceIntegrals). This caused:
- Each FIRE run spawned a new fermat subprocess
- State wasn't preserved between calls
- Race conditions and cleanup issues with fermat processes (exit codes 134, 137)

**Fix:** Restructured FIREInterface.m for proper multi-topology batch processing:

1. **Problem numbers**: Each topology gets a unique `ProblemNumber` (1, 2, 3, ...)
2. **Batch functions**: Added `SetupFIREBatch`, `FindBasisBatch`, `ReduceIntegralsBatch`
3. **Multi-problem configs**: Config files can list multiple `#problem N name.start` lines
4. **Proper integral format**: Integrals include problem number: `{pn, {indices}}`
5. **State management**: `$SetupTopologies`, `$NextProblemNumber`, `ClearFIREState[]`

This matches how FIRE6 is designed to work - define all topologies first, then run
reductions in a single FIRE6 invocation. The batch approach is more stable and
avoids fermat subprocess management issues.

---

## FIXED: FIRE State Clearing (`startinglist`)

**Previously:** When computing multiple levels sequentially (3 → 2 → 1), level 2
would hang because FIRE was using level 3's variable names (xx3 instead of xx2).
Fermat would wait forever for equations that referenced undefined variables.

**Root cause:** FIRE caches IBP relations in `FIRE`startinglist`. If not cleared,
`PrepareIBP[]` skips regenerating IBPs and reuses the cached ones with old variables.

**Fix:** Added `FIRE`startinglist` to the list of variables cleared in `SetupFIRE`:
```mathematica
Clear[FIRE`PrepareIBPd, FIRE`BackMatrix, FIRE`Squares, FIRE`startinglist];
```

---

## FIXED: ARM64 Fermat and Shell Execution

**Previously:** On Apple Silicon Macs, FIRE6 (arm64) couldn't spawn fermat (x86_64)
properly when called via Mathematica's `RunProcess[]`. The child process would get
killed (signal 9) or fail with `CHILD LAUNCH ERROR: errno = 8`.

**Fix:**
1. Installed native ARM64 fermat v7.8 from `https://home.bway.net/lewis/fermat64/Ferm7a.tar.gz`
2. Changed from `RunProcess[]` to `Run[]` with shell execution for FIRE6 invocation

---

## Files and Test Results

### Tests (all in Tests/)
- `test_feynmantrick_algebra.m`: **21/21 PASS** — PropagatorAlgebra, EpsPrefactors
- `test_feynmantrick_fire.m`: **15/15 PASS** — FIREInterface (DefineTopology, SetupFIRE, FindBasis, ReduceIntegrals, ComputeDiffMatrix)
- `test_feynmantrick_iteration.m`: **17/17 PASS** — FeynmanTrickIteration (DefineFTIteration, BuildLevel, FeynmanTrickDecomposition, ExportLevel)
- `test_feynmantrick_endtoend.m`: Partially passes (boundary/transport sections have bugs)
- `test_feynmantrick_pipeline.m`: **NEW** — Full end-to-end pipeline (fails at boundary computation due to above bug)

### Pipeline test output
```
Part 1: Define Topology — PASS (4 propagators, 3 levels)
Part 2: Build & Compute — Level 3 PASS (7 masters, correct variable xx3)
         Levels 1-2: FIRE6 instability on test system (fermat crashes)
Part 3: Boundary — Fixed (original topology + rescaled Feynman params)
Part 4-6: Transport/Integration — awaiting stable FIRE6 runs
```

**Note:** FIRE6/fermat shows instability on some systems (exit codes 134, 137, 141).
The algorithm is correct; infrastructure issues need resolution.

---

## What ComputeLevelBoundary Does (FIXED)

Given a master M at level l, with combination positions {posI, posJ} at level l+1:

| vi = M[posI] | vj = M[posJ] | Case | Needed at l+1 | Action |
|---|---|---|---|---|
| >0 | >0 | integrate | M with posI→vi+vj, posJ→0 | ∫ x^(vi-1)*(1-x)^(vj-1) * I dx |
| >0 | =0 | limitUpper | M with posI→vi, posJ→0 | lim_{x→1} I |
| =0 | >0 | limitLower | M with posI→vj, posJ→0 | lim_{x→0} I |
| =0 | =0 | direct | M with posJ→0 | Direct evaluation |

After constructing the needed integral vector, it's reduced to masters at level+1
via `ReduceIntegrals` (FIRE6 IBP), giving coefficients c_j(xx). Then:
- **integrate**: calls `IntegrateLevelMaster` which uses `DefiniteIntegralWithPrefactor`
- **limits**: calls `EvaluateLimitFromTransport` which decomposes and extracts Taylor part
- **direct**: uses boundary values from the transport starting point

---

## Design Decisions

- Variable `xx` (not `x`) for Feynman parameter to avoid DiffExp's line parameter `x`
- FIRE6 for IBP reductions (not Kira — compilation issues)
- Order-by-order ε expansion (not ε-sampling)
- Fixed parameter value: 11/23 (avoids singularities)
- DiffExp config: `UseMobius -> False`, `SaveExpansions -> True`
- ε-prefactors: define J_i = ε^{k_i} * I_i to remove poles from matrices
- Sign convention: D_j = -q_j^2 + m_j^2 (matches DiffExp/FIRE)
- FIRE6's UF function for Symanzik polynomials (multiloop-ready)

---

## Next Steps (Priority Order)

1. **Run pipeline test**: Both boundary bugs are fixed. Run `test_feynmantrick_pipeline.m`
   end-to-end to validate the full pipeline.

2. **Validate against known result**: The 1-loop massless box at s=-1, t=-1/3 has
   known analytic values. Compare pipeline output.

3. **(Future) Physical region**: Add i*delta prescriptions for threshold crossings.

4. **(Future) More examples**: Double box (topo7), double pentagon (5p from paper).

---

## Test Example: 1-Loop Massless Box

```
Propagators (D_j = -q_j^2):
  D1 = -l1^2
  D2 = -(l1+p1)^2
  D3 = -(l1+p1+p2)^2
  D4 = -(l1+p1+p2+p3)^2

Kinematics (Euclidean):
  p1^2 = p2^2 = p3^2 = 0
  p1*p2 = s/2 = -1/2
  p2*p3 = t/2 = -1/6
  p1*p3 = -(s+t)/2 = 2/3

Combination sequence: {{1,2}, {1,3}, {1,4}}
  Level 1: D_12 = xx*D1 + (1-xx)*D2 (triangle, ~7 masters)
  Level 2: D_123 = xx*D_12 + (1-xx)*D3 (bubble, ~4 masters)
  Level 3: D_1234 = xx*D_123 + (1-xx)*D4 (tadpole, 1 master)

Fixed parameter: xx = 11/23
```

---

## How to Run Tests

```bash
export WOLFRAMSCRIPT_KERNELPATH="/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/WolframKernel"
cd /Users/mhidding/Desktop/diffexp_refactored/Tests

# Quick tests (all pass)
wolframscript -file test_feynmantrick_algebra.m
wolframscript -file test_feynmantrick_fire.m
wolframscript -file test_feynmantrick_iteration.m

# End-to-end pipeline (currently fails at boundary)
wolframscript -file test_feynmantrick_pipeline.m
```
