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
| `FeynmanTrickIteration.m` | Multi-level orchestrator | ✓ Working |
| `BoundaryConditions.m` | Tadpole formula, Symanzik polynomials | **BUG** (see below) |
| `DiffExpIntegration.m` | Transport, integration, pipeline | ✓ Fixed (untested end-to-end) |

---

## Current Bug: `DeepestLevelBoundary` in BoundaryConditions.m

### Problem

The function computes Symanzik polynomials from the DEEPEST-LEVEL topology's
propagators, which are already combined. This is wrong. The output shows:
```
U = {l1}           ← This is just the loop momenta list!
F = {expressions}  ← These are the combined propagators, not Symanzik F!
```

FIRE6's UF function returns `{loopMomenta, propagators, replacements}` unchanged
because:
1. The combined propagators at the deepest level have numerical kinematic values
   already substituted, making UF unable to detect the quadratic structure
2. The UF call may not be reaching FIRE6 at all (context/loading issue)

### Correct Approach

Per the paper (eq 2.16), U_tilde and F_tilde should be computed as:

1. Call `UF[loopMomenta, propagators, replacements]` on the **ORIGINAL** (level 0)
   topology's propagators (symbolic, with proper scalar product rules)
2. This gives `{U(x[1],...,x[n]), F(x[1],...,x[n]), {x[1],...,x[n]}}`
   where x[i] are the standard Feynman parameters
3. Compute the rescaled parameters `x'_j` as functions of the iteration's
   Feynman parameters (following eq 2.11 adapted for the specific combination sequence)
4. Substitute the rescaled parameters into U and F: `U_tilde = U /. x[j] -> x'_j`
5. Set all Feynman parameters to 11/23, giving numerical U_tilde, F_tilde
6. Apply the tadpole formula: `Gamma(v-Ld/2)/Gamma(v) * U_tilde^(...) / F_tilde^(...)`

### Additional Fix Needed

The rescaling formula in eq 2.11 is written for the case where we always combine
positions {1,2} (leftmost two). For a general combination sequence like
`{{1,2}, {1,3}, {1,4}}`, the rescaling needs to track which propagators get
combined at each step and derive the corresponding x' → x mapping.

For a general combination sequence, `D_{1...n} = x'_1*D_1 + ... + x'_n*D_n`
where x'_j can be read off from the recursive combination. The simplest approach:
just compute `D_combined` symbolically with Feynman parameters, then read off the
coefficients of each original D_j to get x'_j.

---

## Current Bug: DiffExp Configuration

When `WorkingPrecision -> 200` is set, DiffExp requires `ChopPrecision` to be
smaller. The TransportLevel function needs to set:
```mathematica
ChopPrecision -> precision - 50  (* or similar *)
```

---

## Files and Test Results

### Tests (all in Tests/)
- `test_feynmantrick_algebra.m`: **21/21 PASS** — PropagatorAlgebra, EpsPrefactors
- `test_feynmantrick_fire.m`: **15/15 PASS** — FIREInterface (DefineTopology, SetupFIRE, FindBasis, ReduceIntegrals, ComputeDiffMatrix)
- `test_feynmantrick_iteration.m`: **17/17 PASS** — FeynmanTrickIteration (DefineFTIteration, BuildLevel, FeynmanTrickDecomposition, ExportLevel)
- `test_feynmantrick_endtoend.m`: Partially passes (boundary/transport sections have bugs)
- `test_feynmantrick_pipeline.m`: **NEW** — Full end-to-end pipeline (fails at boundary computation due to above bug)

### Pipeline test output (what works)
```
Part 1: Define Topology — PASS (4 propagators, 3 levels)
Part 2: Build & Compute — PASS (Level 3: 1 master, Level 2: 4 masters, Level 1: 7 masters)
         Export — PASS (all matrix files created)
Part 3: Boundary — FAIL (UF returns garbage, see bug above)
Part 4-6: Transport/Integration — not reached
```

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

1. **Fix `DeepestLevelBoundary`**: Compute U, F from original topology, apply rescaled
   Feynman parameters, then evaluate tadpole formula. Also fix UF loading/context.

2. **Fix TransportLevel config**: Add `ChopPrecision -> precision - 50` to DiffExp config.

3. **Run pipeline test**: After fixes, `test_feynmantrick_pipeline.m` should work end-to-end.

4. **Validate against known result**: The 1-loop massless box at s=-1, t=-1/3 has
   known analytic values. Compare pipeline output.

5. **(Future) Physical region**: Add i*delta prescriptions for threshold crossings.

6. **(Future) More examples**: Double box (topo7), double pentagon (5p from paper).

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
