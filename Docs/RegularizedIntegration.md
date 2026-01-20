# Regularized Integration

This document describes the `RegularizedIntegration` subpackage for computing definite integrals of DiffExp series expansions using the regularization prescription for Feynman parameter integrals.

## Overview

When integrating series solutions from DiffExp, we may encounter non-integrable singularities at integration boundaries. Near a singularity at x=0, the series decomposes as:

```
f(x, ε) = Σᵢ x^(aᵢ + bᵢε) gᵢ(x, ε)
```

where `gᵢ(x, ε)` is finite (starts at x⁰ or higher). Terms with `aᵢ ≤ -1` are non-integrable at x=0 unless regularized.

## Regularization Formula

For integrals of the form:

```
∫₀ᶜ x^(a+bε) g(x) dx
```

when `a ≤ -1` and `b ≠ 0`, we apply the identity:

```
∫₀ᶜ x^(a+bε) g(x) dx = ∫₀ᶜ [x^(a+1+bε)/(1+a+bε)] × [(2+a+bε)/c × g(x) - (1-x/c) × g'(x)] dx
```

This increases the power by 1. Repeated application resolves non-integrable singularities.

### ε-Pole Generation

When `a = -1` and `b ≠ 0`, the prefactor becomes:

```
1/(1 + a + bε) = 1/(bε) = (1/b) × ε⁻¹
```

This generates an ε⁻¹ pole, which is tracked via an `epsMinPower` offset in the implementation.

## Boundary Limits

Per the regularization prescription, when evaluating limits at a singularity (x → 0):
- **Terms with b ≠ 0**: contribute 0 (they vanish under analytic regularization)
- **Taylor terms (b = 0, a ≥ 0)**: evaluate g(x) at x = 0

## Main Functions

### DefiniteIntegral

```mathematica
DefiniteIntegral[savedData, {a, b}]
```

Computes ∫ f(x) dx from `a` to `b` where x is the main line parameter.

**Arguments:**
- `savedData`: Output from `TransportTo[..., save=True]`
- `{a, b}`: Integration bounds in main line x coordinates

**Returns:** List indexed by `{integralIndex, epsilonOrder}`

**Example:**
```mathematica
(* Transport with saving enabled *)
LoadConfiguration[{
  MatrixDirectory -> "path/to/matrices/",
  UseMobius -> False,  (* Required for integration! *)
  ...
}];

savedData = TransportTo[BCs, <|z -> 1/2|>, 1, True];

(* Compute ∫ f(x) dx from x=0 to x=0.5 *)
result = DefiniteIntegral[savedData, {0, 0.5}];

(* If you need ∫ f(z) dz where main line is z -> z0 + slope*x: *)
(* Multiply by dz/dx = slope *)
resultInZ = slope * result[[integralIdx, epsOrder]];
```

### IndefiniteIntegral

```mathematica
IndefiniteIntegral[savedData]
```

Returns piecewise functions representing the indefinite integral ∫ f(x) dx. The integration constant is fixed such that the integral is 0 at the start of the first segment.

**Returns:** List of piecewise functions indexed by `{integralIndex, epsilonOrder}`

### IntegrateDecomposition

```mathematica
IntegrateDecomposition[decomposition, {xmin, xmax}]
```

Integrates a decomposed series (from `DecomposeSingularity`) over the given interval, applying regularization as needed.

### EvaluateLimitAtSingularity

```mathematica
EvaluateLimitAtSingularity[decomposition, direction]
```

Evaluates the limit of a decomposed series at x=0.
- `direction`: +1 for x→0⁺ or -1 for x→0⁻

### ApplyRegularizationStep

```mathematica
ApplyRegularizationStep[a, b, epsMinPower, gList, c]
```

Applies one step of the regularization formula, increasing the power from `a` to `a+1`.

**Returns:** `{a+1, b, newEpsMinPower, newGList}`

### RegularizeIntegrand

```mathematica
RegularizeIntegrand[a, b, gList, c]
```

Repeatedly applies the regularization formula until `a ≥ 0`.

**Returns:** `{newA, b, epsMinPower, newGList}`

## Requirements

1. **Use `UseMobius -> False`** when transporting. Mobius transforms create nonlinear coordinate transformations that complicate the integration. Linear segment transforms are required.

2. **Transport with `save=True`** to generate the savedData structure needed for integration.

## Coordinate System

All integrals are computed in the main line parameter x:

```
∫ f(x) dx  from x = a to x = b
```

The kinematic invariants are functions of x along the main line. For example, if the main line is `z -> z0 + slope*x`, then:
- At x=0: z = z0
- At x=1: z = z0 + slope
- dz/dx = slope

To convert from an x-integral to a z-integral:
```
∫ f(z) dz = (dz/dx) × ∫ f(x) dx
```

## Connection to Singularity Decomposition

The integration functions use `DecomposeSingularity` (see `Docs/SingularityDecomposition.md`) to extract the singular structure of series near boundaries. Each term `x^(a + bε) g(x,ε)` is then integrated using:
- Direct integration if `a ≥ 0`
- Regularization if `a < 0` and `b ≠ 0`
- Error if `a < 0` and `b = 0` (truly non-integrable)

## References

The regularization formula is similar to the pole extraction in sector decomposition. See:
- Heinrich, G. "Sector Decomposition" (2008)
- Panzer, E. "Algorithms for the symbolic integration of hyperlogarithms" (2014)
