# Singularity Decomposition

Decomposes `IntegrateSystem` output near singularities into canonical form.

## Overview

When evaluating Feynman integrals near a singular point (e.g., a threshold), the result from `IntegrateSystem` contains series with fractional and negative powers of `x`, along with `Logx` terms at various epsilon orders. This module decomposes such output into a canonical form that makes the singular structure explicit:

```
f(x, eps) = sum_i x^(a_i + b_i * eps) * g_i(x, eps)
```

where:
- `a_i` is a rational number (the singular power, can be negative like -3/2, -1, -1/2, etc.)
- `b_i` is a rational number (the epsilon-dependent exponent)
- `g_i(x, eps)` is a series that is **finite at x=0** (starts at x^0 or higher)

## Functions

### DecomposeSingularity

```mathematica
DecomposeSingularity[seriesList]
```

Decomposes a single integral's output (a list of series, one per epsilon order) into canonical form.

**Input:** A list `{f_eps0, f_eps1, f_eps2, ...}` where each `f_epsN` is a `SeriesData` in `x` (or a number).

**Output:** A list of terms `{term1, term2, ...}` where each term is an Association:
```mathematica
<|"a" -> -3/2, "b" -> -3, "g" -> {g0, g1, g2, ...}|>
```
- `"a"`: The leading power (rational)
- `"b"`: The epsilon-dependent exponent (rational)
- `"g"`: List of series `{g_eps0, g_eps1, ...}` that are finite at x=0

### DecomposeSingularityAll

```mathematica
DecomposeSingularityAll[integrateSystemOutput]
```

Applies `DecomposeSingularity` to all integrals in the output of `IntegrateSystem`.

**Input:** The full output from `IntegrateSystem`, a list of lists.

**Output:** A list of decompositions, one per integral.

### PrintDecomposition

```mathematica
PrintDecomposition[decomp]
```

Prints a decomposition in human-readable form.

## Algorithm

The algorithm determines `b` from the `Logx` structure. The key insight is that:

```
x^(b*eps) = exp(b*eps*log(x)) = 1 + b*eps*Logx + (b*eps)^2/2 * Logx^2 + ...
```

So the coefficient of `Logx` at epsilon order n+1 is related to the coefficient at order n by the factor `b`.

**Steps:**
1. Find the leading power `a` (most negative power across all epsilon orders)
2. Find the first non-zero epsilon order and use the ratio of Logx coefficients to determine `b`
3. Factor out `x^(a + b*eps)` to obtain `g(x, eps)`
4. The result `g` should be finite (start at x^0 or higher)

## Configuration

The tolerance for rationalizing the exponents `a` and `b` is controlled by:

```mathematica
UpdateConfiguration[RationalizationTolerance -> 10^-10]  (* default *)
```

This affects:
- When a coefficient is considered "effectively zero"
- The precision for rationalizing `b` to an exact rational

The coefficients of `g` retain their full numerical precision from the original computation.

## Example Usage

```mathematica
(* Load package and configure *)
Get["DiffExp.m"];
LoadConfiguration[{
  DeltaPrescriptions -> {t - 16 + I * delta},
  MatrixDirectory -> "path/to/matrices/",
  EpsilonOrder -> 4
}];

(* Transport to near singularity and integrate *)
PreparedBCs = PrepareBoundaryConditions[bcs, <|t -> -1/x|>];
Results = TransportTo[PreparedBCs, <|t -> 15|>];
intResult = IntegrateSystem[Results, <|t -> 16 + x|>];  (* x=0 is singularity *)

(* Decompose all integrals *)
decompositions = DecomposeSingularityAll[intResult];

(* Print results *)
Do[
  Print["Integral ", i, ":"];
  PrintDecomposition[decompositions[[i]]];
, {i, Length[decompositions]}];

(* Access individual components *)
decomp1 = decompositions[[1]][[1]];  (* First term of first integral *)
a = decomp1["a"];     (* e.g., -3/2 *)
b = decomp1["b"];     (* e.g., -3 *)
g = decomp1["g"];     (* {g_eps0, g_eps1, ...} *)
```

## Example Output

For the equal mass banana integral near t=16:

```
Integral 1: x^(-3/2 - 3*eps) * g(x, eps)
Integral 2: x^(-1/2 - 3*eps) * g(x, eps)
Integral 3: x^0 * g(x, eps)
Integral 4: x^0 * g(x, eps)
```

The `b = -3` corresponds to the dimensional regularization structure of the singularity.
