# Singularity Decomposition Problem

## Goal

Given the output of `IntegrateSystem` near a singular point, decompose it into a canonical form that makes the singular structure explicit.

## Input

The output of `IntegrateSystem[bcs, line]` where `line` is centered at a singularity (e.g., `t = 16 + x` so `x = 0` is the singular point).

This gives, for each integral, a list of `SeriesData` objects - one per epsilon order:
```
{series_at_eps^0, series_at_eps^1, series_at_eps^2, ...}
```

Each `series_at_eps^n` is a `SeriesData` in `x` around `x = 0`, with:
- Possibly fractional powers (e.g., `x^(-3/2)`, `x^(-1/2)`, `x^(1/2)`, ...)
- Coefficients that may contain `Logx` terms

## Desired Output

Decompose each integral's solution into a sum of terms:

```
f(x, ε) = Σᵢ x^(aᵢ + bᵢ·ε) · gᵢ(x, ε)
```

where:
- `aᵢ ∈ ℚ` (rational) - captures the singular power (can be negative like -3/2, -1, -1/2, etc.)
- `bᵢ ∈ ℚ` (rational) - the epsilon-dependent part of the exponent, determined from Logx structure
- `gᵢ(x, ε)` is a series in `x` that is **always finite at x=0**, meaning:
  - `gᵢ` starts at `x^0` or higher powers (no negative powers)
  - `gᵢ` may have fractional powers like `x^(1/2)` if the original series does
  - `gᵢ` has an epsilon expansion: `gᵢ(x, ε) = Σₙ gᵢₙ(x) · εⁿ`

## Key Insight: How b Appears

When we expand `x^(bε)`:
```
x^(bε) = exp(bε·log(x)) = 1 + bε·Logx + (bε)²/2·Logx² + (bε)³/6·Logx³ + ...
```

So the Logx structure at different epsilon orders reveals `b`:
- At ε⁰: no Logx contribution from `x^(bε)`
- At ε¹: coefficient of Logx¹ gets contribution `b · (ε⁰ term)`
- At ε²: coefficient of Logx² gets contribution `b²/2 · (ε⁰ term)`
- etc.

## Algorithm (High-Level)

1. **Identify leading power `a`**: Find the most negative power of `x` across all epsilon orders
2. **Determine `b`**: Use the Logx structure to fit `b` (ratio of coefficients)
3. **Extract `g(x, ε)`**: Factor out `x^(a + bε)` to get the finite part `g`
4. **Compute residual**: Subtract `x^(a + bε) · g(x, ε)` from the original
5. **Repeat**: If residual is non-zero, go back to step 1 with the residual
6. **Terminate**: When residual is effectively zero

## SeriesData Structure Reminder

Mathematica's `SeriesData[var, center, coeffs, nmin, nmax, den]`:
- `coeffs` = list of coefficients
- `nmin` = numerator of minimum power
- `nmax` = numerator of (max power + 1)
- `den` = denominator for fractional powers
- Actual min power = `nmin / den`
- Actual max power = `(nmax - 1) / den`
- Coefficient of `x^(nmin/den + k/den)` is `coeffs[[k + 1]]`

Example: `SeriesData[x, 0, {a, b, c}, -3, 3, 2]` represents:
```
a·x^(-3/2) + b·x^(-1) + c·x^(-1/2) + O(x^1)
```

## Test Setup

Use the equal mass banana example:
1. Load configuration with singularity at `t = 16`
2. Transport boundary conditions to `t = 15` (near singularity)
3. Run `IntegrateSystem[Results15, <|t -> 16 + x|>]` to get expansion around singularity
4. Apply decomposition to the result

## Current Status

- `IntegrateSystem` works correctly and returns series with fractional powers
- Transport to `t = 15` and `IntegrateSystem` around `t = 16` works
- **TODO**: Implement the fitting function `DecomposeSingularity` that performs the decomposition
