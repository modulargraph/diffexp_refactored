# Series Operations, Integration, and Evaluation

This document covers three closely related DiffExp subpackages that form the core computational pipeline for manipulating, integrating, and evaluating power series solutions of differential equations:

1. **SeriesOps.m** (`DiffExp`SeriesOps``) -- Series manipulation and arithmetic
2. **Integration.m** (`DiffExp`Integration``) -- Term-by-term integration of series
3. **Pade.m** (`DiffExp`Pade``) -- Pade approximants and series evaluation

Together, these packages handle the lifecycle of a series from construction through differentiation/integration to final numerical evaluation at a target point.

---

## Mathematical Context

DiffExp solves systems of linear ordinary differential equations of the form:

```
d/dx F(x, eps) = M(x, eps) . F(x, eps)
```

where `F` is a vector of master integrals, `M` is a matrix of rational functions, `x` is the line parameter, and `eps` is the dimensional regulator. The solutions are represented as truncated power series in `x` around regular singular points, with coefficients that are polynomials in a symbolic `Logx` (representing `Log[x]`). This symbolic treatment of logarithms avoids branch-cut issues during intermediate algebraic manipulations.

A typical series element has the structure:

```
SeriesData[x, 0, {c0 + c1*Logx + c2*Logx^2, ...}, nmin, nmax, den]
```

representing a sum of terms `x^(k/den) * (polynomial in Logx)` for `k` from `nmin` to `nmax-1`.

---

## Package 1: SeriesOps.m (`DiffExp`SeriesOps``)

**File:** `/Users/mhidding/Desktop/diffexp_refactored/DiffExp/SeriesOps.m` (178 lines)

**Dependencies:** `DiffExp`Symbols``, `DiffExp`State``, `DiffExp`Utilities``

**Purpose:** Provides the fundamental operations for manipulating Mathematica `SeriesData` objects that arise in the DiffExp computation. All functions are designed to handle edge cases (zero, numeric, non-series inputs) gracefully, and most carry the `Listable` attribute for automatic threading over lists.

### Attributes

The following functions have `Listable` set, meaning they automatically map over list arguments:

```mathematica
SetAttributes[{SApply, SN, SExpand, SN, SSN, SMultiply, SEval, DecreaseSeriesOrderBy,
               SeriesCoefficientMinus, SplitTimes, ApplyAnalyticContinuation}, Listable];
SetAttributes[DiffExpSeries, Listable];
```

---

### Core Series Manipulation

#### SApply

```mathematica
SApply[f_, a_]
```

Applies function `f` to the coefficient list of a `SeriesData` object. This operates on position `[[3]]` of the internal `SeriesData` structure (the list of coefficients) via `MapAt`.

**Behavior by input type:**
- `SApply[f, 0]` returns `0`
- `SApply[f, a_SeriesData]` returns `MapAt[f, a, 3]` (applies `f` to each coefficient)
- `SApply[f, a_]` returns `f[a]` (fallback for non-series)

**Example:**
```mathematica
(* Expand all coefficients of a series *)
SApply[Expand, mySeries]

(* Apply numerical evaluation to coefficients *)
SApply[N, mySeries]
```

---

#### SExpand

```mathematica
SExpand[a_]
```

Expands a series expression and chops small numerical residuals. This is the workhorse simplification function used throughout DiffExp to keep series coefficients clean.

**Implementation:**
- For `SeriesData`: applies `PChop @* Expand` to each coefficient via `SApply`
- For non-series: applies `Expand` directly
- For `0`: returns `0`

`PChop` chops numbers below `10^(-ChopPrecision)`, where `ChopPrecision` defaults to 250. This prevents accumulation of numerical noise in high-precision computations.

**Example:**
```mathematica
(* Clean up a series after arithmetic *)
result = series1 * series2 // SExpand
```

---

#### SN / SSN

```mathematica
SN[a_]
SSN[a_]
```

Apply `N` (numerical evaluation) to series coefficients. Both are defined as `SApply[N, #] &`. These are synonyms provided for convenience.

---

#### SMultiply

```mathematica
SMultiply[a_, b_]
```

Multiplies two series with full expansion and chopping at each step:

```mathematica
SMultiply[a_, b_] := SExpand[a] SExpand[b] // SExpand
```

This ensures both operands are cleaned before multiplication, and the result is cleaned after, preventing coefficient bloat.

---

#### SeriesCoefficientMinus

```mathematica
SeriesCoefficientMinus[a_, k_:1]
```

Extracts the `k`-th coefficient from the **end** (highest power) of a series, multiplied by the appropriate power of `x`. Useful for error estimation based on trailing coefficients.

**Implementation:**
```mathematica
a[[3]][[Max[-k, -Length[a[[3]]]]]] * x^((a[[5]] - k)/a[[6]])
```

For non-`SeriesData` input, returns `0`.

**Parameters:**
- `a` -- A `SeriesData` object
- `k` -- Position from end (default: 1, i.e., last coefficient)

---

### Analytic Continuation and Replacement

#### ApplyAnalyticContinuation

```mathematica
ApplyAnalyticContinuation[s_SeriesData]
```

Applies the current `AnalyticContinuationReplacements` to a series. This is used when crossing branch cuts, where expressions like `Sqrt[f(x)]` need to pick up sign changes or phase factors.

**Implementation:**
1. Converts series to normal expression via `Normal[s]`
2. Applies `AnalyticContinuationReplacements` (stored in `DiffExp`State``)
3. Re-expands back into a series around `x=0`

The truncation order is preserved as `Floor[s[[5]]/s[[6]]]`.

---

#### SafeReplaceSeries11

```mathematica
SafeReplaceSeries11[a_, b_]
```

Performs replacement rules `b` on a series while preserving the series structure and truncation order.

**Behavior by input type:**
- For `SeriesData`: converts to `Normal`, applies replacement, re-expands, then truncates to original order using `O[x]^(a[[5]]/a[[6]])`
- For `List`: maps over elements recursively
- For other: applies replacement directly via `/.`

This is "safe" because naive replacement on `SeriesData` internals can corrupt the series structure.

---

### Logarithm Handling

A key design principle in DiffExp is that `Log[x]` is represented symbolically as `Logx` within series coefficients. This avoids issues with branch cuts and allows algebraic manipulation of logarithmic terms. The following functions manage this representation.

#### MaxLogxPower

```mathematica
MaxLogxPower[ex_]
```

Returns the maximum power of `Logx` appearing in expression `ex`. Searches for patterns of the form `Logx^k` (with `k` defaulting to 1 for bare `Logx`). Returns 0 if no `Logx` is present.

**Implementation:**
```mathematica
Append[GetCases[ex // SExpand, Logx^(k_:1) :> k], 0] // Max
```

---

#### LogxCoeff

```mathematica
LogxCoeff[Ser_, Which_]
```

Extracts the coefficient of `Logx^Which` from a `SeriesData` object, operating on each series coefficient via `SApply`.

**Behavior:**
- If `Which === 0`: replaces `Logx -> 0` in each coefficient (extracts the non-logarithmic part)
- Otherwise: extracts `Coefficient[#, Logx^Which]` from each coefficient

**Example:**
```mathematica
(* Extract the Log[x]^2 coefficient from a series *)
logSquaredPart = LogxCoeff[mySeries, 2]

(* Extract the purely rational part (no logs) *)
rationalPart = LogxCoeff[mySeries, 0]
```

---

#### LogxCoeffNS

```mathematica
LogxCoeffNS[Ser_, Which_]
```

Same as `LogxCoeff` but for non-series (plain) expressions. Does not use `SApply`.

---

#### LogxCoeffList

```mathematica
LogxCoeffList[Ser_]
```

Returns a list of all `Logx` coefficients from power 0 up to `MaxLogxPower[Ser]`:

```mathematica
Table[LogxCoeff[Ser, ord], {ord, 0, MaxLogxPower[Ser]}]
```

The resulting list has length `MaxLogxPower[Ser] + 1`, with element `i` being the coefficient of `Logx^(i-1)`.

---

### Matrix Operations

#### MatrixMultiplySExpand

```mathematica
MatrixMultiplySExpand[MatA_, MatB_]
```

Matrix multiplication where each element-wise product uses `SMultiply` and summation uses `SExpand`. This ensures proper handling of series truncation and numerical noise during matrix operations.

**Implementation:**
1. Checks dimension compatibility (`Dim1[[2]] === Dim2[[1]]`)
2. Computes `(A.B)_{ij} = SExpand[Sum[SMultiply[A_{ik}, B_{kj}], {k}]]`
3. On dimension mismatch, sets `LastErrorContext` and calls `ReportError`

**Note:** A parallel computation version exists in the original code but is commented out in the refactored version.

---

#### MatrixPowerSExpand

```mathematica
MatrixPowerSExpand[a_, n_]
```

Raises a square matrix to the `n`-th power using repeated `MatrixMultiplySExpand`.

**Behavior:**
- `MatrixPowerSExpand[a, 0]` returns `IdentityMatrix[dim]` (after checking squareness)
- `MatrixPowerSExpand[a, n]` uses `Nest[MatrixMultiplySExpand[a, #]&, a, n-1]`

---

### Series Construction

#### DiffExpSeries

```mathematica
DiffExpSeries[Ser_, ord_]
DiffExpSeries[Ser_]
```

Creates a DiffExp series at the configured working precision. Converts `Ser` to numerical precision `FEWorkingPrecision` and wraps in `SeriesAlways`.

- The one-argument form uses `ExpansionOrderVal` as the default order.
- The two-argument form uses the specified `ord`.

---

#### SeriesAlways

```mathematica
SeriesAlways[term_, {a_, b_, c_}, ex_:1]
```

Creates a `SeriesData` object even for constant terms (which Mathematica's `Series` would not wrap in `SeriesData`).

**Behavior:**
- If `term` depends on `a` (the expansion variable): calls `Series[term, {a, b, c}]` with assumption `x > 0`
- If `term` is constant: manually constructs `SeriesData[a, b, {term}, 0, c*ex + 1, ex]`

This ensures a uniform series representation throughout the computation pipeline.

---

#### LeadingCoefficientSeries

```mathematica
LeadingCoefficientSeries[Ser_, AddTo2_:1]
```

Extracts the leading-order term of a series as a new series with minimal truncation. This is used to isolate the dominant behavior near a singular point.

**Implementation details:**
- For x-independent series: returns `Ser + O[x]^(1/AddTo2)`
- For x-dependent series: replaces `x -> x + O[x]^2`, takes the first coefficient, then constructs a new series with the combined fractional power structure
- The `AddTo2` parameter controls the additional fractional power structure (used for matching denominators in fractional-power series)
- Computes the combined denominator from the series' own denominator and `AddTo2` using `FactorInteger` and `Merge[..., Max]`

---

### Series Truncation Information

#### SeriesMinPower

```mathematica
SeriesMinPower[Ser_]
```

Returns the minimum power (starting power) of a `SeriesData` object: `Ser[[4]] * Ser[[6]]`.

Note: In `SeriesData[x, 0, coeffs, nmin, nmax, den]`, the actual minimum power is `nmin/den`, but this function returns `nmin * den` which appears to be a specific convention used in DiffExp's internal logic.

---

#### SeriesMaxPower

```mathematica
SeriesMaxPower[Ser_]
```

Returns the maximum power (truncation order) of a `SeriesData` object: `Ser[[5]] * Ser[[6]]`.

---

#### DecreaseSeriesOrderBy

```mathematica
DecreaseSeriesOrderBy[a_, k_:1]
```

Decreases the truncation order of a series by `k` units (in terms of the series' internal indexing, i.e., `k * den` in actual power).

**Implementation:**
```mathematica
tmp[[5]] = tmp[[5]] - k * tmp[[6]]
```

For non-`SeriesData` input, returns the input unchanged.

---

### Series Differentiation

#### SD

```mathematica
SD[a_, b_]
SD[a_, b__]
SD[a_]
```

Takes the derivative of a series with respect to variable `b`, correctly handling the symbolic `Logx` representation. This is the key function that prevents the Mathematica kernel from introducing explicit `Log[x]` terms during differentiation.

**Algorithm:**
1. Determines `CurrMaxLogPower = MaxLogxPower[a]`
2. For each power `logxord` from 0 to `CurrMaxLogPower`:
   - Extracts the coefficient of `Logx^logxord` via `LogxCoeff`
   - Differentiates the coefficient normally: `D[Tmp, b]`
   - Computes the derivative of `Log[x]^logxord` symbolically, then replaces `Log[x] -> Logx`
   - Combines: `D[coeff, b] * Logx^logxord + coeff * D[Log[x]^logxord, b]|_{Log[x]->Logx}`
3. Sums over all log powers and applies `SExpand`

**Multi-variable form:** `SD[a, b, c, ...]` differentiates sequentially, first with respect to the first variable, then the next, etc.

**Identity form:** `SD[a]` with no derivative variable returns `a` unchanged.

**Mathematical insight:** If `f(x) = g(x) + h(x) * Logx + k(x) * Logx^2`, then:
```
d/dx f = g'(x) + h'(x)*Logx + h(x)/x + k'(x)*Logx^2 + 2*k(x)*Logx/x
```

The `SD` function performs exactly this decomposition, keeping `Logx` symbolic throughout.

---

## Package 2: Integration.m (`DiffExp`Integration``)

**File:** `/Users/mhidding/Desktop/diffexp_refactored/DiffExp/Integration.m` (73 lines)

**Dependencies:** `DiffExp`Symbols``, `DiffExp`State``, `DiffExp`Utilities``, `DiffExp`SeriesOps``

**Purpose:** Provides term-by-term integration of series expressions, handling logarithmic terms via pattern-matched replacement rules. This is the inverse operation to `SD` and is used in the variation of parameters method for solving inhomogeneous ODEs.

### Architecture

The integration works by:
1. Converting the series to a normal expression
2. Replacing `Logx -> Log[x]` (switching from symbolic to standard form)
3. Applying pre-computed replacement rules (`IntReps`) that map each monomial term to its integral
4. Subtracting the integration constant (value at `x=0`)
5. Re-wrapping the result as a series

This approach avoids calling Mathematica's `Integrate` on the full expression, instead using a lookup table of pre-computed integrals for each term type.

---

### Functions

#### DiffExpIntegrate

```mathematica
DiffExpIntegrate[a__]
```

Main entry point for integration. Before delegating to `DiffExpIntegrate1`, it checks whether the maximum `Logx` power in the input has increased beyond the current `IMaxLogOrder`. If so, it calls `UpdateIntReps` to regenerate the replacement rules for the higher log order.

**Implementation:**
```mathematica
DiffExpIntegrate[a__] := Block[{LogOrd},
  LogOrd = Max[Append[GetCases[{a} // SExpand, Logx^(k_:1) | Log[x]^(k_:1) :> k], 1]];
  If[LogOrd > IMaxLogOrder,
    PrintDebug["Encountered ", Log[x]^LogOrd, ". Updating IntReps."][2];
    UpdateIntReps[LogOrd];
  ];
  DiffExpIntegrate1[a]
]
```

**Parameters:**
- `a__` -- One or more series expressions to integrate (passed through to `DiffExpIntegrate1`)

---

#### DiffExpIntegrate1

```mathematica
DiffExpIntegrate1[a_]
DiffExpIntegrate1[exp0_SeriesData, var_]
DiffExpIntegrate1[exp0_ /; NumericQ[exp0], var_]
```

Internal integration function (has `Listable` attribute). The single-argument form defaults the integration variable to `x`.

**For `SeriesData` input:**
1. Expands and converts `Logx -> Log[x]`
2. Replaces `x -> b` (a placeholder variable used in the integration rules)
3. Applies `IntReps` replacement rules to compute the antiderivative
4. Subtracts the constant of integration (the value at `a=0`, where `a` is an auxiliary variable in the rules)
5. Adds back `b * Const` (the linear-in-x part from constants)
6. Substitutes back `a->1, b->x`, converts `Log[x] -> Logx`
7. Wraps in `SeriesAlways` with truncation order `Floor[nmax/den]`

**For numeric input:**
Returns `x * exp0` wrapped in `SeriesAlways` at `ExpansionOrderVal`.

**For unsupported input:**
Calls `ReportError` with a descriptive message.

**Example conceptual flow:**
```
Input:  x^2 + 3*x*Logx
Step 1: x^2 + 3*x*Log[x]           (Logx -> Log[x])
Step 2: b^2 + 3*b*Log[b]            (x -> b)
Step 3: a*b^3/3 + 3*a*(b^2*Log[b]/2 - b^2/4)  (apply IntReps)
Step 4: subtract value at a=0, add constant term
Step 5: substitute a->1, b->x, Log[x]->Logx
Result: x^3/3 + 3*(x^2*Logx/2 - x^2/4) as SeriesData
```

---

#### UpdateIntReps

```mathematica
UpdateIntReps[MaxOrd_]
```

Generates the replacement rules for integrating all monomial types up to `Log[x]^MaxOrd`. Sets `IMaxLogOrder = MaxOrd` and builds the `IntReps` list.

**Rule categories generated (for each `n` from 1 to `MaxOrd`):**

| Pattern | Integral |
|---------|----------|
| `Log[x]^n` | `a * Integrate[Log[x]^n, x]` evaluated at `x->b` |
| `Log[x]^n * x` | `a * Integrate[Log[x]^n * x, x]` evaluated at `x->b` |
| `Log[x]^n * x^m` (m != -1) | `a * Integrate[Log[x]^n * x^m, x]` evaluated at `x->b` |
| `Log[x]^n / x` | `a * Integrate[Log[x]^n / x, x]` evaluated at `x->b` |

**Additional rules (independent of log order):**

| Pattern | Integral |
|---------|----------|
| `x^m` (m != -1) | `a * x^(m+1)/(m+1)` evaluated at `x->b` |
| `x` | `a * x^2/2` evaluated at `x->b` |
| `1/x` | `a * Log[x]` evaluated at `x->b` |

The auxiliary variables `a` and `b` are used to implement definite integration: `a` acts as a switch (set to 0 for the lower limit, 1 for the upper), and `b` is the evaluation point.

The rules are stored in **reversed** order and expanded, so more specific patterns (like `Log[x]^n * x^m`) are tried before less specific ones (like `x^m`).

---

#### IntReps

```mathematica
IntReps
```

The current list of integration replacement rules. Initialized to `{}` and populated by `UpdateIntReps`. This is a module-level variable that persists between calls (caching the rules as long as the log order does not increase).

---

## Package 3: Pade.m (`DiffExp`Pade``)

**File:** `/Users/mhidding/Desktop/diffexp_refactored/DiffExp/Pade.m` (74 lines)

**Dependencies:** `DiffExp`Symbols``, `DiffExp`State``, `DiffExp`Utilities``, `DiffExp`SeriesOps``

**Purpose:** Provides Pade approximant computation and series evaluation functions. Pade approximants are rational function approximations that typically converge better than truncated power series, especially near singularities. This package also handles the final step of evaluating a series at a specific numerical point, including analytic continuation.

### Mathematical Background

A Pade approximant `[M/N]` of a power series `f(x) = sum c_k x^k` is a rational function `P(x)/Q(x)` where `P` has degree `M`, `Q` has degree `N`, and the Taylor expansion of `P/Q` matches `f` to order `M+N`. For a series of order `K`, DiffExp uses the "diagonal" or near-diagonal approximant with:

```
M = N = Floor[(K+1)/2]
```

This symmetric choice often provides the best convergence properties.

---

### Functions

#### GetPade

```mathematica
GetPade[a_]
```

Computes the Pade approximant of a series, handling the `Logx` decomposition.

**Behavior by input type:**
- `GetPade[0]` returns `0`
- `GetPade[a_?NumericQ]` returns `a` unchanged
- `GetPade[a_SeriesData]` computes the Pade approximant as described below

**Algorithm for `SeriesData`:**
1. Sets `$MinPrecision = FEWorkingPrecision` to maintain numerical accuracy
2. Computes `maxPadeOrder = Floor[((nmax - nmin)/den + 1) / 2]` for both numerator and denominator
3. Determines `MaxLogOrder = MaxLogxPower[a]`
4. For each `Logx` power from 0 to `MaxLogOrder`:
   - Extracts the coefficient series via `LogxCoeff[a, ind]`
   - Converts to `Normal` form
   - Applies `Chop[..., 10^(-ChopPrecisionVal)]` to remove numerical noise
   - Calls `PadeApproximant[..., {x, 0, {M, N}}]`
   - If `PadeApproximant` fails (returns unevaluated), falls back to `PChop` of the normal form with a warning
5. Recombines: `Sum[Logx^ind * PadeResult_ind, {ind, 0, MaxLogOrder}]`

**Example:**
```mathematica
(* Compute Pade approximant of a series solution *)
pade = GetPade[mySeries]
(* Result is a rational function in x, possibly with Logx terms *)
```

**Failure handling:** If the Pade computation fails (e.g., due to insufficient non-zero coefficients), a warning is printed and the function returns the `Normal` form instead.

---

#### SEval1

```mathematica
SEval1[a_]
```

Returns the evaluable form of a series -- either the Pade approximant or the normal (polynomial) form, depending on the `UsePade` configuration option.

**Behavior:**
- `SEval1[a_SeriesData]`: returns `GetPade[a]` if `UsePade === True`, else `Normal[a]`
- `SEval1[0]` returns `0`
- `SEval1[a_?NumericQ]` returns `a`

Has `Listable` attribute set.

---

#### SEval2

```mathematica
SEval2[a_, at_]
```

Evaluates an expression (typically the output of `SEval1`) at a specific numerical point `at`, applying analytic continuation replacements.

**Algorithm:**
1. Sets `$MinPrecision = FEWorkingPrecision`
2. Applies `AnalyticContinuationReplacements` with theta-function resolution:
   - If `at >= 0`: sets `theta_p -> 1, theta_m -> 0`
   - If `at < 0`: sets `theta_p -> 0, theta_m -> 1`
3. Replaces `Logx -> Log[x]` (converts symbolic log to standard)
4. Replaces `x -> at` (substitutes the evaluation point)
5. Sets precision to `FEWorkingPrecision`
6. Applies `Expand`

The theta functions (`theta_p`, `theta_m`) implement the sign-dependent branch selection for analytic continuation across branch cuts. When approaching a singularity from the positive side, `theta_p = 1`; from the negative side, `theta_m = 1`.

Has `Listable` attribute set.

**Example:**
```mathematica
(* Evaluate series at x = 0.5 *)
expr = SEval1[mySeries];
result = SEval2[expr, 0.5]
```

---

#### SEval

```mathematica
SEval[a_, at_]
```

Combined evaluation function that performs both `SEval1` and `SEval2` in a single call. This is the primary function for evaluating a series at a point.

**Behavior by input type:**
- `SEval[a_SeriesData, at_]`: applies the full pipeline (Pade/Normal + analytic continuation + substitution)
- `SEval[0, at_]` returns `0`
- `SEval[a_?NumericQ, at_]` returns `a`

**Implementation (for `SeriesData`):**
```mathematica
SEval[a_SeriesData, at_] := Block[{$MinPrecision = FEWorkingPrecision},
  If[UsePade === True, GetPade[a], Normal[a]] /.
    (AnalyticContinuationReplacements /.
      If[at >= 0,
        {theta_p -> 1, theta_m -> 0},
        {theta_p -> 0, theta_m -> 1}
      ]) /. Logx -> Log[x] /. x -> at //
    SetPrecision[#, FEWorkingPrecision]& // Expand
]
```

**Example:**
```mathematica
(* Direct evaluation of a series at a point *)
value = SEval[mySeries, 0.7]
```

---

## Interaction Between Packages

The three packages form a pipeline used during the transport of boundary conditions:

```
                SD (differentiate)
                      |
                      v
  SeriesOps  <-->  Integration  -->  Pade/Evaluation
  (manipulate)    (antiderivative)    (evaluate at endpoint)
```

### Typical Workflow

1. **Series Construction** (`DiffExpSeries`, `SeriesAlways`): Initial boundary conditions and Frobenius solutions are wrapped as series.

2. **Matrix-Series Arithmetic** (`MatrixMultiplySExpand`, `SMultiply`): The connection matrix `M(x)` is multiplied with the current solution vector, producing the right-hand side of the ODE.

3. **Differentiation** (`SD`): Used to verify solutions and compute derivatives needed in the variation of parameters method.

4. **Integration** (`DiffExpIntegrate`): Integrates the inhomogeneous terms in the variation of parameters formula. The `IntReps` rules handle all monomial types including logarithmic terms.

5. **Evaluation** (`SEval`, `GetPade`): At the end of each line segment, the series is evaluated at the endpoint to produce the initial conditions for the next segment.

### The Role of Logx

The symbolic `Logx` is central to the design:

- **SeriesOps** keeps `Logx` symbolic during all manipulations
- **SD** correctly differentiates `Logx^n` terms without introducing `Log[x]`
- **Integration** temporarily converts `Logx -> Log[x]` for pattern matching, then converts back
- **Pade/SEval** converts `Logx -> Log[x]` only at the final evaluation step

This separation ensures that logarithmic branch cuts are handled correctly throughout the computation.

---

## Configuration Dependencies

These packages depend on several configuration values from `DiffExp`State``:

| Variable | Default | Used By | Purpose |
|----------|---------|---------|---------|
| `FEWorkingPrecision` | 500 | All | Numerical precision for computations |
| `ChopPrecision` | 250 | SeriesOps, Pade | Threshold for chopping small numbers |
| `ExpansionOrderVal` | 50 | Integration, SeriesOps | Default series truncation order |
| `UsePade` | False | Pade | Whether to use Pade approximants |
| `IMaxLogOrder` | 1 | Integration | Current maximum log order for IntReps |
| `AnalyticContinuationReplacements` | {} | SeriesOps, Pade | Branch-cut handling rules |

---

## Error Handling

- **MatrixMultiplySExpand**: Stores `{MatA, MatB}` in `LastErrorContext` before aborting on dimension mismatch
- **MatrixPowerSExpand**: Aborts if matrix is not square
- **DiffExpIntegrate1**: Reports error for unsupported argument types
- **GetPade**: Prints warning and falls back to Normal form if Pade computation fails

---

## Usage Examples

### Basic Series Manipulation

```mathematica
(* Create a DiffExp series from an expression *)
ser = DiffExpSeries[1 + x + x^2/2 + x^3/6, 5];

(* Expand and clean *)
clean = SExpand[ser];

(* Extract log structure *)
maxLog = MaxLogxPower[ser];
logCoeffs = LogxCoeffList[ser];

(* Take derivative preserving Logx *)
deriv = SD[ser, x];
```

### Matrix Operations with Series

```mathematica
(* Multiply connection matrix with solution vector *)
matA = {{ser1, ser2}, {ser3, ser4}};
matB = {{ser5}, {ser6}};
product = MatrixMultiplySExpand[matA, matB];

(* Compute matrix power *)
matSquared = MatrixPowerSExpand[matA, 2];
```

### Integration

```mathematica
(* Integrate a series term-by-term *)
integrated = DiffExpIntegrate[mySeries];

(* The integration handles Log[x] terms automatically *)
(* e.g., integrating x^2 * Logx gives the correct antiderivative *)
```

### Pade Evaluation

```mathematica
(* Evaluate series at endpoint of line segment *)
endpointValue = SEval[mySeries, 0.8];

(* Or in two steps for reuse *)
evaluable = SEval1[mySeries];        (* Pade or Normal form *)
val1 = SEval2[evaluable, 0.8];       (* Evaluate at x=0.8 *)
val2 = SEval2[evaluable, -0.3];      (* Evaluate at x=-0.3, uses theta_m *)
```

### Checking Convergence

```mathematica
(* Compare Pade vs Normal evaluation *)
UpdateConfiguration[UsePade -> True];
padeVal = SEval[mySeries, 0.9];

UpdateConfiguration[UsePade -> False];
normalVal = SEval[mySeries, 0.9];

(* Difference indicates convergence quality *)
Print["Pade vs Normal difference: ", Abs[padeVal - normalVal]];
```
