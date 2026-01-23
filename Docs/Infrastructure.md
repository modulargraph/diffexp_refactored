# Infrastructure Subpackages

This document covers three foundational DiffExp subpackages that handle line geometry, Wronskian-based equation reduction, and matrix I/O:

- **LineSegmentation.m** (`DiffExp`LineSegmentation``) -- Line segmentation, validation, and dynamic interval computation
- **Wronskian.m** (`DiffExp`Wronskian``) -- Wronskian computations for combining systems into scalar equations
- **MatrixLoading.m** (`DiffExp`MatrixLoading``) -- Matrix file loading, preparation, and coupling analysis

---

## LineSegmentation.m

**Context:** `DiffExp`LineSegmentation``
**Dependencies:** `DiffExp`Symbols``, `DiffExp`State``, `DiffExp`Utilities``, `DiffExp`SeriesOps``
**Lines:** 262

### Overview

The LineSegmentation subpackage provides functions for relating parametric line segments to one another, finding singularities along lines, estimating series convergence radii, and validating/reparametrising boundary conditions. These are essential when DiffExp decomposes a transport path into multiple segments and must match solutions between them.

A "line" in this context is an `Association` mapping kinematic variables to expressions in the line parameter `x`. For example: `<|s -> 1 + 2*x, t -> 3 - x|>`.

---

### Functions

#### RelateLines

```mathematica
RelateLines[a_Association, b_Association, noerror_: False]
```

Determines the reparametrisation that maps line `a` onto line `b`. Specifically, it finds a function `f` such that `a(f(x)) = b(x)`.

**Algorithm:**
1. Identifies a component of `a` that depends on `x` (using `DependsQ`)
2. Solves the equation `a_k(y) = b_k(x)` for `y` (where `k` is the chosen component)
3. Checks that the solution is unique (errors if multiple solutions are found, indicating non-linearity)
4. Cross-checks the solution against all other components of the association, verifying `a_j(f(x)) - b_j(x) = 0` for all `j`

**Parameters:**
- `a` -- Source line (Association of variable -> expression in `x`)
- `b` -- Target line (Association of variable -> expression in `x`)
- `noerror` (default `False`) -- If `True`, returns `False` on failure instead of calling `ReportError`

**Returns:** The expression `f(x)` (the reparametrisation), or `False` if the lines cannot be related.

**Error Conditions:**
- `"Argument does not depend on the line parameter!"` -- No component of `a` depends on `x`
- `"Multiple solutions encountered while relating lines"` -- Non-linear line segment detected
- `"Could not relate lines."` -- Cross-check failed

---

#### RelateLinesPoint

```mathematica
RelateLinesPoint[a_Association, b_Association, pointb_]
```

Evaluates the line relation at a specific point. This is a convenience wrapper that computes `RelateLines[a, b]` and then substitutes `x -> pointb`.

**Parameters:**
- `a` -- Source line
- `b` -- Target line
- `pointb` -- The point at which to evaluate the relation

**Returns:** The numeric value of the reparametrisation at `pointb`.

---

#### FindMatrixSingularities

```mathematica
FindMatrixSingularities[line_, getcomplex_: False, {fixat_, to_}]
```

Determines all singularity positions along a given line by analyzing the irreducible factors stored in the differential equation matrices.

**Algorithm:**
1. Collects all singular terms by combining `MatricesIrreducibleFactors` (substituted along the line) with denominators of the line components
2. Factors and deduplicates these terms, removing any that do not depend on `x`
3. Solves each term equal to zero for `x`, yielding singularity positions (including complex ones)
4. Projects complex singularities onto the real line: a singularity at position `Re(z) + i*Im(z)` generates projected points at `Re(z) - Im(z)`, `Re(z)`, and `Re(z) + Im(z)` (the imaginary part serves as a convergence radius estimate)
5. Suppresses a projected point if another real singularity already exists in that interval
6. Sorts and deduplicates the final list

**Parameters:**
- `line` -- The line along which to find singularities
- `getcomplex` (default `False`) -- If `True`, also returns the raw complex singularity positions
- `{fixat, to}` -- The integration interval (currently used for context)

**Returns:**
- If `getcomplex` is `False`: a sorted list of projected real singularity positions
- If `getcomplex` is `True`: `{projectedSingularities, rawComplexSingularities}`

---

#### GetLargestTerm

```mathematica
GetLargestTerm[line_]
```

Extracts the largest coefficient from the expanded matrices at a specific power of `x`, for use in error estimation during dynamic segmentation.

**Algorithm:**
- Extracts coefficients of `x^n` from the expanded matrices, where `n = ExpansionOrder - ISafetyExpansionSubtract - (MaxCouplingOrder - 1)`
- Returns the maximum absolute value multiplied by `x^n`

**Parameters:**
- `line` -- The current line

**Returns:** An expression of the form `C * x^n` representing the dominant term at the truncation boundary.

---

#### GetMatricesPrecisionDistance

```mathematica
GetMatricesPrecisionDistance[line_Association]
```

Determines how far along the line parameter `x` the integration can proceed before the series truncation error exceeds the precision goal. This is the core calculation for the dynamic segmentation strategy.

**Algorithm:**
1. Validates that `AccuracyGoal` is numeric
2. Computes `DigitsNeeded = AccuracyGoal + ISafetyDigits`
3. Extracts the largest term from the expanded matrices (via `GetLargestTerm`)
4. Solves the equation `C * x^n = 10^(-DigitsNeeded)` for `x`
5. For odd expansion orders, finds both positive and negative real roots; for even orders, finds one root
6. Returns the minimum absolute value among the roots

**Parameters:**
- `line` -- The current line (must have expanded matrices available)

**Returns:** A positive real number representing the maximum safe step size in `x`.

**Error Conditions:**
- `"Accuracy goal is not given as a number."` -- `AccuracyGoal` configuration is not numeric
- `"Could not determine variation in the expanded matrices."` -- Matrix expansion yields degenerate coefficients

---

#### CheckBoundaryConditionsAndReparametrize

```mathematica
CheckBoundaryConditionsAndReparametrize[bcs3_, line_Association]
```

Validates boundary conditions for consistency and reparametrises them to match the current integration line. Handles both point boundary conditions and asymptotic (line/series) boundary conditions.

**Validation Checks:**
1. Line must be linear in `x` (no powers other than 1 or -1)
2. The BC point/line must fix all kinematic invariants and masses
3. If BCs are at a point, the BC values must not depend on `x`
4. If BCs are on a line, the values must be given as `SeriesData`
5. Dimensions of the BC matrix must match `NumIntegrals`
6. Sufficient epsilon orders must be provided

**Reparametrisation Logic:**
- Computes `lineRelation = RelateLines[line, bcsLine]`
- Evaluates the zero-limit: if `lineRelation -> 0` as `x -> 0`, the BC is at the origin of the current line
  - **Point at origin:** Wraps BC values in trivial series at `x = 0`
  - **Asymptotic line at origin:** Computes the inverse relation, verifies matching orientation, and reparametrises the series BCs via substitution and Log-term expansion
- **Point not at origin:** Simply evaluates the relation to get the `FixAt` parameter

**Parameters:**
- `bcs3` -- Boundary conditions in the form `{point_or_line, matrix_of_values}`
- `line` -- The current integration line

**Returns:** `{reparametrised_bcs, FixAt}` where:
- `reparametrised_bcs` -- `{point, values}` adapted to the current line
- `FixAt` -- The value of `x` at which the boundary conditions are imposed (0 for origin)

---

#### GetMatchingPoint

```mathematica
GetMatchingPoint[line_Association, bcsline_]
```

Determines the value of the line parameter `x` on `line` that corresponds to the boundary condition location `bcsline`.

**Algorithm:**
1. Calls `RelateLines[line, bcsline]`
2. If the result depends on `x`, checks whether it vanishes at `x = 0` (indicating an origin match)
3. If the result is a constant, that constant is the matching point

**Parameters:**
- `line` -- The current line
- `bcsline` -- The boundary condition location (Association)

**Returns:** The value of `x` at the matching point (0 for origin, or a numeric constant).

---

## Wronskian.m

**Context:** `DiffExp`Wronskian``
**Dependencies:** `DiffExp`Symbols``, `DiffExp`State``, `DiffExp`Utilities``, `DiffExp`SeriesOps``
**Lines:** 112

### Overview

The Wronskian subpackage implements the reduction of a first-order matrix ODE system to a higher-order scalar ODE for a single component. This is achieved through repeated differentiation and null-space computation, analogous to constructing a Wronskian matrix. The resulting scalar equation can then be solved more efficiently using Frobenius-type methods.

Given a system `dY/dx = A(x) * Y`, this module derives an n-th order scalar equation for one component `Y_k` by building the matrix of iterated derivatives and finding its null space.

---

### Functions

#### MatrixLogxInverse

```mathematica
MatrixLogxInverse[Mat_]
```

Computes the inverse of a matrix whose entries are polynomials in `Logx` (the log of the expansion variable). Standard matrix inversion does not handle `Logx`-dependent entries correctly when working with truncated series, so this function expands order-by-order.

**Algorithm:**

Given `M = M_0 + M_1 * Logx + M_2 * Logx^2 + ...`, computes `M^{-1}` order-by-order:
1. `M^{-1}_0 = Inverse(M_0)`
2. For `n >= 1`: `M^{-1}_n = -M^{-1}_0 * sum_{j=1}^{n} M_j * M^{-1}_{n-j}`

The maximum power of `Logx` is bounded by `dim(M) - 1`.

**Parameters:**
- `Mat` -- A square matrix with entries that are polynomials in `Logx`

**Returns:** The inverse matrix, also expressed as a polynomial in `Logx`, computed using `SExpand` for series-safe arithmetic.

---

#### NullSpaceTryAgainOnFail

```mathematica
NullSpaceTryAgainOnFail[ex_, r___]
```

A wrapper around Mathematica's `NullSpace` that detects problematic results involving degenerate `SeriesData` objects.

**Algorithm:**
1. Calls `NullSpace[ex, r]`
2. Inspects the result for `SeriesData[x, _, List[], k, k, _]` patterns with `k < 0` (empty series with negative order, indicating numerical breakdown)
3. If detected, issues a warning suggesting the `"HomogeneousSolve" -> "DontExpand"` option and aborts

**Parameters:**
- `ex` -- The matrix to compute the null space of
- `r` -- Additional options passed to `NullSpace`

**Returns:** The null space result if successful; aborts with a warning otherwise.

---

#### CombineDifferentialEquationsHomogeneous

```mathematica
CombineDifferentialEquationsHomogeneous[Amat_, topind_: 1]
```

Derives an n-th order scalar ODE for the component at index `topind` from the matrix system `dY/dx = A * Y`.

**Algorithm:**

Constructs the iterated derivative matrix `M^{(j)}` for `j = 0, 1, ..., n` where `n = dim(A)`:
- `M^{(0)} = I` (identity matrix)
- `M^{(j)} = d/dx M^{(j-1)} + M^{(j-1)} . A`

This gives `Y^{(j)} = M^{(j)} . Y`. The rows at index `topind` form a system:
```
[M^{(0)}[topind], M^{(1)}[topind], ..., M^{(n)}[topind]]^T . Y = [Y_k, Y'_k, ..., Y^{(n)}_k]^T
```

Finding the null space of the transposed matrix of these rows yields the coefficients of the scalar equation.

**Modes (controlled by `"HomogeneousSolve"` option):**
- `"Expand"` (default): Series-expands all matrices first, uses `DivisionFreeRowReduction` method with tolerance for `NullSpace`
- `"DontExpand"`: Works with exact (Together) expressions, standard `NullSpace`

**Parameters:**
- `Amat` -- The system matrix `A(x)` (square, n x n)
- `topind` (default 1) -- The component index for which to derive the scalar equation

**Returns:**
- On success: `{HomogeneousEquation, MtildeMat}` where:
  - `HomogeneousEquation` -- List of coefficients `{c_0, c_1, ..., c_n}` with `c_n = 1` (monic)
  - `MtildeMat` -- The n x n submatrix `M^{(j)}[topind]` for `j = 0..n-1` (the "Mtilde" matrix, invertible)
- On failure (`$Failed`): The null space has dimension > 1, indicating the Mtilde matrix is singular for this pivot

---

#### CombineDifferentialEquationsWithPivotSelection

```mathematica
CombineDifferentialEquationsWithPivotSelection[Amat_]
```

Wrapper that iterates over all possible pivot indices to find one that yields an invertible Mtilde matrix.

**Algorithm:**
1. Tries `CombineDifferentialEquationsHomogeneous[Amat, pivotIndex]` for `pivotIndex = 1, 2, ..., n`
2. Returns on the first successful pivot (uses `Catch`/`Throw` for early exit)
3. If all pivots fail, returns `$NeedsFallback` (triggering the VOPAlt integration strategy)

**Parameters:**
- `Amat` -- The system matrix `A(x)`

**Returns:**
- On success: `{HomogeneousEquation, MtildeMat, pivotIndex}` (the successful pivot index is appended)
- On failure: `$NeedsFallback` (all n pivots yielded singular Mtilde)

---

## MatrixLoading.m

**Context:** `DiffExp`MatrixLoading``
**Dependencies:** `DiffExp`Symbols``, `DiffExp`State``, `DiffExp`Utilities``, `DiffExp`SeriesOps``
**Lines:** 389

### Overview

The MatrixLoading subpackage handles all aspects of loading differential equation matrices from disk, preparing them for integration along specific lines, and analyzing the coupling structure to determine an efficient integration order. It supports two matrix formats:

1. **Order-by-order:** Separate files `d{var}_{order}.m` for each epsilon order (e.g., `ds_0.m`, `ds_1.m`, `dt_2.m`)
2. **Closed-form:** Files `d{var}_d.m` containing the full epsilon-dependent matrix

Additionally, it handles a special canonical matrix `d_1.m` encoding alphabet logarithms.

---

### Functions

#### LoadMatrices

```mathematica
LoadMatrices[Folder_]
```

Loads all partial derivative matrices from the configured `MatrixDirectory`. This is typically called once during initialization.

**Algorithm:**

1. **Variable Detection:**
   - Scans filenames matching `d{var}_{order}.m`
   - Extracts kinematic variable names and determines if matrices are order-by-order or closed-form (detected by `_d.m` suffix)

2. **Canonical Matrix (`d_1.m`):**
   - If present and not using closed-form, loads the first-order canonical matrix
   - Validates it contains only terms of the form `c * Log[expr]`
   - Extracts alphabet logs and assigns them symbolic labels `Log[l_i]`
   - Adds any new variables found in the alphabet to the variable list

3. **Matrix Loading:**
   - Order-by-order: Loads `d{var}_{ord}.m` for each variable and each epsilon order 0 to `EpsilonOrder`
   - Closed-form: Loads `d{var}_d.m` for each variable

4. **Validation:**
   - Checks all loaded matrices have consistent dimensions
   - Detects unsupported function heads (only allows: Association, List, Complex, Integer, Plus, Power, Rational, Symbol, Times)
   - Rejects higher-order roots (denominator of power > 2)

5. **Square Root Handling:**
   - Detects all `Power[a, b]` with `Denominator[b] = 2`
   - Verifies square root arguments are irreducible
   - Assigns `+i*delta` prescriptions to newly encountered square roots
   - Flips signs of square roots with `-i*delta` prescription according to the rule: `Sqrt[-a] -> -I*Sqrt[a]`

6. **Irreducible Factor Extraction:**
   - Collects all factors from denominators and square roots
   - These are used later by `FindMatrixSingularities` for singularity detection

**Parameters:**
- `Folder` -- The matrix directory path (typically from configuration)

**Side Effects (State Variables Set):**
- `ExpansionMatrices` -- Association of `{var, order} -> matrix`
- `ExpansionMatricesClosedForm` -- Association of `var -> matrix`
- `ExpansionMatricesCanonical1` -- The `d_1.m` matrix (or zero matrix)
- `AlphabetLogRules` -- Rules mapping symbolic logs to actual expressions
- `NumIntegrals` -- Size of the system
- `DEqnSquareRoots` -- List of square root arguments in the equations
- `MatricesIrreducibleFactors` -- List of irreducible polynomial factors
- `UsingClosedFormMatrix` -- Boolean flag

---

#### PrepareMatrices

```mathematica
PrepareMatrices[line_Association]
```

Top-level function that prepares all matrices for integration along a given line. Calls `PrepareMatricesFactored` followed by `PrepareMatricesExpanded`. Skips computation if expanded matrices already exist for this line.

**Parameters:**
- `line` -- The line along which to prepare matrices

---

#### PrepareMatricesFrom1

```mathematica
PrepareMatricesFrom1[lineorig_Association, linenew_Association]
```

Prepares factored matrices for `linenew` by reusing the factored matrices from `lineorig`, applying the reparametrisation between the two lines.

**Algorithm:**
1. Computes `ParRelns = RelateLines[lineorig, linenew]`
2. For each epsilon order, transforms: `A_new(x) = A_orig(f(x)) * f'(x)` where `f` is the reparametrisation
3. Applies `FactorOrTogether` for simplification (unless `UseMobius` is enabled)
4. Also transforms alphabet log rules if present

**Parameters:**
- `lineorig` -- Source line (must have factored matrices available)
- `linenew` -- Target line

---

#### PrepareMatricesFrom

```mathematica
PrepareMatricesFrom[lineorig_Association, linenew_Association]
```

Convenience wrapper: calls `PrepareMatricesFrom1` then `PrepareMatricesExpanded`.

---

#### PrepareMatricesFactored

```mathematica
PrepareMatricesFactored[line_Association]
```

Computes the total derivative matrix along a line by summing partial derivatives weighted by their line derivatives.

**Algorithm:**

For order-by-order matrices:
```
A_eps^n(x) = sum_v M_v^n(line(x)) * d(line_v)/dx
```
where the sum runs over all kinematic variables `v`, `M_v^n` is the partial derivative matrix for variable `v` at epsilon order `n`, and `line_v` is the component of the line for variable `v`.

For closed-form matrices:
1. Computes the full epsilon-dependent sum: `A(x, eps) = sum_v M_v(line(x), eps) * d(line_v)/dx`
2. Series-expands in epsilon to extract individual orders
3. Applies `FactorOrTogether` and `PChop` for numerical stability

Also computes `AlphabetLogRulesFactored` by substituting line values into the alphabet log arguments.

**Parameters:**
- `line` -- The integration line

**Side Effects:**
- Sets `DEqnMatricesFactored[line]` -- Association of `epsOrder -> matrix`
- Sets `AlphabetLogRulesFactored[line]` (if alphabet logs present)

---

#### PrepareMatricesExpanded

```mathematica
PrepareMatricesExpanded[line_Association]
```

Series-expands the factored matrices in the line parameter `x` around `x = 0`, up to `ExpansionOrder` terms. This is the form used directly during integration.

**Algorithm:**
1. For each epsilon order, computes `Series[A_n(x), {x, 0, ExpansionOrder}]` at `WorkingPrecision`
2. If alphabet logs are present:
   - Differentiates the factored log rules to get `d/dx Log[alpha_i(x)]`
   - Adds the canonical matrix contribution: `A_expanded[1] += M_canonical . {d/dx Log[alpha_i]}`

**Parameters:**
- `line` -- The integration line (must have factored matrices available)

**Side Effects:**
- Sets `DEqnMatricesExpanded[line]` -- Association of `epsOrder -> SeriesData matrix`
- Sets `AlphabetLogRulesExpanded[line]` (if alphabet logs present)

---

#### ClearMatrices

```mathematica
ClearMatrices[line_]
ClearMatrices[]
```

Removes cached matrix data to free memory.

**Behavior:**
- Single-line form: Removes factored, closed-form factored, and expanded matrices for the given line
- No-argument form: Clears all cached matrices for all lines
- Respects the `"KeepMatrixExpansions"` option: if `True`, does nothing

---

#### InitializeIntegrationSequence

```mathematica
InitializeIntegrationSequence[line_]
```

Analyzes the coupling structure of the differential equation system at leading epsilon order to determine an efficient integration sequence.

**Algorithm:**
1. Constructs a boolean mask of the leading-order matrix `A^{(0)}`: `True` where entries are zero, `False` where non-zero
2. Builds a directed dependency graph: edge `i -> j` means integral `i` couples to integral `j` in the homogeneous equation
3. For each integral, computes its full dependency set via `VertexOutComponent` (transitive closure)
4. Sorts integrals by dependency set size (least coupled first)
5. Groups integrals that share connected components (these must be integrated together as coupled systems)
6. Removes redundant groups (subsets of larger groups)

**Parameters:**
- `line` -- The line (must have factored matrices available)

**Side Effects:**
- Sets `IntegrationSequence` -- Ordered list of integral groups `{{i1}, {i2, i3}, ...}`
- Sets `MaxCouplingOrder` -- Maximum size of any coupled group

**Usage Context:**

The integration sequence allows DiffExp to integrate decoupled integrals independently and in the optimal order: integrals with fewer dependencies are solved first, and their solutions are then used as inhomogeneous terms when solving more coupled integrals.

---

## Relationships Between Subpackages

These three subpackages interact closely during a typical DiffExp transport computation:

1. **MatrixLoading.LoadMatrices** is called once to load all matrix files and extract singularity-relevant factors
2. **LineSegmentation.FindMatrixSingularities** uses `MatricesIrreducibleFactors` (set by `LoadMatrices`) to find singularities along each line
3. **MatrixLoading.PrepareMatrices** prepares matrices for a specific line segment
4. **MatrixLoading.InitializeIntegrationSequence** analyzes the prepared matrices to determine integration order
5. **LineSegmentation.GetMatricesPrecisionDistance** uses the expanded matrices to determine dynamic step sizes
6. **Wronskian.CombineDifferentialEquationsWithPivotSelection** uses the expanded matrices to derive scalar equations for each integral group
7. **LineSegmentation.RelateLines** is used by `PrepareMatricesFrom1` to efficiently reuse matrix computations across successive line segments
8. **LineSegmentation.CheckBoundaryConditionsAndReparametrize** validates and transforms boundary conditions before integration begins

---

## Configuration Options

The following `DiffExpConfiguration` keys are relevant to these subpackages:

| Option | Used By | Description |
|--------|---------|-------------|
| `MatrixDirectory` | `LoadMatrices` | Path to directory containing matrix files |
| `Variables` | `LoadMatrices` | Kinematic variables (auto-detected if empty) |
| `EpsilonOrder` | `LoadMatrices`, `PrepareMatricesFactored` | Maximum epsilon order |
| `AccuracyGoal` | `GetMatricesPrecisionDistance` | Target number of correct digits |
| `WorkingPrecision` | `PrepareMatricesExpanded`, `CheckBoundaryConditionsAndReparametrize` | Internal arithmetic precision |
| `ExpansionOrder` | `GetLargestTerm`, `PrepareMatricesExpanded` | Series truncation order in `x` |
| `DeltaPrescriptions` | `LoadMatrices` | i*delta prescriptions for square roots |
| `HomogeneousSolve` | `CombineDifferentialEquationsHomogeneous` | `"Expand"` or `"DontExpand"` mode |
| `UseMobius` | `PrepareMatricesFrom1` | Skip factoring when Mobius transforms are active |
| `KeepMatrixExpansions` | `ClearMatrices` | Prevent cache clearing |

---

## State Variables

Key state variables managed by these subpackages:

| Variable | Set By | Description |
|----------|--------|-------------|
| `ExpansionMatrices` | `LoadMatrices` | Raw loaded matrices (order-by-order) |
| `ExpansionMatricesClosedForm` | `LoadMatrices` | Raw loaded matrices (closed-form) |
| `ExpansionMatricesCanonical1` | `LoadMatrices` | The `d_1.m` canonical matrix |
| `AlphabetLogRules` | `LoadMatrices` | Symbolic-to-actual log mappings |
| `NumIntegrals` | `LoadMatrices` | System dimension |
| `DEqnSquareRoots` | `LoadMatrices` | Square roots in the equations |
| `MatricesIrreducibleFactors` | `LoadMatrices` | Factors for singularity detection |
| `UsingClosedFormMatrix` | `LoadMatrices` | Boolean: closed-form mode |
| `DEqnMatricesFactored` | `PrepareMatricesFactored` | Factored matrices per line |
| `DEqnMatricesExpanded` | `PrepareMatricesExpanded` | Expanded matrices per line |
| `IntegrationSequence` | `InitializeIntegrationSequence` | Ordered list of integral groups |
| `MaxCouplingOrder` | `InitializeIntegrationSequence` | Largest coupled block size |

---

## File Format Conventions

### Order-by-Order Matrices

Files are named `d{variable}_{epsilonOrder}.m` and contain a Mathematica matrix expression:

```
Tests/Banana_Matrices/
  dmm1_0.m   -- d/d(mm1) at eps^0
  dmm1_1.m   -- d/d(mm1) at eps^1
  ...
  dpsq_4.m   -- d/d(psq) at eps^4
```

### Closed-Form Matrices

Files are named `d{variable}_d.m` and contain the full epsilon-dependent matrix.

### Canonical Matrix

The file `d_1.m` (if present) contains the eps^1 contribution in canonical (dlog) form:
```mathematica
(* Example d_1.m content *)
{{Log[s] + 2*Log[t-4], -Log[s-t]}, {0, 3*Log[t]}}
```

Each entry must be a sum of terms `c * Log[expr]` where `c` is numeric.
