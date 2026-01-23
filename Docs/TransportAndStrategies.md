# Transport, Integration Strategies, and Main Entry Point

This document covers the three core subpackages that form the user-facing API of DiffExp:

1. **DiffExp.m** (`DiffExp``) -- Main entry point and configuration API
2. **Transport.m** (`DiffExp`Transport``) -- Transport and integration logic
3. **IntegrationStrategies.m** (`DiffExp`IntegrationStrategies``) -- Strategy implementations for solving differential equations

Together these packages implement the complete workflow: loading a configuration, preparing boundary conditions, and transporting solutions along parametric lines through kinematic space.

---

## Table of Contents

- [Architecture Overview](#architecture-overview)
- [DiffExp.m -- Main Entry Point](#diffexpm----main-entry-point)
- [Transport.m -- Transport and Integration](#transportm----transport-and-integration)
- [IntegrationStrategies.m -- Strategy Implementations](#integrationstrategiesm----strategy-implementations)
- [Complete Workflow](#complete-workflow)
- [Configuration Reference](#configuration-reference)
- [Error Handling and Diagnostics](#error-handling-and-diagnostics)

---

## Architecture Overview

```
User Code
    |
    v
+------------------+
|   DiffExp.m      |  Configuration API (LoadConfiguration, UpdateConfiguration)
+------------------+
    |
    v
+------------------+
|  Transport.m     |  PrepareBoundaryConditions, IntegrateSystem, TransportTo, ToPiecewise
+------------------+
    |
    v
+----------------------------+
| IntegrationStrategies.m    |  DispatchStrategy -> SolveDefault / SolveVOP / SolveVOPAlt / SolveSimple
+----------------------------+
    |
    v
+---------------------------------------+
| Supporting packages                   |
| Frobenius, Wronskian, SeriesOps,      |
| Integration, Pade, Mobius,            |
| AnalyticContinuation, LineSegmentation|
+---------------------------------------+
```

### Data Flow

```
LoadConfiguration
    |
    +--> Loads matrices from disk
    +--> Sets epsilon order, expansion order, precision
    +--> Parses DeltaPrescriptions
    |
    v
PrepareBoundaryConditions(bcs, line)
    |
    +--> Converts user BCs to internal series form
    +--> Prepares analytic continuation
    |
    v
TransportTo(bcs, line, endpoint)
    |
    +--> Finds singularities on line
    +--> Determines segmentation (pole intervals)
    +--> For each segment:
    |       |
    |       +--> Rescales line (Mobius if enabled)
    |       +--> Expands matrices to series
    |       +--> Calls IntegrateSystem
    |       |       |
    |       |       +--> For each (eps_order, integral_group):
    |       |       |       |
    |       |       |       +--> Computes inhomogeneous term bVec
    |       |       |       +--> Calls DispatchStrategy
    |       |       |       |       |
    |       |       |       |       +--> SolveSimple / SolveDefault / SolveVOP / SolveVOPAlt
    |       |       |       |
    |       |       |       +--> Fixes boundary conditions (linear solve for c_i)
    |       |       |
    |       |       +--> Returns series solutions table
    |       |
    |       +--> Evaluates series at segment boundary (SEval / Pade)
    |       +--> Updates boundary conditions for next segment
    |
    +--> Returns {endpoint_values, [errors], [saved_data]}
```

---

## DiffExp.m -- Main Entry Point

**File:** `/DiffExp/DiffExp.m`
**Context:** `DiffExp``
**Lines:** 197
**Dependencies:** All 16 subpackages

This is the main package file that loads all subpackages and provides the user-facing configuration API. It prints version information on load and re-exports symbols from all subpackages.

### Package Loading

On load, the package:
1. Loads all 16 subpackages via `Get[]` in dependency order
2. Begins the `DiffExp`` package context with all subpackage contexts imported
3. Prints version info (`v1.1 refactored`)

### Functions

#### CurrentConfiguration

```mathematica
CurrentConfiguration[]
```

Returns an Association containing only the user-facing configuration keys:

| Key | Description |
|-----|-------------|
| `AccuracyGoal` | Target accuracy for transport |
| `ChopPrecision` | Digits after decimal before chopping to zero |
| `DeltaPrescriptions` | Singularity prescriptions (branch cuts) |
| `DivisionOrder` | Inverse distance to singularity for segment evaluation |
| `EpsilonOrder` | Maximum order in dimensional regulator |
| `ExpansionOrder` | Maximum power of line parameter in series |
| `IntegrationStrategy` | Strategy for solving DEs ("Default", "VOP", "VOPAlt") |
| `LineParameter` | Symbol used for line parametrization |
| `MatrixDirectory` | Path to matrix files on disk |
| `RadiusOfConvergence` | Rescaling factor for line segments |
| `SegmentationStrategy` | How segments are determined ("Dynamic", "Predivision") |
| `UseMobius` | Whether to use Mobius transformations |
| `UsePade` | Whether to use Pade approximants |
| `Variables` | Kinematic invariants and masses |
| `Verbosity` | Output level (1=info, 2=detailed, 3=debug) |
| `WorkingPrecision` | Internal numerical precision |

**Returns:** An `Association` with the above keys and their current values.

#### LoadConfiguration

```mathematica
LoadConfiguration[optionList_List]
```

Resets all configuration to defaults, then applies the given options. This is the recommended way to initialize DiffExp for a new computation.

**Parameters:**
- `optionList` -- A list of rules (e.g., `{MatrixDirectory -> "path/", EpsilonOrder -> 4, ...}`)

**Behavior:**
1. Resets `DiffExpConfiguration` to `DefaultConfiguration`
2. Calls initialization functions of any registered extensions
3. Delegates to `UpdateConfiguration` for the actual option processing

**Returns:** The result of `CurrentConfiguration[]` after applying options.

#### UpdateConfiguration

```mathematica
UpdateConfiguration[rule1, rule2, ...]
UpdateConfiguration[{rule1, rule2, ...}]
UpdateConfiguration[assoc_Association]
```

Updates the current configuration without resetting to defaults. Accepts rules, a list of rules, or an Association.

**Parameters:**
- Rules or Association with configuration keys and values

**Validation and Side Effects:**

1. **Precision check:** Verifies `ChopPrecision < WorkingPrecision`; aborts with error otherwise
2. **LogFile:** Opens/reopens a log stream if `LogFile` key is provided
3. **CrosscheckLevel:** Updates active crosscheck flags based on level
4. **LineParameter:** Updates the internal `x` symbol; validates it is not in `Variables`
5. **DeltaPrescriptions:** Parses expressions of the form `poly + I*delta` into `{poly, sign}` pairs; validates irreducibility
6. **MatrixDirectory / EpsilonOrder:** Triggers matrix loading from disk via `LoadMatrices`
7. **ExpansionOrder:** Updates `ExpansionOrderVal`
8. **EstimateError:** Validates and normalizes the error estimation setting
9. **Cache clearing:** Resets all cached factored/expanded matrices, analytic continuation data, and integration sequences
10. **Extensions:** Calls registered extension update functions
11. **IntReps:** Calls `UpdateIntReps` to reset integration representations

**Returns:** The result of `CurrentConfiguration[]`.

### Usage Example

```mathematica
(* Load the package *)
Get["DiffExp.m"];

(* Initialize with configuration *)
LoadConfiguration[{
  MatrixDirectory -> "/path/to/matrices/",
  DeltaPrescriptions -> {t - 16 + I * \[Delta]},
  EpsilonOrder -> 4,
  ExpansionOrder -> 50,
  WorkingPrecision -> 500,
  ChopPrecision -> 200,
  UseMobius -> True,
  UsePade -> True,
  Verbosity -> 2
}];

(* Later, update a single option *)
UpdateConfiguration[ExpansionOrder -> 80];
```

---

## Transport.m -- Transport and Integration

**File:** `/DiffExp/Transport.m`
**Context:** `DiffExp`Transport``
**Lines:** 1128
**Dependencies:** Symbols, State, Utilities, SeriesOps, Integration, Pade, Mobius, AnalyticContinuation, LineSegmentation, Frobenius, Wronskian, MatrixLoading, IntegrationStrategies

This package contains the main computational functions that users interact with directly.

### Functions

#### PrepareBoundaryConditions

```mathematica
PrepareBoundaryConditions[bcs_List, line_Association | line_List]
```

Converts user-supplied boundary conditions into the internal format required by `IntegrateSystem` and `TransportTo`.

**Parameters:**
- `bcs` -- A list with one entry per master integral. Each entry can be:
  - A `List` of explicit epsilon-order coefficients: `{c0, c1, c2, ...}`
  - The string `"?"` to indicate unknown/ignored boundary conditions
  - A closed-form expression in `eps` (will be series-expanded)
- `line` -- An Association or List defining the parametric line, e.g., `<|t -> 1/x, s -> 3|>`

**Behavior:**
1. Sorts and validates the line; increases precision if needed
2. Prepares analytic continuation for the line
3. For each integral's BC:
   - **List:** Takes the first `EpsilonOrder + 1` entries as epsilon coefficients
   - **"?":** Marks as ignored (free parameters will be introduced)
   - **Expression:** Series-expands in `eps` to required order; sets sub-leading orders to zero for asymptotic BCs
4. For non-point lines (asymptotic BCs): evaluates coefficients on the line to get series in `x`
5. For point lines: multiplies by unit precision number
6. Replaces `Log[a*x]` with `Log[a] + Logx` for internal representation

**Returns:** `{line, prepared_bcs}` -- A two-element list with the (possibly modified) line and the prepared boundary condition table (dimensions: `NumIntegrals x (EpsilonOrder+1)`).

**Example:**
```mathematica
(* Asymptotic boundary conditions at t -> infinity via t = -1/x *)
bcs = {
  "?",   (* Unknown - will introduce free parameter *)
  "?",
  closedFormExpr,  (* Closed-form in eps *)
  eps^3 * Gamma[eps]^3  (* Another closed-form *)
} // PrepareBoundaryConditions[#, <|t -> -1/x|>] &;
```

#### IntegrateSystem

```mathematica
IntegrateSystem[line_Association | line_List]
IntegrateSystem[bcs_List, line_Association | line_List]
IntegrateSystem[bcs_List, line_Association | line_List, opts_List]
```

Obtains series solutions of the differential equation system along a given line. This is the core solver that computes local Taylor/Laurent series solutions.

**Parameters:**
- `bcs` -- Prepared boundary conditions (from `PrepareBoundaryConditions`), or `"?"` for general solutions with free constants
- `line` -- The parametric line (Association or List of rules)
- `opts` -- Internal options list (e.g., `{"TransportToCall"}` when called from `TransportTo`)

**Algorithm:**

```
1. Validate line fixes all kinematic variables
2. Prepare analytic continuation for the line
3. Expand differential equation matrices to series on the line
4. Initialize integration sequence (topological ordering of coupled integrals)

5. MAIN DOUBLE LOOP:
   For epsord = 0, 1, ..., EpsilonOrder:
     Precompute bVec contributions from lower epsilon orders

     For each intind in IntegrationSequence:
       a) Compute inhomogeneous term bVec:
          bVec = A[0]_{intind,complement} . f_{complement,epsord}    (same eps, other integrals)
                + A[1] . f_{all, epsord-1}                            (one eps order lower)
                + Sum_{l>=2} A[l] . f_{all, epsord-l}                 (lower eps orders)

       b) Dispatch to integration strategy:
          {cIndices, fGeneral, BufferedData} = DispatchStrategy[intind, bVec, line, epsord, BufferedData]

       c) Fix boundary conditions:
          - Extract relevant BCs for current integrals and eps order
          - Set up linear system: fGeneral(x=FixAt) = BC values
          - Solve for constants c_i
          - If underdetermined: introduce free parameters c_{epsord,intind,i}
          - Apply constant replacements to fGeneral

       d) Store results in IntegrationData[{i, epsord}]

6. Cross-check: verify no unexpected multivalued functions remain
7. Apply analytic continuation to final results
```

**Returns:** A table of series solutions with dimensions `NumIntegrals x (EpsilonOrder + 1)`. Each entry is a `SeriesData` object in `x` (possibly containing `Logx` terms).

When called with `"TransportToCall"` option, skips some initialization steps and returns raw (un-continued) series.

**Notes:**
- The function uses a `BufferedData` association to cache expensive computations (Wronskians, FMat/FMatInv) across epsilon orders
- If boundary conditions are insufficient, free parameters `Subscript[c, epsord, intind, i]` are introduced and Pade is disabled
- Cross-checks can be enabled via `CrosscheckLevel` configuration

#### TransportTo

```mathematica
TransportTo[bcs_List, line_Association | line_List, to_:1, save_:False, SampleAtList_List:{}]
```

The main user-facing function. Transports boundary conditions along a parametric line to a target point, handling singularities automatically by segmenting the path.

**Parameters:**
- `bcs` -- Prepared boundary conditions (output of `PrepareBoundaryConditions` or a previous `TransportTo` call)
- `line` -- Parametric line through kinematic space (e.g., `<|t -> x, s -> 3|>` or `<|t -> 10|>` for a point)
- `to` -- Target value of line parameter `x` (default: 1). Ignored if `line` is a point.
- `save` -- If `True`, saves all intermediate segment expansions for later use with `ToPiecewise`
- `SampleAtList` -- Optional: `{{x_values}, directory, [precision]}` to evaluate and export at specific points during transport

**Algorithm (Segment-by-Segment Transport):**

```
1. PREPROCESSING:
   a) Validate inputs; increase precision if needed
   b) If line is a point: construct linear interpolation from BC point to target
   c) Parse indeterminate parameters; disable Pade if present
   d) Check line does not lie on a singularity
   e) Reparametrize BCs to match line

2. SINGULARITY ANALYSIS:
   a) Find all singularities of the DE matrices on the line
   b) For each singularity: compute Mobius-rescaled line segment
   c) Determine pole intervals (regions where singular expansions are valid)

3. SEGMENT COUNTING (Predivision strategy only):
   a) Simulate the transport without integration to count segments
   b) Used for progress reporting and precision estimation

4. MAIN TRANSPORT LOOP:
   currentCenter = FixAt
   CurrLine = rescaled line around currentCenter

   While not Done:
     a) Prepare analytic continuation for current segment
     b) Expand DE matrices on current line segment
     c) Determine valid interval for current expansion

     d) Determine evaluation point:
        - If target is within current interval: evaluate at target, Done=True
        - If a pole interval overlaps: find midpoint of overlap region
        - Otherwise: evaluate at interval boundary

     e) ACCURACY GOAL (Before mode):
        - Estimate matrix expansion error
        - Dynamically adjust ExpansionOrder up/down to meet target accuracy

     f) Call IntegrateSystem on current segment with current BCs

     g) Evaluate series at CurrEvalPoint:
        - Use Pade approximants if enabled
        - Estimate error (Fast mode): compare with lower-order evaluation

     h) ACCURACY GOAL (After mode):
        - If error exceeds target: increase ExpansionOrder, repeat segment

     i) Update BCs for next segment:
        Currbcs = {line /. x -> CurrEvalPoint, CurrEval}

     j) Advance to next segment center:
        - If next is a pole: use pole's Mobius line
        - Otherwise: find next center point

     k) SaveExpansions: store segment data with line relations
     l) SampleAtList: export evaluations at requested points

     m) Handle MultivaluedFail: abort if analytic continuation impossible
```

**Returns:** Depends on options:

| save | EstimateError | Return value |
|------|---------------|--------------|
| False | False | `{endpoint_association, values_table}` |
| False | "Fast" | `{endpoint_association, values_table, errors_table}` |
| True | False | `{{endpoint_association, values_table}, segment_data_list}` |
| True | "Fast" | `{{endpoint_association, values_table, errors_table}, segment_data_list}` |

Where:
- `endpoint_association` -- The kinematic point at the endpoint (e.g., `<|t -> 10, s -> 3|>`)
- `values_table` -- Numerical values, dimensions `NumIntegrals x (EpsilonOrder + 1)`
- `errors_table` -- Error estimates, same dimensions
- `segment_data_list` -- List of `{CurrLine, lineRelation, interval, CurrLineInterval, expansionData}`

**Chaining Transports:**

The output of `TransportTo` can be directly used as input for another transport:

```mathematica
result1 = TransportTo[bcs, <|t -> x, s -> 3|>, 5];
result2 = TransportTo[result1, <|t -> x, s -> 3|>, 10];
```

**Error Estimation:**

When `EstimateError -> "Fast"` is configured:
- At each segment, the series is evaluated at two orders (full and reduced)
- The difference provides an error estimate for that segment
- Errors accumulate additively across segments
- Both endpoint and boundary-fix-point errors are computed; the maximum is used

#### ToPiecewise

```mathematica
ToPiecewise[savedData_List, pade_:False, ord_:Null]
```

Converts saved segment data (from `TransportTo` with `save=True`) into Piecewise pure functions that can be evaluated at any point on the transport line.

**Parameters:**
- `savedData` -- Either the full output of `TransportTo[..., True]` or just the segment data list `savedData[[2]]`
- `pade` -- If `True`, applies Pade approximants to each segment before constructing the piecewise function
- `ord` -- Optional: truncate each segment's series to this order before evaluation

**Behavior:**
1. Extracts segment data (handles both compressed and uncompressed formats)
2. For compressed data stored in files: loads and uncompresses on demand (with memoization)
3. For each integral and epsilon order:
   - Constructs a `Piecewise` expression over all segments
   - Each piece: evaluates the series (Normal or Pade), substitutes `Logx -> Log[x]`, applies theta functions, maps line parameter
   - Domain: the interval covered by each segment on the original line
4. Wraps in a pure function of `x` (the original line parameter)

**Returns:** A table of pure functions, dimensions `NumIntegrals x (EpsilonOrder + 1)`. Each function takes a single argument (the line parameter value) and returns the numerical value of that integral at that epsilon order.

**Example:**
```mathematica
(* Transport with saved expansions *)
result = TransportTo[bcs, <|t -> x|>, 10, True];

(* Convert to piecewise functions *)
fns = ToPiecewise[result];

(* Evaluate integral 3 at epsilon^2 order at t = 5 *)
fns[[3, 3]][5]

(* With Pade for better accuracy *)
fnsPade = ToPiecewise[result, True];
```

### Internal Helper Functions

#### IntervalOverlapQ, IntervalIntersec, IntervalContainsQ

```mathematica
IntervalOverlapQ[intv1_, intv2_]     (* True if intervals overlap *)
IntervalIntersec[intv1_, intv2_]     (* Returns {min, max} of intersection *)
IntervalContainsQ[intv_, point_]     (* True if point in [a,b] *)
```

Utility functions for interval arithmetic used in the segmentation logic.

---

## IntegrationStrategies.m -- Strategy Implementations

**File:** `/DiffExp/IntegrationStrategies.m`
**Context:** `DiffExp`IntegrationStrategies``
**Lines:** 697
**Dependencies:** Symbols, State, Utilities, SeriesOps, Integration, Frobenius, Wronskian

This package implements the different strategies for solving the coupled system of first-order differential equations at each segment and epsilon order.

### Mathematical Background

The system of differential equations has the form:

```
df/dx = A(x, eps) . f + b(x, eps)
```

where `A` is expanded in `eps`:
```
A = A[0] + eps * A[1] + eps^2 * A[2] + ...
```

At each epsilon order, the problem reduces to:
```
df/dx = A[0] . f + bVec(x)
```

where `bVec` contains contributions from lower epsilon orders. The strategies differ in how they find the homogeneous solutions and construct the general solution.

### Shared Helper Functions

These are private functions used by multiple strategies.

#### ComputeMatricesMSupj

```mathematica
ComputeMatricesMSupj[intind_List, line_]
```

Computes the sequence of M^(j) matrices via nested differentiation. These matrices encode the relationship between the first-order system and higher-order scalar equations.

**Definition:**
```
M^(0) = Identity matrix (n x n)
M^(j+1) = dM^(j)/dx + M^(j) . A[0]_{intind,intind}
```

**Parameters:**
- `intind` -- List of integral indices forming a coupled subsystem
- `line` -- Current line segment

**Behavior:**
- If `HomogeneousSolve === "Expand"`: uses series arithmetic (`SExpand`, `SD`)
- Otherwise: uses `Together` with exact differentiation

**Returns:** A list of `n+1` matrices `{M^(0), M^(1), ..., M^(n)}` where `n = Length[intind]`.

**Used by:** SolveVOP, SolveVOPAlt

#### ComputeNthOrderDiffEqns

```mathematica
ComputeNthOrderDiffEqns[MatricesMSupj_List, MyN_Integer]
```

Extracts nth-order scalar differential equations from the M^(j) matrices. Each coupled integral gives rise to a scalar ODE of order at most `n`.

**Algorithm:**
For each integral index `ind`:
1. Form the matrix `MatricesMSupj[[All, ind]]` (all M^(j) restricted to column `ind`)
2. Find null vectors of its transpose (coefficients of the scalar ODE)
3. Select the minimal-length null vector (lowest order ODE)
4. Normalize by dividing by the last coefficient (leading coefficient = 1)

**Parameters:**
- `MatricesMSupj` -- Output of `ComputeMatricesMSupj`
- `MyN` -- Number of integrals (system size)

**Returns:** A list of coefficient vectors, one per integral. Each vector `{a0, a1, ..., an, 1}` represents the ODE: `a0*f + a1*f' + ... + an*f^(n) + f^(n+1) = 0` (or lower order if the null space permits).

**Used by:** SolveVOP, SolveVOPAlt

#### ComputeWronskianMatrix

```mathematica
ComputeWronskianMatrix[homogeneousSolutions_List]
```

Builds the Wronskian matrix from a set of homogeneous solutions.

**Definition:**
```
W[j, i] = d^j/dx^j (solution_i),    j = 0, ..., n-1
```

**Parameters:**
- `homogeneousSolutions` -- List of `n` solutions to a homogeneous ODE

**Returns:** An `n x n` matrix where entry `(j, i)` is the `j`-th derivative of the `i`-th solution.

**Used by:** SolveVOP, SolveVOPAlt, SolveDefault

#### ReduceIndeterminates

```mathematica
ReduceIndeterminates[fGeneral_, vanishingTerms_, numSolutions_Integer]
```

Reduces the number of free constants in a general solution by enforcing consistency constraints order-by-order in the series.

**Algorithm:**
```
1. Extract all constants c_{i,j} appearing in fGeneral
2. If vanishingTerms are not all zero:
   seriesOrder = -1
   While (number of free constants) > numSolutions:
     a) Expand vanishingTerms to current seriesOrder
     b) Extract coefficient equations (including Logx coefficients)
     c) Set up linear system: Cmat . constants = Cb
     d) Solve: particular solution + null space
     e) Express old constants in terms of fewer new constants
     f) Increment seriesOrder

   If too few constants remain: report error
   Apply the constant replacements to fGeneral
3. If vanishingTerms are all zero: just rename constants
```

**Parameters:**
- `fGeneral` -- Current general solution (may have too many free constants)
- `vanishingTerms` -- Expressions that should be zero (e.g., residuals of the DE)
- `numSolutions` -- Target number of free constants

**Returns:** Updated `fGeneral` with reduced constants.

**Used by:** SolveVOP, SolveVOPAlt (inline versions)

#### ComputeGMat

```mathematica
ComputeGMat[FMat_, FMatInv_, bVec_, intind_]
```

Computes the general solution matrix using the fundamental matrix and its inverse.

**Formula:**
```
BMat = (1/n) * [bVec | bVec | ... | bVec]^T    (n copies, transposed)
GMat = FMat . ( Integrate[FMatInv . BMat, x] + diag(c_1, ..., c_n) )
```

This implements variation of parameters: the particular solution is `FMat . Integral[FMatInv . BMat]` and the homogeneous solution is `FMat . diag(c_i)`.

**Parameters:**
- `FMat` -- Fundamental matrix (columns are independent solutions)
- `FMatInv` -- Inverse of fundamental matrix
- `bVec` -- Inhomogeneous term vector
- `intind` -- Integral indices (determines system size)

**Returns:** `{cIndices, GMat}` where:
- `cIndices` -- List of constant symbols `{Subscript[c,1], ..., Subscript[c,n]}`
- `GMat` -- The general solution matrix (each column is a solution)

**Used by:** SolveDefault, SolveVOPAlt

---

### Strategy Functions

#### DispatchStrategy

```mathematica
DispatchStrategy[intind_, bVec_, line_, epsord_, BufferedData_]
```

The main dispatcher that routes to the appropriate strategy based on problem type and configuration.

**Dispatch logic:**

```
1. If Length[intind] == 1 AND A[0]_{intind,intind} == {{0}}:
   -> SolveSimple (trivial case: just integrate)

2. If IntegrationStrategy == "Default":
   -> SolveDefault

3. If IntegrationStrategy == "VOP" or "VariationOfParameters":
   -> SolveVOP

4. If IntegrationStrategy == "VOPAlt":
   -> SolveVOPAlt
```

**Parameters:**
- `intind` -- List of integral indices forming a coupled block
- `bVec` -- Inhomogeneous term vector (from lower eps orders and other integrals)
- `line` -- Current line segment
- `epsord` -- Current epsilon order (0, 1, 2, ...)
- `BufferedData` -- Association for caching expensive computations across eps orders

**Returns:** `{cIndices, fGeneral, BufferedData}` where:
- `cIndices` -- List of undetermined constant symbols
- `fGeneral` -- General solution vector (series in `x`, containing constants)
- `BufferedData` -- Updated buffer (may contain new cached data)

---

#### SolveSimple

```mathematica
SolveSimple[intind_, bVec_, line_, epsord_]
```

Handles the trivial case: a single integral with `A[0] = 0` (no coupling to itself at leading epsilon order).

**Algorithm:**
```
If bVec == {0}:
  f = c_1 + O[x]^(ExpansionOrder + 1)
Else:
  f = Integrate[bVec[[1]], x] + c_1 + O[x]^(ExpansionOrder + 1)
```

**Returns:** `{cIndices, fGeneral}` (no BufferedData needed).

---

#### SolveDefault

```mathematica
SolveDefault[intind_, bVec_, line_, epsord_, BufferedData_]
```

The default strategy using Frobenius solutions and Wronskian inversion.

**Algorithm (at eps^0 -- first call for this integral block):**

```
1. COMBINING DIFFERENTIAL EQUATIONS:
   - Call CombineDifferentialEquationsWithPivotSelection on A[0]_{intind,intind}
   - Selects a pivot integral and combines the first-order system into a
     single nth-order scalar ODE
   - If all pivots fail ($NeedsFallback): fall back to SolveVOPAlt

   Output: {HomogeneousEquation, MtildeMat, selectedPivot}

2. COMPANION MATRIX:
   Construct NMat (companion matrix for the nth-order ODE):
   NMat = [[0,1,0,...],
           [0,0,1,...],
           [...        ],
           [-a0,-a1,...,-a_{n-1}]]

3. FROBENIUS SOLUTIONS:
   Solve the scalar ODE using the Frobenius method
   -> n linearly independent solutions (possibly with Logx terms)

4. WRONSKIAN CONSTRUCTION:
   Build Wronskian matrix: W[j,i] = d^j/dx^j (solution_i)

5. WRONSKIAN INVERSION (selected by InvWronskSolver setting):
   a) "Frobenius" method (for Logx-containing Wronskians):
      - Solve the adjoint (transposed) system
      - Build adjoint Wronskian
      - Multiply and invert the constant part
      - Cross-check: W_inv * W should be identity (up to precision)

   b) "Inverse" method (for purely rational/algebraic Wronskians):
      - Direct matrix inversion with DivisionFreeRowReduction

   c) "InverseLogx" method:
      - Uses MatrixLogxInverse (specialized for Logx structure)

   d) "Auto" (default):
      - Uses "Frobenius" if Logx present, "Inverse" otherwise

6. FUNDAMENTAL MATRIX:
   FMat = MtildeInv . Wronskian      (transforms from scalar ODE to system)
   FMatInv = WronskInv . Mtilde       (inverse transformation)

7. BUFFER:
   Store {FMat, FMatInv} in BufferedData[intind] for reuse at higher eps orders
```

**Algorithm (at eps^k, k > 0 -- using buffered data):**

```
1. Retrieve {FMat, FMatInv} from BufferedData[intind]
2. Compute GMat using ComputeGMat (variation of parameters formula)
3. fGeneral = sum of columns of GMat^T
```

**Cross-checks available** (enabled via CrosscheckLevel):
- `"Wronskians"`: Verifies `dW/dx = NMat . W`
- `"WronskInv"`: Verifies `W_inv . W = Identity`
- `"PeriodMatrix"`: Verifies `dFMat/dx = A[0] . FMat`
- `"GeneralSolutionMatrix"`: Verifies `dGMat/dx = A[0] . GMat + BMat`

**Returns:** `{cIndices, fGeneral, BufferedData}`

---

#### SolveVOP

```mathematica
SolveVOP[intind_, bVec_, line_, epsord_, BufferedData_]
```

Variation of Parameters strategy. Derives nth-order scalar ODEs for each integral independently, solves them, and uses a Cramer's-rule-like formula for the particular solution.

**Algorithm (at eps^0):**

```
1. Compute M^(j) matrices via ComputeMatricesMSupj
2. Extract nth-order ODEs via ComputeNthOrderDiffEqns

3. Determine solving approach:
   - If highest-order ODE has order < n:
     -> Solve ALL integrals' ODEs independently
     -> SolveFrom = all indices
   - If highest-order ODE has order = n:
     -> Solve only the highest-order ODE
     -> SolveFrom = {highest_order_position}
     -> Set up MtildeMatrix for back-transformation

4. For each integral to solve from:
   a) Solve homogeneous ODE via FrobeniusSolutions
   b) Build Wronskian matrix via ComputeWronskianMatrix
   c) Compute 1/det(Wronskian) as a series (inverted once, reused)

5. Buffer all data for reuse
```

**Algorithm (at eps^k, using buffered data):**

```
1. Compute b^(j) vectors (derivatives of bVec transformed by M^(j)):
   b^(0) = bVec
   b^(j+1) = M^(j) . bVec + d/dx b^(j)

2. For each integral to solve from:
   a) Compute inhomogeneous term for the scalar ODE:
      g = ODE_coefficients . b^(j)_values

   b) Apply Cramer's rule variant:
      f = Sum_i [ solution_i * (c_{i} + Integral[
            (1/det(W)) * det(W_with_column_i_replaced) * g
          ]) ]

      where the column replacement puts {0,...,0, g} in column i

3. If SolveFrom has multiple entries (underdetermined):
   - Compute vanishing terms: df/dx - A[0].f - bVec (should be zero)
   - Use ReduceIndeterminates logic inline to eliminate excess constants

4. If SolveFrom has one entry:
   - Compute derivative vector of f
   - Back-transform via MtildeInv to get full system solution
   - Rename constants: c_{myintind,j} -> c_j
```

**Returns:** `{cIndices, fGeneral, BufferedData}`

---

#### SolveVOPAlt

```mathematica
SolveVOPAlt[intind_, bVec_, line_, epsord_, BufferedData_]
```

Alternative Variation of Parameters strategy. Similar to VOP but constructs the fundamental matrix (FMat) directly from homogeneous solutions, and uses `MatrixLogxInverse` for inversion when logarithms are present. This serves as the fallback when SolveDefault fails (all pivots rejected).

**Algorithm (at eps^0):**

```
1. Compute M^(j) matrices and nth-order ODEs (same as VOP)

2. Solve ALL integrals' ODEs (always uses SolveFrom = all indices)

3. Build Wronskian matrices for each integral

4. Construct general homogeneous solution:
   For each integral index myintind:
     f_myintind = Sum_i [ c_{myintind,i} * solution_i ]

5. Reduce indeterminates:
   - Compute vanishing terms: df/dx - A[0].f (should be zero for homogeneous)
   - Eliminate excess constants order-by-order until n remain

6. Construct FMat directly:
   For each remaining constant c_j:
     FMat column j = fGeneral with c_j=1, all other c_k=0

7. Invert: FMatInv = MatrixLogxInverse[FMat]
   (handles the case where FMat contains Logx terms)

8. Buffer {FMat, FMatInv}
```

**Algorithm (at eps^k, using buffered data):**

```
1. Retrieve {FMat, FMatInv} from buffer
2. Compute GMat using ComputeGMat (same as SolveDefault)
3. Cross-check if enabled
4. fGeneral = sum of columns of GMat^T
```

**Returns:** `{cIndices, fGeneral, BufferedData}`

---

### Strategy Comparison

| Feature | SolveSimple | SolveDefault | SolveVOP | SolveVOPAlt |
|---------|-------------|--------------|----------|-------------|
| System size | 1 | Any | Any | Any |
| A[0] requirement | Must be 0 | Any (needs valid pivot) | Any | Any |
| Homogeneous solve | None | Single nth-order ODE | Per-integral ODEs | Per-integral ODEs |
| Wronskian inversion | N/A | Frobenius/Direct/LogxInverse | det(W) inversion | MatrixLogxInverse |
| Handles Logx | N/A | Yes (Frobenius solver) | Yes | Yes |
| Fallback | N/A | Falls back to VOPAlt | None | None |
| Performance | Fastest | Fast (single ODE) | Medium | Slowest (always solves all) |
| When to use | Decoupled, zero A[0] | Default choice | Alternative | Last resort / fallback |

---

## Complete Workflow

### Example: Equal Mass Banana Graph

```mathematica
(* 1. Load package *)
Get["DiffExp.m"];

(* 2. Configure *)
LoadConfiguration[{
  DeltaPrescriptions -> {t - 16 + I * \[Delta]},
  MatrixDirectory -> "Banana_EqualMass_Matrices/",
  EpsilonOrder -> 3,
  ExpansionOrder -> 50,
  WorkingPrecision -> 500,
  ChopPrecision -> 200,
  UseMobius -> True,
  UsePade -> True,
  Verbosity -> 2
}];

(* 3. Prepare asymptotic boundary conditions at t -> infinity *)
bcs = {
  "?",
  "?",
  closedFormExpressionForIntegral3,
  closedFormExpressionForIntegral4
} // PrepareBoundaryConditions[#, <|t -> -1/x|>] &;

(* 4. Transport to finite point t = -1 *)
result1 = TransportTo[bcs, <|t -> -1|>];
(* result1 = {<|t -> -1|>, {{val_1_eps0, val_1_eps1, ...}, ...}} *)

(* 5. Transport along a new line to t = 10, saving expansions *)
result2 = TransportTo[result1, <|t -> x|>, 10, True];
(* result2 = {{<|t -> 10|>, values, errors}, segmentData} *)

(* 6. Convert to piecewise functions for arbitrary evaluation *)
fns = ToPiecewise[result2];

(* 7. Evaluate at any point on the line *)
value = fns[[3, 2]][5]  (* Integral 3, eps^1 order, at t = 5 *)
```

### Example: Using IntegrateSystem Directly

```mathematica
(* Get general series solutions without fixing BCs *)
generalSolutions = IntegrateSystem[<|t -> -1/x, s -> 3|>];

(* Get solutions with BCs fixed *)
fixedSolutions = IntegrateSystem[preparedBCs, <|t -> -1/x, s -> 3|>];
```

### Example: Unequal Mass with Higher Precision

```mathematica
LoadConfiguration[{
  ChopPrecision -> 500,
  DivisionOrder -> 4,
  ExpansionOrder -> 70,
  MatrixDirectory -> "Banana_Matrices/",
  RadiusOfConvergence -> 10,
  UseMobius -> True,
  UsePade -> True,
  WorkingPrecision -> 1000,
  IntegrationStrategy -> "Default",   (* or "VOP", "VOPAlt" *)
  Verbosity -> 2
}];
```

---

## Configuration Reference

### Precision-Related Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `WorkingPrecision` | Integer | 300 | Internal arbitrary-precision digit count |
| `ChopPrecision` | Integer | 100 | Digits below which intermediate values are set to zero |
| `LinearSolveChopPrecision` | Integer | (=ChopPrecision) | Separate precision for linear solves (advanced) |
| `AccuracyGoal` | Integer | -- | Target accuracy digits for transport (enables dynamic order) |

**Constraint:** `ChopPrecision` must be strictly less than `WorkingPrecision`.

### Expansion and Order Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `EpsilonOrder` | Integer | 4 | Highest power of epsilon to compute |
| `ExpansionOrder` | Integer | 50 | Maximum series order in x |
| `DivisionOrder` | Integer | 2 | Segment size = RadiusOfConvergence/DivisionOrder |
| `RadiusOfConvergence` | Number | -- | Rescaling factor for line segments |

### Strategy Options

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `IntegrationStrategy` | "Default", "VOP", "VariationOfParameters", "VOPAlt" | "Default" | Which DE solver to use |
| `SegmentationStrategy` | "Dynamic", "Predivision" | "Predivision" | How to divide the transport path |
| `UseMobius` | True/False | True | Use Mobius transformations for segments |
| `UsePade` | True/False | True | Use Pade approximants for evaluation |

### Error and Diagnostics Options

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `EstimateError` | False, "Fast" | "Fast" | Error estimation mode |
| `CrosscheckLevel` | Integer | 0 | Level of cross-checking (higher = more checks) |
| `CrosscheckFlags` | List | {} | Specific checks to enable |
| `Verbosity` | 1, 2, 3 | 1 | Output verbosity level |
| `LogFile` | String | -- | Path to log file |

### Physics Options

| Option | Type | Description |
|--------|------|-------------|
| `MatrixDirectory` | String | Directory containing `dvar_k.m` matrix files |
| `Variables` | List | Kinematic variables and masses (auto-detected from matrices) |
| `DeltaPrescriptions` | List | Branch cut prescriptions for physical singularities |
| `LineParameter` | Symbol | Symbol for line parametrization (default: `x`) |

---

## Error Handling and Diagnostics

### Common Errors

| Error Message | Cause | Solution |
|---------------|-------|----------|
| "ChopPrecision should be smaller than WorkingPrecision" | Invalid precision settings | Ensure ChopPrecision < WorkingPrecision |
| "Line does not fix all kinematic variables" | Incomplete line specification | Add all variables from `Variables` to the line |
| "Line lies on a singularity" | Line passes through a pole | Choose a different parametrization |
| "Boundary conditions cannot be matched" | Inconsistent BCs | Check BC expressions; try adjusting ChopPrecision |
| "Cross-check failed" | Numerical instability | Increase WorkingPrecision or decrease ExpansionOrder |
| "Error reducing indeterminates" | VOP/VOPAlt instability | Decrease ChopPrecision or increase WorkingPrecision |
| "Analytic continuation failed" | Cannot handle branch cut | Add appropriate DeltaPrescriptions |

### Debugging

- Set `Verbosity -> 3` for detailed debug output
- Access `DiffExp`State`LastErrorContext` after an error for diagnostic data
- Access `DiffExp`State`BenchmarkData` after a computation for timing information
- Enable specific cross-checks: `UpdateConfiguration["CrosscheckFlags" -> {"PeriodMatrix", "GeneralSolution"}]`

### Available Cross-check Flags

| Flag | What it checks |
|------|---------------|
| `"Wronskians"` | dW/dx = NMat . W |
| `"WronskInv"` | W_inv . W = Identity |
| `"PeriodMatrix"` | dFMat/dx = A[0] . FMat |
| `"GeneralSolutionMatrix"` | dGMat/dx = A[0] . GMat + BMat |
| `"GeneralSolution"` | df/dx = A[0] . f + bVec (after fixing BCs) |
| `"SingularityCheck"` | No unexpected multivalued functions |

### Performance Tips

1. **Start with lower precision** for initial runs; increase once results look correct
2. **Use "Predivision" strategy** for predictable segment counts
3. **Use Pade approximants** (`UsePade -> True`) for faster convergence between segments
4. **Use Mobius transformations** (`UseMobius -> True`) to handle singularities gracefully
5. **Buffer reuse:** The `BufferedData` mechanism means higher epsilon orders are much cheaper than eps^0
6. **AccuracyGoal with "Before" validation:** Dynamically adjusts expansion order per segment
7. **SaveExpansions:** Use `save=True` and `ToPiecewise` when you need values at many points
