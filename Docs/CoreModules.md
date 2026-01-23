# Core Modules: Symbols, State, and Utilities

The three foundational subpackages of DiffExp that all other modules depend on. They provide the shared symbol definitions, global configuration and mutable state, and a library of small helper functions used throughout the codebase.

## Dependency Chain

```
Symbols.m  (no dependencies)
    ^
    |
State.m    (imports Symbols)
    ^
    |
Utilities.m (imports Symbols, State)
    ^
    |
[All other DiffExp subpackages]
```

All other DiffExp subpackages (SeriesOps, Integration, Frobenius, Transport, etc.) import these three modules to ensure consistent symbol usage, access to configuration, and availability of helper functions.

---

## Symbols.m (`DiffExp`Symbols``)

**File:** `DiffExp/Symbols.m` (26 lines)

Defines the core symbols used throughout DiffExp. This is the lowest-level subpackage with no dependencies. Every other subpackage imports it to guarantee that the same symbol objects are used consistently across the entire package.

### Exported Symbols

| Symbol | Type | Description |
|--------|------|-------------|
| `\[Epsilon]` | Symbol | The dimensional regulator in dimensional regularization. |
| `eps` | Delayed alias | Alternative ASCII name for `\[Epsilon]`. Defined as `eps := \[Epsilon]`. |
| `Logx` | Symbol | Represents `Log[x]`. Used in series expansions to track logarithmic terms. |
| `\[Theta]p` | Symbol | Represents `HeavisideTheta[x]`, the Heaviside step function for positive argument. |
| `\[Theta]m` | Symbol | Represents `HeavisideTheta[-x]`, the Heaviside step function for negative argument. |
| `x` | Symbol | Internal line parameter symbol. Set to the user-specified `LineParameter` by State.m at initialization. |

### Design Notes

- `eps` is defined with delayed evaluation (`:=`), so it always evaluates to the current `\[Epsilon]` symbol. This allows users to write `eps` in their input while the internal code consistently uses `\[Epsilon]`.
- `x` is initially just a symbol placeholder. When State.m loads, it assigns `DiffExp`Symbols`x = FEC[LineParameter]`, linking it to the user-configured line parameter (defaults to `Global`x`).
- `\[Theta]p` and `\[Theta]m` are used in analytic continuation to encode branch cut prescriptions.

### Usage Example

```mathematica
(* Other subpackages import Symbols to use the shared symbols *)
BeginPackage["DiffExp`MyModule`", {"DiffExp`Symbols`"}];

(* Now \[Epsilon], x, Logx, etc. all refer to the same objects *)
myFunction[expr_] := Coefficient[expr, \[Epsilon], 2];
```

---

## State.m (`DiffExp`State``)

**File:** `DiffExp/State.m` (235 lines)

Handles all configuration management and mutable state for DiffExp. This module defines the default configuration, provides accessor functions for configuration values, and maintains the global state variables that track the progress and intermediate results of computations.

**Imports:** `DiffExp`Symbols``

### Configuration System

#### DefaultConfiguration

An `Association` containing all default configuration values:

```mathematica
DefaultConfiguration = <|
  AccuracyGoal -> "?",
  "AccuracyGoalValidate" -> "Before",
  ChopPrecision -> 250,
  "CrosscheckLevel" -> 0,
  DeltaPrescriptions -> {},
  DivisionOrder -> 3,
  EpsilonOrder -> 4,
  "EstimateError" -> "Fast",
  ExpansionOrder -> 50,
  "IgnoreIndicialCheck" -> False,
  "InvWronskSolver" -> "Auto",
  "KeepMatrixExpansions" -> False,
  "LinearSolveChopPrecision" -> 250,
  LineParameter -> Global`x,
  MatrixDirectory -> "",
  RadiusOfConvergence -> 1,
  RationalizationTolerance -> 10^-10,
  "SaveExpansionsCompress" -> False,
  "HomogeneousSolve" -> "Expand",
  "Parallel" -> False,
  SegmentationStrategy -> "Predivision",
  IntegrationStrategy -> "Default",
  UseMobius -> False,
  UsePade -> False,
  Variables -> {},
  Verbosity -> 1,
  "VerbosityDebug" -> 0,
  WorkingPrecision -> 500
|>
```

#### Configuration Store

| Symbol | Description |
|--------|-------------|
| `DiffExpConfiguration` | The main mutable configuration association. Initialized to `DefaultConfiguration`. |
| `FEC` | Shorthand alias for `DiffExpConfiguration` (delayed: `FEC := DiffExpConfiguration`). Used extensively throughout the codebase. |

### Configuration Keys (Exported Symbols with Usage)

These are the user-facing configuration option names:

| Key | Default | Description |
|-----|---------|-------------|
| `ChopPrecision` | 250 | Number of zeros after the decimal point below which intermediate values are set to 0. |
| `DeltaPrescriptions` | `{}` | List of polynomials containing `+/- I*delta` factors, describing singularities like thresholds or branch points. |
| `DivisionOrder` | 3 | Inverse distance to the nearest singularity for segment evaluation in the predivision strategy. |
| `EpsilonOrder` | 4 | Highest power of `\[Epsilon]` to compute in the result. |
| `ExpansionOrder` | 50 | Maximum power of the line parameter kept in intermediate series expansions. |
| `IntegrationStrategy` | `"Default"` | How the differential equations are solved. Values: `"Default"`, `"VOP"` (Variation of Parameters). |
| `LineParameter` | `Global`x` | The parameter variable used for parsing lines. |
| `LogFile` | (none) | File path for logging session output. |
| `MatrixDirectory` | `""` | Directory containing the partial derivative matrices on disk. |
| `RadiusOfConvergence` | 1 | Rescaling factor for the line parameter on each segment. |
| `RationalizationTolerance` | `10^-10` | Tolerance for rationalizing exponents in singularity decomposition. |
| `SegmentationStrategy` | `"Predivision"` | Which segmentation strategy to use when dividing integration paths. |
| `UseMobius` | `False` | Whether to use Mobius transformations (vs. linear) for line segments. |
| `UsePade` | `False` | Whether to use Pade approximants during boundary condition transport. |
| `Verbosity` | 1 | Level of printed output (0 = silent, higher = more verbose). |

**Note on System symbols:** `AccuracyGoal`, `Variables`, and `WorkingPrecision` are `System`` context symbols. They are used as configuration keys via `System`AccuracyGoal`, `System`Variables`, and `System`WorkingPrecision` to avoid shadowing the built-in symbols.

### Value Accessors

Delayed-evaluation accessors that always return the current configuration value:

| Accessor | Reads From | Description |
|----------|-----------|-------------|
| `ChopPrecisionVal` | `FEC[ChopPrecision]` | Current chop precision. |
| `LinearSolveChopPrecisionVal` | `FEC["LinearSolveChopPrecision"]` | Chop precision for linear solves. |
| `CrosscheckChopPrecision` | (constant 30) | Precision for crosscheck comparisons. |
| `ExternalScalesVal` | `FEC[Variables]` | List of external kinematic variables. |
| `LineParameterVal` | `FEC[LineParameter]` | Current line parameter symbol. |
| `MatrixDirectoryVal` | `FEC[MatrixDirectory]` | Current matrix directory path. |
| `EpsilonOrderVal` | `FEC[EpsilonOrder]` | Current epsilon truncation order. |
| `FEAccuracyGoal` | `FEC[AccuracyGoal]` | Current accuracy goal setting. |
| `FEWorkingPrecision` | `FEC[WorkingPrecision]` | Current working precision (default 500). |
| `DeltaPrescriptionsVal` | `FEC[DeltaPrescriptions]` | Current delta prescriptions list. |
| `UseMobiusVal` | `FEC[UseMobius]` | Whether Mobius transforms are enabled. |
| `RadiusOfConvergenceVal` | `FEC[RadiusOfConvergence]` | Current radius of convergence. |
| `DivisionOrderVal` | `FEC[DivisionOrder]` | Current division order. |
| `ExpansionOrderVal` | (mutable variable) | Current expansion order (may change during computation). |
| `MaxCouplingOrder` | (mutable variable) | Maximum coupling order (default 1). |

### Mutable State Variables

These variables track the state of an ongoing computation:

#### Differential Equation Matrices

| Variable | Initial Value | Description |
|----------|--------------|-------------|
| `DEqnMatricesFactored` | `<||>` | Factored form of loaded DE matrices. |
| `DEqnMatricesFactoredClosedForm` | `<||>` | Closed-form factored matrices. |
| `DEqnMatricesExpanded` | `<||>` | Series-expanded DE matrices. |
| `NumIntegrals` | 0 | Number of master integrals in the system. |
| `UsingClosedFormMatrix` | `False` | Whether the current matrix is in closed form. |
| `DEqnSquareRoots` | `{}` | Square roots appearing in the DE matrices. |

#### Integration State

| Variable | Initial Value | Description |
|----------|--------------|-------------|
| `IntegrationSequence` | (uninitialized) | Order in which variables are integrated. |
| `ExpansionMatrices` | (uninitialized) | Matrices used for series expansion. |
| `ExpansionMatricesCanonical1` | (uninitialized) | Canonical expansion matrices. |
| `ExpansionMatricesClosedForm` | (uninitialized) | Closed-form expansion matrices. |

#### Analytic Continuation

| Variable | Initial Value | Description |
|----------|--------------|-------------|
| `AnalyticContinuationFailed` | `False` | Flag set when analytic continuation encounters an error. |
| `AnalyticContinuationReplacements` | `{}` | List of replacement rules for continuation. |
| `AnalyticContinuationReplacementsAssociation` | `<||>` | Association form of continuation replacements. |
| `CurrentSingularityWasAddedFromSquareRoot` | `False` | Whether the current singularity comes from a square root. |
| `CurrentSingularityHasIDeltaPrescription` | `False` | Whether the current singularity has an i-delta prescription. |
| `MultivaluedFail` | `False` | Flag for multivalued function failure detection. |
| `UserDeltaPrescriptions` | `{}` | User-specified delta prescriptions (parsed from configuration). |

#### Alphabet and Logging

| Variable | Initial Value | Description |
|----------|--------------|-------------|
| `AlphabetLogs` | (uninitialized) | Logarithms forming the symbol alphabet. |
| `AlphabetLogRules` | (uninitialized) | Replacement rules for alphabet logs. |
| `AlphabetLogRulesFactored` | (uninitialized) | Factored form of log rules. |
| `AlphabetLogRulesExpanded` | (uninitialized) | Expanded form of log rules. |
| `MatricesIrreducibleFactors` | (uninitialized) | Irreducible polynomial factors found in the matrices. |

#### Wronskian and Buffered Data

| Variable | Initial Value | Description |
|----------|--------------|-------------|
| `MyWronsk` | (uninitialized) | Storage for computed Wronskian matrices. |
| `MyWronskDetInv` | (uninitialized) | Inverse of the Wronskian determinant. |
| `BufferedData` | (uninitialized) | General-purpose buffer for intermediate computation data. |
| `BenchmarkData` | `<||>` | Timing data for performance benchmarking. |

#### Crosscheck and Debugging

| Variable | Initial Value | Description |
|----------|--------------|-------------|
| `CurrCrosscheckFlags` | `{}` | Currently active crosscheck flags for the computation. |
| `LogStream` | (uninitialized) | Stream object for log file output. |
| `LastErrorContext` | `{}` | Diagnostic context from the last error for debugging. |
| `DiffExpExtensions` | `{}` | List of loaded extensions. |

### Crosscheck Flags

The `CrosscheckFlags` association defines which validation checks are performed and at what crosscheck level:

```mathematica
CrosscheckFlags = {
  "FrobeniusSolutions" -> 1,
  "MatrixDelta" -> 1,
  "Wronskians" -> 1,
  "WronskInv" -> 0,
  "PeriodMatrix" -> 1,
  "GeneralSolutionMatrix" -> 2,
  "GeneralSolution" -> 1,
  "VariationOfParameters" -> 1,
  "SingularityCheck" -> 0
}
```

A check is performed when the `"CrosscheckLevel"` configuration value is >= the flag value. Flags with value 0 are always performed.

### Internal Constants

Tuning parameters for the adaptive series expansion algorithm:

| Constant | Value | Description |
|----------|-------|-------------|
| `ISeriesChangeCoefficient` | 2 | Multiplier for series order adjustments. |
| `IMaxLogOrder` | 1 | Maximum logarithmic order (initialized to default). |
| `IMaxLogOrderDefault` | 1 | Default maximum log order. |
| `ICheckMultivaluedChop` | 5 | Precision for multivalued function detection. |
| `ICrossCheckPrintResultOrder` | 5 | Epsilon order shown in crosscheck output. |
| `ICrossCheckVerifyResultOrder` | 5 | Epsilon order used for crosscheck verification. |
| `ISafetyDigits` | 2 | Extra digits of precision kept as a safety margin. |
| `ISafetyExpansionSubtract` | 5 | Safety subtraction from expansion order. |
| `IExpansionOrdersAveraging` | 3 | Number of orders used for averaging convergence estimates. |
| `IExpansionOrderIncrease` | 10 | Amount to increase expansion order when convergence is slow. |
| `IExpansionOrderDecrease` | 10 | Amount to decrease expansion order when precision is excess. |
| `IExpansionOrderIncrease2` | 25 | Larger expansion order increase for difficult cases. |
| `IDigitsSurplusDecreaseExpansionOrder` | 3 | Surplus digits threshold for decreasing order. |
| `ICurrEvalErrorSeriesDecrease` | `Ceiling[0.7*MaxCouplingOrder] + 2` | Dynamic decrease based on coupling order. |
| `IDecreaseOrderByErrorPrecise` | `MaxCouplingOrder` | Decrease tied to coupling order. |
| `IMinExpansionOrder` | 10 | Floor for the expansion order (never goes below this). |

### Helper Function

#### SquareRootPrescriptionsAdded

```mathematica
SquareRootPrescriptionsAdded[]
```

Returns a list of `{polynomial, 1}` pairs for square roots found in the DE matrices that do not already have user-specified delta prescriptions. Used to automatically add i-delta prescriptions for square-root branch points.

### Usage Example

```mathematica
(* Loading and configuring DiffExp *)
Get["DiffExp.m"];
LoadConfiguration[{
  EpsilonOrder -> 6,
  WorkingPrecision -> 1000,
  ChopPrecision -> 400,
  DeltaPrescriptions -> {t - 16 + I*delta},
  MatrixDirectory -> "path/to/matrices/"
}];

(* Accessing configuration in internal code *)
prec = FEWorkingPrecision;           (* returns 1000 *)
chopPrec = ChopPrecisionVal;         (* returns 400 *)
epsOrd = EpsilonOrderVal;            (* returns 6 *)
```

---

## Utilities.m (`DiffExp`Utilities``)

**File:** `DiffExp/Utilities.m` (124 lines)

A collection of small helper functions used throughout DiffExp. Provides printing/logging utilities, shorthand aliases, numerical chopping functions, geometric predicates for lines and intervals, and expression splitting helpers.

**Imports:** `DiffExp`Symbols``, `DiffExp`State``

### Printing Functions

#### PrintDebug

```mathematica
PrintDebug[args__][lev_]
```

Prints debug messages when the current `"VerbosityDebug"` configuration value is >= `lev`. Messages are indented by `lev - 1` levels and prefixed with `"DiffExp Debug: "`.

**Example:**
```mathematica
PrintDebug["Matrix size: ", n][2];
(* Prints "      DiffExp Debug: Matrix size: 5" if VerbosityDebug >= 2 *)
```

#### PrintInfo

```mathematica
PrintInfo[args__][lev_]
```

Prints informational messages when the current `Verbosity` configuration value is >= `lev`. Messages are indented by `lev - 1` levels and prefixed with `"DiffExp: "`.

**Example:**
```mathematica
PrintInfo["Starting integration..."][1];
(* Prints "DiffExp: Starting integration..." if Verbosity >= 1 *)
```

#### PrintWarning

```mathematica
PrintWarning[args__]
```

Unconditionally prints a warning message prefixed with `"DiffExp Warning: "`.

#### ReportError

```mathematica
ReportError[mes__]
ReportError[mes__, False]
```

Reports an error prefixed with `"DiffExp Error: "`. The first form also calls `Abort[]` to halt execution. The second form (with `False` as the last argument) prints the error but does not abort.

#### PrintMobiusNormalized

```mathematica
PrintMobiusNormalized[a_]
```

Normalizes a rational function of `x` for human-readable display. Divides both numerator and denominator by the largest coefficient so the printed form has coefficients of order 1.

### General Utility Functions

#### AllSameQ

```mathematica
AllSameQ[l_, b_]
```

Returns `True` if all elements of list `l` are identical to `b`. Returns `True` for empty lists.

**Example:**
```mathematica
AllSameQ[{0, 0, 0}, 0]   (* True *)
AllSameQ[{0, 1, 0}, 0]   (* False *)
AllSameQ[{}, 0]           (* True *)
```

#### CA

```mathematica
CA
```

Alias for `ConstantArray`. Shorthand used throughout the codebase.

**Example:**
```mathematica
CA[0, {3, 3}]   (* {{0,0,0},{0,0,0},{0,0,0}} *)
```

#### GetCases

```mathematica
GetCases[expr_, case_]
```

Extracts all subexpressions matching `case` from `expr` at any depth (`Infinity` level). Returns a sorted, deduplicated list.

**Example:**
```mathematica
GetCases[a + b*x + c*x^2, x]   (* {x} *)
GetCases[Sin[x] + Cos[y], _Symbol]   (* {x, y} *)
```

#### DependsQ

```mathematica
DependsQ[a_, b_]
```

Returns `True` if expression `a` contains `b` as a subexpression (at any depth). Implemented via `GetCases`.

**Example:**
```mathematica
DependsQ[x^2 + y, x]   (* True *)
DependsQ[y + z, x]     (* False *)
```

#### ZeroQ

```mathematica
ZeroQ[a_]
```

Returns `True` if `a` is identically zero (using `===`). This is an exact check, not a numerical one.

#### R

```mathematica
R
```

Alias for `ReplaceAll`. Shorthand used throughout the codebase.

**Example:**
```mathematica
(x^2 + x) // R[x -> 3]   (* 12 *)
```

#### FirstOrNull

```mathematica
FirstOrNull[l_]
```

Returns the first element of list `l`, or `Null` if the list is empty. Avoids the error that `First[{}]` would produce.

#### FindPivots

```mathematica
FindPivots[Matrix_]
```

Finds the column position of the first nonzero entry in each row of `Matrix`. Returns a flat list of column indices (rows that are all-zero are excluded via `DeleteCases[..., Null]`).

**Example:**
```mathematica
FindPivots[{{0, 1, 2}, {3, 0, 0}, {0, 0, 0}}]   (* {2, 1} *)
```

#### SplitTimes

```mathematica
SplitTimes[Expr_]
```

Splits a `Times` expression into a list of its factors. If the expression is not a product, wraps it in a singleton list. Has the `Listable` attribute, so it threads over lists automatically.

**Example:**
```mathematica
SplitTimes[a*b*c]   (* {a, b, c} *)
SplitTimes[a + b]   (* {a + b} *)
```

#### SplitSum

```mathematica
SplitSum[Expr_]
```

Splits a `Plus` expression into a list of its summands. If the expression is not a sum, wraps it in a singleton list.

**Example:**
```mathematica
SplitSum[a + b + c]   (* {a, b, c} *)
SplitSum[a*b]         (* {a*b} *)
```

### Chopping Functions

These functions zero out numerical values below a precision threshold. They are essential for high-precision arithmetic where tiny numerical noise must be removed.

#### PChop

```mathematica
PChop
```

A pure function that chops values below `10^(-ChopPrecision)`. Uses the current `ChopPrecision` from the configuration (default: `10^-250`).

**Example:**
```mathematica
{1.0, 10^-300, 10^-200} // PChop   (* {1.0, 0, 10^-200} *)
```

#### LSPChop

```mathematica
LSPChop
```

A pure function that chops values below `10^(-LinearSolveChopPrecision)`. Used after solving linear systems where different precision may be appropriate. Default threshold: `10^-250`.

#### CPChop

```mathematica
CPChop
```

A pure function that chops values below `10^(-CrosscheckChopPrecision)`. Uses a fixed precision of 30 digits, appropriate for crosscheck comparisons where less precision is needed. Default threshold: `10^-30`.

### Line/Point Detection

#### IsPoint

```mathematica
IsPoint[line_]
```

Returns `True` if `line` (an Association of variable -> expression mappings) represents a single point (none of the values depend on `x`).

**Example:**
```mathematica
IsPoint[<|t -> 5, s -> 3|>]       (* True *)
IsPoint[<|t -> 1 + x, s -> 3|>]   (* False *)
```

#### IsLine

```mathematica
IsLine[line_]
```

Returns `True` if `line` depends on `x` (i.e., it is not a point). Complement of `IsPoint`.

### Interval Functions

#### IntervalOverlapQ

```mathematica
IntervalOverlapQ[intv1_, intv2_]
```

Returns `True` if the two intervals (each given as `{min, max}`) have a non-empty intersection. Uses Mathematica's built-in `IntervalIntersection`.

**Example:**
```mathematica
IntervalOverlapQ[{0, 1}, {0.5, 1.5}]   (* True *)
IntervalOverlapQ[{0, 1}, {2, 3}]        (* False *)
```

#### IntervalIntersec

```mathematica
IntervalIntersec[intv1_, intv2_]
```

Returns the intersection of two intervals as a `{min, max}` pair. Assumes the intervals do overlap (use `IntervalOverlapQ` first to verify).

**Example:**
```mathematica
IntervalIntersec[{0, 1}, {0.5, 1.5}]   (* {0.5, 1} *)
```

#### IntervalContainsQ

```mathematica
IntervalContainsQ[intv_, point_]
```

Returns `True` if the interval `intv = {min, max}` contains `point` (inclusive on both ends).

**Example:**
```mathematica
IntervalContainsQ[{0, 1}, 0.5]   (* True *)
IntervalContainsQ[{0, 1}, 1.5]   (* False *)
```

### Line Type Checking

#### ExactLineQ

```mathematica
ExactLineQ[line_Association]
```

Returns `True` if all coefficients of the line's rational functions (numerator and denominator coefficients with respect to `x`) are exact numbers (integers, rationals, or algebraic numbers). Used to determine whether `Factor` can be used instead of the more general `Together`.

**Example:**
```mathematica
ExactLineQ[<|t -> (1 + 2*x)/(3 + 4*x)|>]         (* True *)
ExactLineQ[<|t -> (1.5 + 2*x)/(3 + 4*x)|>]       (* False, 1.5 is inexact *)
```

#### FactorOrTogether

```mathematica
FactorOrTogether[line_Association]
FactorOrTogether[line1_Association, line2_Association]
```

Returns `Factor` if the line(s) have exact coefficients, or `Together` otherwise. `Factor` produces nicer results but only works on exact expressions. The two-argument form requires both lines to be exact.

**Example:**
```mathematica
simplify = FactorOrTogether[myLine];
result = simplify[expr];
```

---

## Cross-Module Interactions

### Symbol Resolution

All modules refer to `\[Epsilon]` and `x` via the `DiffExp`Symbols`` context. This ensures that when a user writes a boundary condition involving `eps`, it resolves to the same symbol object used in the internal computation.

### Configuration Flow

```
User code                     State.m                    Other modules
---------                     -------                    -------------
LoadConfiguration[...]  -->   FEC = <|...|>
                                  |
                                  v
                              ChopPrecisionVal  ------>  PChop (Utilities.m)
                              EpsilonOrderVal   ------>  Frobenius, Transport
                              FEWorkingPrecision -----> SetPrecision calls
                              LineParameterVal  ------>  x (Symbols.m)
```

### Chopping Pipeline

The three chop functions in Utilities.m correspond to three levels of numerical cleaning:

1. **PChop** (default `10^-250`): Standard intermediate computation cleanup. Most matrix operations and series coefficients are chopped at this level.
2. **LSPChop** (default `10^-250`): Applied after `LinearSolve` operations. Can be configured independently if linear system solutions need different treatment.
3. **CPChop** (fixed `10^-30`): Used only during crosscheck comparisons, where the goal is to verify agreement to a modest number of digits rather than full working precision.

### Printing Verbosity Levels

The `PrintInfo` and `PrintDebug` functions use a level-based system:

- **Level 1**: Top-level progress messages (starting integration, loading matrices, etc.)
- **Level 2**: Per-segment or per-step information
- **Level 3+**: Detailed internal algorithm steps

These are controlled by the `Verbosity` (for `PrintInfo`) and `"VerbosityDebug"` (for `PrintDebug`) configuration keys.

---

## Initialization Order

When DiffExp is loaded via `Get["DiffExp.m"]`, the subpackages are loaded in dependency order:

1. **Symbols.m** loads first, creating the symbol objects.
2. **State.m** loads second, importing Symbols. It:
   - Creates the `DefaultConfiguration` association
   - Initializes `DiffExpConfiguration` to the defaults (if not already set)
   - Sets `DiffExp`Symbols`x` to the configured `LineParameter`
   - Initializes all mutable state variables to their default values
3. **Utilities.m** loads third, importing both Symbols and State. It defines helper functions that reference State accessors via delayed evaluation.

This ordering ensures that by the time any computational module (Frobenius, Transport, etc.) loads, all symbols are defined, configuration is available, and utility functions are ready.
