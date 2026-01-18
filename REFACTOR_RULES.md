# DiffExp Refactoring Rules

**CRITICAL: READ THIS FILE AT THE START OF EVERY ITERATION**

This file contains critical rules for refactoring the DiffExp.m Mathematica package. These rules MUST be followed exactly to ensure a successful 1-1 refactor.

## Core Principles

1. **NEVER rename variables or functions** - Keep all names exactly as they are in the original code
2. **NEVER change functionality** - This is a 1-1 refactor; the behavior must be identical
3. **NEVER add or remove Attributes** - Attributes like `SetAttributes[{SExpand,SN,...},Listable]` are critical
4. **NEVER duplicate functions** - Each function must be defined exactly once across all subpackages
5. **NEVER take shortcuts** - Even if code seems redundant, preserve it exactly

## Symbol Namespacing (CRITICAL)

The following variables are used throughout the code and MUST be handled carefully:
- `x` - Line parameter (set to `Global`x` by default, can be changed via LineParameter config)
- `eps` / `\[Epsilon]` - Dimensional regulator
- `Logx` - Symbol for Log[x]
- `\[Theta]p` - HeavisideTheta[x]
- `\[Theta]m` - HeavisideTheta[-x]
- `y` - Temporary variable used in some integration rules

**Key insight from original code:**
```mathematica
LineParameterVal=FEC[LineParameter];
x=FEC[LineParameter];
```
The internal `x` is set to `FEC[LineParameter]` which defaults to `Global`x`.

**Solution: Create a Symbols subpackage that:**
1. Defines these variables
2. Is imported by ALL other subpackages
3. Uses delayed rules (`:=`) so replacements work correctly

## Variable `y` Warning

Check if `y` is used temporarily within functions or carried between functions. This is especially important for integration rules.

## Function Disambiguation

### GetLargestTerm
There may be two versions of this function:
1. One for error checking (looks at largest power term in series)
2. Another version for different purpose
**Action: Verify both usages and ensure correct one is used in each context**

### GetMinSeriesPower
This uses multiplication not division - it gives position in a list. Keep exactly as-is.

## Subpackage Structure

Planned subpackages:
1. **Symbols** - Core variables (x, eps, Logx, theta, etc.)
2. **State** - Configuration and state management (LoadConfiguration, UpdateConfiguration, all state variables)
3. **Utilities** - Small helper functions (SExpand, SApply, PChop, etc.)
4. **SeriesOps** - Series operations (DiffExpSeries, SeriesAlways, LeadingCoefficientSeries, etc.)
5. **Integration** - DiffExpIntegrate, UpdateIntReps, integration rules
6. **Pade** - Pade approximant functions (GetPade, etc.)
7. **Mobius** - Mobius transformation functions (GetMobius, etc.)
8. **LineSegmentation** - Line segmentation strategies
9. **MatrixLoading** - LoadMatrices and matrix preparation functions
10. **Frobenius** - Frobenius solutions (Frobenius1, FrobeniusSolutions)
11. **Wronskian** - Wronskian and F-matrix methods
12. **VOP** - Variation of parameters methods
13. **AnalyticContinuation** - Analytic continuation logic
14. **ErrorEstimates** - Error estimation functions
15. **Transport** - IntegrateSystem and TransportTo (these are complex and interconnected)

## State Variables to Track

### User-configurable (in DefaultConfiguration):
- AccuracyGoal, ChopPrecision, DeltaPrescriptions, DivisionOrder
- EpsilonOrder, ExpansionOrder, IntegrationStrategy, LineParameter
- LogFile, MatrixDirectory, RadiusOfConvergence, SegmentationStrategy
- UseMobius, UsePade, Variables, Verbosity, WorkingPrecision
- And many internal options like "CrosscheckLevel", "EstimateError", etc.

### Internal state variables:
- AnalyticContinuationFailed, AnalyticContinuationReplacements
- AnalyticContinuationReplacementsAssociation
- BenchmarkData, CurrentSingularityWasAddedFromSquareRoot
- CurrentSingularityHasIDeltaPrescription, DEqnSquareRoots
- MultivaluedFail, UserDeltaPrescriptions, UsingClosedFormMatrix
- DEqnMatricesFactored, DEqnMatricesFactoredClosedForm
- DEqnMatricesExpanded, NumIntegrals
- ExpansionOrderVal, MaxCouplingOrder
- IMaxLogOrder, IntReps
- And many more

## Complex Functions Requiring Detailed Plans

1. **LoadMatrices** - Processes matrices and sets state variables
2. **IntegrateSystem** - Checks if called by TransportTo for analytic continuation decisions
3. **TransportTo** - Complex function handling transport and error estimates
4. **PrepareAnalyticContinuation** - Sets analytic continuation variables

## Testing

Primary test: Equal mass banana example
```mathematica
SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
Get[FileNameJoin[{ParentDirectory[scriptDir], "DiffExp.m"}]];
LoadConfiguration[{
  DeltaPrescriptions -> {t - 16 + I * \[Delta]},
  MatrixDirectory -> scriptDir <> "/Banana_EqualMass_Matrices/",
  Verbosity -> 2, UseMobius -> True, UsePade -> True
}];
(* ... rest of test *)
```

Final test: Unequal mass banana example (in Examples folder)

## Integration Rules

The series integration functions (DiffExpIntegrate, IntReps) have rules for powers of Logx. Higher powers update the rules dynamically. These are precomputed for performance. Test thoroughly by comparing with Mathematica's Integrate on monomials of x and Logx.

## File Structure

```
diffexp_refactor/
├── DiffExp.m                    # Main package (loads subpackages)
├── DiffExpSubpackages/
│   ├── Symbols.m
│   ├── State.m
│   ├── Utilities.m
│   ├── SeriesOps.m
│   ├── Integration.m
│   ├── Pade.m
│   ├── Mobius.m
│   ├── LineSegmentation.m
│   ├── MatrixLoading.m
│   ├── Frobenius.m
│   ├── Wronskian.m
│   ├── VOP.m
│   ├── AnalyticContinuation.m
│   ├── ErrorEstimates.m
│   └── Transport.m
├── Tests/
│   ├── test_symbols.m
│   ├── test_integration.m
│   └── test_banana.m
├── Docs/
│   ├── MatrixLoading_LogicFlow.md
│   ├── AnalyticContinuation_LogicFlow.md
│   └── TransportTo_LogicFlow.md
└── REFACTOR_RULES.md            # This file
```

## Order of Operations

1. Create Symbols subpackage and test namespacing
2. Create State subpackage with all configuration/state variables
3. Create Integration subpackage with tests
4. Create remaining utility subpackages
5. Document complex functions (MatrixLoading, AnalyticContinuation, TransportTo)
6. Refactor complex functions based on documentation
7. Create main DiffExp.m that loads all subpackages
8. Test with equal mass banana
9. Test with unequal mass banana

## Completion Criteria

The refactor is complete when:
1. All functions are in appropriate subpackages
2. No function is duplicated
3. All variable names are preserved
4. All Attributes are preserved
5. Equal mass banana example produces identical results
6. Unequal mass banana example produces identical results
