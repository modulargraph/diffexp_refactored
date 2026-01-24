# DiffExp Refactoring Project

## Project Rules

- ONLY modify files within this directory (~/Desktop/diffexp_refactored)
- NEVER access files outside this folder
- Ask before running any destructive commands

## Project Status: REFACTORED (v1.2)

The DiffExp Mathematica package has been successfully refactored from a monolithic file into 16 modular subpackages (~4,850 lines total). All tests pass. The package solves differential equations for Feynman integrals using series expansion methods in dimensional regularization.

## Git Repository

This project now uses git for version control. Always commit your changes with descriptive messages after completing tasks.

## Directory Structure

```
diffexp_refactored/
├── DiffExp.m                          # Root loader (loads DiffExp/DiffExp.m)
├── DiffExp/                           # Refactored modular package (16 subpackages)
│   ├── DiffExp.m                      # Main entry point (197 lines)
│   ├── Symbols.m                      # Core symbols: ε, eps, Logx, θp, θm, x (26 lines)
│   ├── State.m                        # Configuration and state management (235 lines)
│   ├── Utilities.m                    # Helper functions (printing, chopping, etc.) (124 lines)
│   ├── SeriesOps.m                    # Series operations (178 lines)
│   ├── Integration.m                  # Series integration functions (73 lines)
│   ├── Pade.m                         # Padé approximant functions (74 lines)
│   ├── Mobius.m                       # Möbius transformation functions (151 lines)
│   ├── AnalyticContinuation.m         # Analytic continuation handling (107 lines)
│   ├── LineSegmentation.m             # Line segmentation strategies (262 lines)
│   ├── Frobenius.m                    # Frobenius method functions (141 lines)
│   ├── Wronskian.m                    # Wronskian computations (112 lines)
│   ├── MatrixLoading.m                # Matrix loading and processing (389 lines)
│   ├── IntegrationStrategies.m        # Integration strategies package wrapper (38 lines)
│   ├── IntegrationStrategies/          # Strategy implementations
│   │   ├── Helpers.m                   # Shared helpers (ComputeGMat, ComputeWronskianMatrix, etc.)
│   │   ├── Default.m                   # SolveSimple + SolveDefault
│   │   ├── VOP.m                       # SolveVOP + SolveVOPAlt
│   │   ├── Recurrence.m               # Rational + Singular recurrence strategies
│   │   ├── ResonantRecurrence.m        # General singular recurrence (resonant/Jordan) + unified particular solver
│   │   └── Dispatch.m                  # DispatchStrategy routing
│   ├── Transport.m                    # TransportTo, IntegrateSystem, ToPiecewise (1128 lines)
│   ├── SingularityDecomposition.m     # Singularity analysis (232 lines)
│   └── RegularizedIntegration.m       # Regularized integration (724 lines)
├── Reference/                         # Original code and examples for reference
│   ├── DiffExp_original.m             # Original monolithic DiffExp.m
│   └── Examples/                      # Example scripts
│       ├── Banana_example.m
│       ├── MultiplePolylogarithms_example.m
│       ├── FivePointNonPlanar_example.m
│       └── Hypergeometric2F1_example.m
├── Tests/                             # Test scripts and test data
│   ├── Banana_EqualMass_Matrices/     # Equal mass matrices (5 files)
│   ├── Banana_Matrices/               # Unequal mass matrices (25 files)
│   ├── Hypergeometric2F1_Matrices/    # Hypergeometric function matrices (5 files)
│   ├── test_package_loading.m         # Package loading test (19/19 pass)
│   ├── test_symbols_namespacing.m     # Symbol namespacing test (all pass)
│   ├── test_unequal_mass_full.m       # Full transport test: Default vs Recurrence comparison
│   ├── test_banana_refactored.m       # Equal mass banana test
│   ├── test_decomposition.m           # Singularity decomposition tests
│   ├── test_regularized_integration.m # Regularized integration tests
│   ├── test_hypergeometric2f1.m       # Hypergeometric function tests
│   ├── test_resonant_2f1.m           # Resonant 2F1 (c=2, eigenvalues {-2,0})
│   ├── test_resonant_banana.m        # Resonant equal-mass banana (Jordan+resonance)
│   ├── test_multiple_polylogarithms.m # Multiple polylogarithms test
│   ├── test_five_point_nonplanar.m    # Five-point nonplanar test
│   ├── test_topiecewise.m             # Piecewise function conversion test
│   ├── test_feynmantrick_algebra.m    # FeynmanTrick propagator algebra tests
│   ├── test_feynmantrick_fire.m       # FeynmanTrick FIRE6 interface tests
│   ├── test_feynmantrick_iteration.m  # FeynmanTrick iteration/combination tests
│   ├── test_feynmantrick_endtoend.m   # FeynmanTrick end-to-end pipeline tests
│   └── test_feynmantrick_pipeline.m   # Full Feynman trick integration pipeline test
├── Docs/                              # Documentation
│   ├── CoreModules.md                 # Symbols.m, State.m, Utilities.m
│   ├── SeriesAndEvaluation.md         # SeriesOps.m, Integration.m, Pade.m
│   ├── Transformations.md            # Mobius.m, AnalyticContinuation.m, Frobenius.m
│   ├── Infrastructure.md             # LineSegmentation.m, Wronskian.m, MatrixLoading.m
│   ├── TransportAndStrategies.md      # Transport.m, IntegrationStrategies.m, DiffExp.m
│   ├── SingularityDecomposition.md    # Singularity decomposition algorithm
│   ├── RegularizedIntegration.md      # Regularized integration formulas
│   └── Kira.md                        # Kira IBP reduction tool guide
├── FeynmanTrick/                      # Feynman trick integration package
│   ├── FeynmanTrick.m                 # Main loader (loads submodules)
│   ├── PropagatorAlgebra.m            # Symbolic propagator manipulation
│   ├── FIREInterface.m                # FIRE6 IBP reduction interface
│   ├── MatrixExport.m                 # Matrix export to DiffExp format
│   ├── EpsPrefactors.m                # ε-prefactor computation for poles
│   ├── FeynmanTrickIteration.m        # Multi-level iteration logic
│   ├── BoundaryConditions.m           # Deepest-level boundary (generalized tadpole)
│   ├── DiffExpIntegration.m           # DiffExp transport + integration bridge
│   └── IMPLEMENTATION_STATUS.md       # Current pipeline status and known bugs
├── Papers/                            # Reference papers (gitignored)
│   └── FeynmanTrick/                  # Hidding & Usovitsch (2022)
│       └── main.tex                   # "Feynman parameter integration through differential equations"
└── CLAUDE.md                          # This file
```

## Key Refactoring Decisions

### Namespace Management
- **Symbols.m**: Defines core symbols (x, eps, Logx, θp, θm) that are imported by all other packages
- **System symbols**: AccuracyGoal, Variables, WorkingPrecision use `System\`` context explicitly to avoid shadowing
- **Configuration keys**: Other config options (ChopPrecision, UseMobius, etc.) are in `DiffExp\`State\`` context
- **Expected warning**: `Global\`x::shdw` appears because LineParameter defaults to Global\`x - this is normal

### State Management
- All global state is in State.m
- Configuration via LoadConfiguration[], UpdateConfiguration[], CurrentConfiguration[]
- FEC is shorthand for DiffExpConfiguration (the main config association)

### Integration Strategies
- **Default**: Frobenius/Wronskian method with automatic pivot selection. Falls back to VOPAlt if all pivots fail.
- **VOP**: Variation of Parameters method. Derives nth-order scalar ODEs per integral.
- **VOPAlt**: Alternative VOP. Constructs FMat directly from homogeneous solutions. Used as fallback.
- **Recurrence** (UseRationalRecurrence -> True): Direct recurrence on series coefficients. Includes:
  - Rational recurrence for non-singular points
  - Singular recurrence for non-resonant regular singular points
  - General singular recurrence (ResonantRecurrence.m) for resonant/non-diagonalizable cases with unified particular solution solver

### Important Rules
1. Stay focused on the task at hand - don't fix or refactor unrelated code
2. Each function should be defined in exactly one place
3. Use `DiffExp\`Symbols\`x` for the line parameter variable throughout

## FeynmanTrick Package

Implements the method from **Papers/FeynmanTrick/main.tex** (Hidding & Usovitsch, "Feynman parameter integration through differential equations", JHEP 2022). When working on the FeynmanTrick code, read this paper for the full mathematical context.

### Key Formulas (from paper)

**Feynman's trick** (eq. 2.1): Combines two propagators into one parametric integral:
```
1/(D_i^v_i * D_j^v_j) = Gamma(v_i+v_j)/(Gamma(v_i)*Gamma(v_j)) *
    ∫₀¹ x^(v_i-1)*(1-x)^(v_j-1) / [x*D_i + (1-x)*D_j]^(v_i+v_j) dx
```

**Rescaled Feynman parameters** (eq. 2.11):
```
x_1' = prod_{i=1}^{n-1} x_i
x_j' = (1-x_{j-1}) * prod_{i=j}^{n-1} x_i   for j=2,...,n-1
x_n' = (1-x_{n-1})
```

**Generalized tadpole** (eq. 2.16, multiloop-ready):
```
I_v = Gamma(v - L*d/2) / Gamma(v) * U^(v-(L+1)*d/2) / F^(v-L*d/2)
```
where U, F are Symanzik polynomials (computed via FIRE6's UF function), L is loop number, d=4-2ε.

**Differential equation** (eq. 2.9): The combined integral satisfies:
```
d/dx I^(k)_{v1+v2,...} = A(x) * I^(k)_{masters}
```
where A(x) is computed via IBP reductions (FIRE6) of propagator derivatives.

### Architecture

The pipeline works bottom-up through levels:
```
Level N (tadpole) → boundary via generalized tadpole formula
  ↓ DiffExp transport (piecewise series in [0,1])
  ↓ Integration: ∫ x^(v1-1)*(1-x)^(v2-1) * r(x) * f(x) dx
Level N-1 → boundary values
  ↓ ... repeat ...
Level 0 (original topology) → final result
```

### Key Design Decisions
- Variable `xx` (not `x`) for Feynman parameter to avoid DiffExp conflict
- FIRE6's UF function for Symanzik polynomials (multiloop-ready from start)
- Order-by-order ε expansion (not ε-sampling)
- ε-prefactors: define J_i = ε^{k_i} * I_i to remove poles from matrices
- DiffExp config: `UseMobius -> False`, `SaveExpansions -> True` for integration
- Fixed parameter value: 11/23 (avoids singularities)

### Current Test Results (FeynmanTrick)
- test_feynmantrick_algebra.m: **21/21 PASS**
- test_feynmantrick_fire.m: **15/15 PASS**
- test_feynmantrick_iteration.m: **17/17 PASS**
- test_feynmantrick_endtoend.m: Partially passes (boundary bugs)
- test_feynmantrick_pipeline.m: Fails at boundary (DeepestLevelBoundary bug)

### Known Issues
See `FeynmanTrick/IMPLEMENTATION_STATUS.md` for detailed status.
- **FIXED**: `ComputeLevelBoundary` rewritten with correct recursion mapping
- **BUG**: `DeepestLevelBoundary` computes Symanzik polys from wrong topology (should use level-0 original, not combined)
- **BUG**: `TransportLevel` missing `ChopPrecision` in DiffExp config

## Running Tests

```bash
export WOLFRAMSCRIPT_KERNELPATH="/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/WolframKernel"

# Package loading test (quick)
cd Tests && wolframscript -file test_package_loading.m

# Symbol namespacing test (quick)
cd Tests && wolframscript -file test_symbols_namespacing.m

# Full transport test (Default vs Recurrence comparison)
cd Tests && wolframscript -file test_unequal_mass_full.m

# Singularity decomposition test
cd Tests && wolframscript -file test_decomposition.m

# Regularized integration test
cd Tests && wolframscript -file test_regularized_integration.m

# Hypergeometric function test
cd Tests && wolframscript -file test_hypergeometric2f1.m

# Multiple polylogarithms test
cd Tests && wolframscript -file test_multiple_polylogarithms.m

# Five-point nonplanar test (requires internet for arXiv files)
cd Tests && wolframscript -file test_five_point_nonplanar.m

# Piecewise conversion test
cd Tests && wolframscript -file test_topiecewise.m

# FeynmanTrick tests
cd Tests && wolframscript -file test_feynmantrick_algebra.m
cd Tests && wolframscript -file test_feynmantrick_fire.m
cd Tests && wolframscript -file test_feynmantrick_iteration.m
cd Tests && wolframscript -file test_feynmantrick_endtoend.m
cd Tests && wolframscript -file test_feynmantrick_pipeline.m
```

## Running Tests in tmux

For long-running tests, use the tmux session `diffexp-tests` so output can be monitored independently:

```bash
# Start (or attach to) the test session
tmux new-session -d -s diffexp-tests 2>/dev/null || true

# Run a test in the tmux session (example)
tmux send-keys -t diffexp-tests "cd /Users/mhidding/Desktop/diffexp_refactored/Tests && export WOLFRAMSCRIPT_KERNELPATH=\"/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/WolframKernel\" && wolframscript -file test_recurrence_speedup.m 2>&1 | tee /tmp/diffexp_test_output.log" Enter

# Check progress manually
tmux attach -t diffexp-tests        # attach to watch live
# or
wc -l /tmp/diffexp_test_output.log  # check output size before reading
tail -50 /tmp/diffexp_test_output.log  # see recent output
```

**Important:** Before reading test output files, always check their size first with `wc -l` or `ls -lh` to avoid loading excessively large output into context.

## Current Test Results
- test_package_loading.m: **19/19 PASS**
- test_symbols_namespacing.m: **ALL PASS**
- test_unequal_mass_full.m: **3/3 PASS** (Default vs Recurrence: identical results, 2.4x speedup)
- test_hypergeometric2f1.m: **PASS** (50+ digits accuracy)
- test_resonant_2f1.m: **PASS** (52 digits, resonant eigenvalues {-2, 0})
- test_resonant_banana.m: **2/2 PASS** (matches Default to 10^-31)

## Documentation

Detailed documentation for all subpackages is in the `Docs/` folder:
- **CoreModules.md**: Symbols, State (all config keys/defaults), Utilities (all helpers)
- **SeriesAndEvaluation.md**: SeriesOps (SD, SExpand, etc.), Integration (IntReps), Padé (GetPade, SEval)
- **Transformations.md**: Möbius transforms, analytic continuation (θ functions), Frobenius method
- **Infrastructure.md**: Line segmentation, Wronskian (pivot selection), matrix loading (file formats)
- **TransportAndStrategies.md**: TransportTo algorithm, IntegrateSystem, all integration strategies
- **SingularityDecomposition.md**: Decomposing series near singularities into x^(a+bε)·g(x,ε)
- **RegularizedIntegration.md**: Regularization prescription for non-integrable singularities

## Task Instructions

When given a task:
1. Read this prompt to understand the project state
2. Use git to track changes (commit after completing work)
3. Run relevant tests to verify changes don't break anything

## Constraints
- ONLY modify files within ~/Desktop/diffexp_refactored
- NEVER access files outside this folder
- Ask before running any destructive commands
- Use wolframscript (not Mathematica notebooks) for testing
