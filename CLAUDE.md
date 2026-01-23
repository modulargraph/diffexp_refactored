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
│   ├── test_package_loading.m         # Package loading test (18/18 pass)
│   ├── test_symbols_namespacing.m     # Symbol namespacing test (all pass)
│   ├── test_unequal_mass_banana.m     # Simple configuration test (5/5 pass)
│   ├── test_unequal_mass_full.m       # Full transport test (passes)
│   ├── test_banana_refactored.m       # Equal mass banana test
│   ├── test_decomposition.m           # Singularity decomposition tests
│   ├── test_regularized_integration.m # Regularized integration tests
│   ├── test_hypergeometric2f1.m       # Hypergeometric function tests
│   ├── test_multiple_polylogarithms.m # Multiple polylogarithms test
│   ├── test_five_point_nonplanar.m    # Five-point nonplanar test
│   └── test_topiecewise.m             # Piecewise function conversion test
├── Docs/                              # Documentation
│   ├── CoreModules.md                 # Symbols.m, State.m, Utilities.m
│   ├── SeriesAndEvaluation.md         # SeriesOps.m, Integration.m, Pade.m
│   ├── Transformations.md            # Mobius.m, AnalyticContinuation.m, Frobenius.m
│   ├── Infrastructure.md             # LineSegmentation.m, Wronskian.m, MatrixLoading.m
│   ├── TransportAndStrategies.md      # Transport.m, IntegrationStrategies.m, DiffExp.m
│   ├── SingularityDecomposition.md    # Singularity decomposition algorithm
│   ├── RegularizedIntegration.md      # Regularized integration formulas
│   └── Kira.md                        # Kira IBP reduction tool guide
├── Paper/                             # Academic paper
│   ├── main.tex                       # Paper source
│   ├── main.pdf                       # Compiled paper
│   ├── refs.bib                       # Bibliography
│   └── Plots/                         # Generated plots
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

### Important Rules
1. Stay focused on the task at hand - don't fix or refactor unrelated code
2. Each function should be defined in exactly one place
3. Use `DiffExp\`Symbols\`x` for the line parameter variable throughout

## Running Tests

```bash
export WOLFRAMSCRIPT_KERNELPATH="/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/WolframKernel"

# Package loading test (quick)
cd Tests && wolframscript -file test_package_loading.m

# Symbol namespacing test (quick)
cd Tests && wolframscript -file test_symbols_namespacing.m

# Simple configuration test (quick)
cd Tests && wolframscript -file test_unequal_mass_banana.m

# Full transport test
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
- test_package_loading.m: **18/18 PASS**
- test_symbols_namespacing.m: **ALL PASS**
- test_unequal_mass_banana.m: **5/5 PASS**
- test_unequal_mass_full.m: **PASS** (completes transport to mm1=2, mm2=3/2, mm3=4/3, mm4=1)

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
