# DiffExp Refactoring Project

## Project Rules

- ONLY modify files within this directory (~/Desktop/diffexp_refactor)
- NEVER access files outside this folder
- Ask before running any destructive commands

## Project Status: REFACTORED (v1.2)

The DiffExp Mathematica package has been successfully refactored from a monolithic file into modular subpackages. All tests pass including the unequal mass banana example.

## Git Repository

This project now uses git for version control. Always commit your changes with descriptive messages after completing tasks.

## Directory Structure

```
diffexp_refactor/
├── DiffExp.m                    # Root loader (loads DiffExp/DiffExp.m)
├── DiffExp/                     # Refactored modular package
│   ├── DiffExp.m                # Main entry point (loads all subpackages)
│   ├── Symbols.m                # Core symbols (x, eps, Logx, θp, θm)
│   ├── State.m                  # Configuration and state management
│   ├── Utilities.m              # Helper functions (printing, chopping, etc.)
│   ├── SeriesOps.m              # Series operations
│   ├── Integration.m            # Series integration functions
│   ├── Pade.m                   # Pade approximant functions
│   ├── Mobius.m                 # Mobius transformation functions
│   ├── AnalyticContinuation.m   # Analytic continuation handling
│   ├── LineSegmentation.m       # Line segmentation strategies
│   ├── Frobenius.m              # Frobenius method functions
│   ├── Wronskian.m              # Wronskian computations
│   ├── MatrixLoading.m          # Matrix loading and processing
│   └── Transport.m              # TransportTo and IntegrateSystem
├── Reference/                   # Original code and examples for reference
│   ├── DiffExp_original.m       # Original monolithic DiffExp.m
│   └── Examples/                # Example scripts (moved from root)
│       ├── Banana_example.m
│       ├── MultiplePolylogarithms_example.m
│       └── FivePointNonPlanar_example.m
├── Tests/                       # Test scripts and test data
│   ├── Banana_EqualMass_Matrices/   # Equal mass matrices (moved from Examples)
│   ├── Banana_Matrices/             # Unequal mass matrices (moved from Examples)
│   ├── test_package_loading.m       # Package loading test (18/18 pass)
│   ├── test_symbols_namespacing.m   # Symbol namespacing test (all pass)
│   ├── test_unequal_mass_banana.m   # Simple configuration test (5/5 pass)
│   ├── test_unequal_mass_full.m     # Full transport test (passes)
│   ├── test_banana_refactored.m     # Equal mass banana test
│   ├── test_multiple_polylogarithms.m  # Multiple polylogarithms test
│   └── test_five_point_nonplanar.m     # Five-point nonplanar test
└── CLAUDE.md                    # This file
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

### Important Rules (preserved from original prompt)
1. NEVER rename variables or functions from the original code
2. NEVER change functionality - this is a 1-1 refactor
3. Keep all Attributes on functions exactly as they were
4. Each function should be defined in exactly one place
5. Use `DiffExp\`Symbols\`x` for the line parameter variable throughout

## Running Tests

```bash
export WOLFRAMSCRIPT_KERNELPATH="/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/WolframKernel"

# Package loading test (quick)
cd Tests && wolframscript -file test_package_loading.m

# Symbol namespacing test (quick)
cd Tests && wolframscript -file test_symbols_namespacing.m

# Simple configuration test (quick)
cd Tests && wolframscript -file test_unequal_mass_banana.m

# Full transport test (~2 minutes)
cd Tests && wolframscript -file test_unequal_mass_full.m

# Multiple polylogarithms test (requires computation)
cd Tests && wolframscript -file test_multiple_polylogarithms.m

# Five-point nonplanar test (requires internet for arXiv files)
cd Tests && wolframscript -file test_five_point_nonplanar.m
```

## Current Test Results
- test_package_loading.m: **18/18 PASS**
- test_symbols_namespacing.m: **ALL PASS**
- test_unequal_mass_banana.m: **5/5 PASS**
- test_unequal_mass_full.m: **PASS** (completes transport to mm1=2, mm2=3/2, mm3=4/3, mm4=1)

## Task Instructions

When given a task:
1. Read this prompt to understand the project state
2. Use git to track changes (commit after completing work)
3. Run relevant tests to verify changes don't break anything

## Constraints
- ONLY modify files within ~/Desktop/diffexp_refactor
- NEVER access files outside this folder
- Ask before running any destructive commands
- Use wolframscript (not Mathematica notebooks) for testing
