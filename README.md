# DiffExp Refactored

DiffExp is a Mathematica package for integrating Feynman integrals in terms of series expansions, using differential equations.

This repository is an experimental, AI-assisted refactoring and upgrade project based on the original DiffExp code. The original reference repository remains at [gitlab.com/hiddingm/diffexp](https://gitlab.com/hiddingm/diffexp).

## Version

This refactored repository is version 2.0.

## What Is New

- Modular package layout under `DiffExp/` instead of a single generated package file.
- A stricter recurrence-based transport path enabled by `UseRationalRecurrence -> True`.
- Recursive finite-width coefficient solvers intended to avoid long truncated products.
- Support for ordinary points, regular singular points, resonant/logarithmic sectors, fractional powers, singular transport, infinity transport, and inhomogeneous multi-sector source terms.
- Focused regression tests comparing recurrence transport against the standard solver on equal-mass and unequal-mass banana examples.

The recurrence methods are still an active experimental upgrade. The older Wronskian/Frobenius strategies remain available, but recurrence failures should be reported directly rather than silently hidden by fallback behavior when `UseRationalRecurrence -> True`.

## Loading

From Mathematica or Wolfram Engine:

```mathematica
Get["/path/to/diffexp_refactored/DiffExp.m"]
```

The package prints its version on load.

## Tests

Focused recurrence and singularity tests can be run from the repository root:

```sh
wolframscript -file Tests/test_package_loading.m
wolframscript -file Tests/test_recurrence_no_fallback.m
wolframscript -file Tests/test_local_series_edge_cases.m
wolframscript -file Tests/test_resonant_2f1.m
wolframscript -file Tests/test_singular_recurrence.m
wolframscript -file Tests/test_unequal_mass_banana_parity_speed.m
```

Additional implementation notes are in `Docs/`.
