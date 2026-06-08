# Repository Instructions

## Scope

- Keep changes focused on the requested behavior. Avoid unrelated refactors and metadata churn.
- Do not revert or overwrite unrelated local changes. If the worktree is dirty, stage only the files needed for the task.
- Use `wolframscript` or `WolframKernel` command-line runs for verification; do not rely on notebooks for tests.

## Package Notes

- `DiffExp.m` is the package loader. Implementation modules live under `DiffExp/`.
- The internal line parameter is ``DiffExp`Symbols`x``; use that symbol consistently in package code.
- Configuration state is centralized in `DiffExp/State.m`. Prefer the existing configuration helpers over ad hoc global state.
- Some option keys intentionally use built-in Wolfram contexts, for example `AccuracyGoal`, `Variables`, and `WorkingPrecision`.
- A ``Global`x::shdw`` warning can appear during test setup because the package exports `x`; this warning is expected.

## Recurrence Solver Notes

- `UseRationalRecurrence -> True` is the strict recursive finite-width solver path. Do not hide recurrence failures by silently falling back to the older Wronskian/Frobenius/default solver.
- The recurrence layer should work with finite-width systems of the form `q(z, eps) theta[g] == C(z, eps) g`, using coefficient recurrences rather than long truncated products.
- Preserve support for ordinary rational recurrence, regular-singular recurrence, resonant/log sectors, fractional powers, inhomogeneous multi-sector sources, and singularity/infinity transport.
- When handling singular points, keep the recurrence path active there as well. If the recurrence cannot support a case, return a diagnostic failure rather than taking an old integration path.
- Avoid introducing `SeriesData`-dependent logic in the recurrence core; prefer explicit coefficient extraction and finite block solves.

## Focused Tests

For recurrence or singularity work, run the relevant subset from the repository root:

```sh
wolframscript -file Tests/test_package_loading.m
wolframscript -file Tests/test_recurrence_no_fallback.m
wolframscript -file Tests/test_local_series_edge_cases.m
wolframscript -file Tests/test_resonant_2f1.m
wolframscript -file Tests/test_singular_recurrence.m
wolframscript -file Tests/test_unequal_mass_banana_parity_speed.m
```

The docs in `Docs/` describe the main package areas: core modules, series/evaluation, transformations, infrastructure, transport strategies, singularity decomposition, and regularized integration.
