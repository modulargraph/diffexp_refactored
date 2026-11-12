# IBP Solver core

Core headers from https://github.com/modulargraph/ibp-solver, commit `8d803554c33ac691877ef72c3dd338929c862abc`, extended for DiffExp's recursion. Licensed under GPL-3.0; see LICENSE.

This source copy lets DiffExp build offline. DiffExp uses the `InputGeometry`, `ParametricProgram`, `ArithmeticTrace`, `Field`, `Solver` and `BatchSolver` interfaces through its adapter (`include/diffexp/ibp_solver_provider.hpp`). The adapter compiles the parameterized equations once, eliminates the selected equations at the sample points of each reconstruction prime, and independently checks full reductions at validation points. GPU batch replay remains in the standalone package.

Changes against the upstream commit:

- `geometry.hpp`: an integral holds up to 32 propagator powers (`IBP_SOLVER_INTEGRAL_SLOTS`); sector masks are `size_t`; no caps on loops, denominators or seed powers beyond what an `int8_t` power holds.
- `parametric.hpp`: an optional per-sector seed bound for `ParametricProgram::compile` (an empty bound keeps the upstream behaviour); the origin (seed and operator) of every equation in `ParametricProgram::origins`; coefficients interned and equations stored as (column, form) pairs that solvers read directly.
- `solver.hpp`: eliminations continued from the pivots of a subsystem (`solve` with `kept`, `release`, `restart`); reductions of columns other than the program's targets; a row under reduction held in a dense vector with an occupancy bitmap; whole systems eliminated with complete rows; a pivot-fill budget derived from available memory instead of a fixed cap.
- `batch.hpp` (new): `BatchSolver` eliminates one selection at several sample points at once (eight by default), with one pivot search and row layout per batch and unreduced accumulation of products. With AVX-512 IFMA available at compile time it uses a vector kernel; otherwise a portable one. A lane whose pivot vanishes falls back to the scalar solver.
- `trace.hpp`: traces read the interned equations in place. The adapter learns a trace only where replaying it is cheaper than eliminating.

Update the headers together from an audited upstream commit and re-apply these changes.
