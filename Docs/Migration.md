# Migrating from DiffExp 1

DiffExp 2 keeps the line-transport workflow but replaces the strategy stack
and collapsed epsilon-series representation.  The central change is that local
``x^(a+b eps) Log[x]^p`` behavior remains exact solver data.

## Loading

Old:

```mathematica
Get["/path/to/DiffExp/DiffExp.m"];
```

Current DiffExp 2 snapshot:

```mathematica
Get["/path/to/DiffExp2/DiffExp2/DiffExp2.m"];
```

The repository-root ``DiffExp.m`` still loads the legacy modular package in
the prototype snapshot.  A clean release must replace it with a DiffExp 2
loader or remove the ambiguity; the documentation does not pretend this has
already happened.

## Workflow mapping

| DiffExp 1 concept | DiffExp 2 replacement |
| --- | --- |
| ``LoadConfiguration`` | ``DiffExp2`Config`LoadConfiguration`` |
| matrix directory with epsilon slices | ``DiffExp2`API`LoadSystem`` with an exact matrix or ``d<var>_full.m`` |
| ``PrepareBoundaryConditions`` | finite coefficient array at a regular anchor; a friendly closed-form wrapper is still missing |
| ``TransportTo`` | ``DiffExp2`API`TransportEndpoint`` |
| explicit line segmentation | ``SegmentLine`` followed by ``TransportLine`` |
| saved segment expansions | ``transport["Charts"]`` with one exact ``LocalSolution`` per chart |
| ``DefiniteIntegral`` / prefactor variants | ``DiffExp2`API`LineIntegral`` |
| singular endpoint limit | ``DiffExp2`API`EndpointLimitValues`` |
| ``DecomposeSingularity`` | ``DiffExp2`SectorSeries`SectorDecomposition`` |
| ordinary epsilon ``SeriesData`` | ``DiffExp2`EpsSeries`` with an honest Laurent window |
| Wronskian/VOP/recurrence strategy selection | one strict recurrence solver, with C++ or Wolfram execution backend |
| ``ToPiecewise`` / plot helpers | inspect and evaluate ``result["Charts"]``; convenience wrappers are not yet ported |

## Configuration changes

The following concepts remain:

- working/chop precision;
- expansion and epsilon orders;
- division order and radius-of-convergence rescaling;
- delta prescriptions;
- Pade evaluation;
- error estimation and verbosity.

The following old strategy controls are deliberately rejected:

- ``IntegrationStrategy``;
- ``UseRationalRecurrence``;
- ``InvWronskSolver``;
- ``HomogeneousSolve``;
- ``IgnoreIndicialCheck``;
- ``UseMobius``.

There is one recurrence mathematics path.  ``RecurrenceBackend`` chooses only
where its finite-width coefficient work executes.  Requesting C++ is strict
and never silently falls back to Wolfram.

``SegmentationStrategy`` currently accepts only ``"Predivision"``.  Dynamic
segmentation and Mobius chart maps are not part of the current solver.

## Boundary values

DiffExp 1 accepted several convenient closed-form and wildcard boundary
formats.  The current direct DiffExp 2 API expects a finite rectangular array
at a regular anchor:

```text
boundary[[master, k+1]] = coefficient of eps^k.
```

Example:

```mathematica
boundary = Transpose@Table[
  {
    SeriesCoefficient[f1[eps], {eps, 0, k}],
    SeriesCoefficient[f2[eps], {eps, 0, k}]
  },
  {k, 0, epsilonOrder}
];
```

A release-friendly wrapper that accepts closed expressions and exposes their
normalization is still needed.  Do not emulate missing data by inserting zero
coefficients: DiffExp 2's ``EpsSeries`` window exists specifically to prevent
that ambiguity.

## Matrix format

DiffExp 1 commonly loaded files such as ``dx_0.m``, ``dx_1.m``, and so on.
DiffExp 2 requires the exact epsilon-rational matrix:

```mathematica
DiffExp2`API`LoadSystem[<|
  "FullMatrixFile" -> "/path/to/dx_full.m",
  "Variable" -> Global`x
|>]
```

This is required for exact indicial exponents, resonance decisions, Laurent
window budgeting, and analytic-regulator compatibility.

## Singularities and analytic continuation

Keep explicit prescriptions:

```mathematica
"DeltaPrescriptions" -> {{factor1, 1}, {factor2, -1}}
```

DiffExp 2 derives the chart side from every odd-multiplicity vanishing factor
and fails on conflicts.  Interior singular crossings transform exact
``{a,b,p}`` sectors, including log-chain mixing.

The current input limitation is narrower and stricter than a generic
“algebraic not supported” statement: a non-integer power whose base depends on
the line variable, such as ``Sqrt[x]`` in the matrix, is rejected.  Rationalize
the system with an external basis/variable transformation before loading it.

## Result access

Old code often returned a list of epsilon coefficients and optionally saved
global segment state.  DiffExp 2 returns an association:

```mathematica
result["Value"]
result["Final"]
result["Charts"]
result["ErrorEstimate"]
```

Read values with ``ESCoefficient`` and retain the window metadata until final
presentation.  To inspect exact singular behavior:

```mathematica
DiffExp2`SectorSeries`SectorDecomposition[result["Final"]]
```

## Integration

``LineIntegral`` takes a coefficient vector with one rational function per
master and combines the scalar observable before endpoint regularization.  If
an IBP coefficient introduces additional poles, pass their factors through
``"ExtraSingularFactors"`` so segmentation sees them before multiplication.

This order replaces several DiffExp 1 regularized-integration and fitted-sector
helpers.  It is also why separate divergent master integrals should not be
computed and added manually when only their combination is finite.

## Feynman-trick migration

Use ``Scripts/run_ft_stepwise2.m``, not the legacy
``Scripts/run_ft_stepwise.m``.  The DiffExp 2 runner consumes the exact
in-memory differential matrix, transports exact local sectors, supports
separate FIRE preparation caches and ladder checkpoints, and refuses
incomplete epsilon windows.

The current CLI remains more mature than the public Wolfram facade.  Existing
research notebooks should invoke a named example or extract a stable public
driver before release, rather than calling private runner helpers.

## Features not yet restored as friendly wrappers

- closed-form boundary preparation and ``"?"`` wildcard handling;
- high-level full-line and arbitrary line-chain convenience functions;
- ``ToPiecewise`` and plot helpers;
- a top-level namespace with short names;
- a one-call Feynman-trick API;
- automatic rationalizing transformations for square roots in the input
  basis.

The lower-level data required for line segments, plots, and exact local powers
is already present.  The examples show how to use it without claiming the
missing wrappers exist.
