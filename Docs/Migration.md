# Migrating from DiffExp 1

DiffExp 2 keeps the line-transport workflow but replaces the strategy stack
and collapsed epsilon-series representation.  The central change is that local
``x^(a+b eps) Log[x]^p`` behavior remains exact solver data.

## Loading

Old:

```mathematica
Get["/path/to/DiffExp/DiffExp.m"];
```

DiffExp 2:

```mathematica
Get["/path/to/diffexp2/DiffExp2.m"];
```

During development, repository-root ``DiffExp.m`` still loads the legacy
modular package. Load ``DiffExp2.m`` explicitly; the clean release removes
that ambiguity.

## Workflow mapping

| DiffExp 1 concept | DiffExp 2 replacement |
| --- | --- |
| ``LoadConfiguration`` | ``DiffExp2`LoadConfiguration`` |
| matrix directory with epsilon slices | ``DiffExp2`LoadSystem`` with an exact matrix or ``d<var>_full.m`` |
| ``PrepareBoundaryConditions`` | ``DiffExp2`PrepareBoundary`` for finite regular-anchor expressions |
| ``TransportTo`` | ``DiffExp2`TransportEndpoint`` |
| explicit line segmentation | ``PlanLine`` followed by ``TransportLine`` |
| saved segment expansions | ``LineSegments[result]`` with one exact ``LocalSolution`` per chart |
| ``DefiniteIntegral`` / prefactor variants | ``DiffExp2`IntegrateLine`` |
| singular endpoint limit | ``DiffExp2`EndpointLimit`` |
| ``DecomposeSingularity`` | ``DiffExp2`LocalBehavior`` / ``ExactSectors`` |
| ordinary epsilon ``SeriesData`` | honest ``EpsilonWindow`` and coefficient accessors |
| Wronskian/VOP/recurrence strategy selection | one strict recurrence solver, with C++ or Wolfram execution backend |
| ``ToPiecewise`` / plot helpers | ``PiecewiseSolution`` and ``LineSegments`` |

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
formats. The direct DiffExp 2 API accepts a finite rectangular array at a
regular anchor:

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

For finite expressions, ``PrepareBoundary[{f1[eps],f2[eps]}]`` constructs the
same array at the configured epsilon order. Pole-normalized boundaries still
need an explicit ``LocalSolution`` so their negative lower window is not
discarded. Do not emulate missing data by inserting zero coefficients.

## Matrix format

DiffExp 1 commonly loaded files such as ``dx_0.m``, ``dx_1.m``, and so on.
DiffExp 2 requires the exact epsilon-rational matrix:

```mathematica
DiffExp2`LoadSystem[
  "/path/to/dx_full.m", "Variable" -> Global`x
]
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

Read values with ``EpsilonCoefficient`` and retain the window metadata until final
presentation.  To inspect exact singular behavior:

```mathematica
DiffExp2`ExactSectors[result]
```

## Integration

``IntegrateLine`` takes a coefficient vector with one rational function per
master and combines the scalar observable before endpoint regularization.  If
an IBP coefficient introduces additional poles, pass their factors through
``"ExtraSingularFactors"`` so segmentation sees them before multiplication.

This order replaces several DiffExp 1 regularized-integration and fitted-sector
helpers.  It is also why separate divergent master integrals should not be
computed and added manually when only their combination is finite.

## Feynman-trick migration

Use ``FeynmanTrick`RunIntegrationPipeline`` (or its underlying
``Scripts/run_ft_stepwise2.m`` runner), not the legacy
``Scripts/run_ft_stepwise.m``. The DiffExp 2 runner consumes the exact
in-memory differential matrix, transports exact local sectors, supports
separate FIRE preparation caches and ladder checkpoints, and refuses
incomplete epsilon windows.

The public facade accepts either one named registry example or an exact raw
family/unprepared topology with selected output integrals or `All`. It
preserves the runner's prepared FIRE cache and atomic ladder checkpoints, and
binds custom-family requests and dynamic `All` resolutions into those cache
identities.

## Features not yet restored as friendly wrappers

- ``"?"`` wildcard boundary handling and pole-normalized closed-form
  boundaries (finite regular-anchor expressions are supported);
- arbitrary multi-leg line-chain convenience functions;
- a dedicated plotting-style helper (piecewise evaluation is supported);
- automatic rationalizing transformations for square roots in the input
  basis.

The custom-family facade rejects FIRE-added numerator slots during `All`
discovery, numerators at merge positions, and custom analytic-prescription or
kinematic-assumption fields until their recursion and branch semantics are
exposed safely. Explicit-target preparation may retain FIRE auxiliary slots,
but they are not part of the public output-selection contract.

The lower-level data required for line segments, plots, and exact local powers
is already present.  The examples show how to use it without claiming the
missing wrappers exist.
