# Analytic Continuation and Exact Local Sectors

DiffExp 2 treats analytic continuation as part of transport geometry.  A chart
knows the singular factors that vanish at its center, the side prescribed for
each odd-multiplicity factor, and the exact local powers carried by its
solution.  Crossing is therefore a transformation of tagged sectors, not a
late numerical replacement in a collapsed epsilon series.

## Singularities and chart radii

```mathematica
singularities = DiffExp2`Transport`FindSingularities[sys];
```

The result includes exact roots of the matrix and extra singular factors.
True complex roots constrain chart radii through their complex-plane distance.
Their classic real projections, `Re[z]` and `Re[z] +/- Im[z]`, can also appear
as regular planning waypoints.  A projected waypoint is not marked singular
unless it is itself an exact real root.

Line geometry uses an exact epsilon-zero alphabet, not the moving roots of a
finite regulator.  Before finding roots, the planner removes each factor's
overall epsilon valuation and takes its first nonzero epsilon coefficient.
For example, `eps (x + eps)` projects to `x`; `1 + eps x` projects to a
constant and therefore has no finite root.  This projection never modifies
the differential matrix: `LoadSystem` retains the full factors as
`SingularFactorsExact` and exposes the planner factors as `SingularFactors`.
The same projection is applied to `DeltaPrescriptions` before attaching branch
metadata to a chart, so an exact prescription for `x + eps` is not lost when
the limiting chart lies at `x = 0`.

This is a geometry and branch-metadata projection, not an assertion that every
moving matrix pole has a local Frobenius expansion in the current solver.  If
such a pole degenerates exactly onto a chart center, denominator clearing
returns the documented loud `E3` failure rather than expanding in `t/eps` or
silently treating the point as ordinary.

This distinction matters for stable matching: projection points can force a
better chart layout without fabricating Frobenius behavior on the real line.

## Delta prescriptions

Configure a physical side before planning:

```mathematica
x = Global`x;

DiffExp2`LoadConfiguration[
  "DeltaPrescriptions" -> {
    {x, 1},
    {x - 1/2, -1},
    {1 - x, 1}
  }
];
```

Each pair is `{polynomial, sign}`, with `sign` equal to `1` or `-1`.  The
equivalent ``polynomial +/- I Global`\[Delta]`` input form is also parsed.  The
configuration layer:

- pins non-system symbols into the Wolfram global context;
- canonicalizes an overall sign of the polynomial together with its
  prescription sign;
- requires an irreducible nonzero polynomial;
- rejects conflicting prescriptions for the same canonical factor.

At a chart, the transport layer combines every vanishing prescribed factor,
its multiplicity, and its leading-coefficient orientation into one imaginary
side.  Conflicting odd-multiplicity requirements are a loud error.

The branch convention is regression-pinned.  For example, transporting
`f' = eps f/(x-1/2)` from `1/4` to `3/4` with
`{x-1/2,1}` produces the rightward factor `Exp[-I Pi eps]`; reversing the line
produces its inverse.  A material multivalued sector cannot cross without a
derivable prescription.

## Exact local representation

A `LocalSolution` contains a list of sectors.  Each sector has exact tags

```text
<|"a" -> a, "b" -> b, "p" -> p, "Coeffs" -> ...|>
```

and represents

```text
t^(a + b eps) (eps Log[t])^p/p!
```

times a regular vector Taylor series.  `t` is the affine chart coordinate.
Resonant log chains occupy distinct values of `p`; sectors with different `b`
remain distinct even when they collide at `eps = 0`.

Read this structure directly:

```mathematica
behavior = DiffExp2`LocalBehavior[localSolution];
sectors = DiffExp2`ExactSectors[localSolution];
tags = KeyTake[#, {"a", "b", "p", "Exponent"}] & /@ sectors;
```

This is the exact replacement for the original `DecomposeSingularity`
fitting workflow.  There is no Prony fit or averaged exponent between the
solver and endpoint integration.

## Evaluating either side of a chart

```mathematica
evaluation = DiffExp2`EvaluateLocal[
  localSolution,
  tValue,
  "ImSign" -> Automatic,
  "UsePade" -> False
];
```

For positive `t`, ordinary real powers and logs are used.  For negative `t`,
the exact power phase and the triangular mixing of log-chain sectors are
determined by the chart's imaginary sign.  Evaluation fails if the point is
outside the local radius or if a multivalued negative point lacks a branch
side.

The convention agrees literally with Mathematica on the upper rim.  Writing
`t=-r`, `r>0`,

```text
+i0:  Log[t] = Log[r] + i pi,   t^alpha = r^alpha exp(+i pi alpha)
-i0:  Log[t] = Log[r] - i pi,   t^alpha = r^alpha exp(-i pi alpha).
```

Thus the lower-rim value is Mathematica's principal negative-axis value with
`Log[t] -> Log[t] - 2 pi i` and
`t^alpha -> exp(-2 pi i alpha) t^alpha`.  The same logarithm is expanded
exactly when `alpha` contains `b eps`; it is not inferred from a numerical
fit.

Do not evaluate a multivalued sector at `t=0`.  Use
`EndpointLimit`, which implements the dimensional-regularization endpoint
rule and checks for uncancelled divergences.

## Crossing an interior singularity

`TransportLine` handles a supported interior singularity in three stages:

1. approach on the incoming side at a point certified inside both adjacent
   convergence disks;
2. match to the exact singular local basis;
3. apply the crossing operator exactly once and continue on the far side.

Pole-normalized Frobenius columns can make a finite epsilon window look
artificially incomplete if the physical boundary is applied to them
directly. For an interior crossing, DiffExp 2 instead normalizes the singular
basis to the identity at the incoming match point, evaluates that transfer
frame on the outgoing side, removes only certified polar cancellation
remnants, and then applies the user's boundary window. The private basis
guard is checked against the requested complete order before marching
continues; missing higher physical boundary coefficients are not fabricated
or requested.

The crossing operator contains both the phase of `a+b eps` and the unipotent
mixing generated by shifting `Log[t]`.  Treating a `p>0` sector as a scalar
phase would be incorrect.

The same sign requirement applies when a singular endpoint is approached on
its negative local arm and when a branch-sensitive local solution is
integrated wholly on the negative arm.  DiffExp 2 does not silently replace
a missing sign by `+i0`; these cases fail loudly.  A sign is unnecessary only
when all material sectors are single-valued integer powers (`a` integer,
`b=0`, `p=0`).

## Analytically regularized endpoint integration

For `b != 0`, the lower endpoint of a sector is defined by analytic
continuation in epsilon.  Integrals of pole sectors can therefore deepen the
Laurent window.  DiffExp 2 represents that change explicitly; division by a
denominator beginning at epsilon order one shifts both edges of an
`EpsSeries` window.

For `b == 0`, a genuine power or logarithmic divergence is rejected unless it
cancels after the requested master combination is formed. `IntegrateLine` and
`EndpointLimit` combine components before applying this gate.

One convention is intentionally separate from contour continuation:
integer-`a`, `b=0` poles and log terms at an *interior* integration point use
the real principal-value/Hadamard finite part.  This is the Euclidean
Feynman-trick convention validated by the box examples and is independent of
the delta sign.  A future physical contour mode would instead add the
corresponding `-i pi sigma` residue term; it must be an explicit option, not a
silent change of the current default.  Noninteger-`a` interior sectors are
multivalued and already use the prescribed contour phase.

## Verified limitation: roots in the input basis

The current `LoadSystem` implementation rejects an entry containing
`Power[base, exponent]` when:

- `base` depends on the line variable; and
- `exponent` is not an integer.

Thus `Sqrt[x]`, `Sqrt[1-x]`, and similar algebraic dependence in the input
differential matrix are out of scope.  This is the precise source-verified
meaning behind the shorthand “roots in the basis are not yet supported.”

The restriction is not broader than that:

- algebraic constants in coefficients are handled by Wolfram recurrence paths
  and some numeric C++ paths;
- exact algebraic singular locations are accepted by the planner;
- algebraic values of the indicial constant `a` are represented by sector
  tags;
- the C++ symbolic-regulator coefficient field is narrower: it accepts exact
  rational functions over `Q`, not algebraic symbolic coefficients.

Do not rationalize away an actual square root to bypass the check.  Use a
rationalizing basis/variable transformation outside DiffExp 2, then load the
exact transformed system.

## Current public-interface gaps

- `PrepareBoundary` accepts finite epsilon expansions at a regular anchor,
  while `PrepareLaurentBoundary` constructs typed regular-anchor point data
  with an honest Laurent window. Neither reproduces the old
  `PrepareBoundaryConditions` mode that
  accepts an asymptotic/singular-start expression depending on the line
  coordinate and converts its `Log[x]` and fractional powers into a tagged
  `LocalSolution`. Starting a transport at a singular point is currently
  unsupported. Ordinary Euclidean-point boundaries and all subsequent
  transport crossings are branch-aware as described above.
- The public system consumes an already pulled-back one-variable matrix and
  prescriptions in that variable.  It does not yet accept a multivariate
  physical line association and pull physical prescription factors back
  automatically as DiffExp 1 did.
- `AbortOnAnalyticContinuationFail` remains a compatibility schema key, but
  DiffExp 2 currently has no flagged-incomplete soft mode: branch ambiguity
  is always loud.  Passing `False` must not be interpreted as permission to
  choose `+i0` silently; implementing or removing that compatibility key is
  still pending.
- Automatic prescriptions derived from square roots in the input matrix do not
  exist because that matrix class is rejected.
- There is no dedicated branch-plot styling function; use `LineSegments` or
  `PiecewiseSolution` as shown in the direct example.
- Weighted `EndpointLimit` accepts epsilon-free endpoint scalars. General
  rational endpoint coefficients are supported by the lower-level
  `MultiplyRational` plus combination path, but do not yet have a dedicated
  public wrapper.
