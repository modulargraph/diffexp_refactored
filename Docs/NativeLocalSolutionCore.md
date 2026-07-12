# Native local-solution core

This note records the first native data/evaluation seam for the persistent
C++ solver.  It is a map of the current Wolfram contracts, not a claim that
the complete persistent solver already exists.

## Source contracts mapped

The authoritative Wolfram representation is `DiffExp2/SectorSeries.m`:

- a sector is
  `t^(a+b eps) (eps Log[t])^p/p! Sum_n c[k,n,component] t^n`;
- `a` and `b` are exact tags and `p` is a nonnegative integer;
- all sectors share one honest epsilon window and Taylor width;
- prescriptions are exact factor records.  Odd multiplicities require
  `sigma = Sign * LeadingCoeffSign`, even multiplicities impose no
  constraint, and conflicting odd requirements are fatal;
- on the negative arm,
  `Log[t] = Log[|t|] + i pi sigma` and
  `t^a = |t|^a exp(i pi sigma a)`.  Expanding `exp(b eps Log[t])`
  supplies the complete epsilon-dependent phase;
- coefficients below `Min` are structural zero.  Coefficients above
  `CompleteMax` are unknown and must never be zero-padded into a proof.

`DiffExp2/Transport.m` adds the matching obligations:

1. evaluate the incoming solution and every local basis column at the same
   certified overlap point and on the same sheet;
2. epsilon-saturate a degenerate basis before solving (column valuation
   shifts, certified null combinations divided by epsilon, then optional
   epsilon-independent right-frame normalization);
3. solve the square Laurent system with full row/column pivoting and honest
   quotient windows;
4. perform at most two residual corrections;
5. prove the final residual against the original, untrimmed matrix and
   right-hand side.  A correction may not pass merely by shortening the
   complete upper window;
6. apply the crossing operator exactly once when changing the local
   coordinate to the far-side positive arm.  Direct negative-arm evaluation
   and the crossing operator must not both contribute the phase.

`DiffExp2/Solve.m` evaluates both `f` and `theta f`, substitutes the original
theta-form matrix, uses the union of lower supports and the intersection of
complete upper supports, and compares a componentwise scaled residual.
Formal analytic regulators are currently checked at several exact rational
specializations rather than numerically erasing the regulator field.

## Native types

`cpp/include/diffexp2/series_types.hpp` owns `EpsilonWindow`.  It is intended
to be shared by recurrence assembly, matching, endpoint limits and line/tile
integration.

`cpp/include/diffexp2/local_solution.hpp` adds:

- `ExactScalarDescriptor`: immutable canonical exact text, exact domain,
  symbol names, structural zero/integer/sign facts, and an optional Acb
  specialization.  Structural decisions use the facts, never the ball;
- `LocalSector<Scalar>`: tags plus a flat
  `[epsilon][Taylor][component]` coefficient tensor, matching the ordering of
  `AssembledResult` in `recurrence.hpp`;
- `LocalSolution<Scalar>`: chart geometry, frames, sectors, prescriptions,
  error provenance and a stable checkpoint identity;
- `EpsilonVector` and `EpsilonMatrix`: numeric framed values used at the
  evaluation/matching/residual boundary;
- `Magnitude` and `ErrorEnvelope`: nonnegative FLINT magnitude bounds with an
  explicit `None`, `Advisory` or `Certified` guarantee.

The exact descriptor intentionally outlives a numeric specialization.  For
example a checkpoint can retain `b=(1+rho)/(2-rho)` exactly, while a numeric
transport installs a certified Acb value at one regulator sample.  Unknown
integer/zero facts are treated conservatively as branch-sensitive.

## Implemented primitives

`evaluate_local_solution` performs direct Taylor evaluation in Acb and
returns both `f` and `theta f`.  It supports:

- arbitrary honest Laurent frame bases, including negative orders;
- fractional `a`, epsilon-dependent `b`, and explicit log sectors;
- exact `+i0` and `-i0` conventions without asking `acb_pow` to choose a
  negative-axis branch (Acb sees the positive modulus; the phase is explicit);
- exact topology of the point (`-1`, `0`, `+1`) separate from its enclosure;
- a reduced Taylor order for the existing full-versus-reduced probe;
- a geometric top-column tail estimate marked **advisory**.

The theta transform is performed on the stored sectors before evaluation:

```text
theta [c[k,n] t^(a+b eps+n) L_p]
  = ((a+n)c[k,n] + b c[k-1,n]) L_p
    + c[k-1,n] L_(p-1),
```

where `L_p=(eps Log[t])^p/p!`.  This retains the exact epsilon/log shifts and
does not differentiate a collapsed numerical value.

`certify_theta_residual` accepts an already evaluated Laurent matrix and an
optional source.  Acb upper and lower magnitude bounds yield three outcomes:

- `Pass`: the upper residual bound is below `tol * scaleLower`;
- `Fail`: the lower residual bound is above `tol * scaleUpper`;
- `Inconclusive`: the enclosures overlap the threshold.

This is a genuine enclosure for the **stored truncation**.  The generic
full-local residual request remains inconclusive; in particular it never uses
the geometric top-column estimate as a proof.

`cpp/include/diffexp2/tail_majorant.hpp` provides one deliberately narrow
full-function exception.  For a retained ordinary homogeneous chart it first
proves all of the following from the prepared operator and run:

- `a=b=0`, no log sector, no source, identity assembly, and singleton zero
  indicial blocks;
- every retained epsilon shift in `q` and `N` is exactly zero (Rational and
  Acb coefficient fields are both supported);
- `N(0)=0`, so `q(t) theta f=N(t)f` is an ordinary regular ODE;
- the stored Taylor tensor forward-encloses the recurrence from the parsed
  initial frame through its claimed complete order.

For an exact rational witness radius `R` strictly inside the chart, it proves
`q` is separated from zero using
`lower(|q0|)-sum_(j>0) upper(|qj|) R^j > 0`.  An infinity-norm ODE bound and
Gronwall then bound the solution on `|t|=R`; Cauchy's coefficient estimate
certifies the unseen value, theta-value, and integrated line tails at strict
subradius points.  The resulting `ErrorEnvelope` is marked `Certified` and
retains the operator identity, local checkpoint identity, witness disk, and
analytic-prescription provenance.  A line result is promoted from
`StoredTruncation` to `FullLocalWithCertifiedTail` only after that envelope is
installed.

A disk on which the triangle bound cannot separate `q` is
`Inconclusive`.  Singular/logarithmic or regulator-dependent powers,
epsilon-coupled operators, nonidentity assemblies, sourced recurrences, and
epsilon rows outside the retained model are explicitly `Unsupported`.  These
outcomes are load-bearing: none falls back to the advisory last-column model.

Evaluation at `t=0` is deliberately strict.  Only known nonnegative integer
`a`, identically zero `b`, and `p=0` are admitted.  Negative integer powers
are redirected to endpoint-limit logic; this follows the written
SectorSeries contract and closes a permissive edge in the current Wolfram
implementation.

## Shared integration API assumptions

The native endpoint/integration layer can consume `LocalSolution<Scalar>`
directly.  Its scalar Laurent arithmetic may use a separate
`EpsilonFrame<Scalar>`, but it must use the same `EpsilonWindow` and the same
exact tag descriptors.  In particular:

- an endpoint decision uses `ExactScalarDescriptor::is_zero/is_integer/sign`,
  never `specialization.contains_zero()`;
- divergent cancellation groups by exact canonical absolute power and log
  depth;
- `b != 0` analytic-regularization drops and Laurent deepening are exact
  tag operations;
- negative-arm line/tile integration derives the sign from the shared
  prescriptions and uses the same explicit phase as local evaluation;
- any returned error bound carries an explicit guarantee level and source
  epsilon frame.

## Gaps before a persistent native session

The following are explicit missing pieces, not silent fallbacks:

1. **Session and serialization.**  There is no persistent handle table or
   JSON/checkpoint codec yet for charts, exact descriptors, LocalSolutions,
   SCC plans or native match state.
2. **Prepared operator evaluation.**  The residual primitive currently takes
   `B(t,eps)` after evaluation.  A persistent chart must retain its cleared
   operator and evaluate that frame natively at each probe/match point.
3. **Unresolved symbolic evaluation.**  `SymbolicRational` has exact field
   arithmetic but no substitution/evaluation API.  Numeric local evaluation
   therefore requires an explicit Acb specialization.  A future session
   should materialize all coefficients once per regulator sample, not call a
   parser for every coefficient.  Also, the present global FLINT multivariate
   context can represent only one live symbol field safely per session.
4. **Exact algebraic field.**  Algebraic tags retain canonical Wolfram text
   and a certified ball, but C++ has no exact number-field/`Root` arithmetic.
   Wolfram must currently supply structural facts.  The symbolic C++ domain
   remains rational functions over `Q`, not algebraic extensions.
5. **Pade.**  Direct Taylor evaluation is implemented.  Native Pade and its
   loud per-series fallback record are not.  No request is silently changed
   from Pade to direct evaluation.
6. **Crossing object transform.**  The scalar branch rule is implemented;
   the exact sector-level crossing transform (phase plus triangular log-chain
   mixing and `(-1)^n`) still lives in Wolfram.
7. **Laurent matching.**  Epsilon saturation, full-pivot Laurent elimination,
   bounded refinement and the original-system match residual proof still live
   in Wolfram.  These should be moved together; porting only the linear solve
   would omit the banana-critical lattice normalization.
8. **General full error proof.**  The regular homogeneous, epsilon-decoupled
   strict-subradius slice now has a rigorous Gronwall/Cauchy tail envelope.
   Singular/log sectors, inhomogeneous source tails, epsilon coupling,
   nonidentity assembly, and promotion of the generic theta-residual API still
   need corresponding majorants.  Outside the implemented slice, the
   geometric Taylor tail remains advisory.
9. **Endpoint and integration.**  Endpoint limits, analytic-regularized
   primitives, interior pairing and line/tile accumulation are a separate
   native layer built on these shared types.
10. **Allocation/performance pass.**  The first slice favors transparent
    semantics.  A persistent engine should reuse evaluation scratch buffers,
    precompute factorial/exponential convolutions, and batch all basis columns
    at one point rather than allocate one result per column.

## Recommended integration order

1. Add persistent session handles and codecs for exact tags, prepared charts,
   SCC order and LocalSolutions.
2. Assemble recurrence output directly into `LocalSolution<ComplexBall>` (or
   exact coefficient objects followed by one bulk specialization).
3. Add prepared-operator evaluation and batch local value/theta evaluation.
4. Port the complete Laurent saturation/matching/refinement unit and retain
   its original-system residual certificate.
5. Add the exact crossing transform and checkpoint the sheet/state transition.
6. Connect native endpoint/line/tile integration through the shared types.
7. Replace the advisory tail model with a rigorous majorant, then enable
   full-solution certificates.

The narrow C++ smokes cover direct value/theta evaluation, an Acb residual
certificate, both fractional-power rims, epsilon-exponent convolution,
explicit-log epsilon shifting, loud missing/conflicting prescriptions, and a
Rational/Acb regular-tail vertical slice whose point and line envelopes are
checked against `exp(t)`.
