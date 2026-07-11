# Algebraic Square-Root Matrix Support Plan

Status: design only. No implementation is authorized by this document.

## 1. Recommendation

The safest first implementation of differential matrices containing square
roots is an exact rational companion lift performed before the system reaches
the existing recurrence solver. It is not a new recurrence fallback and does
not weaken the rational finite-width contract.

For one radical

```text
r^2 = P(x)
A(x, eps) = A0(x, eps) + r A1(x, eps),
```

introduce the auxiliary master vector `Z = r Y`. The original rank-`d`
system becomes the rational rank-`2d` system

```text
d/dx (Y) = ( A0      A1                 ) (Y)
     (Z)   ( P A1    A0 + P'/(2 P) I_d ) (Z).
```

The public physical masters are the first `d` components. A typed,
prescription-aware boundary operation initializes the auxiliary block as

```text
Z(x0) = r_sigma(x0) Y(x0),
```

where `r_sigma` is the root on the sheet fixed by the effective `i delta`
prescription.

This route is recommended because the lifted matrix is rational in the line
variable and epsilon. `Indicial.m`, `Solve.m`, `SectorSeries.m`, the
denominator-cleared recurrence, and the C++ kernel can therefore retain their
current strict behavior. In particular, unsupported cases remain loud
failures; there is no silent fallback to a truncated product or an older
solver.

## 2. Why simply accepting `Sqrt` is insufficient

The current rejection in `DiffExp2/API.m` is only the first visible barrier.
Removing it alone would violate assumptions throughout the implementation:

- `Indicial.m` computes exact valuations and Laurent coefficients of rational
  functions in the chart coordinate.
- `Solve.m` clears a polynomial denominator and prepares finite Taylor lags.
- `Transport.m` derives the singular alphabet from rational denominator
  factors; it would otherwise miss branch points at radicand zeros and poles.
- `SectorSeries.m` supports fractional powers in *solutions*, but its
  multiplier algebra intentionally accepts rational functions.
- the C++ exact symbolic domain represents rational functions over `Q`; it
  does not represent an algebraic function generator depending on the chart
  variable.

A native square-root recurrence would consequently be a cross-module change,
not a scalar encoder change. The companion lift confines the new algebra to a
preprocessing and API layer and reuses the proven rational core.

## 3. Exact quadratic lift

Let `r^2 = P` and suppose every matrix entry has been reduced exactly in the
quadratic function field to `A0 + r A1`. Since

```text
r' / r = P' / (2 P),
```

the equations for `Y` and `Z = r Y` are

```text
Y' = A0 Y + A1 Z,
Z' = P A1 Y + (A0 + P'/(2 P) I_d) Z.
```

Define the lifted matrix `B` by the block matrix above and, on a chosen sheet,

```text
C_sigma(x) = Stack[I_d, r_sigma(x) I_d].
```

The construction certificate is

```text
B C_sigma = C_sigma' + C_sigma A_sigma
```

exactly modulo `r^2 - P`. This identity must be checked during lift
construction. It proves that every solution of the chosen original sheet
embeds into the rational system.

The sheet constraint

```text
Z - r_sigma Y = 0
```

is then invariant under exact transport. It is also an essential independent
runtime certificate; see section 6.

### 3.1 Exact field reduction

The loader must replace each distinct radical by an inert generator before
performing algebra. `PowerExpand` is forbidden.

For a single generator, an entry rational in `r` can be reduced by polynomial
inversion modulo `r^2 - P`. In the common linear-denominator case,

```text
(N0 + r N1)/(D0 + r D1)
  = ((N0 D0 - P N1 D1) + r (N1 D0 - N0 D1))
      /(D0^2 - P D1^2).
```

The denominator norm must be retained exactly. A zero norm or a noninvertible
denominator is a load error, not a numerical special case.

Radicands should be normalized by their square classes over the base rational
function field. Repeated and dependent forms must not create duplicate
generators. Relations such as `Sqrt[x]` versus `Sqrt[4 x]`, or `Sqrt[P]`
versus `Sqrt[-P]`, must be established with exact field algebra plus an
anchor-branch check, never syntactic simplification alone.

## 4. Multiquadratic lift

For independent generators `r_i^2 = P_i`, use one block for every subset of
generators. For a subset `S`, define

```text
r_S = Product[r_i, i in S]
W_S = r_S Y
A = Sum[A_T r_T, T subset of {1,...,s}].
```

Then

```text
W_S' = 1/2 Sum[P_i'/P_i, i in S] W_S
       + Sum[
           Product[P_i, i in Intersection[S,T]] A_T W_(S symmetric-difference T),
           T
         ].
```

All coefficients are rational. The boundary lift is

```text
W_S(x0) = r_(S,sigma)(x0) Y(x0),
```

and `W_empty = Y` is the physical output.

The solver dimension is `2^s d`. Square-class normalization is therefore a
correctness and performance requirement. The initial production cap should be
one independent radical. Degree four, corresponding to two independent
radicals, should be enabled only after sparse-block benchmarks and memory
measurements.

## 5. Branch and prescription semantics

Old DiffExp's square-root support did more than expand a matrix. It detected
radicands, auto-added `+i delta` prescriptions when the user supplied none,
gave user prescriptions precedence, and normalized literal matrix roots to
the selected sheet. See `DiffExp/MatrixLoading.m` and
`DiffExp/AnalyticContinuation.m`.

The companion lift still needs the first three behaviors. Literal matrix
flips are replaced by the typed boundary lift: choosing the sign of
`r_sigma(x0)` chooses the invariant sheet subspace of the rational system.

Required rules are:

1. Factor the square-free numerator and denominator of each radicand.
2. Treat every irreducible factor of odd valuation as a branch factor,
   including an odd valuation at infinity.
3. Add a system-scoped `+i delta` prescription for an unprescribed radical
   factor.
4. Let an explicit user prescription override the automatic one.
5. Preserve prescription provenance (`"User"` or
   `"AlgebraicRootAuto"`) in the system and plan.
6. Merge prescriptions using Config's sign-aware canonical orientation;
   sign-blind `Union` is forbidden.
7. Compute the anchor root through a typed branch evaluator. A generic
   Mathematica `Sqrt` at a negative real anchor would silently select the
   upper rim and is forbidden at this boundary.
8. Include normalized radicands, the anchor sheet signature, and effective
   prescriptions in solve, checkpoint, and preparation-cache hashes.

The intended real-axis convention matches the existing continuation contract:
Mathematica's principal value supplies the `+i delta` value, while a
`-i delta` prescription applies the lower-rim correction on the negative real
axis. The exact rule must be shared with the analytic-regularization branch
implementation rather than duplicated in the algebraic loader.

For a simple zero of `P`, the lifted block contains the residue `1/2` from
`P'/(2P)`. Existing fractional-sector crossing then transports the auxiliary
root factor with the same `ImSign` as the physical chart. No new crossing
formula is required.

## 6. Sheet constraint and matching

The rational lift contains all conjugate sheets. The chosen boundary lies in
one `d`-dimensional invariant subspace of the `2^s d`-dimensional rational
system. A residual check against the lifted ODE cannot detect a small
admixture of a different sheet because that admixture is also a valid lifted
solution.

Every match, regular value handoff, singular crossing, and final evaluation
must therefore certify

```text
W_S - r_(S,continued) W_empty = 0
```

for every nonempty subset `S`, at a deterministic interior probe where the
roots are finite. The continued root values must come from the same branch
state used by `ChartImSign` and `ApplyCrossing`.

The certificate must use the established uncertainty-aware residual scale.
Ambiguity or a residual above tolerance is a loud error. The implementation
must not repair a failed constraint by silently setting `W_S = r_S W_empty`.

Internal chart solutions must retain all lifted components until every match
and crossing is complete. Public projection to the first `d` components can
happen only afterward.

## 7. Compatibility with the existing recurrence

After the lift, the solver again receives a rational theta-form system

```text
D(t, eps) theta U = N(t, eps) U + S.
```

For a target sector `a + b eps`, a residue Jordan block
`(a_i + b_i eps) I + N_i`, Taylor order `n`, and log level `p`, the existing
block equation is

```text
d0 ( ((a+n-a_i) + (b-b_i) eps) I - N_i ) u_(n,p)
  = R_(n,p) - d0 eps u_(n,p+1).
```

Consequences:

- CASE T, CASE P, and CASE R are unchanged.
- epsilon-affine indicial roots remain a hard requirement.
- rational or half-integer shifts introduced by the algebraic lift are exact
  `a` tags already supported by SectorSeries.
- Jordan chains, true resonance, log sectors, pseudo-resonant epsilon shifts,
  and epsilon-window budgeting remain unchanged.
- inhomogeneous algebraic sources lift by the same subset multiplication rule
  as the matrix.
- an irregular lifted system or non-affine lifted spectrum must raise the
  existing exact error. It must not fall back to a dense sampled solve.

The C++ recurrence request schema need not change for this first slice because
the prepared lifted coefficients are rational/numeric exactly as today.

## 8. Supported classes and limitations

### 8.1 Recommended first class

- one epsilon-independent square root `r^2 = P(x)`;
- exact square-free rational or polynomial `P`;
- finite regular boundary anchor with `P(x0) != 0, Infinity`;
- matrix entries polynomial or affine in `r`, with base-field denominators;
- finite real transport lines, including a supported regular-singular branch
  crossing or branch endpoint;
- exact analytic-regulator dependence remains rational in the declared
  regulators.

This class includes roots analytic on every chart as well as simple
square-root branch points. It requires no global rational parametrization of
the curve.

### 8.2 Conjugate-sheet poles

General rational functions of a radical introduce norm denominators. For
example,

```text
1/(1+r) = (1-r)/(1-P).
```

The selected sheet can be regular at a zero of `1-P` while the conjugate sheet
is singular. The full rational lift then oversegments, may forbid an otherwise
regular anchor, or may be irregular because of a sheet the user did not
select.

This is a real conservative limitation of restriction of scalars. The first
slice should reject or explicitly classify sheet-selective norm poles in the
transport domain. A later implementation may cross a certified apparent norm
pole, but must not claim that the lift supports every sheet-specific
algebraic system.

### 8.3 Moving branch points

Radicands depending on epsilon or an unspecialized kinematic regulator move
their branch locus with the formal parameter. That conflicts with the fixed
exact singular alphabet used for segmentation and indicial analysis. Reject
these in the first implementation.

### 8.4 Algebraic constants and the C++ symbolic field

Fixed algebraic chart centers and coefficients can use the existing numeric
Acb route. Exact symbolic analytic regulators multiplied by algebraic
constants are not representable by the current C++ `Q(vars)` scalar domain;
such cases require the Wolfram backend or a later number-field coefficient
extension. This is distinct from algebraic *function* support and should not
be hidden by numerical specialization.

### 8.5 Infinity

For `s = 1/x`, an odd valuation of a radicand at infinity is a ramification
point. The companion matrix remains rational after pullback, so the local
recurrence needs no new mathematics. The public finite-line planner does not
currently own an infinity chart, however. Infinity support should first be
tested through an explicit transformed system and exposed only with a typed
infinity-plan API.

## 9. Alternatives and later optimizations

### 9.1 Local analytic-root recurrence

At a regular chart, write

```text
P(t) = Sum[p_n t^n]
r(t) = Sum[r_n t^n],  r_0 != 0.
```

The root coefficients obey

```text
r_n = (p_n - Sum[r_j r_(n-j), j=1..n-1])/(2 r_0).
```

This permits a direct rank-`d` dense-lag recurrence on a chosen sheet. It may
eventually be a useful optimization for regular value transport, but it is
not the first implementation because it introduces chart-local branch state,
order-growing convolution width, and a new C++ preparation path.

### 9.2 Puiseux uniformizers

At a ramification point of index `e`, choose `t = u^e` and solve the pulled
system in `u`. A local term

```text
u^(alpha + beta eps + n)
```

must be split by `n mod e` into physical-`t` sectors, with

```text
a_t = (alpha + residue-class)/e,
b_t = beta/e,
Log[u]^p = Log[t]^p/e^p.
```

For a general algebraic curve, selecting and normalizing a local place,
handling multiple Puiseux branches, and proving a convergence radius are
substantial new obligations. This path should be pursued only if profiling
shows that companion-lift dimension growth is unacceptable.

### 9.3 Rationalizing substitutions

`Sqrt[x-a]` is rationalized by `x = a + u^2`, and some quadratic radicands
admit a conic parametrization. A general square-free cubic or quartic defines
a genus-one curve; higher degrees are generally hyperelliptic. Automatic
global rationalization is therefore not a general solution.

Even in genus zero it introduces non-affine chart maps, inverse-branch
selection, transformed singular geometry, and integration Jacobians. Treat a
rationalizing substitution as a later typed optimization or a user-supplied
transform, not the default implementation.

### 9.4 General algebraic extensions

For a separable minimal polynomial `F(x,y) = 0`, one may use the companion
states

```text
Y, y Y, ..., y^(m-1) Y
```

and reduce products modulo `F`, with `y' = -F_x/F_y`. This again gives a
rational system over the base field. Near discriminant points a power basis
can introduce severe apparent poles, so normalization and integral bases are
required before this becomes production-ready.

## 10. API and record implications

An algebraic system should preserve the original physical dimension while
carrying an internal lifted solver system. A suitable record is

```text
DiffExp2.AlgebraicSystem/v1
  OriginalMatrix
  Variable
  Dimension
  SolverSystem
  SolverDimension
  AlgebraicStructure
    Kind
    Generators
    SquareClassBasis
    CoefficientBlocks
    BranchFactors
    AutoPrescriptions
  PhysicalComponentIndices
  ConstructionCertificate
```

Plans and results additionally need

```text
AlgebraicSignature
AnchorBranchValues
EffectiveDeltaPrescriptions
SheetConstraintDiagnostics
```

Public operations should behave as follows:

- `TransportLine` lifts the boundary, transports internally, and projects
  physical values only after completing each internal crossing/match.
- `ExactSectors` and `LocalBehavior` expose the physical component slice.
- `EndpointLimit` applies to the physical slice.
- `IntegrateLine` pads the physical coefficient vector with zeros for all
  auxiliary components and integrates the full internal chart object.
- retained/precomputed charts must preserve their lifted local solution for
  future transport and integration even if the public segment also carries a
  physical projection.
- checkpoint metadata must distinguish field generators and branch sheets.

A lifted sector may have exact zero leading Taylor columns in the physical
slice while an auxiliary component starts at its base exponent. Public
projection may retain that honest zero prefix. It may shift the sector's `a`
tag only after proving the entire discarded physical prefix is exactly zero.

## 11. Failure modes

New failures should be structured and loud:

- unsupported nonquadratic or undeclared algebraic dependence;
- inexact radicand or coefficient input;
- epsilon-dependent/moving branch locus;
- branch anchor at a radicand zero or pole;
- noninvertible denominator in the quotient field;
- dependent generators not reducible to a certified square-class basis;
- requested field degree above the configured production cap;
- missing or conflicting effective prescriptions;
- sheet constraint residual above tolerance;
- conjugate-sheet norm pole at an unsupported anchor or path point;
- lifted irregular singularity;
- lifted non-affine epsilon spectrum;
- C++ exact symbolic coefficient outside `Q(regulators)`.

None of these may trigger a Wolfram/C++ fallback, root sampling, or branch
guess based on floating-point sign.

## 12. Staged implementation

### Stage A: exact single-root construction

1. Add an isolated algebraic-lift module above `API.m`.
2. Detect and normalize one square root.
3. Reduce every matrix entry to `A0 + r A1` exactly.
4. Construct the rational companion matrix.
5. Verify the intertwining certificate exactly.
6. Add unit tests only; do not expose it publicly yet.

### Stage B: one-root public transport

1. Restrict to the class in section 8.1.
2. Implement system-scoped automatic prescriptions with user precedence.
3. Add the typed branch-aware boundary lift.
4. Add physical result projection.
5. Add sheet-constraint certificates at every handoff.
6. Verify both Wolfram and C++ recurrence backends.

### Stage C: public-operation closure

1. Endpoint limits and singular branch endpoints.
2. Line integration and precomputed-chart reuse.
3. Piecewise evaluation and exact physical sector inspection.
4. Checkpoint/cache fingerprints.
5. Feynman-trick runner compatibility.

### Stage D: degree-four multiquadratic support

1. Exact square-class basis reduction.
2. Sparse subset-block construction.
3. All subset sheet constraints.
4. Performance and memory gates before enabling by default.

### Stage E: wider algebraic scope

1. Norm-denominator and conjugate-only-pole classification.
2. General rational functions of radicals where the full lift remains
   regular singular.
3. General separable algebraic companion systems.
4. Only then consider native Puiseux or structured C++ optimizations.

## 13. Acceptance fixtures

The minimum focused suite should include:

1. `Y' = Sqrt[x] Y`, transported entirely in the positive region.
2. The same equation from `x=1` to `x=-1` across the branch point:
   - `+i delta` oracle `Exp[(2/3) (-I - 1)]`;
   - `-i delta` oracle `Exp[(2/3) ( I - 1)]`.
3. The singular endpoint `x=0`, with `Y(0)=Exp[-2/3]` and `Z(0)=0`.
4. `Y' = Y/Sqrt[x]`.
5. A denominator such as `1/(1+Sqrt[x])`, demonstrating a conjugate-only
   norm pole and the required diagnostic.
6. An epsilon-affine singular system with a radical regular term.
7. A true-resonant/Jordan system with a commuting radical factor.
8. A pseudo-resonant different-`b` collision with a radical factor.
9. An irregular case such as `Y' = x^(-3/2) Y`, which must fail loudly.
10. Exact analytic-regulator parity between Wolfram and C++ where the lifted
    coefficients remain in `Q(regulators)`.
11. Two roots, for example `A = Sqrt[x] + Sqrt[1-x]`.
12. Dependent roots such as `Sqrt[x]` and `Sqrt[4 x]`.
13. Coincident branch factors with consistent and conflicting prescriptions.
14. Odd and even radicand valuations at infinity after an explicit `s=1/x`
    pullback.
15. Public `ExactSectors`, `EndpointLimit`, and `IntegrateLine` projection
    tests.

Property sweeps should generate small exact `P`, `A0`, and `A1` and verify:

- the exact intertwining identity;
- absence of line-variable radicals in the lifted matrix;
- original-sheet versus lifted numerical parity;
- the sheet constraint at every retained chart;
- Wolfram/C++ parity;
- strict preservation of epsilon windows and existing recurrence failures.

## 14. Prototype evidence

The companion lift was tested without modifying package code on

```text
Y' = Sqrt[x] Y,  Y(1) = 1,
B = {{0, 1}, {x, 1/(2 x)}}.
```

At working precision 120 and expansion order 60, transport from `1` to `-1`
agreed with both prescribed-sheet oracles to approximately `1.25*10^-33`.
The Wolfram and C++ recurrence backends agreed, and a clean `HEAD` archive
independently reproduced the Wolfram result. This is evidence that the
existing fractional-sector crossing machinery can carry the companion root
state correctly; it does not replace the staged acceptance suite above.
