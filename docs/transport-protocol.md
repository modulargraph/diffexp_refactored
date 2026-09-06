# Native differential-equation transport

`diffexp transport request.json` (or `transport -` for stdin) integrates a
matrix differential equation along an explicitly supplied path. It is the
same interface used by the Mathematica wrapper and has no example-name input.
The input and output each contain one JSON object; progress goes to stderr.

Run the two small complete examples:

```sh
build/diffexp transport examples/transport/ordinary.json
build/diffexp transport examples/transport/partial-boundary.json
```

The first solves `dy/dt = y/(1-t)`, from `t=0` to `t=1/2`, with `y(0)=1`.
Its endpoint is 2. The second fixes only the first component of a coupled
regular-singular system by `y1=x+O(x^2)`; Frobenius matching determines both
components and transports them to `(1,1)`.

## Request

The schema is `DiffExp.Transport/v1`. Matrix positions and epsilon powers are
zero-based. The retained epsilon window is `0..epsilon_order`.

| Field | Meaning |
| --- | --- |
| `dimension` | Number of components |
| `paths` | Object mapping each kinematic variable to an expression in `x`; native transport runs from `x=0` to `x=1` |
| `entries` | Nonzero differential-equation entries |
| `boundary` | Component-by-epsilon matrix of numerical expression strings at the path start |
| `boundary_errors` | Optional nonnegative real absolute uncertainties with the same shape |
| `basis_prefactors` | Optional algebraic normalization expression string per component; convert continued roots to the principal endpoint basis |
| `asymptotic` | Partial power/log constraints at `x=0`, replacing `boundary` |
| `initial_only` | With `asymptotic`, return a regular numerical seed at a small positive `x` |
| `taylor_order` | Retained local series order; default 50 |
| `working_bits` | Binary working precision; default 384 |
| `accuracy_goal` | Requested decimal digits under the reported error estimate; 0 disables this check |
| `recurrence` | `auto` (default) selects adaptive spectral or local transport; `spectral` requires adaptive collocation; `taylor` selects local finite-lag/Taylor; `series` forces Taylor convolution |
| `division_order` | Local step control; default 4 |
| `save_segments` | Export local retained Taylor coefficients; default false, with a 64 MiB estimated output budget |

An ordinary entry is, for example,
`{"row":0,"column":1,"epsilon":2,"variable":"t","expression":"1/(1-t)"}`.
It contributes `epsilon^2/(1-t)` to the `dt` connection. Native pullback
multiplies by the derivative of `paths.t`.

An entry with `"variable":"form"` explicitly
supplies an already pulled-back coefficient of `dx` when no path variable named
`form` exists. A path variable of that name keeps its ordinary derivative
semantics; rename it to use the explicit marker. Supplied forms accept an optional exact rational
`coefficient`; its expression is evaluated along the request paths but is not
multiplied by another path derivative. This preserves a sparse constant-matrix
sum of shared scalar one-forms without combining radical expressions per matrix
position. Supplied forms at a singular initial point require explicit asymptotic
boundary conditions; finite-boundary dlog initialization remains available.

A logarithmic entry uses `"variable":"dlog"`, its letter in `expression`,
and an optional exact rational string `coefficient`. The native code takes
the logarithmic derivative along the supplied path. Compatible epsilon-linear
systems first attempt an exact diagonal root rescaling to a rational connection.
Compatible systems use shared rational-product recurrences with O(N) arithmetic
cost at fixed polynomial degree, precision and epsilon depth. The plan is reused
across charts. Noncanonical systems, inconsistent root gauges, preparation-budget
limits, low orders relative to polynomial degree, and ill-conditioned charts
retain the previous series solver. The selector uses the supplied mathematics,
not a family name.

The same default applies through the Mathematica process wrapper. Saved segments
continue to use the series solver so their exported coefficients remain in the
physical basis. Singular starts retain the existing boundary treatment. Tail
estimates are transported back through the gauge and fed into the existing
accuracy checks; they remain estimates, not certified infinite-tail bounds.
Before local preparation, `auto` may attempt [adaptive spectral transport](adaptive-spectral.md)
for ordinary epsilon-linear systems with at least 16 components, positive
accuracy goal and no saved segments. Ordinary derivative entries are shared
internally when eligible. Spectral transport continues independent polynomial
square roots and honors endpoint `basis_prefactors`. It does not require an
exact rationalizing gauge. Explicit `spectral` also permits smaller systems.
The `spectral` output object records whether it was attempted/accepted, tested
resolutions and the rejection reason. `spectral_charts` is one on success.

The output `recurrence` object reports accepted finite-lag and series chart counts,
rejected candidates in adaptive checking, and any preparation fallback reason.

Exact expressions support integers, rationals, arithmetic, `I` and independent
polynomial `Sqrt[...]` radicands. Paths must fix all matrix variables. The native
parameter `x`, `I`, and internal root symbols `r0`, `r1`, etc. are reserved.
Nested algebraic towers and nonpolynomial radicands are rejected. Decimal
boundary input may carry Mathematica precision (single-backtick) or absolute
accuracy (double-backtick) annotations, including uncertain zeros; these become input
uncertainty, not exact rational data. Responses are decimal strings and never
require evaluation as code.

## Algebraic basis conventions

Roots in the differential equation are continued along the entire path. If the
integrals are defined with algebraic factors evaluated on the principal sheet at
each physical point, supply those factors in `basis_prefactors`. For a component
`F_i = g_i J_i`, the returned value is multiplied by
`g_i(principal endpoint roots) / g_i(continued endpoint roots)`. The initial
boundary must use the principal starting basis. This is explicit mathematical
input; the solver never infers factors from reference values or a family name.
Use `"1"` for components with no algebraic normalization.

For example, a component with prefactor `"2+Sqrt[t]"` needs this conversion
when the path winds around `t=0`. Roots appearing only in the prefactors are
also continued. Endpoint factors must have finite values and a nonzero
continued divisor. Exact endpoint coordinates are substituted before numerical
root evaluation to avoid artificial branch-cut uncertainty from cancellation.
Asymptotic initialization currently requires rational input, including its
prefactors; use algebraic prefactors on the subsequent regular transport.

For a completed ordinary transport with prefactors supplied, `values` and
`errors` use the principal endpoint basis and the response says
`basis_convention: "principal_endpoint"`.
Optional saved Taylor segments remain in their continued basis, marked by
`segments_basis_convention: "continued"`; both endpoint root arrays are
reported. Without prefactors, the returned values retain the continued basis.

## Partial asymptotic boundaries

Each `asymptotic.constraints` entry specifies `row`, `epsilon`, exact rational
`power`, integer `log_degree`, and numerical `value`. Repeated entries at the
same monomial are summed. Each `asymptotic.cutoffs` entry gives a `row` and a
power below which unlisted power/log coefficients vanish. Explicit terms above
a cutoff are still constraints. Omit a row to leave it initially unknown.

For example, `y1=x+O(x^2)` is a coefficient 1 at row 0, epsilon 0, power `"1"`,
log degree 0, together with cutoff `"2"` at row 0. Exact linear elimination in
the Frobenius frame must determine all constants; an underdetermined or
inconsistent boundary returns an error. General asymptotic matching currently
requires a rational connection. No Wolfram kernel participates in the solve.

## Response and accuracy

The response schema is `DiffExp.TransportResult/v1`. Each component/epsilon
entry has decimal-string `real_midpoint`, `imaginary_midpoint`, and `radius`.
The reported radius bounds retained ball arithmetic with injected local
truncation estimates; `errors` reports the corresponding estimates separately.
`parameter` is `"1"` at the requested endpoint, or the seed parameter for an
`initial_only` request. The response includes chart count and preparation,
numerical and total seconds. Initial-only responses report preparation and
total time.

`AccuracyGoal` uses an absolute-plus-relative criterion
`estimated_error <= 10^(-goal) (1+abs(value))`. The engine can shorten individual
steps within finite retry budgets; it does not rerun the complete calculation
or silently relax the goal. Insufficient input precision or excessive propagated
uncertainty yields an explicit error. The last retained terms estimate the
omitted tail; they do not certify it. Consequently `omitted_tails_certified`
is false. See [validation](validation.md) for the independent comparisons.

The Mathematica wrapper internally renames reserved native variable names, so a
kinematic `x` can coexist with a separately configured `LineParameter -> t`.
