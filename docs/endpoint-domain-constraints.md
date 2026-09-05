# Physical endpoint domain constraints

The full Henn trial `ft-henn-all108-verified-columns-acceptance.json` reused
all six exact systems and independently verified its saved69-column lower
endpoint, then rejected its first endpoint operation because an unregulated
fixed logarithm exists in the unrestricted fundamental matrix.

A read-only diagnostic using that same saved endpoint finds69 affected lower
endpoint operations (24 one-sided limits,45 beta integrals). Every offending
coefficient involves only fundamental columns28,30,32 (zero-based). The
sectors are fixed x^-1 terms for integrals and fixed x^0 log(x) terms for
limits. The first operation has logarithmic row
`3 eps c28 - eps c30 + (2 eps+2)c32/3`.
This diagnoses a domain issue; it does not prove that the physical boundary
annihilates these rows. That remains a required numerical check.

For an observable expansion P F c, fixed divergent terms impose linear
conditions on the matching constants c. `AffineFrobeniusSeries::dr_domain`
returns the admissible expansion and separate exact coefficient rows indexed
by the output row, x power and logarithmic degree. Contributions from distinct
fundamental columns are combined within each condition. The existing strict
DR endpoint and integration APIs still reject fixed divergences without a
physical boundary.

The recursive evaluator carries each condition as an auxiliary functional.
Lower and upper conditions have distinct output coordinates. They are
transported and contracted with the same immutable leaf source as the
requested integrals; no new independent copy of a boundary is introduced.
Every retained Laurent coefficient of each condition must have finite absolute
bound at most the configured tolerance (default1e-20). Only then are the
auxiliary rows removed from the returned values and parent expression. The
value route checks the corresponding rows on its matched physical constants.
This is a numerical consistency check through the computed epsilon window,
not a theorem of exact vanishing or a bound on omitted coefficients. General
recursive results remain explicitly uncertified for omitted tails.

Observable projection combines contributions directly into exact coefficient
coordinates instead of first allocating an expanded list of duplicate terms.
The retained-coordinate budget and a separate finite multiplication budget
still apply. This avoids rejecting a small final expansion merely because its
uncombined convolution was large.

`project_endpoint_domain` computes only fixed sectors with nonpositive powers:
these contain every constant and domain condition needed by the DR endpoint
functional. It requests only the necessary Taylor coefficients of each
observable row. Beta integrals and matching frames continue to use the full
requested series. The selected endpoint expansion carries no coherent-frontier
or Wronskian metadata. The endpoint-specific projection is covered by exact
coefficient regressions and the physical-boundary tests described below.

The read-only Henn projection preflight reuses all saved columns without a
recurrence. It deliberately does not repeat the independent ODE verifier and
cannot publish an integral or an exact certificate. The production evaluator
still independently verifies the cached endpoint before numerical use.

The standalone stage-export helper currently refuses stages with undisclosed
physical-domain conditions; use the full evaluator for those stages. Stages
with no such conditions retain their existing API and cache representation.

The analytic regressions exercise constant physical solutions embedded in
larger systems with fixed power or logarithmic modes, all three numerical
methods, lower/upper/beta operations, and physical boundaries that genuinely
excite the forbidden modes. Equal and opposite divergences at separate
endpoints must also fail, rather than cancel as a principal value.

Section3.2 of [Hidding and Usovitsch](https://arxiv.org/pdf/2206.14790)
describes DR endpoint extraction from physical solutions and reports the
absence of unregulated divergent sectors empirically for their examples.
It does not imply that every unrestricted homogeneous solution has that
property. The domain formulation above is an implementation inference from
linear ODE structure, tested independently of the Henn reference table.
