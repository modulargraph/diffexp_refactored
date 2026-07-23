# Example recovery log

This is the short causal ledger for the opt-in end-to-end examples. Each entry
records the observed failure, the controlled change, and the result so that a
later investigation does not repeat an already-disproven retry.

## Massless pentagon

| ID | Observed issue | Controlled change | Result |
|---|---|---|---|
| P-01 | One global matching-digit setting made a lower consumer stricter while leaving its producer unchanged. | Bind matching digits per producing level and bind checkpoints to the producer's value. | Producer retries now tighten only the preceding level. |
| P-02 | Tightening one producer exposed the same uncertainty at the next boundary. | Enforce a transitive two-digit safety ladder through all earlier producers. | Levels 2, 3, and 4 stabilize at 12, 14, and 16 digits for the eight-digit public profile. |
| P-03 | A singular rational chart had a passing midpoint and uncertainty only in its solved weights. | Retry locally with the retained exact Rational shadow and certify the physical residual. | The level-2 upper singular match passes without another upstream retry. |
| P-04 | Level-1 match 12, a regular receiver, had basis=pass, incoming=pass, weights=inconclusive. | Re-solve once with zero-radius midpoint data as a proposal, retaining the untouched full-ball forward residual as publication authority. | Match 12 passes. A negative unit test proves a midpoint candidate is rejected when the full-ball residual remains inconclusive. |
| P-05 | After P-04, level-1 match 24 honestly reached only epsilon power 2 while power 3 was required. | Increase only the private level-1 matching halo from 2 to 3. | Fresh transport completed in 92.67 seconds with finite value `-0.0254117799`; no discovery retry was used. |

Current pinned production profile: working precision 300, Taylor order 25,
matching digits 8, private halos `<|1 -> 3, 2 -> 1|>`, and producer digits
`<|2 -> 12, 3 -> 14, 4 -> 16|>`.
