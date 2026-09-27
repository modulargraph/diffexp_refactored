# Verified epsilon windows for homogeneous uncertainty maps

[UncertaintyWindow.lean](lean/UncertaintyWindow.lean) proves the finite exact-rational contract behind shortening a homogeneous uncertainty map. The module passed Lean 4.34.0 kernel checking with no `sorry`, custom axioms, or native proof evaluation. It modifies no C++ or earlier Lean modules.

## Exact statement

For a retained epsilon width `W`, output degree `k < W`, arbitrary rational map coefficients `T[e]`, and arbitrary rational incoming uncertainty coefficients `δ[d]`, define the finite convolution by summing

```text
if e <= k then T[e] * δ[k-e] else 0
```

for `e = 0, …, W-1`. The guard is inside the definition and precedes subtraction; out-of-range terms are zero rather than interpreting saturating natural subtraction as a legitimate index.

Assume **every** incoming uncertainty coefficient below `L` is exactly zero: `∀ d < L, δ[d] = 0`. Then summing only `e < W-L` produces exactly the same coefficient for every `k < W`.

When `L < W`, the retained map degrees are `0, …, W-1-L`. When `L ≥ W`, the retained map is empty. The implementation-safe formulation is a count `W-L` using a checked/saturating difference, not an unconditional unsigned expression `W-1-L`. The proof handles both `L=0` and `L≥W`.

The key inequality is explicit: if `k<W`, `e≤k`, and `e≥W-L`, then `k-e<L`. Thus every omitted admitted term multiplies a proved zero incoming coefficient. `sumN_shorten` proves a finite sum may omit all these zero terms; this is not merely an assertion that a truncated candidate is valid.

## Checked theorems

| Theorem | Result |
| --- | --- |
| `omitted_input_below_low` | The integer dependency implication, with its subtraction guard. |
| `truncate_uncertainty_map` | Full and shortened finite convolutions agree at every retained output degree. |
| `unknown_map_independent` | Two maps equal only below `W-L` give equal retained output coefficients; their higher coefficients may differ arbitrarily. |
| `truncate_each_column` | A finite matrix row can use its own bound `L_j` for each column, provided **all** lower incoming uncertainties in each column vanish. |
| `affine_decomposition` | Finite matrix-series convolution satisfies `T(m+δ)+f = (Tm+f)+Tδ`. |
| `affine_with_truncated_uncertainty` | Combines the full midpoint/forcing calculation with shortened homogeneous uncertainty columns, preserving every retained coefficient. |

`matrixConvolution` is a finite sum over input columns for one arbitrary output row. Applying the universally quantified theorem to each row gives the finite-vector result. The coefficient-level affine forcing is arbitrary; in particular it is added only once to the midpoint calculation. The homogeneous map does not receive a second copy of the forcing.

The shortened definition physically sums only `W-L` terms. `unknown_map_independent` separately establishes that no value assigned to an unknown higher coefficient can affect the retained exact result under the hypotheses. This is stronger than filling unknown map coefficients with zero and testing a particular example.

## Production boundary

The application-time partial-map marker must enforce the proved predicate on **all incoming uncertainty coefficients**, not on the midpoint or on a sample of coefficients. A numerically small coefficient is not an exact structural zero. If the predicate fails, this proof does not authorize use of the shortened map.

The model uses nonnegative epsilon degrees after any required normalization, arbitrary exact `Rat` coefficients, a common retained output width, and the same map `T` in the full and shortened expressions. It assumes no reduction in arithmetic precision. It neither chooses precision nor proves an ACB enclosure, rounding-error identity, tail bound, endpoint contract, compiled-loop refinement, or cache-marker correctness. If the actual map has negative epsilon powers, a separate valuation/shift argument is required before applying this theorem.

The equality concerns exact retained values. Reordering finite-precision arithmetic can change interval widths even when exact values agree. Numerical acceptance and enclosure guarantees remain the responsibility of the unchanged arithmetic/accuracy policy and its implementation validation.

## Replay

From the package root, using the pinned Lean 4.34.0 toolchain:

```sh
LEAN_BIN=/path/to/lean scripts/check_recurrence_proofs.sh
```

This checks all packaged modules and writes intermediate files to a temporary
directory that is removed afterward.

Exported theorems report only Lean's standard `propext`, `Classical.choice`,
and `Quot.sound`; the integer implication uses only `propext` and `Quot.sound`.
There are no custom axioms or incomplete proofs. The toolchain and standard
foundations remain the trust base; see
[Lean's axiom documentation](https://lean-lang.org/doc/reference/latest/Axioms/).
