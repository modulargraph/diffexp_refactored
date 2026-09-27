# Verified contracts for avoiding redundant fallback solves

Two new Lean modules close a specific gap in the earlier review and establish a conditional early-accept policy. Both were kernel-checked on 27 September 2026 with pinned Lean 4.34.0. They contain no `sorry`, custom axioms, or native-evaluation shortcut. They do **not** verify the production C++ or its ACB arithmetic.

## Actual finite polynomial quotient, rather than an abstract tail

[FiniteQuotient.lean](lean/FiniteQuotient.lean) now connects the earlier general triangular recurrence to an actual finite list of denominator monomials. The earlier report's arbitrary-tail limitation is superseded **for this rational-coefficient model**, not for complex fields or production compilation.

The denominator representation is

`Q(t,ε) = q + Σ term.coefficient * t^term.a * ε^term.b`,

where `q ≠ 0`, exponents are natural numbers, and every listed monomial is nonconstant. Repeated monomials are summed, so support need not be canonicalized. A positive epsilon width `width` maps a flattened index `n` to `(n / width, n % width)`. This is a window in **shifted** epsilon order: actual negative Laurent orders are not discarded; the caller must prove the lower bound and apply its shift consistently.

`offConvolution` explicitly traverses the finite monomial list. A monomial contributes its coefficient times `W[(row-a)*width+(column-b)]` if both indices are in bounds, and contributes zero otherwise. It does not contain an assumed abstract convolution oracle. Every admitted predecessor is proved smaller than `n` using the earlier flattened-order theorem. `tail_is_convolution` proves that the bounded predecessor accessor used by `solve` implements precisely this convolution.

The checked results are:

- `quotient_equation`: the constructed stream satisfies `q*W[n] + offConvolution(Q,n,W) = Y[n]` for every stored cell, hence every retained coefficient of `QW=Y` in this representation.
- `quotient_unique_prefix`: any stream satisfying those equations only below `cutoff` equals the constructed stream below `cutoff`. Setting `cutoff = TaylorRows * width` gives a finite rectangular Taylor/epsilon window. No equation above the cutoff is assumed.
- `source_prefix_independent`: changing unknown source coefficients above the cutoff cannot change the solution prefix.

This is a useful proof boundary for replacing a second exact quotient implementation: once a backend is connected to the represented finite support and recurrence, agreement follows from the universal theorem. It is not merely a comparison on a few sample coefficients.

Remaining boundaries are explicit: the exact compiler must produce the represented polynomial, the C++ loops must implement the checked recurrence and indexing, epsilon valuation normalization must be correct, the coefficient field must match, and machine arithmetic must enclose the retained exact solution. The module proves no ODE integration step, complex-field extension, analytic radius, omitted Taylor/epsilon tail, or endpoint condition. Finite-window source independence also does not prove independence when the epsilon **width** changes; that requires a coordinate-preserving embedding theorem.

## A second enclosure can tighten a result without being required for correctness

[EnclosurePolicy.lean](lean/EnclosurePolicy.lean) models a scalar retained exact value `x : Rat` and rational-endpoint intervals. Its fixed acceptance rule is absolute interval width at most a specified tolerance.

The substantive interval proofs establish that:

1. Intersecting two intervals enclosing the **same** retained value still encloses it.
2. Intersection raises the lower endpoint and lowers the upper endpoint, so its width cannot exceed the candidate's width.
3. An accepted candidate therefore remains accepted after intersection.
4. A lazy selector returns the candidate whenever it is accepted; its result is independent of every possible fallback result on that branch.
5. The full selector is sound if the candidate encloses `x`, and a rejected candidate's fallback encloses `x` and meets the same acceptance policy. Crucially, the fallback contract is required **only when the candidate is rejected**.

`second_enclosure_optional` combines the enclosure and width results: both stopping with an already accepted candidate and intersecting it with a second sound enclosure meet the stated correctness contract. It does not assert equal outputs, equal error bars, or equal acceptance rates. Small exact witnesses `[0,4]` and `[1,2]`, both enclosing `x=1` at tolerance 4, are reduced with `decide +kernel`; this does not invoke native execution as a proof oracle.

For the production policy, this supports the following design **conditional on an implementation enclosure argument**:

```text
compute compact candidate
check unchanged acceptance policy
if accepted:
    return candidate
otherwise:
    run certified fallback/recovery
    check the same acceptance policy before returning
```

An unconditional second solve/intersection can improve tightness or diagnose implementation defects. It is not logically necessary for a contract that already follows from the first candidate's enclosure and acceptance. Conversely, these proofs do not justify accepting a rejected candidate, relaxing tolerances, assuming a residual is a full error bound, or removing a fallback that supplies a missing enclosure or tail certificate.

The checked width policy is deliberately precise. Relative-accuracy-bit tests can depend on the midpoint; the intersection-width theorem must not be advertised as monotonicity of an arbitrary relative-accuracy predicate. The early-return principle applies when the actual candidate has already met the unchanged required predicate, but a direct formal connection to the production predicate remains to be proved. Likewise, real/complex ACB balls are not represented by rational scalar intervals in this module.

## Exact replay and results

Run `scripts/check_recurrence_proofs.sh` with Lean 4.34.0. An existing toolchain
can be selected with `LEAN_BIN=/absolute/path/to/lean`. The script uses an isolated
temporary module directory and leaves no compiled proof files in the source tree.
No mathlib dependency is needed.

The modules report their transitive dependencies with `#print axioms`.
Reported dependencies are limited to Lean's standard `propext`,
`Classical.choice`, and `Quot.sound`; `intersection_encloses` has none.
The examples use kernel reduction. See Lean's
[axiom documentation](https://lean-lang.org/doc/reference/latest/Axioms/)
for the meaning of this audit.
