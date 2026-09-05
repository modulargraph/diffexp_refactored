# Native equation generation and reduction

`ibp.hpp` compiles quadratic propagators into exact integration-by-parts (IBP)
identities without Wolfram or a subprocess. Its finite exact reducer is useful
for small systems and for checking source-identity witnesses. `fire.hpp` now
provides native start/configuration export, bounded invocation of the external
FIRE C++ executable, and data-only import of exact reductions. No Mathematica
preparation or evaluation is required.

`level_preparation.hpp` batches unresolved requests and derivative demands until
a selected ordered basis closes. `level_cache.hpp` durably stores that exact
closure and rechecks its mathematical witnesses on reuse. These are preparation
operations; they do not claim large-family completion or a numerical FT result.

## Mathematical representation

The loop scalar products are `l_i.l_j` and `l_i.p_j`; external scalar products
come from an explicit symmetric Gram matrix. Each propagator is an affine form
in these independent loop scalar products, with exact rational-function
coefficients. The sign convention is `D=m^2-q^2`, matching the Symanzik and gamma
boundary modules. Feynman merging acts on entire quadratic forms.

Physical denominator slots are distinct from auxiliary numerator slots.
Explicit numerator forms keep their order, which is necessary for Henn's
published scalar targets. Exact rank selection fills missing auxiliary slots
deterministically, and one exact inverse rewrites scalar products in that basis.
Dependent physical propagators require partial fractions and fail explicitly;
they cannot be silently discarded during rank completion.

For each loop derivative `l_i` and loop/external vector `v`, the compiler
rewrites `v . dD_k/dl_i = c_k0 + sum_j c_kj D_j`. The identity is

```
0 = delta_(v,l_i) d I(a)
    - sum_k a_k [ c_k0 I(a+e_k) + sum_j c_kj I(a+e_k-e_j) ].
```

Parameter derivatives use the same shifted-integral representation. This API
differentiates propagator parameters at fixed dimension and fixed external
Gram matrix. Differentiating kinematic invariants needs its own vector-field
operator; changing the Gram matrix through this API is rejected.

## Finite work and exact evidence

Seed ranges have explicit bounds on extra denominator powers, total numerator
degree, seed count and equation count. Generated equations retain **every**
term, including integrals outside the seed range. Uneliminated terms remain
unresolved; the reducer never declares them a complete master basis.

Sparse exact elimination uses one strict integral ordering. Every row retains
its exact linear combination of original equations. A returned reduction can
be checked by reconstructing `input - remainder` from that witness. The
selected-basis compiler verifies linear independence, expresses each parameter
derivative in that basis, then rechecks every resulting equation against the
original identities. Insufficient identities or insufficient basis elements
produce an explicit closure failure.

This establishes the closure needed to transport the selected integrals. It is
not a proof of the global master count or of coverage of an unrequested
observable. The level-preparation loop includes the requested observable span and its
derivative demands. Henn must supply the full 257-target union before its
all-component pipeline can be accepted.

## C++ example

```cpp
#include <diffexp/ibp.hpp>
#include <diffexp/families.hpp>
using namespace diffexp;

ExactField field({"x", "eps"});
Exact x(field, "x");
auto example = feynman::example_family("bubble");
auto family = ibp::quadratic_family(example.momenta, x, example.physical_count);
ibp::Generator generator(ibp::PropagatorBasis(ibp::merge(family, 0, 1, x)),
                         Exact(field, "2-2*eps"));
ibp::ExactReducer reducer(x, 1000);
ibp::for_each_seed(1, 2, {2, 2, 100}, [&](const ibp::Integral& seed) {
  for (auto& identity : generator.relations(seed))
    reducer.insert(std::move(identity));
});
auto system = ibp::differential_system(generator, reducer, {{2, 0}}, 0, x);
// system.matrix[0][0] == -(1+eps)*(1-2*x)/(1+x-x*x)
```

The complete bubble runner uses native generation/reduction and now certifies
epsilon coefficients 0–4, including the solver's omitted Taylor tails.
Independent rigorous quadrature verifies the same window. The specialized
sunrise runner executes both FT levels at `d=2-2*eps`; its finite part agrees
with the certified Bessel oracle by about `1.17e-24`. Its scalar boundaries are
certified, while its full transport/endpoint remainder is not yet bounded.

## Native FIRE and level closure

A provider receives a shared batch of typed integral vectors. Reduction rows
become complete exact relations; omitted requested rows fail explicitly, and
zero reductions remain actual relations rather than missing data. Candidate
masters and derivative demands grow until exact coordinate maps exist for the
selected basis and every requested target. Limits cover passes, masters,
demands, source equations and total work time. A finite failure reports its
reason instead of treating unresolved integrals as a successful result.

Every differential row and observable row retains a source-identity witness.
The resulting proof is **conditional on the imported IBP identities**. Native
verification reconstructs those linear combinations; it does not independently
re-prove every identity emitted by an external provider or assert a global
master count.

The `prepare FAMILY` command exposes exact graph preparation, with optional
`--fire PATH`, `--cache DIRECTORY`, and the Henn-only `--henn-basis FILE`.
A persistent cache directory is recommended when upstream work must survive
system temporary-directory cleanup. Both verified Banana4 preparations
closed 15-, 7- and 3-master levels before reaching the scalar leaf; numerical
FT evaluation remains pending.

## Henn observables and numerator pullbacks

`henn_observables.hpp` reads `dlogBasisXB.txt` as data. It imports all 108
observables, preserves the 11-slot index vectors, and retains all 257 distinct
source targets, including 174 with negative powers. Exact substitution at X0
leaves 241 nonzero targets. The full source union is retained for basis discovery.

X0 is `{s12,s23,s34,s45,s15}={3,-1,1,1,-1}`. Coefficients use the exact algebraic
field `Q(eps)[eps5]/(eps5^2+3)` with `eps5=+I*Sqrt[3]`. The per-integral factor
`(-1)^Sum(indices)` converts the paper's `q^2` convention to native `-q^2`.
`henn_boundary.hpp` applies the final factor `eps^4 Exp[2 eps EulerGamma]`
during numerical observable reconstruction. It derives the scalar epsilon
lookahead from the exact canonical coefficients and retains the full negative
window for a separate pole audit. The `ft --henn-basis` command uses this map;
full numerical Henn reconstruction has not yet passed its acceptance gate.

The first canonical component is `+2` times the native scalar integral with
indices `{1,1,0,1,1,0,1,0,0,0,0}`. Tests reconstruct every original canonical
expression independently after scalar assignments. Nonlinear products,
unsupported functions, bad arities and positive genuine-numerator slots fail.

`merge_pullback.hpp` expands negative powers into exact source-integral rows
in the next complete propagator basis. It preserves cancellations, beta
normalization, and the distinction between integration, endpoint and direct
operations. All Henn targets are covered by the pullback regression. A deepest
single-propagator Gaussian reducer supports bounded arbitrary numerator
polynomials. Neither operation by itself completes Henn's numerical recursion.

## Exact cache boundary

`cached_level::prepare` returns the ordinary `level::Result` together with
cache-hit and content-identity information. Its scientific key covers the exact
Gram matrix, scalar-product ordering, ordered denominators and physical count,
dimension, parameter, exact-field symbols, ordered requests, normalization,
branch, and explicit proven-zero-sector assumptions. Runtime budgets and the
FIRE executable path are excluded.

A hit restores exact coefficients into the caller's field and regenerates
parameter derivatives. It verifies basis independence and every stored source
witness for both matrix and target rows. This permits reuse with changed
receiving epsilon/Taylor/precision demands or an unavailable FIRE executable.
Corruption or failed mathematical reconstruction aborts reuse explicitly.
The durable metadata uses certificate type `exact`, scope `exact_closure`,
and states its dependence on imported identities. Requested numerical precision
is never promoted to a numerical certificate.
