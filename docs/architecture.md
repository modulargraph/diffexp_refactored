# Architecture

DiffExp 2.1 separates exact equation preparation, finite series construction,
numerical continuation and scientific acceptance. A cache hit means that an
artifact's mathematical checks passed; it does not turn a requested precision
into an error certificate.

| Layer | Main headers | Responsibility |
| --- | --- | --- |
| Exact algebra | `exact.hpp`, `univariate_rational.hpp` | Owned FLINT fields and exact rational functions |
| Differential systems | `system.hpp`, `frobenius.hpp`, `affine_frobenius.hpp` | Ordinary and regular-singular series recurrences |
| Family preparation | `families.hpp`, `ibp.hpp`, `fire.hpp`, `fire_modular.hpp` | Quadratic families, exact IBP witnesses and external C++ reduction |
| Recursive integration | `merge_pullback.hpp`, `recursion_graph.hpp`, `recursion_pipeline.hpp` | Feynman merging, beta/endpoint functionals and recursive evaluation |
| Numerical transport | `adjoint_transport.hpp`, `factored_transport.hpp`, `affine_matching.hpp` | Continuation and matching of retained series |
| Durable work | `artifact_store.hpp`, `level_cache.hpp`, `cached_affine.hpp` | Bounded, verified reuse of exact and numerical artifacts |
| User interfaces | `src/main.cpp`, `Mathematica/DiffExp.wl` | CLI and one-process JSON calls |

Feynman recursion distinguishes beta integrals, endpoint limits and direct
terms when a merged propagator is absent. Numerator polynomials remain explicit
through merge pullbacks and the Gaussian scalar leaf. The scheduler separates
scientific identity from resource demand, so a local refinement can extend the
deficient artifact without repeating its already adequate dependencies.

Endpoint series are saved incrementally by column and checked independently
against the retained differential equation and normalization conditions when
restored. Rational denominators are cleared to derive finite-lag recurrences
where useful. Unsupported spectra, insufficient exact closure and exhausted
budgets fail explicitly.

Observable maps compose through a shared scalar boundary before its uncertainty
is applied. This retains correlations and cancellations between recursion
levels, which was crucial to the four-loop banana calculations. Adjoint and
factored methods offer different costs without changing the intended integral.
Conditioning checks can subdivide a chart or use the original recurrence.

Fixed divergent endpoint sectors produce explicit conditions on the physical
boundary. Those conditions are checked after contraction with the shared
source. This does not constitute a proof of omitted Taylor/epsilon tails.
See [endpoint domains](endpoint-domain-constraints.md),
[conditioning](adjoint-conditioning-fallback.md),
[native IBP](native-ibp.md), and [validation](validation.md).

The C++ core is independent of Wolfram. Data imports from historical `.m` files
parse a restricted data grammar without evaluating Wolfram code. The optional
prepared-input compatibility runtime contains inherited C++ regression
machinery and is disabled by default. Neither the public wrapper nor the
native recursion pipeline needs a LibraryLink interface.


## Public input interfaces

The native `transport` command accepts actual matrix entries, kinematic paths
and ordinary or partial asymptotic boundaries through the
[transport JSON protocol](transport-protocol.md). Its canonical implementation
shares both letter expansions and identical letter/source convolutions across
rows. Noncanonical systems can apply inherited uncertainty through a local
homogeneous map to reduce repeated interval overestimation. These choices leave
series order, working precision and error acceptance controls explicit.

The `prepare` and `ft` commands accept complete
[Feynman-family configurations](feynman-families.md). Propagator geometry,
requested integer powers, merge ladders, matching anchors and causal input are
separate configuration data. A label never chooses a calculation. Prepared
graphs retain the actual family definition so default Euclidean prescriptions
are checked against its Symanzik polynomials.

The Mathematica package is an input adapter: it reads the original matrix files,
normalizes boundary expressions and launches the executable using `RunProcess`.
It does not own a parallel numerical integration engine. See the
[compatibility guide](mathematica.md) for supported original functions and limits.
