# Spectral transport inside Feynman-trick recursion

DiffExp 2.1 can accelerate ordinary transport inside the generic FT recursion. It identifies structure from the exact connection, including epsilon-dependent rational denominators; no family-name switch or supplied numerical boundary is used. Endpoint Frobenius expansions, recursive child demands, and the common child-source representation are retained.

## Use from Mathematica or C++

The automatic choice is the default for `DiffExpFeynmanTrick` as well as the command line:

```wl
family = DiffExpFamilyTemplate["pentagon_massive"];
result = DiffExpFeynmanTrick[family, {
  "--fire", "/path/to/FIRE7", "--cache", "/path/to/cache"
}];
```

For an explicit comparison, pass `"--ft-transport", "taylor"` or `"--ft-transport", "spectral"` in that same argument list. The latter attempts larger coupled blocks and higher node counts before falling back. It can be slower than `auto`.

The equivalent family configuration is:

```json
"numerical": {
  "transport": "auto",
  "transport_digits": 40,
  "endpoint_order": 32,
  "ordinary_order": 80,
  "working_bits": 384,
  "leaf_digits": 28
}
```

`transport_digits` controls the local spectral convergence target. If omitted, it is `leaf_digits + 12`. This is not a guarantee of that many digits in the final integral: endpoint truncation and conditioning also matter. The CLI override is `--ft-transport-digits N`. `--method` continues to choose the adjoint, factored, or value formulation; `--ft-transport` selects the ordinary solver within the adjoint/factored formulations. The value formulation retains its existing solver.

Existing ordinary Taylor checkpoints are resumed using their original algorithm. Completed spectral arms have separate versioned, checksummed checkpoints that bind the equations, boundary balls, path, precision, and solver settings. A matching corrupt file is an error. `--no-numerical-cache` disables both numerical caches while retaining exact reductions.

## Mathematical design

The transported observable rows satisfy

\[
 g_i' = f_i-\sum_j g_j A_{ji}(x,\epsilon).
\]

The native epsilon gauge first ensures nonnegative epsilon powers in the connection. Where possible, FLINT identifies an exact diagonal integrating factor

\[
 h_i'/h_i=-A_{ii}(x,0),\qquad h_i=\prod_a P_a(x)^{c_a}.
\]

The identification requires an exact reconstructed identity; only integer and half-integer exponents are used by this backend. Setting `g_i=h_i u_i` removes that diagonal, scales an edge by `h_j/h_i`, and scales the forcing by `1/h_i`. Each factor is normalized at the start of a path leg. Square roots are continued along the original contour, and the physical basis is restored at the end. Original poles remain in the geometry checks even if the transformed equation cancels them.

The epsilon-zero dependency graph is decomposed into strongly connected blocks. Acyclic couplings are evaluated in dependency order; positive-epsilon couplings use earlier epsilon layers. Thus the method does not require a globally epsilon-linear form or extra epsilon shifts that increase child demands.

At each Chebyshev resolution, a block operator is factored once and reused for every epsilon coefficient and observable row. Scalar blocks whose diagonal was removed share one integration inverse, including across path legs. Rational coefficients are compiled into polynomial arrays once and sampled numerically, avoiding repeated expression interpretation.

For `m` nodes, block sizes `b`, `K` epsilon coefficients, `r` observable rows, and `E` coupling entries, dense block factorization costs roughly `sum((m b)^3)` per resolution. Reusing these factors costs roughly `K r sum((m b)^2)`. General epsilon convolutions add `O(E r m K²)` work, less when only a few epsilon powers occur. The common scalar integration inverse avoids refactoring the same operator for every component. These are arithmetic-operation estimates; FLINT precision and coefficient complexity also affect runtime. This does not replace the existing finite-lag Taylor recurrence with a claim of universally linear complexity.

## Selection and limitations

Automatic selection currently limits epsilon-zero blocks to scalar components, uses bounded node counts and nearby-pole forecasts, and falls back to the original complete-arm recurrence when necessary. No partially accepted spectral arm is fed into that fallback. This avoids increased interval wrapping observed in the mixed-segment experiment. Larger blocks remain available through explicit spectral selection.

A contour-preserving exponential reparameterization was tested to cluster nodes near endpoints. Its mathematical tests pass, but it did not improve these FT benchmarks, so it is disabled by default and remains a C++ experimental option.

Acceptance compares at least three resolutions, includes arithmetic uncertainty, checks polynomial degree and original singularities, and adds a guarded estimate of the omitted spectral tail. This is an estimate, not a certified infinite-tail bound. General FT results still report `omitted_tails_certified: false`.

The focused tests cover rational and square-root gauges, monodromy, epsilon-zero cycles and triangular feedback, negative Laurent boundaries and forcing, epsilon-dependent denominators, canceled singularities, insufficient precision, polynomial aliasing, and completed-arm cache integrity. The independent integral checks and matched timings are recorded below. All 21
affected/interface regression tests passed. A complete massive-pentagon run
through the Mathematica Association interface also selected all six spectral
arms and passed its independent reference check (6.21 seconds process time in
that separate validation run). [Test and wrapper records](validation/ft-spectral-tests.json)
include the exact timings; this was not a paired Mathematica speed comparison.

## Matched laptop timings

September 6, 2026, Apple M4, macOS ARM64, FLINT 3.4.0, Release build. Each entry is the median of three sequential runs per method, alternating their order. Both methods reuse the same exact IBP reductions and recompute all numerical stages (`--no-numerical-cache`). These are warm-exact-preparation comparisons, not cold FIRE timings. There were no concurrent builds or numerical jobs.

| Configuration | Complete process: Taylor → auto | Ordinary FT transport: Taylor → auto | Result |
|---|---:|---:|---|
| Massive pentagon, epsilon through 2 | **5.78 → 3.85 s** | **3.02 → 1.05 s** | 1.50× overall; 2.88× transport |
| Massive pentagon, higher precision, epsilon through 4 | **32.94 → 21.43 s** | **14.82 → 4.15 s** | 1.54× overall; 3.57× transport |
| Box-triangle | 5.05 → 5.10 s | 3.21 → 3.20 s | No meaningful gain |
| Unequal three-loop banana | 19.39 → 19.61 s | 8.35 → 8.52 s | Existing recurrence retained |
| Sunrise | 0.832 → 0.838 s | 0.224 → 0.227 s | Existing recurrence retained |

Standard settings are endpoint order 32, ordinary order 80, 384 working bits, 28 leaf digits and a 40-digit local spectral target. The higher-precision pentagon uses endpoint order 64, ordinary order 160, 512 working bits, 60 leaf digits and a 72-digit local target. No settings were reduced for the accelerated runs. Preparation itself takes about 0.003–0.024 seconds here because exact reductions are already cached.

Both pentagon configurations use spectral transport on all six arms. Box-triangle accepts one arm and falls back on the others; the saved work is offset by selection overhead. Banana has genuine coupled epsilon-zero blocks, and the conservative selector retains local transport. Sunrise is so small that the speculative path does not provide a useful end-to-end gain.

All **30 full runs pass independent references**, including forbidden-pole checks: pentagon coefficients 0–2 against Feynman-parameter pins at `1e-18`, box-triangle against its Mellin–Barnes/original-IBP reference at `1e-20`, and banana/sunrise finite parts against certified Bessel quadrature at `1e-20`. Higher pentagon epsilon coefficients are compared between solvers, not independently pinned by those references.

All **81 matched coefficient comparisons** pass a normalized `1e-25` bound including the printed arithmetic radii. The largest bound is `5.58e-34` for box-triangle; the high-precision pentagon agrees within `1.83e-75`. This agreement does not certify omitted endpoint tails. See [every configuration, timing, coefficient and reference check](validation/ft-spectral-results.json).

The high-precision pentagon's remaining roughly 17 seconds are outside ordinary transport. Accelerating that transport further cannot remove endpoint construction and the other recursive work. No Henn reconstruction or heavy four-loop banana reconstruction was run for this comparison.

## Experiments retained as evidence

Removing epsilon-zero diagonals reduced a real four-component box-triangle stage from 64/76 factorizations across its arms to 6/8 after inverse sharing. Compiling rational coefficient polynomials once was another useful improvement.

Two attempted extensions did not justify automatic deployment: exponential endpoint clustering failed to improve complete timings, and feeding accepted spectral segments into local recurrence increased interval wrapping and made box-triangle slower. Automatic transport therefore keeps its original contour and falls back on the complete original arm. A coupled-block solver is implemented and tested, but the current banana blocks were not a performance win. The [exact structure audit](../tools/performance/ft_spectral/README.md) records those distinctions so future work can start from the measured bottlenecks.

## Reproducing a matched comparison

Prepare each configuration once using `prepare`, then supply its configuration
and exact cache to the sequential runner. It alternates execution order, uses
fresh numerical calculations for both methods, records every coefficient and
reference check, and rejects unexpected exact-cache misses:

```sh
python3 tools/performance/benchmark_ft_transport.py \
  --executable /path/to/build/diffexp \
  --checker /path/to/build/diffexp_check_ft_laptop \
  --fire /path/to/FIRE7 --output /tmp/ft-comparison \
  --case pentagon pentagon_massive examples/feynman/pentagon_massive.json /path/to/exact-cache
```

`--case` takes a display label, independent-reference identifier, configuration
and cache. Repeat it for more cases. Each solver process has a configurable
limit of at most 600 seconds; reference checks are capped at 120 seconds.
Compilation and reference evaluation are outside the solver timings.
