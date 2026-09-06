# Finite-lag transport as the automatic default

Compatible canonical differential equations now use the rational finite-lag recurrence automatically, including through the Mathematica process wrapper. The solver derives a diagonal root gauge from the supplied connection, keeps the original dimension, and prepares shared rational products once. Distinct coefficient polynomials are translated once per chart and shared over source components. No example-name dispatch is involved.

This changes ordinary differential-equation transport. The separate Feynman-trick and singular-start algorithms retain their existing implementations. The retained recurrence has O(N) arithmetic cost at fixed polynomial degree, epsilon depth and working precision; fallback charts still use Taylor convolutions.

## Full-path one-loop check

The original 13-component planar one-loop example was used instead of a full `zzz` run. The higher-order comparison used epsilon orders 0–4, order 200, 499 working bits and a 40-digit accuracy goal. Both methods completed the same 116-chart path. Timings include exact preparation, chart translation, root continuation, tail estimates and adaptive checks.

| Method | Preparation | Numerical transport | Total |
|---|---:|---:|---:|
| Previous series solver | 0.180 s | 75.978 s | 76.158 s |
| Automatic default | 0.362 s | 65.465 s | 65.828 s |

The default was **1.157× faster** in this sequential full-path comparison. It accepted 41 finite-lag charts and used the original solver for 75 charts. Across all 65 endpoint coefficients, the maximum normalized L1 midpoint difference was **3.973e-98** and the maximum normalized output radius was **1.344e-87**. Normalization divides by max(1, the reference midpoint L1 norm). Both runs passed their final accuracy checks. This is a single paired observation, not a statistical benchmark or a full-family speed guarantee.

The ordinary order-40 configuration selects the previous solver because the maximum rational polynomial degree is 48. Its endpoint coefficients agree exactly with the forced-series run. The final default took 10.40 s, versus 10.09 s for the series baseline (including about 0.18 s extra preparation). This avoids using a high-degree recurrence where its setup and finite-lag work are not justified by the expansion order.

Machine: Apple M4, macOS 26.2, FLINT 3.4.0, C++20 `-O3 -DNDEBUG`. The benchmark processes ran sequentially, without concurrent builds. Each process had a finite timeout; the large reconstruction examples were not started.

The production default was also checked on the first high-precision `zzz` chart
(order 230, 931 bits). It selected finite-lag transport, passed the local tail
check and took **5.59 s versus 10.05 s** for the original chart, excluding exact
preparation. Its normalized midpoint difference was **1.014e-143**. The production CLI also accepted that segment with the original 128-digit goal
after reparameterizing it onto [0,1], using one finite-lag chart and no fallback.
No full high-precision `zzz` run was repeated.

All six focused regressions passed: transport, the new finite-lag test, frontend,
CLI, Mathematica smoke and Mathematica compatibility (25.13 s total).

## Selection and accuracy safeguards

The automatic path requires a consistent exact root gauge, bounded preparation, and an expansion order at least four times the largest polynomial degree. It falls back when a denominator/gauge is noninvertible or the scaled denominator feedback norm `sum |q_j h^j| / |q_0|` exceeds 3/4. It also checks interval amplification. If a candidate fails the adaptive tail or roundoff budget, the original solver is tried at the same step before step refinement.

These safeguards matter: an initial unguarded full-path trial failed local refinement, and an intermediate variant failed its final propagated-error check. Those variants were rejected. Sharing polynomial translations and checking denominator feedback before solving resolved the tested path. A higher-precision baseline at order 160 and 700 bits hit its 150-second experimental cap; it is not used in the timing comparison above.

All incoming boundary uncertainty is retained. Roots continue on the current sheet. The existing last-four-term tail estimate is applied in the transformed basis and divided by the endpoint gauge before physical-basis error propagation. It remains an estimate, not a certificate for the infinite remainder. `omitted_tails_certified` remains false.

Saved Taylor segments use the original solver so their coefficients retain the physical-basis convention. Noncanonical systems and singular starts also retain their existing treatment. Set `"recurrence":"series"` in a native transport request to force the previous solver; the default is `"auto"`. Output reports the accepted chart counts and preparation fallback reason.

## Reproduce

After building the current executable:

```sh
python3 tools/performance/run_planar_default.py --executable build/diffexp --output /tmp/diffexp-default-comparison
```

The compressed exported one-loop input is included. The runner performs sequential order-40 and order-200 comparisons with a 180-second cap per process, verifies endpoint agreement and preserves input/output JSON and logs. No Mathematica process or data download is required.

See [recorded results](validation/finite-lag-default.json), [the initial zzz experiment](algebraic-zzz-experiment.md), and [the transport protocol](transport-protocol.md). The focused regression covers an analytic algebraic solution, translated polynomials, both root sheets, carried uncertainty, an entire complex path, saved physical coefficients, an inconsistent gauge and the denominator-feedback fallback.
