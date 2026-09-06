# Linear-order recurrence experiment for high-precision zzz

An exact diagonal change of basis removes the square roots from the actual high-precision `zzz` differential equation without increasing its 86-component dimension. A standalone finite-lag recurrence is **1.80–1.86 times faster at order 230**, and **3.13 times faster at order 460**, on the first ordinary chart. The validated method is now integrated as the automatic default for compatible ordinary charts, with degree, conditioning and accuracy fallbacks. The measurements below describe the original isolated experiment; see [the production follow-up](finite-lag-default.md).

## Timings

Apple M4, macOS 26.2, FLINT 3.4.0; C++20 `-O3 -DNDEBUG`. Sequential processes, 931-bit complex ball arithmetic, epsilon orders 0–4 (430 values), original PH1→PH2 input and boundary. The endpoint is the production first proposed chart endpoint, approximately 0.000223678484193389. No full-path run was performed.

| Taylor order | Current chart | Experimental chart | Speedup |
|---:|---:|---:|---:|
| 115 | 2.79 s | 2.57 s | 1.09× |
| 230 | 10.03 s | 5.40 s | 1.86× |
| 460 | 35.88 s | 11.47 s | 3.13× |

Order 115 is only a scaling probe, not a 128-digit result. Two additional order-230 measurements gave 9.85/5.47 seconds and 10.63/5.74 seconds; the latter ran the candidate first. Doubling order 230→460 increased candidate time 2.12× and current time 3.58×.

These are kernel wall times, excluding physical gauge conversion and output serialization. Original input compilation takes 0.66–0.74 s; compilation plus exact algebraic decomposition takes about 1.01 s, and preparation of the gauge/product recurrences adds 1.89–2.16 s. Preparation is outside the table. The current chart includes its existing error-estimation work; the prototype does not yet implement full adaptive acceptance or a certified truncation error bound. These measurements cannot establish a whole-path speedup.

## Why the recurrence is linear in order

Write the three roots as r_a with r_a²=R_a(x). Exact decomposition of every active logarithmic derivative gives one rational function times a root-product mask (or zero). For every matrix coupling i←j, the constraints g_i XOR g_j=mask are consistent. Rescale each original integral as Z_i=(product of roots in g_i) I_i. All transformed couplings, including the diagonal derivative of the gauge, become rational functions of x. Exact reconstruction verifies the algebraic decomposition.

The prototype shares rational products across target rows. For f=(p/q)Z_j it computes

```
q_0 f_n = sum_{a=0}^{deg p} p_a Z_{j,n-a}
          - sum_{a=1}^{deg q} q_a f_{n-a}
(n+1) Z_{i,n+1} = sum of weighted f_n
```

Negative-index coefficients are zero, and epsilon shifts are included. There are 1,357 product groups, 3,839 weighted entries and maximum polynomial degree 55. A ring buffer retains only the needed previous product coefficients. At fixed system, degree, epsilon depth and working precision, the arithmetic operation count is O(N). This is not a claim of linear bit complexity when precision also increases.

An initial attempt using the existing row-wise common-denominator compiler was rejected by its cost heuristic. Sharing individual rational products avoids that row-denominator growth. No quadratic fallback is used in the successful experiment.

## Accuracy checks and limits

An independent exact rational calculation checks all 120 Taylor coefficients of a small square-root system, including its transformed companion, through order 14 and epsilon order 3. Three boundary realizations are enclosed when incoming uncertainty is carried through the recurrence. The test passes.

On `zzz`, the maximum normalized complex midpoint difference uses the L1 norm divided by max(1, reference L1 norm), across all 430 coefficients:

- Candidate order 230 vs candidate 460: 1.2532e-143.
- Current order 230 vs current 460: 9.2069e-145.
- Candidate order 230 vs current 460: 1.2532e-143.
- Candidate order 460 vs current 460: 2.6848e-148.

This supports convergence beyond the requested 128 digits on this chart. At order 230, 391/430 retained output balls overlap; at 460, all 430 overlap. Truncating in the rescaled basis and then dividing by the endpoint root gives different omitted tails from truncating in the original basis. The retained arithmetic balls do not certify those tails. Incoming boundary uncertainty is preserved; convergence of midpoints is not a proof of an absolute error bound.

The subsequent production integration validates a complete smaller-family path and adds conditioning and adaptive-error fallbacks; see the production follow-up linked above. The gauge exists for this family and path; it is not assumed to exist for every algebraic system.

## Reproduce

The standalone sources and compressed exported scientific input are in [tools/performance/algebraic_zzz](../tools/performance/algebraic_zzz). No Mathematica session is required. From the repository root:

```sh
python3 tools/performance/algebraic_zzz/run.py --output /tmp/zzz-linear-experiment
```

The runner builds sequentially, runs exact checks and decomposition, then orders 115, 230 and 460. Each numerical benchmark has a 180-second timeout. Use `--prefix` for a dependency prefix other than `/opt/homebrew`; `--orders 8` performs a quick smoke check. Raw endpoint arrays and logs are retained in the output directory. Set `CANDIDATE_FIRST=1` to reverse kernel timing order. Recorded measurements and convergence checks are in [the validation report](validation/algebraic-zzz-experiment.json).
