# Arithmetic conditioning of retained adjoint charts

Compact denominator-cleared recurrences can amplify interval uncertainty through cancellation even when their exact recurrence and floating-point centers are correct. The box-triangle diagnostic identified this before uncertain leaf data entered the calculation; see `ft-box-triangle-enclosure-diagnostic.md`.

The adjoint chart now checks each coefficient's radius magnitude divided by `max(1,abs(midpoint))`. It computes the original rational recurrence on the **same input boundary, center, step, epsilon window and Taylor order** when either:

- the polynomial coefficient is nonfinite;
- normalized radius exceeds the fixed reserve `2^(-working_bits/2)`;
- normalized radius exceeds 256 times the input coefficient's normalized radius, with rounding floor `2^(-working_bits+16)`.

The fixed reserve is independent of chart count, so a sequence of individually modest losses cannot consume all working precision unnoticed. This is a backend selection rule, not a requested-accuracy certificate. If both finite candidates are available, each real and imaginary interval is intersected: both algorithms enclose the same retained polynomial. An empty intersection is an implementation error. A nonfinite candidate is replaced by the finite one; two nonfinite candidates are rejected. Reference compilation is lazy, shared per contour leg, and reused by independent row batches. Geometry and finite budgets are unchanged.

`AdjointOptions::conditioning_stats` optionally points to counters for polynomial charts, rational cross-checks, and lazy rational compilations. No diagnostic output is emitted by core transport.

At 384 bits and Taylor order 80, the analytic adversarial system `y'=y/(1+x)^20`, `h=1/4`, with input radius `2^-300`, gives polynomial radius about `1.94e-26`; the guarded result has radius about `5.34e-86`. Its retained polynomial is independently evaluated with exact rational coefficients in `test_conditioned_adjoint.cpp`. Tests cover the fixed reserve, incompatible intersections, nonfinite candidates, lazy reference compilation, and row batching.

A cache-only seven-master banana chart with three observable rows, order 80 and 384 bits retained the polynomial fast path: one polynomial chart, zero cross-checks, 22 polynomial rows, zero degree/cost fallbacks, and maximum denominator degree 10. Five ordinary charts took 3.09 seconds with the original recurrence and 0.51 seconds with the polynomial recurrence in this local comparison. The inputs and chart points are identical, and all coefficient differences remain below `1e-35`.

No precision or Taylor-order increase is performed. Omitted FT Taylor tails remain uncertified; these intersections improve arithmetic enclosures of the retained polynomial only.

The full box-triangle replay passed at the unchanged N32/N80, 384-bit settings in 71.44 seconds. All epsilon coefficients from -4 through 0 satisfy the 1e-20 independent numerical comparison; forbidden coefficients -7 through -5 are exact zero. The finite-part difference ball has radius 8.76e-34, replacing the previous arithmetic radius near 1e177. The reference has an estimated, nonrigorous 1e-30 absolute error. See [validation](validation.md) for the release acceptance scope.
