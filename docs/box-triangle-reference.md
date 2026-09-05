# Independent box-triangle numerical reference

The native target is the six-line scalar `{0,1,1,1,1,1,1,0,0}` in the original planar double-box basis, with auxiliary slots `-(l1+p123)^2` and `-(l2+p1)^2`. All denominators use `D=-q^2`, the measure is the product of `d^d l/(i*pi^(d/2))`, `d=4-2eps`, `s=-1`, `t=-1/3`; no exponential Euler-gamma normalization is included.

The original unmerged native IBP reduction was generated with one FIRE worker, 180-second timeout and 6 GiB memory cap. Its equations, requests, tables, receipt and reproducer are in `box-triangle-original-ibp/`; its human-readable rows are in `box-triangle-original-ibp.txt`. This is separate from the recursive FT graph cache. The scalar reduces to two massless sunsets, one gamma-function bubble-triangle, a five-line box with an insertion and a five-line diagonal box.

The diagonal Mellin–Barnes representation is [Smirnov–Veretin, hep-ph/9907385, equation (21)](https://arxiv.org/pdf/hep-ph/9907385). Its native polynomials are checked exactly in the test: `U=(a+b)*(d+e)+c*(a+b+d+e)` and `F=c*(b*d+a*e/3)`. The insertion representation follows by integrating the bubble and applying the same paper's Mellin–Barnes identity to the one-loop box. Both remaining contour integrals have single poles at `z=n` and double poles at `z=-2eps+n`. Their paired residue series converge at `t/s=1/3`.

`box_triangle_oracle.hpp` evaluates these independent analytic representations at complex epsilon points. It removes the maximum fourth-order pole and extracts the first five Taylor coefficients by discrete Cauchy sums. The two settings `(radius,samples,residue terms)=(1/16,48,128)` and `(1/20,64,160)` disagree by at most `1.2e-41`. A conservative *estimated* absolute reference error of `1e-30` is used. No theorem bounding residue truncation or discrete-Cauchy aliasing is claimed. The displayed Arb radii cover arithmetic only.

The resulting coefficients are approximately:

| epsilon order | coefficient |
|---|---:|
| -4 | 2.25 |
| -3 | 2.346284806949595738549299160780554731224 |
| -2 | 4.924446627146843202796905296608428289181 |
| -1 | 11.94391135250396922949858412436307897759 |
| 0 | 9.504808900347392415026700059310948833574 |

The tests check the exact native diagonal geometry, the leading pole and the independent resolution comparison. This reference does not use the native FT numerical output. It must remain labelled a numerical reference estimate even when a subsequent FT comparison passes.
