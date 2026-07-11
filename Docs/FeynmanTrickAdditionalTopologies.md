# Additional Euclidean FeynmanTrick topologies

## DiffExp2 box-bubble validation

The existing massless two-loop `box_bubble` example was run from scratch at
the Euclidean point `s = -1`, `t = -1/3` with WP300, expansion order 40,
division order 3, and boundary extra order 16.  FIRE preparation and all four
DiffExp2 ladder levels completed in about 55 seconds.  The level-0 result was

```text
 eps^-3   0.4999999999999719409794842591023
 eps^-2   0.4227843350986460144608829749893
 eps^-1   0.3562795605799888368173675114533
 eps^0   -4.877139662454516944362515776161
```

Every coefficient agrees with `Tests/refs/oracle_dumps/m0_oracle_boxbubble.log`
to about `1.4e-12` absolute or better.  The first run wrote a reusable
prepared FIRE snapshot and synchronous per-arm/per-level ladder checkpoints;
no transport was accepted without an honest complete epsilon window.

## Kite convention

`kite` denotes the fully massive equal-mass two-loop five-propagator kite in
`D = 2 - 2 eps` at `p^2 = -1`:

```text
1-l1^2, 1-(l1-p)^2, 1-l2^2, 1-(l2-p)^2, 1-(l1-l2)^2.
```

This fixes the otherwise ambiguous name “kite”; it is not the distinct
three-mass/two-massless four-dimensional convention.

The completed DiffExp2 run used WP300, expansion order 40, division order 3,
boundary extra order 16, and measured per-level epsilon halos `{0,4,7,7}`.
The first two attempts stopped honestly when the level 3→2 and level 2→1
handoffs were short by three and one coefficients respectively.  The final
level-0 result is

```text
 eps^-1  -9.8559958318033201e-21
 eps^0    0.223983919107444028404077222049
```

An independent direct `NIntegrate` of the exported four-dimensional
Feynman-parameter representation gave `0.223983917019259...` with a
conservative `1.8e-7` error estimate, agreeing with the FT result by
`2.1e-9`.

At the 16-master level-1 endpoint this run exposed a factorial matching
implementation: determinant valuation used a Leibniz sum over all `16!`
permutations.  It is now computed in polynomial time from traces and Newton
identities over the honest truncated epsilon-polynomial ring.  The transport
suite pins exact parity with the old formula for a random dimension-4 case
and a dimension-16 valuation case; the real kite level-1 transport then
completed without permutation or memory warnings.

## Four-loop banana convention

`banana4` is the equal-mass four-loop/five-line banana in `D = 2 - 2 eps`
at `p^2 = -1`.  Its physical propagators are the four massive loop lines and
the massive closing line carrying `-l1-l2-l3-l4+p`.  FIRE completes the five
physical denominators with irreducible-numerator slots for the remaining
loop scalar products; those auxiliaries are restricted from positive powers.

At `eps = 0`, Wick rotation of the package normalization
`d^2 l/(i pi)` gives `d^2 k/pi` for each loop, without an additional sign.
Writing `Q^2 = -p^2 = 1`, the physical integral is therefore

```text
I4 = Integral Product[a=1..4, d^2 k_a/pi]
     / (Product[a=1..4, (k_a^2+1)] ((Q-Sum[k_a])^2+1)).
```

The two-dimensional Fourier transform

```text
Integral d^2 k/(2 pi)^2 Exp[i k.r]/(k^2+1) = K0(r)/(2 pi)
```

reduces the four convolutions to the one-dimensional Bessel moment

```text
I4 = 16 Integral_0^Infinity r J0(r) K0(r)^5 dr
   = 39.6555268342976525299928230469335811564460602187101868181...
```

The factor is `2^4`: `d^2 k/pi = 4 pi d^2 k/(2 pi)^2`, so each
loop convolution contributes `2 K0`; the final `K0/(2 pi)` and angular
integral supply the remaining factors.  Equivalently, the normalized
Feynman-parameter representation is

```text
I4 = Integral_Delta4 d^5 x delta(1-Sum[x_i])/(U+P),
U = Sum_i Product_(j!=i) x_j,   P = Product_i x_i.
```

Here `Gamma(5-4 D/2)=Gamma(1)`, the power of `U` is zero, and
`p^2=-1` gives `F=U+P`.  As a normalization cross-check, the same
coordinate-space formula for the established four-line banana gives
`8.268104535868968731543015345479988868728618483845...`, matching the
stored equal-mass checkpoint.

The standalone oracle uses direct and logarithmic-coordinate quadratures,
plus analytic bounds on both truncated tails.  Its default target is 50
significant digits:

```sh
wolframscript -file Scripts/banana4_bessel_oracle.m
BANANA4_BESSEL_DIGITS=70 wolframscript -file Scripts/banana4_bessel_oracle.m
```

This value is an independent `D=2`, `eps^0` oracle.  It is not a claim that
the full `banana4` DiffExp2 ladder has completed.
