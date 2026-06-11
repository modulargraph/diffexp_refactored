#!/usr/bin/env python3
"""Exact eps-expansion of the box needed-integral spec {2,0,1,1}(t):
I(t) = Gamma(2+eps) Int_0^1 dz0 dz1 z0 z1^2 F^(-2-eps),
F = z0 z1 (A + B z1 + C z0 z1),  A=(1-t)/3, B=(4t-1)/3, C=-t.
=> I = Gamma(2+eps) Int dz0 dz1 z0^(-1-eps) z1^(-eps) G^(-2-eps),
G = A + B z1 + C z0 z1.
z0-endpoint subtraction:
Int_0^1 z0^(-1-eps) G^(-2-eps) dz0
  = -(1/eps) G0^(-2-eps) + Int_0^1 z0^(-1-eps) [G^(-2-eps) - G0^(-2-eps)] dz0,
G0 = A + B z1 (z0->0).  Everything else is regular; expand in eps with
log-moments and integrate numerically at high precision."""
import mpmath as mp
from functools import lru_cache

mp.mp.dps = 30

def spec2011(t, orders=2):
    A = (1 - t) / mp.mpf(3)
    B = (4 * t - 1) / mp.mpf(3)
    C = -t

    def G(z0, z1):
        return A + B * z1 + C * z0 * z1

    def G0(z1):
        return A + B * z1

    # I/Gamma(2+eps) = Int dz1 z1^(-eps) [ -(1/eps) G0^(-2-eps)
    #     + Int dz0 z0^(-1-eps) (G^(-2-eps) - G0^(-2-eps)) ]
    # expand all in eps: x^(-eps) = sum (-eps log x)^k/k!,
    # H^(-2-eps) = H^-2 sum (-eps log H)^k/k!
    # term1: -(1/eps) * z1^(-eps) G0^(-2-eps)
    #   coefficient of eps^m needs k-sum with prefactor -(1/eps)
    # build I as Laurent in eps: orders eps^-1 .. eps^orders
    # moments via 1-d (and 2-d) quadrature
    def m1(k):  # Int dz1 [ -log(z1) - log(G0) ]^k-structure: full log of z1^(-eps)G0^(-2-eps)...
        # define L(z1) = -log z1 - log G0 ... but the powers differ (1 vs 2+):
        # z1^(-eps) G0^(-2-eps) = G0^-2 Exp[-eps(log z1 + log G0)]
        f = lambda z1: G0(z1) ** -2 * (-(mp.log(z1) + mp.log(G0(z1)))) ** k
        return mp.quad(f, [0, 1])

    def m2(k):  # double integral of the subtracted regular piece
        def f(z0, z1):
            g, g0 = G(z0, z1), G0(z1)
            base = g ** -2 - g0 ** -2 if k == 0 else None
            # general: z0^(-1) [ G^(-2) Lg^k_G - G0^(-2) Lg0^k ] with the
            # z0^(-eps) z1^(-eps) logs included:
            Lfull = -(mp.log(z0) + mp.log(z1) + mp.log(g))
            L0 = -(mp.log(z0) + mp.log(z1) + mp.log(g0))
            return (g ** -2 * Lfull ** k - g0 ** -2 * L0 ** k) / z0
        return mp.quad(lambda z0: mp.quad(lambda z1: f(z0, z1), [0, 1]),
                       [0, mp.mpf(1) / 2, 1])

    # series assembly: I/Gamma = -(1/eps) sum_k (eps^k/k!) m1(k)  +  sum_k (eps^k/k!) m2(k)
    # Laurent coefficients of I/Gamma:
    coeffs = {}
    for p in range(-1, orders + 1):
        c = mp.mpf(0)
        # from term1: -(1/eps) eps^k/k! m1(k) at eps^p -> k = p+1
        k = p + 1
        if k >= 0:
            c += -m1(k) / mp.factorial(k)
        # from term2: eps^k/k! m2(k) at eps^p -> k = p (p>=0)
        if p >= 0:
            c += m2(p) / mp.factorial(p)
        coeffs[p] = c

    # multiply by Gamma(2+eps) expansion
    eps = None
    gam = mp.taylor(lambda e: mp.gamma(2 + e), 0, orders + 2)
    out = {}
    for p in range(-1, orders + 1):
        v = mp.mpf(0)
        for q in range(-1, p + 1):
            g_idx = p - q
            if g_idx < len(gam):
                v += coeffs.get(q, mp.mpf(0)) * gam[g_idx]
        out[p] = v
    return out

for tt in [mp.mpf(1) / 20, mp.mpf(2) / 5, mp.mpf(19) / 20]:
    res = spec2011(tt)
    print(f"t = {tt}:")
    for p in sorted(res):
        print(f"   eps^{p:>2}: {mp.nstr(res[p], 16)}")
