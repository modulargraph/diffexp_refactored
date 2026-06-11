#!/usr/bin/env python3
"""Independent evaluation of the exported box family specs (levels 1-3,
0- and 1-variable cases) via sympy/mpmath: value(eps) =
prefactor(eps) * Integral(U^UPower * F^FPower * Remainder, z in simplex),
expanded as a Laurent series in eps. Cross-checks both pySecDec and the
FT STEPWISE rows."""
import json
import re
import sys

import mpmath as mp
import sympy as sp

mp.mp.dps = 30
EPS = sp.symbols("eps")
ORDER_MAX = 2  # expand through eps^2


def laurent_of_prefactor(pref_str):
    expr = sp.sympify(pref_str.replace("gamma", "Gamma"), locals={"eps": EPS, "Gamma": sp.gamma})
    s = sp.series(expr, EPS, 0, ORDER_MAX + 3)
    return s.removeO()


def spec_value(spec):
    variables = spec["Variables"]
    if len(variables) > 1:
        return None
    F = sp.sympify(spec["F"], locals={"eps": EPS})
    U = sp.sympify(spec["U"], locals={"eps": EPS})
    fpow = sp.sympify(spec["FPower"], locals={"eps": EPS})
    upow = sp.sympify(spec["UPower"], locals={"eps": EPS})
    rem = sp.sympify(spec["Remainder"], locals={"eps": EPS})
    if U != 1:
        # U^upow with U==1 in all current specs; bail loudly otherwise
        return None
    # fpow = alpha + beta*eps (linear in eps for these families)
    alpha = fpow.subs(EPS, 0)
    beta = sp.expand(fpow - alpha).coeff(EPS)
    assert sp.expand(fpow - alpha - beta * EPS) == 0

    if not variables:
        F0 = sp.nsimplify(F)
        base = sp.N(F0, 40)
        mom = [mp.mpf(str(base)) ** mp.mpf(str(sp.N(alpha, 40))) *
               mp.log(mp.mpf(str(base))) ** k * mp.mpf(str(sp.N(rem, 40)))
               for k in range(ORDER_MAX + 4)]
    else:
        z = sp.symbols(variables[0])
        Ffun = sp.lambdify(z, F, "mpmath")
        remfun = sp.lambdify(z, rem, "mpmath") if rem.free_symbols else (lambda t, c=mp.mpf(str(sp.N(rem, 40))): c)
        alpha_m = mp.mpf(str(sp.N(alpha, 40)))
        mom = []
        for k in range(ORDER_MAX + 4):
            def integrand(t, k=k):
                Fv = Ffun(t)
                return (Fv ** alpha_m) * mp.log(Fv) ** k * remfun(t)
            mom.append(mp.quad(integrand, [0, mp.mpf(1) / 2, 1]))

    # sum_k moments * (beta*eps)^k / k!
    beta_m = mp.mpf(str(sp.N(beta, 40)))
    integral_series = sum(
        sp.Float(str(mom[k] * beta_m ** k / mp.factorial(k)), 30) * EPS ** k
        for k in range(ORDER_MAX + 4)
    )
    total = sp.expand(laurent_of_prefactor(spec["Prefactor"]) * integral_series)
    out = {}
    for p in range(-3, ORDER_MAX + 1):
        c = total.coeff(EPS, p)
        if c != 0:
            out[p] = complex(sp.N(c, 25))
    return out


def main():
    specs = json.loads(open(sys.argv[1]).read())
    for spec in specs:
        if spec["Name"].endswith("_needed"):
            continue
        val = spec_value(spec)
        if val is None:
            print(f"{spec['Name']}: skipped (needs sector decomposition)")
            continue
        print(f"{spec['Name']} level={spec['Level']} master={spec['Master']}:")
        for p in sorted(val):
            v = val[p]
            print(f"   eps^{p:>3}: {v.real:+.12f}" + (f" {v.imag:+.2e}i" if abs(v.imag) > 1e-15 else ""))


if __name__ == "__main__":
    main()
