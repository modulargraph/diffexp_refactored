#!/usr/bin/env python3
"""Four-way consistency at real eps = -1/20 (everything converges naively):
A) I_box: Gamma(2+eps) Int over 3-simplex of F_box^(-2-eps)  (3-d quad)
B) Int_0^1 dt spec2011(t):  Gamma(2+eps) Int z0 z1^2 F(t,z)^(-2-eps) (3-d quad)
C) pin Laurent series evaluated at eps  (pySecDec/analytic box)
D) FT dump-total Laurent evaluated at eps
A must equal B exactly (Feynman-trick identity).  C approximates A
(truncated series).  D tells whose side the FT machinery is on."""
import mpmath as mp

mp.mp.dps = 20
eps = mp.mpf(-1) / 20

# --- A: box via the exporter's simplex map (vars z0,z1,z2) ---
# F_box = -1/3 z1 z2 (-1 + z0 + z2 + z0(-4+3 z1) z2), U = 1,
# remainder from unitSimplexMap for n=4: rescaled vars + jacobian.
# Use the exporter spec directly: box_L0_1x1x1x1:
#   F(z0,z1,z2) as above, prefactor Gamma(2+eps), FPower = -2-eps,
#   Remainder = jacobian * prod(rescaled^(powers-1)) -- from the spec JSON.
import json
spec = None
for s in json.loads(open("/tmp/pysecdec_box_specs.json").read()):
    if s["Name"] == "box_L0_1x1x1x1":
        spec = s
        break
import sympy as sp
z0, z1, z2 = sp.symbols("z0 z1 z2")
Fb = sp.sympify(spec["F"])
remb = sp.sympify(spec["Remainder"])
Fb_f = sp.lambdify((z0, z1, z2), Fb, "mpmath")
remb_f = sp.lambdify((z0, z1, z2), remb, "mpmath")

def integrandA(a, b, c):
    Fv = Fb_f(a, b, c)
    return remb_f(a, b, c) * Fv ** (-2 - eps)

A = mp.gamma(2 + eps) * mp.quad(
    lambda a: mp.quad(lambda b: mp.quad(lambda c: integrandA(a, b, c), [0, 1]), [0, 1]),
    [0, 1])
print("A (box 3d quad)      :", mp.nstr(A, 12))

# --- B: Int dt of the needed spec ---
def integrandB(t, a, b):
    A_ = (1 - t) / mp.mpf(3)
    B_ = (4 * t - 1) / mp.mpf(3)
    C_ = -t
    F = a * b * (A_ + B_ * b + C_ * a * b)
    return a * b * b * F ** (-2 - eps)

B = mp.gamma(2 + eps) * mp.quad(
    lambda t: mp.quad(lambda a: mp.quad(lambda b: integrandB(t, a, b), [0, 1]), [0, 1]),
    [0, 1])
print("B (int dt spec 3d)   :", mp.nstr(B, 12))

# --- C: pin Laurent at eps ---
pin = [(-2, mp.mpf("12")), (-1, mp.mpf("-0.334914246809752")),
       (0, mp.mpf("-41.28416739757452")), (1, mp.mpf("-70.1654426808142"))]
C = sum(c * eps ** p for p, c in pin)
print("C (pin Laurent)      :", mp.nstr(C, 12), " (truncated at eps^1)")

# --- D: FT dump-total Laurent at eps ---
ft = [(-2, mp.mpf("12")), (-1, mp.mpf("-0.33491424680973618")),
      (0, mp.mpf("-47.97558570044525")), (1, mp.mpf("-64.13156740051324")),
      (2, mp.mpf("-20.16318084574723")), (3, mp.mpf("63.14432733013388")),
      (4, mp.mpf("159.7546282483413")), (5, mp.mpf("148.8565881592727"))]
D = sum(c * eps ** p for p, c in ft)
print("D (FT Laurent)       :", mp.nstr(D, 12), " (through eps^5)")
EOF