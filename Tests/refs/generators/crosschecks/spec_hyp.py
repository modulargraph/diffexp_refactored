#!/usr/bin/env python3
"""spec_{2011}(t; eps) in closed-inner-form:
spec = -Gamma(2+eps)/eps * Int_0^1 dz1 z1^(-eps) P^(-2-eps)
           * hyp2f1(2+eps, -eps, 1-eps, -Q/P),
P = A + B z1, Q = C z1;  A=(1-t)/3, B=(4t-1)/3, C=-t.
1) Laurent of spec at three t's via an eps-ladder (replaces noisy refs).
2) Int_0^1 dt spec(t, -1/20) for the four-way test leg B.
"""
import mpmath as mp

mp.mp.dps = 25

def spec(t, eps):
    A = (1 - t) / mp.mpf(3)
    B = (4 * t - 1) / mp.mpf(3)
    C = -t

    def f(z1):
        P = A + B * z1
        Q = C * z1
        return z1 ** (-eps) * P ** (-2 - eps) * mp.hyp2f1(2 + eps, -eps, 1 - eps, -Q / P)

    quad = mp.quad(f, [0, mp.mpf(1)/2, mp.mpf(9)/10, mp.mpf(99)/100,
                       mp.mpf(999)/1000, 1])
    return -mp.gamma(2 + eps) / eps * quad

# --- pointwise Laurent via eps-ladder (spec = s_-1/eps + s_0 + s_1 eps + ...) ---
def laurent(t, hmax=3):
    # use eps = +/- h, +/- h/2 ... fit s_-1, s_0, s_1, s_2
    hs = [mp.mpf(-1) / 40, mp.mpf(-1) / 80, mp.mpf(-1) / 160, mp.mpf(-1) / 320]
    vals = [spec(t, h) for h in hs]
    # solve Vandermonde in powers eps^-1..eps^2
    Amat = mp.matrix([[h ** p for p in range(-1, 3)] for h in hs])
    rhs = mp.matrix(vals)
    sol = mp.lu_solve(Amat, rhs)
    return list(sol)

for tt in [mp.mpf(1) / 20, mp.mpf(2) / 5, mp.mpf(19) / 20]:
    s = laurent(tt)
    print(f"t={mp.nstr(tt,4)}: eps^-1={mp.nstr(s[0],12)} eps^0={mp.nstr(s[1],12)} eps^1={mp.nstr(s[2],10)}")
    print(f"   check -3/(t(1-t)) = {mp.nstr(-3/(tt*(1-tt)),12)}")

# --- leg B: Int dt at eps = -1/20 ---
epsB = mp.mpf(-1) / 20
B = mp.quad(lambda t: spec(t, epsB), [0, mp.mpf(1) / 4, mp.mpf(1) / 2, 1])
print("legB Int dt spec(-1/20) :", mp.nstr(B, 14))
pin = 12 * epsB ** -2 - mp.mpf("0.334914246809752") / epsB - mp.mpf("41.28416739757452") - mp.mpf("70.1654426808142") * epsB
print("legC pin Laurent(-1/20) :", mp.nstr(pin, 14))
ftv = (12 * epsB ** -2 - mp.mpf("0.33491424680973618") / epsB - mp.mpf("47.97558570044525")
       - mp.mpf("64.13156740051324") * epsB - mp.mpf("20.16318084574723") * epsB ** 2
       + mp.mpf("63.14432733013388") * epsB ** 3 + mp.mpf("159.7546282483413") * epsB ** 4
       + mp.mpf("148.8565881592727") * epsB ** 5)
print("legD FT Laurent(-1/20)  :", mp.nstr(ftv, 14))
