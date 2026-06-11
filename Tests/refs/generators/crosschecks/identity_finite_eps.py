#!/usr/bin/env python3
"""Test FIRE's reduction identity at finite eps = -1/20, t = 1/20:
LHS = spec_{2011}(t)  (2-d quad, converges at negative eps)
RHS = c1(t,eps) M1(t,eps) + c2(t,eps) M2(t,eps)  (closed forms)
c1 = (18-6d)/(t-4t^2), c2 = 18(d-3)/(1-5t+4t^2), d = 4-2eps
M1 = Gamma(eps) t^(-eps) B(1-eps,1-eps)            [F = t z(1-z)]
M2 = Gamma(eps) ((1-t)/3)^(-eps) B(1-eps,1-eps)    [F = (1-t)/3 z(1-z)]
"""
import mpmath as mp

mp.mp.dps = 25
pass

for eps in [mp.mpf(-1)/20, mp.mpf(-1)/100, mp.mpf(-1)/300]:
  d = 4 - 2 * eps
  for t in [mp.mpf(2) / 5]:
      A = (1 - t) / mp.mpf(3)
      B_ = (4 * t - 1) / mp.mpf(3)
      C = -t

      def f(z0, z1):
          G = A + B_ * z1 + C * z0 * z1
          return z0 ** (-1 - eps) * z1 ** (-eps) * G ** (-2 - eps)

      inner = lambda z0: mp.quad(lambda z1: f(z0, z1), [0, mp.mpf(1)/2, mp.mpf(9)/10, 1])
      J = mp.quad(inner, [0, mp.mpf(1)/2, mp.mpf(9)/10, 1])
      lhs = mp.gamma(2 + eps) * J

      Bfun = mp.beta(1 - eps, 1 - eps)
      M1 = mp.gamma(eps) * t ** (-eps) * Bfun
      M2 = mp.gamma(eps) * ((1 - t) / 3) ** (-eps) * Bfun
      c1 = (18 - 6 * d) / (t - 4 * t ** 2)
      c2 = 18 * (d - 3) / (1 - 5 * t + 4 * t ** 2)
      rhs = c1 * M1 + c2 * M2

      print(f"eps = {mp.nstr(eps,6)}, t = {mp.nstr(t,4)}:")
      print("  LHS spec quad :", mp.nstr(lhs, 14))
      print("  RHS c1M1+c2M2 :", mp.nstr(rhs, 14))
      print("  ratio RHS/LHS :", mp.nstr(rhs / lhs, 12))
