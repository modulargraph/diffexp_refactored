#!/usr/bin/env python3
"""One-command GO/NO-GO verdict for a DiffExp2 box stepwise run.

Usage: python3 Scripts/box_verdict.py <stepwise.log>

Exit 0 = GO (all rows matched), exit 1 = NO-GO (any mismatch / missing row).

Why not compare_stepwise_log.py: its --refs route covers only eps^-1 (the
pySecDec summary constant term carries no *eps^0 suffix, so TERM_RE skips
it) and the known-garbage box_L0_1x1x1x1 entry forces exit 1 on a perfect
run (Tests/refs/MANIFEST.md trust notes).  This script pins every level:

- L3/L2/L1 boundary rows against the 50-digit old-core oracle values
  (Tests/refs/oracle_logs/l2_box.log STEPWISE), themselves verified exact
  against independent sympy/mpmath quadrature (Tests/PINS.md section 3.2)
  and against the closed forms Gamma[1+eps] B(1-eps,1-eps) z^(-eps).
- L0 against the pySecDec loop_package pin (Tests/PINS.md section 2.1,
  box1l_pin_result.log): {12, -0.334914246809752, -41.28416739757452}.
- The FINAL row's Finite field (the runner's eps^0 coefficient).

Tolerances: relative 1e-9 on oracle rows (validated DiffExp2 examples match
oracles to >= 1e-13; the old box defect was +6.69 absolute — any real
defect is orders of magnitude outside this band); L0 eps^0 absolute 5e-9
(pin Cuhre error 4.7e-10; old core agreed to 2.8e-11 absolute).
Imaginary parts must be below 1e-10 (Euclidean region: results are real).
"""
import json
import re
import sys

# (level, master) -> {order: reference}; sources in the docstring
ORACLE = {
    (3, (1, 0, 0, 0)): {-1: -0.048051509995660579307224614971576250675266995559435,
                        0: -0.16617490584461899870664997372983976074672382345969},
    (2, (1, 0, 0, 1)): {-1: 1.0,
                        0: 2.7182390927653145607641470084444778498842087780169},
    (2, (1, 0, 0, 0)): {-1: -0.11933919618640585189446864469466589956439549601381,
                        0: -0.30414426986534937364544566985822745886852079202134},
    (1, (1, 0, 1, 0)): {-1: 1.0,
                        0: 2.1603832782292462860752534045687487022469901996850},
    (1, (1, 0, 0, 1)): {-1: 1.0,
                        0: 3.1719841899077262113657764887832286293696442261781},
}
# L0: the pySecDec loop_package pin (sign +1 for 4 propagators)
PIN_L0 = {-2: 12.0, -1: -0.334914246809752, 0: -41.28416739757452}
RTOL = 1e-9          # oracle rows + L0 poles
ATOL_L0_EPS0 = 5e-9  # pin precision floor at the finite part
IM_TOL = 1e-10


def num(c):
    if isinstance(c, (int, float)):
        return complex(c)
    return complex(c["Re"], c["Im"])


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    rows, final = [], None
    for line in open(sys.argv[1]):
        if line.startswith("STEPWISE "):
            payload = re.sub(r"(\d)\.e", r"\1.0e", line[len("STEPWISE "):])
            rows.append(json.loads(payload))
        elif line.startswith("FINAL "):
            payload = re.sub(r"(\d)\.e", r"\1.0e", line[len("FINAL "):])
            final = json.loads(payload)

    failures = 0

    def check(label, got, ref, tol_abs):
        nonlocal failures
        diff = abs(got.real - ref)
        ok = diff <= tol_abs and abs(got.imag) <= IM_TOL
        print(f"  {label}: {got.real:+.15g}  ref {ref:+.15g}  "
              f"diff {diff:.2e}{'' if abs(got.imag) <= IM_TOL else f'  IM {got.imag:.2e}'}"
              f"  {'PASS' if ok else 'FAIL'}")
        if not ok:
            failures += 1

    seen = set()
    for row in rows:
        if row.get("Example") != "box":
            continue
        level = int(round(row["Level"]))
        master = tuple(int(round(v)) for v in row["Master"])
        coeffs = {int(round(p)): num(c) for p, c in row["Coefficients"]}
        if (level, master) in ORACLE:
            seen.add((level, master))
            print(f"L{level} {list(master)} (oracle):")
            for p, ref in ORACLE[(level, master)].items():
                if p not in coeffs:
                    print(f"  eps^{p}: MISSING")
                    failures += 1
                    continue
                check(f"eps^{p}", coeffs[p], ref, RTOL * max(1.0, abs(ref)))
        elif level == 0 and master == (1, 1, 1, 1):
            seen.add((level, master))
            print(f"L0 {list(master)} (pySecDec pin):")
            for p, ref in PIN_L0.items():
                if p not in coeffs:
                    print(f"  eps^{p}: MISSING")
                    failures += 1
                    continue
                tol = ATOL_L0_EPS0 if p == 0 else RTOL * max(1.0, abs(ref))
                check(f"eps^{p}", coeffs[p], ref, tol)

    missing = (set(ORACLE) | {(0, (1, 1, 1, 1))}) - seen
    for level, master in sorted(missing):
        print(f"L{level} {list(master)}: ROW MISSING")
        failures += 1

    if final is None:
        print("FINAL row: MISSING")
        failures += 1
    else:
        print("FINAL Finite (eps^0):")
        check("Finite", num(final["Finite"]), PIN_L0[0], ATOL_L0_EPS0)

    print(f"\nbox verdict: {'GO (all pins matched)' if failures == 0 else f'NO-GO ({failures} failures)'}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
