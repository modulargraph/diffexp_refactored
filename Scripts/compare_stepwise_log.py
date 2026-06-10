#!/usr/bin/env python3
"""Compare STEPWISE boundary rows from a stepwise-check log against
pySecDec references. Usage: compare_stepwise_log.py <logfile> [level]"""
import json
import re
import sys

REFS = {
    # banana level 1 at z0 = 11/23 (pySecDec, family exporter)
    (1, (2, 0, 1, 1)): {0: 5.4025802965},
    (1, (1, 0, 1, 0)): {-2: -2.0037878788, -1: 2.2439088015, 0: -13.5042252469},
    (1, (1, 0, 0, 1)): {-2: -2.0037878788, -1: 2.2439088015, 0: -13.5042252469},
    (1, (1, 0, 1, 1)): {-3: 1/3, -2: -0.1144868285, -1: 1.3904392399, 0: -5.8738356912},
    (1, (1, 0, 1, 2)): {-2: 0.5, -1: -0.1717302427, 0: -0.9025092631},
    (1, (1, 0, 2, 1)): {-2: 0.5, -1: -0.1717302427, 0: -0.9025092631},
    # no direct pySecDec reference for the numerator master {1,-1,1,1}
    # banana level 0 (direct Feynman-parameter integration)
    (0, (1, 1, 1, 1)): {0: 8.26810451329511583109184},
}

TOL = 5e-6


def main():
    path = sys.argv[1]
    only_level = int(sys.argv[2]) if len(sys.argv) > 2 else None
    rows = []
    for line in open(path):
        if line.startswith("STEPWISE "):
            payload = line[len("STEPWISE "):]
            # Mathematica RawJSON can emit zero-precision zeros like 0.e-63,
            # which is invalid JSON (a digit is required after the dot).
            payload = re.sub(r"(\d)\.e", r"\1.0e", payload)
            rows.append(json.loads(payload))
    failures = 0
    for row in rows:
        level = int(round(row["Level"]))
        if only_level is not None and level != only_level:
            continue
        master = tuple(int(round(v)) for v in row["Master"])[:4]
        coeffs = {}
        for p, c in row["Coefficients"]:
            val = c if isinstance(c, (int, float)) else complex(c["Re"], c["Im"])
            coeffs[int(round(p))] = val
        ref = REFS.get((level, master))
        print(f"L{level} {list(master)}:")
        for p in sorted(coeffs):
            v = coeffs[p]
            vr = v.real if isinstance(v, complex) else v
            vi = v.imag if isinstance(v, complex) else 0.0
            line = f"   eps^{p:>3}: {vr:+.10f}" + (f" {vi:+.2e}i" if abs(vi) > 1e-12 else "")
            if ref is not None and p in ref:
                ok = abs(vr - ref[p]) < TOL * max(1.0, abs(ref[p])) and abs(vi) < 1e-8
                line += f"   ref {ref[p]:+.10f}  {'PASS' if ok else 'FAIL'}"
                if not ok:
                    failures += 1
            print(line)
        if ref is not None:
            extra = [p for p in ref if p not in coeffs]
            if extra:
                failures += 1
                print(f"   MISSING reference orders: {extra}")
    print(f"\nfailures: {failures}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
