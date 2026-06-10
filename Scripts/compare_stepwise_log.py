#!/usr/bin/env python3
"""Compare STEPWISE boundary rows from a stepwise-check log against
pySecDec references.

Usage: compare_stepwise_log.py <logfile> [level] [--refs results.json ...]

References come from two sources:
  - hardcoded banana values below (family exporter, z0 = 11/23);
  - any number of --refs files: results.json as written by
    Scripts/pysecdec_family_driver.py (keys (level, master) are read from
    each entry's spec, values parsed from the pySecDec summary line).
"""
import argparse
import json
import re
import sys

REFS = {
    # banana level 1 at z0 = 11/23 (pySecDec, family exporter)
    ("banana", 1, (2, 0, 1, 1)): {0: 5.4025802965},
    ("banana", 1, (1, 0, 1, 0)): {-2: -2.0037878788, -1: 2.2439088015, 0: -13.5042252469},
    ("banana", 1, (1, 0, 0, 1)): {-2: -2.0037878788, -1: 2.2439088015, 0: -13.5042252469},
    ("banana", 1, (1, 0, 1, 1)): {-3: 1/3, -2: -0.1144868285, -1: 1.3904392399, 0: -5.8738356912},
    ("banana", 1, (1, 0, 1, 2)): {-2: 0.5, -1: -0.1717302427, 0: -0.9025092631},
    ("banana", 1, (1, 0, 2, 1)): {-2: 0.5, -1: -0.1717302427, 0: -0.9025092631},
    # no direct pySecDec reference for the numerator master {1,-1,1,1}
    # banana level 0 (direct Feynman-parameter integration)
    ("banana", 0, (1, 1, 1, 1)): {0: 8.26810451329511583109184},
}

TOL = 5e-6

# one Laurent term of a pySecDec summary line, real or complex payload:
#   + (1.234e+00 +/- 2e-09)*eps^-2     or
#   + ((1.2e+00,3.4e-05) +/- (2e-09,1e-09))*eps^0
TERM_RE = re.compile(
    r"([+-])\s*\(\s*(\(?[^()]*?\)?)\s*\+/-\s*(\(?[^()]*?\)?)\s*\)\s*\*\s*eps\^\(?(-?\d+)\)?"
)


def parse_value(tok: str) -> complex:
    tok = tok.strip()
    if tok.startswith("(") and tok.endswith(")"):
        re_s, im_s = tok[1:-1].split(",")
        return complex(float(re_s), float(im_s))
    return complex(float(tok), 0.0)


def parse_summary(summary: str) -> dict:
    out = {}
    for sign, val, _err, power in TERM_RE.findall(summary):
        v = parse_value(val)
        if sign == "-":
            v = -v
        out[int(power)] = v
    return out


def load_refs_file(path: str) -> dict:
    refs = {}
    for entry in json.loads(open(path).read()):
        spec = entry.get("spec", {})
        if "Level" not in spec or "Master" not in spec:
            continue
        # skip the duplicated "_needed" exports of the same master
        if entry.get("source_name", "").endswith("_needed"):
            continue
        coeffs = parse_summary(entry.get("summary", ""))
        if not coeffs:
            continue
        key = (spec.get("Example", ""), int(spec["Level"]),
               tuple(int(v) for v in spec["Master"]))
        refs[key] = {p: c.real for p, c in coeffs.items() if abs(c.imag) < 1e-6 * max(1.0, abs(c.real))}
    return refs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("logfile")
    ap.add_argument("level", nargs="?", type=int, default=None)
    ap.add_argument("--refs", action="append", default=[],
                    help="pysecdec_family_driver results.json (repeatable)")
    args = ap.parse_args()

    refs = dict(REFS)
    for rf in args.refs:
        refs.update(load_refs_file(rf))

    rows = []
    for line in open(args.logfile):
        if line.startswith("STEPWISE "):
            payload = line[len("STEPWISE "):]
            # Mathematica RawJSON can emit zero-precision zeros like 0.e-63,
            # which is invalid JSON (a digit is required after the dot).
            payload = re.sub(r"(\d)\.e", r"\1.0e", payload)
            rows.append(json.loads(payload))
    failures = 0
    for row in rows:
        level = int(round(row["Level"]))
        if args.level is not None and level != args.level:
            continue
        master = tuple(int(round(v)) for v in row["Master"])
        example = row.get("Example", "")
        coeffs = {}
        for p, c in row["Coefficients"]:
            val = c if isinstance(c, (int, float)) else complex(c["Re"], c["Im"])
            coeffs[int(round(p))] = val
        ref = refs.get((example, level, master),
                       refs.get((example, level, master[:4])))
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
