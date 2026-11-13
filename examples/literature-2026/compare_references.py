#!/usr/bin/env python3
"""Compare DiffExp results with the numbers printed in the source papers.

A family takes part when its ``families/<name>.meta.json`` has a
``comparison`` object:

    "comparison": {
      "quantity": "mi[PBF,3] of the paper",          # label only
      "terms": [{"integral": 0, "coefficient": "-604019161/42119136"}],
      "eps_shift": 4,                  # multiply by eps^4 (pure-basis normalisation)
      "exp_eps": "2*gammaE",           # multiply by exp(c*eps); may use gammaE, log(...)
      "reference": {"0": ["-0.0729166...", "0"], "1": ["0.5653...", "-0.4581..."]},
      "conjugate": false,              # compare with the complex conjugate (other i0 side)
      "notes": "sign: (-1)^7 from 1/(q^2-m^2) propagators is in the coefficient"
    }

``terms`` index the family's ``integrals`` rows; coefficients are exact
rationals (all sign conventions folded in).  ``reference`` maps an epsilon power
of the converted quantity to [real, imaginary] decimal strings.  The script
prints, per order, the number of agreeing digits, and stores a one-line summary
as ``check`` in the metadata when called with ``--write``.

Arithmetic uses 60-digit decimals; only the standard library is needed.
"""

import argparse
import json
import re
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path

getcontext().prec = 60
HERE = Path(__file__).resolve().parent
GAMMA_E = Decimal("0.577215664901532860606512090082402431042159335939923598805767")


def parse_ball(text):
    match = re.match(r"\[\s*([^\s+]*)\s*(?:\+/-\s*([^\]]+))?\]", text.strip())
    mid = match.group(1) or "0"
    return Decimal(mid), Decimal(match.group(2) or "0")


def to_decimal(value):
    value = str(value).strip()
    if "/" in value:
        fraction = Fraction(value)
        return Decimal(fraction.numerator) / Decimal(fraction.denominator)
    return Decimal(value)


def eval_constant(expression):
    """Evaluate a small real expression with gammaE, log(...), integers, rationals."""
    expression = expression.replace("gammaE", "G")
    tokens = re.findall(r"\d+/\d+|\d+\.\d*|\d+|[A-Za-z_]+|[()+\-*/]", expression)
    position = 0

    def peek():
        return tokens[position] if position < len(tokens) else None

    def take():
        nonlocal position
        position += 1
        return tokens[position - 1]

    def atom():
        token = take()
        if token == "(":
            value = expr()
            take()
            return value
        if token == "-":
            return -atom()
        if token == "G":
            return GAMMA_E
        if token == "log":
            take()
            value = expr()
            take()
            return value.ln()
        return to_decimal(token)

    def term():
        value = atom()
        while peek() in ("*", "/"):
            value = value * atom() if take() == "*" else value / atom()
        return value

    def expr():
        value = term()
        while peek() in ("+", "-"):
            value = value + term() if take() == "+" else value - term()
        return value

    return expr()


def series_mul(a, b, top):
    """Truncated product of {power: (re, im)} series."""
    out = {}
    for pa, (ra, ia) in a.items():
        for pb, (rb, ib) in b.items():
            power = pa + pb
            if power > top:
                continue
            r, i = out.get(power, (Decimal(0), Decimal(0)))
            out[power] = (r + ra * rb - ia * ib, i + ra * ib + ia * rb)
    return out


def digits_agree(value, reference):
    difference = abs(complex(float(value[0] - reference[0]), float(value[1] - reference[1])))
    scale = abs(complex(float(reference[0]), float(reference[1])))
    if scale == 0:
        return None if difference == 0 else -1
    if difference == 0:
        return 99
    from math import log10
    return max(0.0, -log10(difference / scale))


def compare(name, write):
    meta_path = HERE / "families" / f"{name}.meta.json"
    meta = json.loads(meta_path.read_text())
    specs = meta.get("comparison")
    result_path = HERE / "results" / f"{name}.json"
    if not specs or not result_path.exists():
        return None
    result = json.loads(result_path.read_text())
    response = result.get("response") or {}
    if response.get("status") != "completed":
        return None
    summaries, worst_all = [], None
    for spec in specs if isinstance(specs, list) else [specs]:
        worst, summary = compare_one(name, spec, response)
        summaries.append(summary)
        if worst is not None:
            worst_all = worst if worst_all is None else min(worst_all, worst)
    if write:
        meta["check"] = "vs paper: " + "; ".join(summaries)
        meta_path.write_text(json.dumps(meta, indent=1) + "\n")
    return worst_all


def compare_one(name, spec, response):
    low = response["epsilon_low"]
    rows = response["coefficients"]
    top_reference = max(int(p) for p in spec["reference"])
    shift = int(spec.get("eps_shift", 0))
    top = top_reference - shift  # highest DiffExp order needed

    combined, radius = {}, Decimal(0)
    for term in spec["terms"]:
        coefficient = to_decimal(term["coefficient"])
        for offset, entry in enumerate(rows[term["integral"]]):
            power = low + offset
            if power > top:
                continue
            (re_mid, re_rad), (im_mid, im_rad) = parse_ball(entry["real"]), parse_ball(entry["imaginary"])
            r, i = combined.get(power, (Decimal(0), Decimal(0)))
            combined[power] = (r + coefficient * re_mid, i + coefficient * im_mid)
            radius = max(radius, abs(coefficient) * max(re_rad, im_rad))
    available_top = low + len(rows[0]) - 1
    # exp(c eps) series
    c = eval_constant(spec["exp_eps"]) if spec.get("exp_eps") else Decimal(0)
    exp_series, term_value = {}, Decimal(1)
    for n in range(0, top - low + 1):
        exp_series[n] = (term_value, Decimal(0))
        term_value = term_value * c / (n + 1)
    converted = series_mul(combined, exp_series, top)
    converted = {p + shift: v for p, v in converted.items()}

    lines, worst = [], None
    for power_text, (re_text, im_text) in sorted(spec["reference"].items(), key=lambda kv: int(kv[0])):
        power = int(power_text)
        if power - shift > available_top:
            lines.append(f"  eps^{power}: not computed (DiffExp stops at eps^{available_top + shift})")
            continue
        reference = (to_decimal(re_text), to_decimal(im_text))
        value = converted.get(power, (Decimal(0), Decimal(0)))
        if spec.get("conjugate"):
            value = (value[0], -value[1])
        agree = digits_agree(value, reference)
        shown = "exact zero" if agree is None else f"{agree:5.1f} digits"
        lines.append(f"  eps^{power}: DiffExp {float(value[0]):+.15e} {float(value[1]):+.15e}i"
                     f" | paper {float(reference[0]):+.15e} {float(reference[1]):+.15e}i | {shown}")
        if agree is not None:
            worst = agree if worst is None else min(worst, agree)
    compared = sum(1 for line in lines if "digits" in line or "exact zero" in line)
    label = spec.get("quantity", "reference")
    summary = (f"{label} {worst:.0f} digits ({compared} orders)" if worst is not None
               else f"{label}: no comparable orders")
    print(f"{name}: {summary}")
    print("\n".join(lines))
    return worst, summary


def compare_causal():
    """Compare every completed side choice written by causal_explore.py."""
    for path in sorted((HERE / "causal" / "results").glob("*.json")):
        record = json.loads(path.read_text())
        if record.get("status") != "completed":
            continue
        meta = json.loads((HERE / "families" / f"{record['family']}.meta.json").read_text())
        specs = meta.get("comparison") or []
        sides = ", ".join(f"level {p['level']} x={p['x'][:8]} side {p['side']:+d}" for p in record["points"])
        print(f"== {record['family']} with {sides or 'no explicit sides'}")
        for spec in specs if isinstance(specs, list) else [specs]:
            compare_one(record["family"], spec, record["response"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("families", nargs="*")
    parser.add_argument("--write", action="store_true", help="store the summary as 'check' in the metadata")
    parser.add_argument("--causal", action="store_true", help="compare the side choices in causal/results")
    args = parser.parse_args()
    if args.causal:
        compare_causal()
        return
    names = args.families or sorted(p.name[: -len(".meta.json")] for p in (HERE / "families").glob("*.meta.json"))
    for name in names:
        compare(name, args.write)


if __name__ == "__main__":
    main()
