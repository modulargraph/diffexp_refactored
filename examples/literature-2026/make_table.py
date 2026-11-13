#!/usr/bin/env python3
"""Summarize results/<name>.json and families/<name>.meta.json as a Markdown table.

    python3 make_table.py            # prints the table
    python3 make_table.py --write    # replaces the table between markers in README.md
"""

import argparse
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
START, END = "<!-- results-table:start -->", "<!-- results-table:end -->"


def ball_midpoint(ball):
    """'[0.4603523 +/- 6.66e-39]' -> ('0.4603523', '6.66e-39'); '[+/- 2e-30]' -> ('0', '2e-30')."""
    match = re.match(r"\[\s*([^\s+]*)\s*(?:\+/-\s*([^\]]+))?\]", ball.strip())
    if not match:
        return ball, None
    return (match.group(1) or "0"), match.group(2)


def short(number, digits=12):
    mantissa = re.match(r"(-?)(\d*)\.?(\d*)(e[-+]?\d+)?$", number)
    if not mantissa:
        return number
    sign, whole, frac, exp = mantissa.groups()
    if whole not in ("", "0"):
        text = whole + ("." + frac[: max(0, digits - len(whole))] if frac else "")
    else:
        text = "0." + frac[:digits] if frac else "0"
    return sign + text + (exp or "")


def finite_parts(response):
    """Epsilon^0 coefficient (real part, and imaginary part if nonzero) of each integral."""
    values = []
    low = response.get("epsilon_low", 0)
    for row in response.get("coefficients", []):
        index = -low
        if not 0 <= index < len(row):
            values.append("n/a")
            continue
        real, real_radius = ball_midpoint(row[index]["real"])
        imag, imag_radius = ball_midpoint(row[index]["imaginary"])
        text = short(real)
        if imag not in ("0", "") and imag_radius is not None:
            try:
                if abs(float(imag)) > 10 * float(imag_radius):
                    text += f" {'-' if imag.startswith('-') else '+'} {short(imag.lstrip('-'))} i"
            except ValueError:
                pass
        values.append(text)
    return values


def lowest_pole(response):
    """Lowest epsilon power whose coefficient is not consistent with zero."""
    low = response.get("epsilon_low", 0)
    rows = response.get("coefficients", [])

    def magnitude(entry):
        mid, radius = ball_midpoint(entry["real"])
        try:
            return abs(float(mid)), float(radius) if radius else 0.0
        except ValueError:
            return float("inf"), 0.0

    # Cancelled pole coefficients come back as tiny nonzero balls (~1e-36 for a
    # 38-digit run); call them zero relative to the integral's own scale.
    scale = max((magnitude(e)[0] for row in rows for e in row), default=0.0)
    for offset in range(len(rows[0]) if rows else 0):
        for row in rows:
            value, radius = magnitude(row[offset])
            if value > max(10 * radius, 1e-25 * scale):
                return low + offset
    return None


def paper_links(text):
    """Link every arXiv number in the metadata's 'arxiv' field."""
    numbers = re.findall(r"\d{4}\.\d{4,5}", text or "")
    return ", ".join(f"[{n}](https://arxiv.org/abs/{n})" for n in dict.fromkeys(numbers)) or (text or "")


def brief(text, limit=110):
    text = re.sub(r"\s*\[Euclidean fallback point\]", "", text or "")
    return text if len(text) <= limit else text[: limit - 1].rsplit(" ", 1)[0] + " …"


def seconds(value):
    return "" if value is None else (f"{value:.1f}" if value < 100 else f"{value:.0f}")


def rows():
    out = []
    for meta_path in sorted((HERE / "families").glob("*.meta.json")):
        name = meta_path.name[: -len(".meta.json")]
        meta = json.loads(meta_path.read_text())
        result_path = HERE / "results" / f"{name}.json"
        result = json.loads(result_path.read_text()) if result_path.exists() else {}
        response = result.get("response") or {}
        timings = response.get("timings", {})
        status = result.get("status", "not run")
        if status == "completed":
            status = "✅ works"
        elif status == "timeout":
            status = f"⏱ timeout ({result.get('wall_seconds', 0):.0f} s)"
        elif status in ("error", "crashed"):
            message = (result.get("error") or "crashed").split(":")[0][:60]
            status = f"❌ {message}"
        pole = lowest_pole(response) if response.get("coefficients") else None
        out.append({
            "Family": f"`{name}`",
            "Paper": paper_links(meta.get("arxiv")),
            "What": brief(meta.get("description", "")),
            "Point": meta.get("point_kind", ""),
            "Status": status,
            "Prep [s]": seconds(timings.get("preparation_seconds")),
            "Numerics [s]": seconds(timings.get("numerical_seconds")),
            "Total [s]": seconds(result.get("wall_seconds")) if result else "",
            "Leading ε power": "" if pole is None else str(pole),
            "ε⁰ coefficients": "; ".join(finite_parts(response)) if response.get("coefficients") else "",
            "Check": meta.get("check", ""),
        })
    return out


def markdown(table):
    if not table:
        return "_no families_\n"
    keys = list(table[0])
    lines = ["| " + " | ".join(keys) + " |", "|" + "|".join("---" for _ in keys) + "|"]
    for row in table:
        lines.append("| " + " | ".join(str(row[k]).replace("|", "\\|") for k in keys) + " |")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    text = markdown(rows())
    if not args.write:
        print(text)
        return
    readme = HERE / "README.md"
    content = readme.read_text()
    begin, end = content.index(START) + len(START), content.index(END)
    readme.write_text(content[:begin] + "\n" + text + content[end:])


if __name__ == "__main__":
    main()
