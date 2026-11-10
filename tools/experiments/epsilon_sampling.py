#!/usr/bin/env python3
"""Laurent coefficients of a Feynman-trick family from numeric-epsilon samples.

Each sample runs `diffexp ft CONFIG --epsilon-order 0 --epsilon-shift Q`, which
evaluates the integral in d0-2Q dimensions; the samples at Q=+-step, +-2 step,
... are fitted by the polynomial through eps^poles * I(eps). This is the
experiment recorded in docs/recursion-cost-experiments.md: it reproduces the
expanded recursion and needs about one sample per 1.4 digits. Samples are
independent (--jobs). Only the standard library is used.

    epsilon_sampling.py examples/feynman/box.json --poles 2 --samples 10
"""
import argparse, concurrent.futures, json, re, subprocess, sys, tempfile
from fractions import Fraction

def midpoint(text):
    match = re.match(r"\[([^ \]]+) \+/-", text.strip())
    if match is None:
        return Fraction(0) if text.strip().startswith("[+/-") else Fraction(text.strip())
    return Fraction(match.group(1))

def sample(arguments, shift):
    command = [arguments.diffexp, "ft", arguments.config, "--cache", arguments.cache, "--no-numerical-cache",
               "--json", "--epsilon-order", "0", "--epsilon-shift", str(shift)] + arguments.extra
    done = subprocess.run(command, capture_output=True, text=True, env=dict(__import__("os").environ, OMP_NUM_THREADS="1"))
    report = json.loads(done.stdout.strip().splitlines()[-1])
    if "coefficients" not in report:
        raise RuntimeError(f"sample {shift}: {report.get('message')}")
    value = report["coefficients"][0][-report["epsilon_low"]]
    return shift, (midpoint(value["real"]), midpoint(value["imaginary"])), report["timings"]["numerical_seconds"]

def fit(points, poles):
    """Coefficients of the polynomial through (x, x^poles*y): exact Newton form."""
    xs = [x for x, _ in points]
    table = [y * x ** poles for x, y in points]
    for level in range(1, len(xs)):
        for i in range(len(xs) - 1, level - 1, -1):
            table[i] = (table[i] - table[i - 1]) / (xs[i] - xs[i - level])
    coefficients = [Fraction(0)] * len(xs)
    for i in range(len(xs) - 1, -1, -1):  # Horner on the Newton basis
        shifted = [Fraction(0)] + coefficients[:-1]
        coefficients = [s - xs[i] * c for s, c in zip(shifted, coefficients)]
        coefficients[0] += table[i]
    return coefficients

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config")
    parser.add_argument("--diffexp", default="build/diffexp")
    parser.add_argument("--cache", default=None, help="exact cache directory (default: a temporary one)")
    parser.add_argument("--step", default="1/400", help="spacing of the symmetric nodes")
    parser.add_argument("--samples", type=int, default=8)
    parser.add_argument("--poles", type=int, default=0, help="order of the leading epsilon pole")
    parser.add_argument("--orders", type=int, default=2, help="coefficients to print beyond the pole terms")
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--save", default=None, help="write the sample values to this JSON file")
    parser.add_argument("extra", nargs="*", help="further diffexp ft options, after --")
    arguments = parser.parse_args()
    holder = None
    if arguments.cache is None:
        holder = tempfile.TemporaryDirectory(); arguments.cache = holder.name
    step = Fraction(arguments.step)
    nodes = [sign * k * step for k in range(1, arguments.samples // 2 + 1) for sign in (1, -1)]
    first = sample(arguments, nodes[0])  # fills the exact cache once
    with concurrent.futures.ThreadPoolExecutor(arguments.jobs) as pool:
        results = [first] + list(pool.map(lambda q: sample(arguments, q), nodes[1:]))
    if arguments.save:
        with open(arguments.save, "w") as out:
            json.dump([[str(q), str(v[0]), str(v[1]), t] for q, v, t in results], out)
    print(f"{len(results)} samples, numerical seconds: {sum(t for *_, t in results):.1f} in total")
    for part, name in ((0, "real"), (1, "imaginary")):
        points = [(q, value[part]) for q, value, _ in results]
        if not any(y for _, y in points):
            continue
        final, previous = fit(points, arguments.poles), fit(points[:-2], arguments.poles)
        for index in range(min(len(previous), arguments.poles + arguments.orders)):
            print(f"{name:9s} eps^{index - arguments.poles:+d}: {float(final[index]):.17g}   (change from the last two samples {abs(float(final[index] - previous[index])):.2e})")

if __name__ == "__main__":
    sys.exit(main())
