#!/usr/bin/env python3
"""Independent finite D=2 Bessel oracle for an unequal-mass banana.

The script intentionally imports no project code.  It evaluates the same
one-dimensional moment in r and in the independently parameterized r=t^2
chart, reports their difference and a rigorous massive upper-tail bound.
"""

from __future__ import annotations

import argparse
import json
import sys
import time

import mpmath as mp


DEFAULT_MASS_SQUARED = ("2", "3/2", "4/3", "5/4", "1")


def parse_mpf(value: str) -> mp.mpf:
    """Parse a decimal or a simple rational without using eval."""
    text = value.strip()
    if text.count("/") == 1:
        numerator, denominator = text.split("/", 1)
        result = mp.mpf(numerator) / mp.mpf(denominator)
    elif "/" not in text:
        result = mp.mpf(text)
    else:
        raise ValueError(f"invalid number: {value!r}")
    if not mp.isfinite(result):
        raise ValueError(f"number must be finite: {value!r}")
    return result


def upper_tail_bound(cutoff: mp.mpf, masses: tuple[mp.mpf, ...], loops: int) -> mp.mpf:
    """Bound the omitted eps=0 tail using K0(x)<sqrt(pi/(2x)) exp(-x)."""
    line_count = len(masses)
    mass_sum = mp.fsum(masses)
    mass_product = mp.fprod(masses)
    power = mp.mpf(2) - mp.mpf(line_count) / 2
    constant = (
        mp.power(2, loops)
        * mp.power(mp.pi / 2, mp.mpf(line_count) / 2)
        / mp.sqrt(mass_product)
    )
    return constant * mp.power(mass_sum, -power) * mp.gammainc(
        power, mass_sum * cutoff, mp.inf
    )


def choose_cutoff(
    masses: tuple[mp.mpf, ...], loops: int, digits: int, radial_scale: mp.mpf
) -> tuple[mp.mpf, mp.mpf]:
    mass_sum = mp.fsum(masses)
    cutoff = max(
        8 * radial_scale,
        mp.mpf(digits + 8) * mp.log(10) / mass_sum,
    )
    target = mp.power(10, -(digits + 5))
    bound = upper_tail_bound(cutoff, masses, loops)
    while bound > target:
        cutoff *= mp.mpf("1.25")
        bound = upper_tail_bound(cutoff, masses, loops)
    return cutoff, bound


def radial_breaks(cutoff: mp.mpf, scale: mp.mpf) -> list[mp.mpf]:
    multipliers = (
        "0",
        "1e-6",
        "1e-4",
        "1e-3",
        "1e-2",
        "0.1",
        "0.5",
        "1",
        "2",
        "4",
        "8",
    )
    points = [mp.mpf(value) * scale for value in multipliers]
    next_point = 16 * scale
    while next_point < cutoff:
        points.append(next_point)
        next_point *= 2
    points.append(cutoff)
    return sorted(set(point for point in points if 0 <= point <= cutoff))


def integrate_piecewise(function, points: list[mp.mpf]) -> mp.mpf:
    return mp.fsum(
        mp.quad(function, [left, right])
        for left, right in zip(points, points[1:])
    )


def number_string(value: mp.mpf | mp.mpc, digits: int) -> str:
    return mp.nstr(value, n=digits, strip_zeros=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--digits", type=int, default=30, help="requested agreement digits (default: 30)"
    )
    parser.add_argument(
        "--guard-digits",
        type=int,
        default=15,
        help="extra mpmath working digits (default: 15)",
    )
    parser.add_argument(
        "--q", default="1", help="positive Euclidean momentum magnitude sqrt(-p^2)"
    )
    parser.add_argument(
        "--cutoff", default=None, help="optional finite upper radial cutoff"
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--mass-squared",
        nargs="+",
        metavar="M2",
        help="positive squared masses (decimals or simple rationals)",
    )
    group.add_argument(
        "--masses",
        nargs="+",
        metavar="M",
        help="positive physical masses (decimals or simple rationals)",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.digits < 10:
        raise ValueError("--digits must be at least 10")
    if args.guard_digits < 5:
        raise ValueError("--guard-digits must be at least 5")

    mp.mp.dps = args.digits + args.guard_digits
    q = parse_mpf(args.q)
    if q < 0:
        raise ValueError("--q must be nonnegative")

    if args.masses is not None:
        masses = tuple(parse_mpf(value) for value in args.masses)
        mass_squared = tuple(mass * mass for mass in masses)
    else:
        raw_mass_squared = args.mass_squared or DEFAULT_MASS_SQUARED
        mass_squared = tuple(parse_mpf(value) for value in raw_mass_squared)
        if any(value <= 0 for value in mass_squared):
            raise ValueError("all squared masses must be positive")
        masses = tuple(mp.sqrt(value) for value in mass_squared)

    if len(masses) < 2 or any(value <= 0 for value in masses):
        raise ValueError("at least two positive masses are required")

    loops = len(masses) - 1
    radial_scale = 1 / max((*masses, q))
    if args.cutoff is None:
        cutoff, tail_bound = choose_cutoff(
            masses, loops, args.digits, radial_scale
        )
    else:
        cutoff = parse_mpf(args.cutoff)
        if cutoff <= 0:
            raise ValueError("--cutoff must be positive")
        tail_bound = upper_tail_bound(cutoff, masses, loops)

    prefactor = mp.power(2, loops)

    def density_r(radius: mp.mpf) -> mp.mpf:
        if radius == 0:
            return mp.mpf(0)
        return (
            prefactor
            * radius
            * mp.besselj(0, q * radius)
            * mp.fprod(mp.besselk(0, mass * radius) for mass in masses)
        )

    def density_t(parameter: mp.mpf) -> mp.mpf:
        if parameter == 0:
            return mp.mpf(0)
        return 2 * parameter * density_r(parameter * parameter)

    r_points = radial_breaks(cutoff, radial_scale)
    t_points = [mp.sqrt(point) for point in r_points]

    started = time.perf_counter()
    direct_value = integrate_piecewise(density_r, r_points)
    direct_seconds = time.perf_counter() - started

    started = time.perf_counter()
    squared_value = integrate_piecewise(density_t, t_points)
    squared_seconds = time.perf_counter() - started

    difference = abs(direct_value - squared_value)
    tolerance = mp.power(10, -args.digits)
    accepted = difference + tail_bound < tolerance
    shown_digits = args.digits + 5
    record = {
        "loops": loops,
        "line_count": len(masses),
        "dimension": 2,
        "q": number_string(q, shown_digits),
        "mass_squared": [number_string(value, shown_digits) for value in mass_squared],
        "masses": [number_string(value, shown_digits) for value in masses],
        "requested_digits": args.digits,
        "working_precision": mp.mp.dps,
        "cutoff": number_string(cutoff, shown_digits),
        "tail_bound": number_string(tail_bound, shown_digits),
        "direct_r": number_string(direct_value, shown_digits),
        "squared_chart": number_string(squared_value, shown_digits),
        "absolute_difference": number_string(difference, shown_digits),
        "direct_seconds": round(direct_seconds, 6),
        "squared_chart_seconds": round(squared_seconds, 6),
        "total_seconds": round(direct_seconds + squared_seconds, 6),
        "status": "PASS" if accepted else "FAIL",
    }
    print("BANANA_BESSEL_ORACLE " + json.dumps(record, separators=(",", ":")))
    return 0 if accepted else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, ZeroDivisionError) as error:
        print(f"banana_bessel_oracle.py: {error}", file=sys.stderr)
        sys.exit(2)
