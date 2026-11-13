#!/usr/bin/env python3
"""Give a physical-point family a Euclidean continuation direction.

    python3 tools/add_continuation.py PHYSICAL.json EUCLIDEAN.json OUT.json

The direction V is the Gram matrix and masses of EUCLIDEAN.json (same
propagators, a point where U and F have positive monomials). DiffExp then
starts at P + t V and continues to P along a path in the upper half plane
(docs/feynman-families.md, "Continuation from a Euclidean point").
Numerator (ISP) slots get mass 0 in the direction.

A mass that is positive at the physical point keeps that value in the
direction (and so do the Gram entries tied to it, such as on-shell p^2 = m^2).
Otherwise the line would cross m^2 = 0 at x = 1 + m_P^2/(t m_V^2), right
next to the physical point when m_P^2 is small, and the transport could not
reach it. Pass --keep-direction-masses to use the Euclidean file's masses.
"""
import json
import sys
from fractions import Fraction

keep = "--keep-direction-masses" in sys.argv
physical, euclidean, out = [a for a in sys.argv[1:] if not a.startswith("--")][:3]
p = json.load(open(physical))
e = json.load(open(euclidean))
assert len(p["propagators"]) == len(e["propagators"]), "different propagator counts"
for a, b in zip(p["propagators"], e["propagators"]):
    assert a["loop_coefficients"] == b["loop_coefficients"] and a["external_coefficients"] == b["external_coefficients"], \
        "the two files must have the same propagators"
count = p.get("physical_count", len(p["propagators"]))
masses = [line["mass_squared"] if i < count else "0" for i, line in enumerate(e["propagators"])]
gram = [row[:] for row in e["external_gram"]]
if not keep:
    for i, (a, b) in enumerate(zip(p["propagators"][:count], e["propagators"][:count])):
        mp, mv = Fraction(a["mass_squared"]), Fraction(b["mass_squared"])
        if mp > 0 and mv != mp:
            for r in range(len(gram)):
                for c in range(len(gram)):
                    if Fraction(p["external_gram"][r][c]) == mp and Fraction(gram[r][c]) == mv:
                        gram[r][c] = a["mass_squared"]
            masses = [a["mass_squared"] if Fraction(m) == mv and Fraction(q["mass_squared"]) == mp else m
                      for m, q in zip(masses, p["propagators"])]
p.pop("causal", None)
p["continuation"] = {"euclidean_direction": {"external_gram": gram, "mass_squared": masses}}
json.dump(p, open(out, "w"), indent=1)
print(out)
