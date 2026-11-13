#!/usr/bin/env python3
"""Give a physical-point family a Euclidean continuation direction.

    python3 tools/add_continuation.py PHYSICAL.json EUCLIDEAN.json OUT.json

The direction V is the Gram matrix and masses of EUCLIDEAN.json (same
propagators, a point where U and F have positive monomials). DiffExp then
starts at P + t V and continues to P along a path in the upper half plane
(docs/feynman-families.md, "Continuation from a Euclidean point").
Numerator (ISP) slots get mass 0 in the direction.
"""
import json
import sys

physical, euclidean, out = sys.argv[1:4]
p = json.load(open(physical))
e = json.load(open(euclidean))
assert len(p["propagators"]) == len(e["propagators"]), "different propagator counts"
for a, b in zip(p["propagators"], e["propagators"]):
    assert a["loop_coefficients"] == b["loop_coefficients"] and a["external_coefficients"] == b["external_coefficients"], \
        "the two files must have the same propagators"
count = p.get("physical_count", len(p["propagators"]))
masses = [line["mass_squared"] if i < count else "0" for i, line in enumerate(e["propagators"])]
p.pop("causal", None)
p["continuation"] = {"euclidean_direction": {"external_gram": e["external_gram"], "mass_squared": masses}}
json.dump(p, open(out, "w"), indent=1)
print(out)
