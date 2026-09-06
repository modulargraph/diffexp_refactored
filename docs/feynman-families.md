# Configuring Feynman-trick families

`diffexp ft family.json --json` accepts an explicit family definition. The
calculation uses its propagators, scalar products and requested powers; its
`name` is a label. No family registration or source change is required.
Every supplied example has an editable JSON configuration in
[`examples/feynman`](../examples/feynman/).

For example, this two-propagator family has masses squared 2, external momentum
squared -3 and dimension 2-2 epsilon:

```json
{
  "schema": "DiffExp.FeynmanFamily/v1",
  "name": "my-two-propagator-family",
  "loops": 1,
  "external_gram": [["-3"]],
  "propagators": [
    {"loop_coefficients": ["1"], "external_coefficients": ["0"], "mass_squared": "2"},
    {"loop_coefficients": ["1"], "external_coefficients": ["-1"], "mass_squared": "2"}
  ],
  "physical_count": 2,
  "dimension_at_epsilon_zero": 2,
  "integrals": [[1, 1]],
  "epsilon_order": 0,
  "numerical": {
    "endpoint_order": 32,
    "ordinary_order": 80,
    "working_bits": 384,
    "leaf_digits": 28,
    "method": "adjoint"
  }
}
```

Each propagator is `(sum a_i k_i + sum b_j p_j)^2 - mass_squared` in the
Minkowski convention. `external_gram[i][j]` is `p_i . p_j`. The dimension is
`dimension_at_epsilon_zero - 2 epsilon`. Integers and rational strings are
accepted for exact data; floating-point JSON coefficients are rejected.

The first `physical_count` lines are physical denominators; any remaining
lines specify auxiliary denominator coordinates. Missing scalar-product
coordinates are completed by the native IBP machinery. Each `integrals` row
contains integer powers in this order; omitted auxiliary powers are padded
with zero. If `integrals` is omitted, the target has unit physical powers and
zero auxiliary powers. Negative powers represent numerators. All requested
integrals share preparation and recursive source work.

```sh
build/diffexp family-template sunrise > my-family.json
# Edit masses, scalar products, propagators, powers and controls.
build/diffexp prepare my-family.json --fire /path/to/FIRE7 --cache /path/to/cache --json
build/diffexp ft my-family.json --fire /path/to/FIRE7 --cache /path/to/cache --json
```

Use `-` as the configuration path to read JSON from stdin. Template names
remain a convenience, and older named CLI invocations still work. The primary
interface is the complete configuration. `prepare` constructs and verifies the
recursion without reporting an integral value. Cold IBP preparation requires
FIRE when the internal scalar/reduction cases do not suffice. Previously
verified preparation can be reused from the selected cache.

For runs that do not need numerical restart files, add `--no-numerical-cache`.
This avoids writing endpoint-series and ordinary-transport checkpoints while
retaining the exact reductions under `--cache`. In-memory sharing remains active.
The response reports `persistent_numerical_cache: false`; numerical stages are
recomputed on a later invocation. This is useful when disk space is limited.
See [the laptop examples](ft-laptop-examples.md) for cold runs with this option.

Optional `preparation` fields are `anchors`, `merges`, `total_seconds`,
`level_seconds`, `fire_seconds`, and `max_sources_per_level`. A merge is a pair
of zero-based positions in the current physical propagator list, after earlier
merges. Anchors are exact rationals. CLI resource options override JSON fields.
Optional numerical fields also include `contour_height` and `overlap`.

For Euclidean families the native code verifies positivity of the actual
Symanzik polynomials before selecting its contour prescription. Other physical
regions require an explicit `causal` object:

```json
"causal": {
  "f_rim": -1,
  "level_signs": [1, -1],
  "provenance": "Description of the continuation prescription for these kinematics"
}
```

There must be one nonzero sign per merge level. A prescription is mathematical
input; the parser does not establish that an arbitrary supplied contour is the
correct physical continuation. Unknown fields and inconsistent dimensions are
rejected. The experimental full Henn reconstruction remains frozen; including
its configuration does not claim a completed all-component calculation. The
optional `--henn-basis` importer applies the published Henn canonical projection
only after verifying the actual X0 geometry and denominator order; changing a
label cannot authorize that projection for another integral family.

The response reports separate preparation, numerical and total seconds, cache
build/reuse counts, requested coefficient windows and `omitted_tails_certified`.
General recursion still has the numerical limitations described in
[validation.md](validation.md).

In Mathematica the same configuration is an ordinary Association:

```wolfram
family = DiffExpFamilyTemplate["sunrise"];
family["name"] = "my-family";
family["external_gram"] = {{"-3"}};
family = ReplacePart[family,
  {Key["propagators"], 1, Key["mass_squared"]} -> "2"];
result = DiffExpFeynmanTrick[family,
  {"--fire", "/path/to/FIRE7", "--cache", "/path/to/cache"}];
result["timings"]
```

A JSON filename can be passed in place of the Association. Editing a label
alone does not change the integral; editing the geometry does.
