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
build/diffexp prepare my-family.json --ibp-provider ibp-solver --cache /path/to/cache --json
build/diffexp ft my-family.json --ibp-provider ibp-solver --cache /path/to/cache --json
```

Use `-` as the configuration path to read JSON from stdin. Template names
remain a convenience, and older named CLI invocations still work. The primary
interface is the complete configuration. `prepare` constructs and verifies the
recursion without reporting an integral value. Cold IBP preparation defaults to the bundled IBP Solver core when no explicit
FIRE executable is supplied. It reconstructs the full parameter-dependent
coefficients; no family-name special case is used. `--fire /path/to/FIRE7` and
`--fire-prime /path/to/FIRE7p` select the existing providers. The native provider
has bounded seed and reconstruction budgets; see [scope and controls](ibp-solver.md). Verified
preparation is reused from the selected cache.

For runs that do not need numerical restart files, add `--no-numerical-cache`.
This avoids writing the internal-basis choice and ordinary-transport
checkpoints while retaining the exact reductions under `--cache`. Endpoint
series are never cached. In-memory sharing remains active.
The response reports `persistent_numerical_cache: false`; numerical stages are
recomputed on a later invocation. This is useful when disk space is limited.
See [timings](timings.md) for cold runs with this option.

Optional `preparation` fields are `anchors`, `merges`, `total_seconds`,
`level_seconds`, `fire_seconds`, `max_sources_per_level`, `ibp_provider`,
`ibp_dots`, and `ibp_numerators`. Provider choices are `auto`, `ibp-solver`,
`fire`, and `fire-modular`. Seed controls are integers from 0 through 8. A merge is a pair
of zero-based positions in the current physical propagator list, after earlier
merges. Without `merges` each level merges the pair with the fewest master integrals
([merge selection](performance.md#what-the-defaults-do)); explicit `merges`,
`anchors`, `causal.level_signs` or `causal.points` keep the ladder as written. Anchors are exact rationals. CLI resource options override JSON fields.
Without `anchors` a level whose endpoint qualifies takes its boundary at x = 0 or 1 instead of at an interior anchor
([endpoint boundary](performance.md#what-the-defaults-do), `--boundary anchor|endpoint|auto`); the report's `ladder` lists the
masters and the boundary of each level.
Optional numerical fields also include `contour_height`, `overlap`, `overlap_ratio`
and `endpoint_order` (the last two are chosen together when omitted), and
`epsilon_strategy` (`"expand"` or `"sample"`), `samples` and `sample_step` for
[epsilon sampling](performance.md#epsilon-sampling).

## The Feynman i0

Denominators are `D = m^2 - q^2` and the propagator is `1/(D - i0)`, so the
second Symanzik polynomial enters as `F - i0`. Nothing has to be stated for
that: the contour of every level is derived from the family.

- A family whose `U` and `F` have only positive monomials (the Euclidean
  region) has no cut on the integration domain. This is checked on the exact
  polynomials; every real singular point of a level is then passed on the
  same side.
- Otherwise each real singular point `x*` of each level is treated on its
  own. With `F(alpha, x)` the second Symanzik polynomial of the level's merged
  family, the contour passes on the side

  ```
  sign(Im x) = -sign(dF/dx)   at the Landau pinch of F at x*,
  ```

  the side on which `F` acquires a negative imaginary part where it
  vanishes. A point without a pinch (a pseudo-threshold, or a pole of the
  master basis only) is free: the physical branch is analytic there and both
  sides give the same value. Neighbouring points share a detour only if they
  ask for the same side; between two that ask for opposite sides the path
  returns to the real axis (an epsilon sample, whose path is lifted off the
  axis, crosses it there).
- The scalar leaf integrates `F(x)^(-a + L eps)` in closed form over `[0,1]`.
  On `F < 0` its logarithm is `log|F| - i pi`, and a simple zero of `F`
  inside the interval is passed off the axis on the side
  `Im x = -sign(F'(x*))`, inside a disk that holds no other singular point of
  the integrand; the logarithm is continued along that path and compared with
  the rim value on return to the axis.

How a point is decided, in ball arithmetic on its isolating interval:

1. Every monomial coefficient of `F` in the Feynman parameters is nonzero at
   `x*` with one common sign: `F` does not vanish on the open simplex near
   `x*`. Free.
2. Every monomial coefficient of `dF/dx` is identically zero or nonzero at
   `x*`, the nonzero ones with one common sign `s`: `Im F` has the sign of
   `s Im x` on the whole real simplex. The side is `-s`. This holds for any
   number of loops and for massless lines (it is the rule "`s -> s + i0`"
   when `F` depends on an invariant monotonically).
3. One-loop levels, where `F = alpha^T M(x) alpha`: a pinch on a face `S` of
   the simplex is a positive null vector of the principal submatrix `M_SS`.
   A face with `det M_SS(x*) != 0` has none. A face whose determinant
   vanishes for every `x` (soft and collinear regions of massless lines) is
   skipped when the rank at `x*` equals the generic rank: its null vectors
   are present for every `x` and `dF/dx` vanishes on them. On a remaining
   face, entries of `dM_SS/dx` of one sign give the side; on an edge with a
   massless line, where `F = a_j (2 b a_k + c a_j)`, the zero of the second
   factor enters through the massless vertex with the sign of `b'`; otherwise
   the null vector is read from the adjugate and is either of mixed sign (no
   pinch) or positive, with the sign of `alpha^T M' alpha`.

Statement 2 gives a side without asking whether the point has a pinch at all;
a point reported as oriented may therefore be free, which is harmless. The
converse is what matters: a point is called free only when a pinch is
excluded. A side is never guessed. If none of these statements decides a
point that the path has to pass, the run stops:

```
the causal side of the real singular point x = 0.0912... of level 1 could not
be derived from the level's Symanzik F (...); supply an explicit `causal`
object that states its side, "causal": {"points": [{"level": 1, "x":
"0.0912801957751", "side": 1}]}, with side 1 (above the real axis) or -1
```

What is not decided:

- Levels with two or more loops beyond statements 1 and 2. Their faces with
  vanishing `U` need the Symanzik polynomials of contracted subgraphs. A
  two-loop level above a threshold in which `F` is not monotone (the equal
  mass sunrise at `p^2 > 9 m^2`) is refused.
- One-loop pinches that sit on the boundary of a face, which happens when a
  massless line meets a threshold of its neighbours in a way statement 2 does
  not cover, and points where two pinches ask for opposite sides.
- A repeated zero of the leaf's `F` inside `(0,1)` (a family exactly at a
  threshold): the contour is pinched and the run stops with
  `F has a repeated zero at x = ... the contour is pinched there`.

A detour of finite size is the infinitesimal one only if no singular point
of the level lies between it and the real axis. When points are treated one
by one, the complex roots of the level's denominators are therefore kept out
of every detour: a triangle shrinks, a group whose rectangle would hold one
is split, and the lifted path of an epsilon sample is lowered below them.
(The contour of a Euclidean family, and one stated by `level_signs` alone,
is laid out without this check.) The JSON report has `causal.mode`
(`euclidean`, `derived` or `supplied`) and, when derived, the oriented and
the free points of each level.

Physical-region examples: `bubble_timelike.json`, `bubble_above_threshold.json`,
`triangle_physical.json` (its first level passes its two cuts on opposite
sides) and `box_physical.json` (massless box at `s > 0`, `t < 0`).

### Explicit contour

An optional `causal` object overrides the derivation:

```json
"causal": {
  "f_rim": -1,
  "level_signs": [1, -1],
  "provenance": "Description of the continuation prescription for these kinematics"
}
```

With `level_signs` there must be one nonzero sign per merge level, and every
real singular point of a level is passed on that level's side; the scalar
leaf uses `f_rim` (and derives the side of a zero of its own `F` from it).
`provenance` is then required. Such a prescription is mathematical input; the
parser does not establish that a supplied contour is the correct physical
continuation, and one side per level is not enough for every family.
Without `level_signs` the object only sets `f_rim` (`1` is the complex
conjugate prescription `F + i0`) and the sides are derived as above.

`points` states the side of single real singular points and takes precedence
over both the level sides and the derivation; it is how a point that could
not be derived is supplied:

```json
"causal": {
  "points": [{"level": 1, "x": "0.0912801957751", "side": 1}]
}
```

`level` counts merge levels from one, `x` (a rational or a decimal) locates
the point to within `1e-6`, and `side` is the sign of `Im x` on the detour.
An entry that names no real singular point of its level is an error that
lists the level's points. Unknown fields and inconsistent dimensions are
rejected.

A requested anchor that turns out to be singular in the master basis
the reduction selected is moved to the simplest regular rational within 3/20
that no other real singular point of the level separates from it; the progress
output names the replacement. The full Henn reconstruction remains
experimental; including its configuration does not claim a completed
all-component calculation. The
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
  {"--ibp-provider", "ibp-solver", "--cache", "/path/to/cache"}];
result["timings"]
```

A JSON filename can be passed in place of the Association. Editing a label
alone does not change the integral; editing the geometry does.

Ordinary transport inside the FT recursion uses Taylor charts
(`"transport": "taylor"` or `"auto"`).
