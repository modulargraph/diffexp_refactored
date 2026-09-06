# FT connection structure audit

This symbolic diagnostic reconstructs the actual homogeneous connection used by `recursion::Evaluator::plan`: an exact closure matrix followed by the existing diagonal epsilon gauge. It does not compute endpoints, forcing rows, numerical transport, or integral values. Indices in the JSON are one-based. The adjoint equation transposes/reverses dependency direction; SCC sizes are unchanged.

Run with Wolfram Kernel available:

```
python3 audit.py --output audit-results.json pentagon_massive=/path/to/cache box_triangle=/path/to/cache banana_unequal=/path/to/cache sunrise=/path/to/cache
```

The supplied results use retained September 6 laptop and AMFlow-comparison exact caches. The extra `sunrise_fixture` record is the existing product-recurrence fixture, separately labelled. Source paths and ordered bases are recorded. There are no Henn inputs. Wolfram is used only for this optional diagnostic, not by the proposed numerical backend.

| Family | Physical propagators at closure | Dimension | Epsilon-zero SCC sizes | Nonzero offdiagonal A0 entries |
|---|---:|---:|---|---:|
| pentagon_massive | 2 | 3 | 1,1,1 | 2 |
| pentagon_massive | 3 | 7 | seven scalar | 10 |
| pentagon_massive | 4 | 15 | fifteen scalar | 34 |
| box_triangle | 2 | 2 | 1,1 | 1 |
| box_triangle | 3 | 3 | 1,1,1 | 2 |
| box_triangle | 4 | 4 | 1,1,1,1 | 0 |
| box_triangle | 5 | 4 | 1,1,1,1 | 3 |
| banana_unequal | 2 | 3 | 2,1 | 4 |
| banana_unequal | 3 | 7 | 5,1,1 | 15 |
| sunrise | 2 | 3 | 1,1,1 | 3 |

The comparison bubble uses one scalar-leaf level and builds no connection system; it offers no ordinary-stage spectral target. All native epsilon-gauge shifts are zero except box_triangle at physical_count=5, whose shifts are [1,1,1,0].

## Concrete transforms

Write physical columns as `y=diag(g_i) z`; the new diagonal is `A_ii-g_i'/g_i`. Every box_triangle A0 diagonal is a rational logarithmic derivative. The physical_count=4 connection has the same diagonal `1/x+1/(x-1)` in all four components and no A0 offdiagonals: the common gauge `g=x(x-1)` removes A0 completely. This is the cheapest full canonicalization target.

Other box_triangle rational gauges, in stored basis order (arbitrary nonzero constants omitted):

- physical_count=2: `x(77x-81)^3/(7x-27)^5`, `1/(259x-270)`.
- physical_count=3: `x(x-1)^3/(7x-9)^5`, `1/(47x-45)`, `(x-1)^2/[x^2(25x-27)]`.
- physical_count=5, after existing epsilon gauge: `x`, `x-1`, `1`, `1`.

All pentagon A0 diagonals are logarithmic derivatives with integer or half-integer residues; removing them generally requires polynomial square roots. The JSON lists exact polynomial factors and exponents, so `g_i=product(factor^exponent)`. Such gauges require continued branches on ordinary arms. Banana's 5-component SCC contains three diagonal entries that cannot be expressed as constant combinations of irreducible polynomial logarithmic derivatives; diagonal rational/algebraic gauges are insufficient.

After diagonal removal, all pentagon and box_triangle matrices admit an additional epsilon shear making every nonzero entry strictly positive in epsilon (maximum additional shift two and one, respectively). The JSON contains exact shift vectors. Neither banana nor sunrise admits this shear, even after removing all A0 diagonal entries. Sunrise's higher-epsilon feedback obstructs positive shearing despite scalar A0 SCCs.

For production, factoring the epsilon-zero node-by-SCC operator once and reusing it across epsilon orders avoids the extra boundary/child demands caused by these shears. Rational gauges can further simplify selected scalar operators. Apply a gauge only on the regular ordinary arms and transform forcing and boundary rows consistently; the existing singular endpoint construction should remain in its native basis.

This is a structure audit, not a numerical convergence or independent integral validation. The log-derivative test matches exact coefficients over the rational polynomial field; `log_derivative=false` rules out that form, not arbitrary non-diagonal transformations.
