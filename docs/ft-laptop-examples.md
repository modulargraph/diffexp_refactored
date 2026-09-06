# Four FT examples that fit on a laptop

All four examples completed and passed independent numerical checks. These are full recursive Feynman-trick calculations from propagator configurations, including exact reductions and numerical boundary generation. No supplied numerical boundary is required.

## Cold timings

Apple M4, macOS ARM64, Release C++20 build, FLINT 3.4.0. Runs were sequential, with one FIRE worker and one simplifier worker, fresh exact caches, endpoint order 32, ordinary order 80 and 384 working bits. The outer process-group limit was 600 seconds per integral, including preparation. No precision or order was reduced to make a run fit.

| Family | Loops / physical propagators | Coefficients checked | Cold wall time | Preparation | Numerical calculation |
|---|---:|---|---:|---:|---:|
| Massive pentagon, five different masses | 1 / 5 | epsilon 0, 1, 2 | **7.85 s** | 1.88 s | 5.94 s |
| Massless box-triangle | 2 / 6 | epsilon -4 through 0 | **18.88 s** | 13.32 s | 5.54 s |
| Unequal-mass banana | 3 / 4 | finite part | **23.86 s** | 3.71 s | 20.12 s |
| Massless planar double-box | 2 / 7 | epsilon -4 through 0 | **7m 45s** | 26.93 s | 7m 18s |

Combined solver wall time was **8m 35s**. Independent reference checks took about 1.41 seconds in total and are excluded from the solver table. Compilation and the initial disk-space interruption are also excluded. These are single cold observations, not statistical averages or promises for arbitrary kinematics. “Cold” means an empty package preparation cache, not flushed operating-system caches.

The configurations are [massive pentagon](../examples/feynman/pentagon_massive.json), [box-triangle](../examples/feynman/box_triangle.json), [unequal three-loop banana](../examples/feynman/banana_unequal.json), and [planar double-box](../examples/feynman/double_box_planar.json). The runner changes their display labels to `laptop_...` and passes the complete JSON to `ft`; it does not use a family-name solver shortcut. No geometry is changed.

The massive pentagon has squared masses `(1, 3/2, 4/3, 5/4, 6/5)` and the supplied off-shell Euclidean Gram matrix. The banana has squared masses `(2, 3/2, 4/3, 1)`, external momentum squared `-1`, and dimension `2-2 epsilon`. The box-triangle and double-box use `s=-1`, `t=-1/3` in `4-2 epsilon` dimensions.

## Independent checks

- **Massive pentagon:** all three coefficients pass the `1e-18` absolute tolerance against independent Feynman-parameter integration values. These are numerical reference values, not rigorous reference enclosures. The leading coefficient is approximately `0.0181337866863019576423`.
- **Box-triangle:** all five Laurent coefficients pass `1e-20`, using the independent original-IBP/Mellin–Barnes reference described in [the reference note](box-triangle-reference.md). Its estimated reference error is `1e-30`, not a theorem bounding the remainder. The finite part is approximately `9.50480890034739241503`.
- **Unequal banana:** the finite part passes `1e-20` against certified coordinate-space Bessel quadrature. It is approximately `5.83402729266214946741`.
- **Planar double-box:** all five coefficients pass `1e-20` against the independent analytic Smirnov Laurent coefficients. The finite part is approximately `-247.003678220592567009`.

Forbidden lower poles also pass `1e-20`: one for the massive pentagon, three for the box-triangle, five for the banana and eight for the double-box. The reference checker verifies the configuration geometry and unit-power target before choosing a reference; changing only a display label does not affect it.

The recursive solver still reports `omitted_tails_certified: false`. Reference agreement validates these observations but does not turn the solver's retained arithmetic balls into certified full-integral enclosures. See [the raw results](validation/ft-laptop-examples.json) for every coefficient, reference comparison, pole check and separate timing.

## Avoiding large checkpoint files

The initial attempt was interrupted by a nearly full filesystem. Unused build outputs and package download caches were cleaned, and these successful runs used **`--no-numerical-cache`**. This keeps exact IBP reductions on disk and shares numerical work in memory, but skips the large endpoint-series and ordinary-transport checkpoint files. Later runs recompute the numerical stages. The mathematical recurrences and accuracy settings are unchanged.

All successful runs built their exact systems from scratch: 3, 4, 2 and 5 systems respectively, with zero previously cached systems reused. This timing table therefore does not benefit from saved numerical answers or completed-arm replay. The FT pipeline uses its own adjoint and polynomial recurrences.

## Reproduce

Build the executable and independent checker with tests enabled:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --target diffexp diffexp_check_ft_laptop -j2
python3 tools/performance/run_ft_laptop.py --fire /path/to/FIRE7 --output /tmp/diffexp-ft-laptop-new
```

Use an empty output directory. The runner records each actual configuration, process result, solver output and independent check. Select a subset with `--cases pentagon_massive banana_unequal`. Each solver process has a ten-minute cap and each reference process has a two-minute cap; timed-out process groups are terminated.

A direct generic invocation is:

```sh
build/diffexp ft examples/feynman/double_box_planar.json \
  --fire /path/to/FIRE7 --cache /tmp/double-box-exact \
  --no-numerical-cache --json
```

The ten-minute outer timeout is provided by the Python runner, not by `--total-seconds`, which only limits exact preparation. From Mathematica the same invocation is available through `DiffExpFeynmanTrick[configuration, {"--fire", "/path/to/FIRE7", "--cache", "/tmp/double-box-exact", "--no-numerical-cache"}]`.
