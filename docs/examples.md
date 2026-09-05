# Examples

Each Feynman-trick example has a complete editable JSON file in
`examples/feynman/`. Pass that file to `diffexp ft` or import it as a Mathematica
Association. Change its propagators, kinematics and requested powers to evaluate
a new family; no name registration is needed. See
[family configurations](feynman-families.md).

For an independent comparison, the opt-in acceptance runner reports coefficients,
reference differences and forbidden-pole checks as JSON:

```sh
build/diffexp_ft_examples sunrise /path/to/cache /path/to/FIRE7
build/diffexp_ft_examples banana4 /path/to/cache /path/to/FIRE7
build/diffexp_ft_examples banana4_unequal /path/to/cache /path/to/FIRE7
```

Its optional arguments are `FAMILY CACHE FIRE ENDPOINT_ORDER ORDINARY_ORDER BITS
GRAPH_SECONDS METHOD`, with defaults 32/80/384/600/adjoint. Use `-` for FIRE to
require previously prepared exact caches. These calculations can take minutes
or longer and do not run automatically. The full Henn FT reconstruction is
frozen and is excluded from release acceptance.

## Original DiffExp examples

The three download scripts fetch upstream matrices and verify their SHA-256
hashes. They require curl and shasum or sha256sum. Reference boundary files are
included; downloaded upstream data are ignored by Git.

```sh
bash scripts/fetch_original_banana_data.sh
bash scripts/fetch_henn_nonplanar_data.sh
bash scripts/fetch_planar_one_mass_data.sh
cmake -S . -B build -DBUILD_TESTING=ON -DDIFFEXP_ORIGINAL_EXAMPLES_DIR="$PWD/examples/original"
cmake --build build -j2
ctest --test-dir build -L original --output-on-failure -j1
```

The nine registered groups cover Henn nonplanar, planar `1loop`, `zmz`, `mzz`,
`zzz`, equal-banana singular continuation, the historical equal-banana endpoint
32, both unequal-banana routes, and short/weight-20 MPL. MPL requires no external
data and is also included in the default suite. These are transport tests from
supplied boundaries, separate from FT reconstruction of a boundary.

Individual commands include:

```sh
build/diffexp mpl --short
build/diffexp henn-nonplanar examples/original/Data/HennNonplanar
build/diffexp planar examples/original/Data/PlanarOneMass 1loop
build/diffexp banana-equal examples/original --to 32
build/diffexp banana-unequal examples/original --route mass-first
build/diffexp singular-endpoint
```

`singular-endpoint` emits JSON with the retained Frobenius sectors, evaluated
endpoint and plot samples. The remaining original-example commands print
comparison reports and use a nonzero exit code on failure.
