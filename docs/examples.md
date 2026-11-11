# Examples

Each Feynman-trick example has a complete editable JSON file in
`examples/feynman/`. Pass that file to `diffexp ft` or import it as a Mathematica
Association. Change its propagators, kinematics and requested powers to evaluate
a new family; no name registration is needed. See
[family configurations](feynman-families.md).

For an independent comparison, the opt-in acceptance runner reports coefficients,
reference differences and forbidden-pole checks as JSON:

```sh
build/diffexp_ft_examples sunrise /path/to/cache
build/diffexp_ft_examples banana4 /path/to/cache
build/diffexp_ft_examples banana4_unequal /path/to/cache
```

Its optional arguments are `FAMILY CACHE FIRE ENDPOINT_ORDER ORDINARY_ORDER BITS
GRAPH_SECONDS METHOD`, with defaults 32/80/384/600/adjoint. Without a FIRE
path (or with `-` in its place) the built-in IBP solver prepares the levels,
with the settings of `diffexp ft`. The four-loop bananas take minutes and none
of these runs automatically.

`examples/feynman/henn_double_pentagon_x0.json` is experimental and has no
acceptance reference. `diffexp prepare` completes its exact recursion in about
three minutes. The epsilon-expanded `diffexp ft` then runs for hours; the
practical route is sampling,

```sh
build/diffexp ft examples/feynman/henn_double_pentagon_x0.json --cache /path/to/cache \
  --epsilon-strategy sample --samples 12 --json
```

([details and timing](performance.md)).

`examples/feynman/higgs_jet_family_f.json` holds the seven propagators of the
non-planar Higgs+jet family F of arXiv:1907.13156 (four massive lines) at the
Euclidean point s=-1, t=-2, u=-3, m^2=1, also without a reference. It prepares
in about a minute and a half and one epsilon sample
(`--epsilon-shift 1/400 --epsilon-order 0 --leaf-digits 12`) takes under two
minutes on four threads; two different merge ladders give the same sample
value to 28 digits.

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
