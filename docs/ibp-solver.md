# IBP Solver integration

DiffExp includes the reduction core from the separate [IBP Solver package](https://github.com/modulargraph/ibp-solver). `prepare` and `ft` use it when no explicit FIRE executable is supplied. The same provider works through Mathematica's ordinary `RunProcess` wrapper.

```sh
build/diffexp ft examples/feynman/sunrise.json --ibp-provider ibp-solver \
  --cache /path/to/cache --no-numerical-cache --json

# Existing providers remain available.
build/diffexp ft examples/feynman/sunrise.json --fire /path/to/FIRE7 --json
build/diffexp ft examples/feynman/sunrise.json --fire-prime /path/to/FIRE7p --json
```

`--ibp-provider auto|ibp-solver|fire|fire-modular` selects the provider. `auto` honors `--fire-prime` first, then `--fire`, otherwise selects IBP Solver. Explicit provider selection takes precedence. `backend-info` reports availability and the default. Scalar boundary stages need no IBP reduction and can report zero fresh probes.

In Mathematica:

```wolfram
family = DiffExpFamilyTemplate["sunrise"];
result = DiffExpFeynmanTrick[family,
  {"--ibp-provider", "ibp-solver", "--cache", "/path/to/cache"}];
result["timings"]
result["ibp_statistics"]
```

Alternatively set `preparation.ibp_provider` to `"ibp-solver"` in the family JSON or Association. `preparation.ibp_dots` and `preparation.ibp_numerators` set the initial seed bounds; CLI `--ibp-dots` and `--ibp-numerators` override them. Defaults are one additional dot and numerator degree two. The adapter raises these bounds when requested powers require it, up to the upstream cap of eight. Seed bounds define which equations are generated; all terms produced by those equations are retained. Unsupported or exhausted requests report an error, without silently substituting another provider.

## Mathematics and checks

At each recursive stage, DiffExp already owns the exact affine denominator basis and its inverse scalar-product map. The adapter converts its differentiated contractions to finite-field coefficients and calls the upstream `Geometry`, `generate`, `Field` and `Solver` interfaces. It preserves the supplied denominator order and auxiliary coordinates. It never dispatches on a family name or substitutes a fixed merge anchor for the variable being integrated.

The modular reconstruction driver fits rational functions in all active variables, including the current merge parameter and dimension/epsilon. Reconstruction starts at the 61-bit prime 2^61−1, then uses distinct descending 61-bit primes. Two further primes, each at three fresh points, validate the result. Applicable native IBP identities are checked exactly; the level preparer separately checks exact derivative and target closure relative to the imported relations. This is not a proof of globally minimal masters or an unconditional exact certificate for every reconstructed coefficient.

Completed reconstructions and sample records are stored under `--cache/ibp-solver`. Provider/version, seed controls, scientific inputs and prime policy identify sample records. Completion reload checks independently retained validation samples. The existing exact level cache can also reuse previously verified FIRE-produced closures; `systems_reused` and `ibp_statistics.fresh_probes` distinguish reuse from new solver work.

The adapter supports at most four loops, sixteen scalar products and twelve physical denominators, with finite time, seed-state, equation-term and pivot-fill budgets. These are upper bounds, not a promise that every family of that size fits the chosen budgets. Large existing families remain runnable with their previous explicit FIRE provider. The frozen Henn calculation was not run.

## Measured complete runs

Apple M4, native release build, separate cold caches and no persistent numerical checkpoint files. Both providers receive the same family configurations and requested epsilon coefficients. Sunrise figures are medians of three interleaved runs; box and box-triangle are single observations. These are complete process times, not isolated finite-field kernel speeds.

| Family | IBP Solver | FIRE7 | IBP preparation / FIRE preparation |
|---|---:|---:|---:|
| Sunrise | 1.21 s | 0.99 s | 0.326 / 0.219 s |
| Box | 0.41 s | 0.72 s | 0.110 / 0.422 s |
| Box-triangle | 37.13 s | 17.96 s | 32.09 / 12.91 s |

The two providers can choose different valid bases, which can also change numerical transport cost. The integration improves the box run; it does not yet beat FIRE on all complete recursive workloads. The multivariate sampling strategy currently rebuilds and solves the bounded finite-field system at each point. It does not yet exploit upstream arithmetic traces across varying merge parameters, or GPU batch replay. Consequently the standalone package's warmed trace/GPU speedups do not transfer directly to these complete FT timings.

The provider tests also check ten ordinary double-box targets against exact FIRE reductions after basis conversion at two dimensions and two 61-bit primes. Complete native FT, cache reuse, Mathematica and installed C++ consumer tests pass.

The full comparison records, coefficient differences and executable hash are in [integration.json](benchmarks/ibp-solver/integration.json). Reproduce them with:

```sh
python3 docs/benchmarks/ibp-solver/run.py --fire /path/to/FIRE7 --output /new/result/directory
```

Reported Arb radii still omit the general recursion's unbounded series tails. Agreement between providers is a numerical regression check, not a new full-integral error certificate.

## C++ package and source provenance

Include `diffexp/ibp_solver_provider.hpp`, construct `diffexp::ibp_solver::Session` with the stage basis, dimension, exact field and a cache directory, then supply it through `diffexp::level::Provider`. `Sampler` is available for finite-field reductions. Qualify the existing `diffexp::ibp` namespace explicitly when also using the upstream `::ibp` types.

The unmodified core headers are pinned at upstream commit `574c5bd3a0140a141515facafb44b563727efd4f` under `third_party/ibp-solver`, with the upstream license. Builds and installed CMake consumers need no network download or separate IBP executable. The standalone repository retains its own CLI, reconstruction experiments and Metal batch implementation.
