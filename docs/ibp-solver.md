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

At each recursive stage, DiffExp owns the exact affine denominator basis and its inverse scalar-product map. The adapter compiles the IBP equation structure once per requested batch and seed policy using upstream `InputGeometry` and `ParametricProgram`. The dimension and changing contraction coefficients are explicit inputs. A full reduction learns the relevant equations; `ArithmeticTrace` replays their arithmetic at subsequent parameter points. Integral and row identities remain fixed even if a coefficient vanishes at a sample. The supplied denominator order and auxiliary coordinates are preserved; there is no family-name dispatch or fixed merge anchor substituted for the integration variable.

The integral ordering compares total extra denominator and numerator powers, eliminating numerator integrals first at equal degree. This general basis heuristic removes a numerical integration penalty seen in Sunrise. It does not claim to find an optimal basis for every family.

Each prime has its own trace. Zero pivots and sample-specific cancellations trigger full reduction and relearning. The first changed sample is also compared with full elimination. Reconstruction fits rational functions in every active variable, including the current merge parameter and dimension/epsilon. After degree discovery, later primes fit only the observed nonzero monomial support, using fewer samples. That support is a hypothesis that must pass independent validation.

Reconstruction starts at the 61-bit prime 2^61−1 and uses distinct descending 61-bit primes. Two further primes, each at three fresh points, validate the result **using full elimination, without the reconstruction trace**. Applicable native IBP identities are checked exactly; the level preparer separately checks exact derivative and target closure relative to the imported relations. These checks do not prove globally minimal masters or provide an unconditional exact certificate for every reconstructed coefficient.

`ibp_statistics` reports compiled templates, equations actually generated, full solves, successful trace replays and fallbacks. It separates full-elimination, trace-learning and trace-replay time. The first changed replay is also fully checked, so replay and full-solve counts can overlap. Generation and elimination totals include surrounding adapter work; the detailed times are components of those totals, not additional costs. Fresh native rows remain in memory instead of being serialized and immediately reparsed, while durable sample checkpoints are retained.

Completed reconstructions and sample records are stored under `--cache/ibp-solver`. Provider/version, seed controls, scientific inputs and prime policy identify sample records. Completion reload checks independently retained validation samples. The existing exact level cache can also reuse previously verified FIRE-produced closures; `systems_reused` and `ibp_statistics.fresh_probes` distinguish reuse from new solver work.

The adapter supports at most four loops, sixteen scalar products and twelve physical denominators, with finite time, seed-state, equation-term and pivot-fill budgets. These are upper bounds, not a promise that every family of that size fits the chosen budgets. Large existing families remain runnable with their previous explicit FIRE provider. The frozen Henn calculation was not run.

## Measured complete runs

Apple M4, native release build, separate cold caches and no persistent numerical checkpoint files. Both providers receive the same family configurations and requested epsilon coefficients. All figures below are medians of three interleaved runs. Complete times include process startup, preparation and numerical integration.

| Family | IBP Solver complete | FIRE7 complete | Complete speedup | IBP preparation / FIRE preparation |
|---|---:|---:|---:|---:|
| Sunrise | 0.89 s | 0.99 s | 1.12× | 0.121 / 0.218 s |
| Box | 0.35 s | 0.67 s | 1.95× | 0.047 / 0.420 s |
| Box-triangle | 9.27 s | 17.95 s | 1.94× | 4.260 / 12.943 s |

The first adapter rebuilt and fully eliminated a system at every point. Its box-triangle run took 37.13 seconds, including 32.09 seconds of IBP preparation; the [earlier measurements](benchmarks/ibp-solver/integration-before-parametric.json) are retained. Parametric equation reuse and guarded replay remove that repeated work. Full independent verification remains a substantial part of the improved box-triangle preparation cost.

The standalone package's 10.7× result measures fixed-kinematics rational reconstruction in the dimension for ten double-box targets. It is a different workload and is not an end-to-end FT speedup claim. The FT adapter currently uses CPU parametric replay; standalone Metal batch timings do not apply to these runs.

The provider tests also check ten ordinary double-box targets against exact FIRE reductions after basis conversion at two dimensions and two 61-bit primes. Complete native FT, cache reuse, Mathematica and installed C++ consumer tests pass.

The full comparison records, coefficient differences and executable hash are in [integration.json](benchmarks/ibp-solver/integration.json). Reproduce them with:

```sh
python3 docs/benchmarks/ibp-solver/run.py --fire /path/to/FIRE7 --output /new/result/directory --repeats 3
```

Reported Arb radii still omit the general recursion's unbounded series tails. Agreement between providers is a numerical regression check, not a new full-integral error certificate.

## C++ package and source provenance

Include `diffexp/ibp_solver_provider.hpp`, construct `diffexp::ibp_solver::Session` with the stage basis, dimension, exact field and a cache directory, then supply it through `diffexp::level::Provider`. `Sampler` is available for finite-field reductions. Qualify the existing `diffexp::ibp` namespace explicitly when also using the upstream `::ibp` types.

The unmodified core headers are pinned at upstream commit `8d803554c33ac691877ef72c3dd338929c862abc` under `third_party/ibp-solver`, with the upstream license. Builds and installed CMake consumers need no network download or separate IBP executable. The standalone repository retains its own CLI, reconstruction experiments and Metal batch implementation.
