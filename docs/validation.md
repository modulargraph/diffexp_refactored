# Validation and limitations

The release retains the mathematical implementation tested under the DiffExp3
working name. Packaging changes rename the public C++ namespace, include paths,
CMake package and executable; the CLI adds JSON output and stdin requests.
Scientific cache identifiers remain unchanged.

The publication check passes **91/91 native, original-example, external-reduction
and interface tests**, plus **52/52 optional compatibility-runtime tests**.
A clean source-only build, installed external C++ consumer, and installed
Mathematica wrapper with a cached sunrise calculation also pass. See the
[publication report](validation/publication.json).

## Numerical coverage

These are recorded end-to-end acceptance results from development. Renaming the
package does not constitute a new independent numerical run of each expensive
family. The publication regression and wrapper checks are recorded separately
in `validation/publication.json`.

| Example | Accepted coefficient window | Scope |
| --- | --- | --- |
| Specialized bubble | epsilon 0..4 | Certified enclosures, including omitted tails |
| Sunrise; equal/unequal three-loop bananas | finite part | Independent numerical reference comparisons |
| Equal/unequal four-loop bananas | finite part | Independent comparisons and retained forbidden-pole checks |
| Box; massless pentagon | epsilon -2..0 | Numerical reference comparisons |
| Massive pentagon | epsilon 0..2 | Numerical reference values; tolerance 1e-18 |
| Box-bubble | epsilon -3..0 | Numerical reference comparisons |
| Box-triangle; planar double-box | epsilon -4..0 | Numerical reference comparisons and lower-pole checks |
| Kite | finite part | Independent quadrature with an estimated error; tolerance 5e-7 |
| Henn canonical component 1 | epsilon 0..4 | Numerical reconstruction against the published boundary |
| Henn all 108 components | Not accepted | Experimental; frozen by the user on 2026-09-05 |
| Original DiffExp examples | Nine registered groups | Transport/reference tests, including both unequal-banana routes |

For the equal-mass four-loop banana, the recorded finite part is
`39.655526834297652529992823046933581156446060`.
For the unequal-mass case it is
`28.6471675002335701630167055405033591695048469`.
Their differences from the independent references are below `8e-36`, at endpoint
order 32, ordinary order 80, 384 working bits and 28 leaf digits. Both passed all
ten retained forbidden-pole checks at absolute tolerance `1e-20`. The
[equal-mass report](validation/ft-banana4-full-conditioned-acceptance.json) and
[unequal-mass report](validation/ft-banana4-unequal-continuation-acceptance.json)
preserve the comparison data. The calculations took approximately 20 and 46
minutes respectively on the development machine, including reference work.

## What the error statements mean

General FT results enclose arithmetic on retained series, but currently exclude
omitted endpoint and ordinary-transport tails. Numerical agreement and small
arithmetic radii are not a rigorous full-integral error certificate. Several
independent reference tables also have estimated or rounding-only errors.
The CLI makes the omitted-tail status explicit in JSON.

Exact reduction checks establish closure for the selected observables and
basis. They do not prove a global master count. Modular reconstruction uses
fresh-prime checks and applicable exact identities; those checks are not a
universal proof of every possible IBP relation.

The full Henn FT trial was stopped by user request after approximately 81
minutes. It had prepared and verified endpoint work and reached the second
recursion level, but had not completed the physical-domain or 108-component
reference comparison. Its checkpoints and research logs remain preserved in
the development checkout. Henn transport from supplied published boundaries is
a separate, passing test. No full Henn reconstruction runs in the default suite.

## Reproduce the release checks

```sh
cmake -S . -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j2
ctest --test-dir build --output-on-failure -j2
```

See [examples](examples.md) to enable all original data comparisons. Optional
`DIFFEXP_FIRE_EXECUTABLE`, `DIFFEXP_FIRE_PRIME_EXECUTABLE` and
`DIFFEXP_HENN_BASIS_FILE` CMake paths enable external-reduction and ancillary
import checks. `DIFFEXP_WOLFRAMSCRIPT` enables the real Mathematica process-wrapper
test; an explicit `DIFFEXP_WOLFRAM_KERNEL` can select the licensed installation.
None of these optional dependencies is needed for the default native build.

The optional inherited C++ compatibility runtime can be tested separately with
`-DDIFFEXP_BUILD_KERNEL_RUNTIME=ON` and `ctest -R '^kernel\.'`. It is not used by
the Mathematica wrapper.
