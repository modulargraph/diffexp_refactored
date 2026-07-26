# Original DiffExp example recovery

This ledger tracks reproductions of the examples shipped with
[DiffExp 1](https://gitlab.com/hiddingm/diffexp). The goal is numerical
agreement with the original or published reference, plus a like-for-like
transport timing where the old notebook saved one.

The original repository revision used for the inventory is
`784c8229bf92369a03f011a48e161522c8c54bbd`.

| Original example | System | DiffExp 2 status |
| --- | --- | --- |
| `5pNonPlanar.nb` | Henn et al. nonplanar five-point canonical system, 108 masters and 31 letters | passing |
| `5pPlanar1Mass.nb` | one-loop and `zzz`/`zmz`/`mzz` planar one-mass canonical systems | passing |
| `Banana.nb` | banana differential equations | pending direct notebook parity |
| `MultiplePolylogarithms.nb` | multiple-polylogarithm demonstration | passing |

## Henn nonplanar five-point system

Run:

```sh
Scripts/fetch_henn_nonplanar_data.sh
wolframscript -file Examples/OriginalDiffExp/HennNonplanarCanonical.wl
```

The fetcher pins the three ancillary files from arXiv:1812.11160v2 by
SHA-256. They are not silently downloaded by tests.

The example transports the published boundary at

```text
(v1,v2,v3,v4,v5) = (3,-1,1,1,-1)
```

to

```text
(4,-113/47,281/149,349/257,-863/541).
```

It uses the canonical system

```text
d f = eps Sum_i C_i dlog(W_i) f
```

directly. Five letters contain the Gram-determinant square root. This does
not pass an irrational matrix through the rational `LoadSystem` solver;
`LoadCanonicalSystem` retains the constant matrices and algebraic letters
separately.

On the July 2026 development machine, with working precision 50, expansion
order 25, Padé matching, and ε orders 0 through 4:

| Implementation | Transport time | Maximum absolute reference error |
| --- | ---: | ---: |
| DiffExp 1 saved notebook output | 49.765221 s | approximately `1.7e-17` saved segment estimate |
| DiffExp 2 canonical transport | about 5.5 s | `1.8e-12` |

Matrix extraction, system validation, and data import add about 0.38 s to
the DiffExp 2 run. The transport is about 9.0 times faster than the saved
DiffExp 1 order-25 timing. Times are machine-dependent; correctness is the
primary contract.

The starting point lies on six linearly vanishing letters. Their combined
residue has rank 7. The boundary is a finite regular solution, so the
epsilon-triangular singular recurrence is

```text
(n+1) f[k,n+1] =
    R f[k-1,n+1] + Sum[m=0..n] B[m] f[k-1,n-m].
```

In particular, the residue acts on the previous epsilon coefficient. It
must not be moved to the left-hand side as an indicial solve for the current
coefficient.

## Planar one-mass canonical systems

Run one family or all four:

```sh
Scripts/fetch_planar_one_mass_data.sh
PLANAR_ONE_MASS_FAMILY=1loop \
  wolframscript -file Examples/OriginalDiffExp/PlanarOneMassCanonical.wl
PLANAR_ONE_MASS_FAMILY=all \
  wolframscript -file Examples/OriginalDiffExp/PlanarOneMassCanonical.wl
```

The fetcher pins the nine required files from arXiv:2005.04195v1 by
SHA-256. The example transports the published PH1 boundary to the PH6
reference for the 13-, 75-, 74-, and 86-master systems, through epsilon
order 4.

The original path crosses algebraic singular surfaces with a `+i0`
prescription. DiffExp 2 uses an endpoint-fixed, all-positive complex detour.
Slightly different positive direction magnitudes prevent exact cancellations
from leaving a singularity on the real line parameter without changing the
published homotopy.

Uniform charts are not valid here. One active singularity can lie only about
`5e-4` from the parameter interval, so a visually fine 200-chart grid still
places endpoints outside some Taylor convergence disks. Increasing the
expansion order then makes the answer worse. `CanonicalLineChartGeometry`
instead:

1. ignores alphabet letters whose constant matrices vanish;
2. eliminates the three square roots algebraically;
3. locates all finite complex zeros, poles, and branch points; and
4. keeps every chart endpoint within one third of the nearest singularity.

On the July 2026 development machine, at working precision 50 and expansion
order 25:

| Family | Masters | Certified charts | DiffExp 1 saved time | DiffExp 2 comparable time | Speedup | Maximum error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `1loop` | 13 | 107 | 85.556919 s | 24.19 s | 3.54x | `3.6e-12` |
| `zmz` | 75 | 121 | 1391.141062 s | 130.50 s | 10.66x | `2.2e-11` |
| `mzz` | 74 | 122 | 992.459149 s | 119.38 s | 8.31x | `3.2e-11` |
| `zzz` | 86 | 119 | 1670.664508 s | 145.29 s | 11.50x | `2.5e-11` |

The DiffExp 2 comparable time includes chart construction and the complete
transport call, but excludes ancillary import/system validation. Times are
machine-dependent. Correctness against every published endpoint coefficient
is the primary contract.

The canonical recurrence now evaluates the Taylor convolution one matrix
times a coefficient block at a time instead of issuing a separate matrix-
vector product for every `(matrix power, solution power)` pair. Sparse
constant matrices and their sparse linear combinations remain sparse.

## Multiple polylogarithms

Run:

```sh
wolframscript -file \
  Examples/OriginalDiffExp/MultiplePolylogarithms.wl
```

The original notebook evaluates a weight-\(w\) Goncharov polylogarithm by
constructing a triangular rational ODE. The DiffExp 2 port multiplies that
connection by an auxiliary epsilon and reads the coefficient of
`eps^w`. This is an exact weight grading:

```text
d f = eps Sum_i E_(i,i+1) dlog(t-a_i) f,
G(a_1,...,a_w;t) = coefficient of eps^w in f_1(t).
```

It also aligns ordinary GPL logarithms with DiffExp 2's exact
`(eps Log[t])^p/p!` sector normalization. The original shuffle rules are
retained to regularize trailing-zero base-point logarithms. A lower-half-
plane contour implements `t-a-I0`, reproducing Mathematica's principal
`+I Pi` continuation after a positive real letter is crossed.

All seven saved notebook evaluations pass. On the July 2026 development
machine:

| Case | Expansion order | DiffExp 1 | DiffExp 2 | Absolute error |
| --- | ---: | ---: | ---: | ---: |
| `G[1,0,1;4]` | 50 | 0.169930 s | 0.211 s | `7.1e-26` |
| `G[1,-10,0;4]` | 50 | 0.504031 s | 0.575 s | `2.7e-27` |
| complex-letter weight 4 | 50 | 0.141181 s | 0.053 s | `7.7e-31` |
| `G[1,0,1;4]` | 75 | 0.306044 s | 0.371 s | `1.4e-37` |
| `G[1,-10,0;4]` | 75 | 0.657636 s | 0.960 s | `1.3e-37` |
| complex-letter weight 4 | 75 | 0.138740 s | 0.090 s | `8.8e-41` |
| `G[1,...,20;21]` | 100 | 38.701777 s | 64.51 s | `7.2e-42` |

The weight-20 case uses 106 affine charts at a one-half clearance ratio.
It is slower than the original notebook's specialized Möbius path, but
retains roughly 41 correct digits and completes in about a minute.

## Optional timing gates

The full 108-master example is deliberately opt-in:

```sh
DE2_RUN_HENN_NONPLANAR_TIMING=1 Scripts/run_release_tests.sh
```

or directly:

```sh
Scripts/run_henn_nonplanar_timing_regression.sh
```

The default hard deadline is 120 seconds and the default transport ceiling
is 30 seconds. Override them with
`HENN_NONPLANAR_DEADLINE_SECONDS` and
`HENN_NONPLANAR_MAX_TRANSPORT_SECONDS`.

All four planar systems are also opt-in:

```sh
DE2_RUN_PLANAR_ONE_MASS_TIMING=1 Scripts/run_release_tests.sh
```

or directly:

```sh
Scripts/run_planar_one_mass_timing_regression.sh
```

The default hard deadline is 900 seconds. Per-family ceilings default to
60 seconds for `1loop`, 300 seconds for `zmz` and `mzz`, and 350 seconds for
`zzz`; override them with the corresponding
`PLANAR_ONE_MASS_MAX_<FAMILY>_SECONDS` variables.

The complete MPL notebook is opt-in:

```sh
DE2_RUN_MPL_TIMING=1 Scripts/run_release_tests.sh
```

or directly:

```sh
Scripts/run_mpl_timing_regression.sh
```

The default hard deadline is 180 seconds and the weight-20 ceiling is 120
seconds. Override them with `MPL_DEADLINE_SECONDS` and
`MPL_MAX_WEIGHT20_SECONDS`.
