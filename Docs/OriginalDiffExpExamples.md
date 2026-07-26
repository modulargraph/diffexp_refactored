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
| `5pPlanar1Mass.nb` | one-loop and `zzz`/`zmz`/`mzz` planar one-mass canonical systems | next |
| `Banana.nb` | banana differential equations | pending direct notebook parity |
| `MultiplePolylogarithms.nb` | multiple-polylogarithm demonstration | pending |

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
| DiffExp 2 canonical transport | 8.27 s | `1.6e-12` |

Matrix extraction, system validation, and data import add about 0.38 s to
the DiffExp 2 run. The transport is about 6.01 times faster than the saved
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

## Optional timing gate

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
