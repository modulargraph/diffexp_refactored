# Verified Results

This page separates mathematical values, numerical settings, and timed scope.
It is a curated release summary, not a benchmark claim for every machine.
Cold Feynman-trick runs include FIRE preparation; warm runs begin from a
prepared cache.  Those scopes are never mixed in one timing row.

## Small Euclidean fixtures

The repository contains parity-oracle logs for equal-mass examples in
``D=2-2 eps`` at ``p^2=-1``:

| Example | Reported epsilon-zero coefficient | Provenance |
| --- | ---: | --- |
| one-loop bubble | ``0.8608178819280080777765623653643473221`` | stored old-core parity oracle |
| two-loop sunrise | ``2.2367927002126465108229117827758723767`` | stored old-core parity oracle; displayed imaginary part about ``3.1e-21`` |
| three-loop equal-mass banana | ``8.2681045358689687593219901952256396`` | stored completed ladder log; displayed imaginary part about ``9.4e-17`` |

The corresponding records are
[bubble/sunrise](../Tests/refs/oracle_logs/l2_bubsun.log) and
[banana](../Tests/refs/oracle_logs/l2_banana.log).  The bubble and sunrise
files are parity artifacts rather than surviving independent pySecDec runs;
that distinction is recorded in
[their manifest](../Tests/refs/oracle_logs/MANIFEST.md).

## Unequal-mass three-loop banana

At

```text
m_i^2 = {2, 3/2, 4/3, 1},  p^2 = -1,  D = 2 - 2 eps,
```

the independent two-dimensional Bessel oracle is

```text
5.83402729266214946740741989567969814964058746213209...
```

Two compiled routes were compared on the development Mac:

| Route | Settings | Timed scope | Seconds | Scalar | Oracle agreement |
| --- | --- | --- | ---: | ---: | ---: |
| direct differential system (15 masters) | WP250, EO50, epsilon through 4, division 3 | exact slice audit, reconstruction, planning, transport | 46.770 | ``5.8340272926621494708226551820`` | 18.2 digits |
| Feynman-trick ladder | WP500, EO70, finite coefficient, extra order 10, halos ``0,4,7`` | analytic deepest boundary and full warm ladder | 390.286 | ``5.8340272926621494674057011919652947`` | 20.8 digits |

The direct timing excludes FIRE matrix generation and its frozen seed
generation.  The Feynman-trick timing excludes FIRE/IBP preparation.  The
direct route also returns four positive epsilon orders while the ladder row
targets the finite coefficient, so the table establishes agreement rather
than an identical-workload race.  Details and reproduction commands are in
[C++ Recurrence Backend](CppBackend.md#unequal-mass-three-loop-banana-comparison).

## Fully massive pentagon

For the one-loop five-propagator pentagon in ``D=4-2 eps`` with squared
internal masses ``{1,3/2,4/3,5/4,6/5}``, ``p_i^2=-1``, and
``p_i.p_j=1/4`` for ``i != j``, the completed ladder returned

```text
eps^0    0.018133786686301957642296020548326099321491282408573
eps^1    0.0076131154161440535648564784592130217117538189027361
eps^2    0.0052144755784776811411112368938047226922782076161285
```

All three coefficients differ by at most ``2.3e-21`` from independently
integrated, 20-digit Feynman-parameter pins.  The final WP300, expansion-order
50 numerical pass resumed at level 1 from the same campaign's boundary
checkpoint and warm FIRE preparation and took ``1900.513`` seconds (about
31 minutes 41 seconds).  This is deliberately reported as a resumed warm
timing, not as a cold end-to-end benchmark.  The independent evaluator is
[Scripts/pentagon_massive_oracle.m](../Scripts/pentagon_massive_oracle.m).

## Box-bubble

For the massless two-loop ``box_bubble`` fixture at ``s=-1``, ``t=-1/3`` in
``D=4-2 eps``, a from-scratch run at WP300, expansion order 40, division order
3, and boundary extra order 16 completed FIRE preparation and four DiffExp 2
levels in about 55 seconds on the development machine:

```text
eps^-3    0.4999999999999719409794842591023
eps^-2    0.4227843350986460144608829749893
eps^-1    0.3562795605799888368173675114533
eps^0    -4.877139662454516944362515776161
```

Every displayed coefficient agrees with the stored oracle to about
``1.4e-12`` absolute or better. The full development audit records the same
coefficients and normalization; this page keeps the curated release pin.

## Massive kite

The ``kite`` fixture is the fully massive equal-mass five-propagator kite in
``D=2-2 eps`` at ``p^2=-1``.  With WP300, expansion order 40, division order 3,
boundary extra order 16, and halos ``0,4,7,7``, the completed result was:

```text
eps^-1   -9.8559958318033201e-21
eps^0     0.223983919107444028404077222049
```

An independent direct Feynman-parameter integration gave
``0.223983917019259`` with a conservative ``1.8e-7`` error estimate.  This is
a useful normalization check, but it does not certify every displayed digit.

## Four-loop banana status

For the equal-mass five-line banana in ``D=2`` at ``p^2=-1``, the independent
configuration-space Bessel moment is

```text
16 Integral_0^Infinity r J0(r) K0(r)^5 dr
  = 39.6555268342976525299928230469335811564460602187101868181...
```

The standalone evaluator is
[Scripts/banana4_bessel_oracle.m](../Scripts/banana4_bessel_oracle.m). This
number is an independent
epsilon-zero oracle, not evidence that the complete DiffExp 2 four-loop ladder
finished.

The investigation is intentionally paused as of 2026-07-26 so development can
focus on other families. No solver value has passed this oracle. The canonical
handoff, including the last tested hypothesis and resume protocol, is
[`Banana4Status.md`](Banana4Status.md).

The unequal-mass fixture with squared masses
``{2,3/2,4/3,5/4,1}`` is present as ``banana4_unequal``. No completed,
source-controlled ladder result is currently present, so the
release example is explicitly experimental and this page does not publish a
solver value for it.

## Recurrence-kernel speedup

On the committed seven-master ``banana_L1`` chart fixture:

| Case | Wolfram recurrence | C++ recurrence | Speedup |
| --- | ---: | ---: | ---: |
| WP100, EO50, epsilon order 5 | 81.268 s | 4.380 s | 18.6x |
| WP500, EO50, epsilon order 11 | about 246.7 s | 11.492 s | about 21.5x |

These rows time the recurrence-oriented chart workload, not FIRE or an entire
Feynman-trick computation.  Small systems can be dominated by Wolfram
preparation, JSON/LibraryLink transfer, and validation.

## Original Henn nonplanar five-point example

The 108-master, 31-letter canonical system from arXiv:1812.11160 is
reproduced through epsilon order 4. At working precision 50 and expansion
order 25, the maximum absolute error against the published endpoint table is
`1.6e-12`.

The saved DiffExp 1 notebook reports 49.765221 seconds for its order-25
transport. DiffExp 2 takes about 5.5 seconds for canonical transport on the July
2026 development machine, plus about 0.38 seconds for matrix extraction,
validation, and ancillary-data import. See
[Original DiffExp example recovery](OriginalDiffExpExamples.md) for the
pinned data fetcher, settings, and opt-in timing gate.

## Original planar one-mass five-point example

All four canonical systems from arXiv:2005.04195 are reproduced from the
published PH1 boundary to PH6 through epsilon order 4. At working precision
50 and expansion order 25, their maximum absolute errors range from
`3.6e-12` to `3.2e-11`.

Clearance-certified algebraic charts are essential: uniform charts can
cross a Taylor convergence disk even when the grid looks fine, and raising
the series order then amplifies the error. Including chart planning and the
complete transport call, DiffExp 2 takes about 24, 130, 119, and 145 seconds
for `1loop`, `zmz`, `mzz`, and `zzz`, respectively, on the July 2026
development machine. The corresponding saved DiffExp 1 times are about 86,
1391, 992, and 1671 seconds. See
[Original DiffExp example recovery](OriginalDiffExpExamples.md) for exact
settings, hash-pinned ancillary data, and the opt-in timing gate.

## Reproducibility checklist

A release result should record:

- propagators and whether each mass entry is a mass or squared mass;
- loop/external momenta, invariant replacements, and the dimension convention;
- master normalization and any epsilon prefactor;
- backend and thread count;
- working precision, expansion order, epsilon order, division order, radius,
  boundary lookahead, and per-level halos;
- whether the run was cold, warm-prepared, or resumed;
- the independent oracle and its own uncertainty;
- full Laurent coefficients, not only a field named ``Finite``.

Development campaign logs and large oracle dumps remain on the development
branch. This release retains only the curated provenance above and the focused
fixtures needed by its tests; see [Release Manifest](ReleaseManifest.md).
