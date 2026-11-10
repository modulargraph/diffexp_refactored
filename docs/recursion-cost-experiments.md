# What limits the Feynman-trick recursion

Measurements made for 2.2.0 on a 4-core x86-64 Linux container (GCC 13, FLINT
3.4.0), prompted by the two-loop double pentagon, whose numerical recursion
runs for hours. Changing the working precision is no remedy. Epsilon sampling
does not pay for the smaller families but does for the double pentagon, where
one sample now completes in under four minutes and samples run in parallel;
the cost model at the end explains both. Timings taken while
another job shared the machine are marked "shared".

## Threads (adopted)

| Case | Serial | Four threads |
| --- | ---: | ---: |
| Double pentagon, level-1 endpoint series (69 columns) | 583 s | 163 s |
| Double pentagon, level-3 lower transport leg | 335 s | 107 s |
| Kite transport (three threads, shared) | 9.2 s | 6.0 s |
| Banana, massive pentagon and smaller | unchanged | unchanged |

## Working precision

Lower precision and order help a little. Relative differences are against the
default run (384 bits, ordinary order 80, endpoint order 32).

| Settings | Kite | Planar double box |
| --- | --- | --- |
| 384 bits / 80 / 32 | 10.9 s | 30.2 s |
| 256 bits / 50 / 24 | 7.2 s, 4e-31 | 22.6 s, 1e-27 |
| 192 bits / 40 / 20 | 6.5 s, 3e-25 | 20.6 s, 2e-21 |

Higher precision does not remove the slow legs of the double pentagon. Its
level-3 upper leg took 784 s at 384 bits (the lower leg 107 s); at 640 bits it
was 57 % through after 600 s (shared), the same pace. The slow legs are not
caused by a lost precision reserve.

The `factored` method is no alternative: on the double box its transport took
68 s against 32 s for `adjoint` (shared).

## Epsilon sampling

`diffexp ft --epsilon-shift Q` (experimental) expands the prepared recursion
around epsilon = Q instead of zero. At a generic Q no level has epsilon poles,
so the order-zero result is the integral in d0-2Q dimensions and every level
needs only order zero of its child: the widening epsilon windows disappear.
The scalar leaf is re-expanded from its Laurent series at zero.
`tools/experiments/epsilon_sampling.py` runs symmetric samples and fits the
Laurent coefficients.

It is correct. Box (poles to epsilon^-2), ten samples at spacing 1/200:
12, -0.33491424680974, -41.2841673975516 against the expanded
12, -0.33491424680973618, -41.2841673975465. Four-loop equal-mass banana,
samples at +-k/400, against the reference 39.655526834297652...:

| Samples | Error of the finite part |
| ---: | ---: |
| 2 | 2.0e-2 |
| 4 | 1.5e-5 |
| 6 | 1.8e-8 |
| 8 | 3.3e-11 |

For the smaller families it is not faster. One banana4 sample costs 65 s (shared) against 279 s for
the whole expanded run: the transport drops from 238 s to 20 s, but the exact
endpoint series (21 s) and the matching do not shrink. The series in epsilon
has a small radius of convergence (the first-order coefficient is -439), so
each pair of samples buys about 2.7 digits and 30 digits need about 22
samples: roughly 1200 s of CPU against 279 s. On four cores that is parity.
For the box, ten samples take 5.6 s against 0.1 s.

The reason is the cost model below: the expanded transport is already linear
in the number of epsilon orders, so one sample costs about one order, and
more samples are needed than orders are carried. Samples are independent, so
with one core per sample the wall time is that of one sample.

For the double pentagon the balance is different, because its inner levels
carry epsilon through 12 to 16 and their legs slow down far more than
linearly. The first complete sample (epsilon = 1/100, four threads, 384 bits,
leaf digits 12) took 1360 s; its value is
-74166278.2233783249 - 17089677.0664946935 i.

| Leg | Expanded recursion | First sample |
| --- | ---: | ---: |
| Level 4, lower / upper | 26.8 / 19.4 s | 6.0 / 5.0 s |
| Level 3, lower / upper | 107 / 784 s | 24.7 / 70.0 s |
| Level 2, lower | 1790 s | 61.2 s |
| Level 2, upper | hours (not finished) | 163.5 s |
| Level 1, lower / upper | not reached | 21.3 / 60.9 s |

## Making a sample cheap

Three changes followed, each measured on four threads (numerical seconds;
the banana sample is at epsilon = 1/50 except in the last row, 1/400).

| | Four-loop banana, expanded | Four-loop banana, one sample | Double pentagon, one sample |
| --- | ---: | ---: | ---: |
| Start | 168.4 | 63.9 | 1360 |
| Levels specialised to the numeric dimension | 168.4 | 28.0 | 296 |
| Probe propagation of uncertainty | 53.9 | 14.4 | 227 |
| Taylor tail 2^-128 per chart for samples | 53.9 | 10.2 | 181 |

**Specialising.** Order zero of a recursion expanded around a generic
epsilon = Q depends only on the equations at that dimension. Every level
above the scalar leaf therefore substitutes epsilon = Q outright, and its
endpoint series, matching and transport handle rational functions of x with
no epsilon polynomials. The endpoint work of a double-pentagon sample fell
from about 945 s to about 35 s. One rule needed care: once epsilon is a
number, the dimensionally regulated endpoint sectors x^(a+b*epsilon) are no
longer marked by a slope. At a generic Q they are exactly the non-integer
powers, and the endpoint series labels them so. Values are unchanged to
every printed digit.

**Probes.** At a numeric dimension the connection has pairs of nearly
coincident poles with large opposite exponents (the moving poles of the
master basis, 1/Q apart in exponent). A chart's Taylor recurrence then
amplifies interval radii by 2^100 and more, although the linear map itself
is harmless, and the transport fell back to "midpoint plus chart map times
uncertainty": one expansion per master on nearly every chart. The same
happens at the one-row top level of the expanded four-loop banana, which
spent 118 of its 168 s there. The recursion now transports midpoints from the
first such chart of a leg and follows the size of the discarded radii along
two extra homogeneous rows (`--uncertainty probe`, the default;
`--uncertainty map` restores the map). The radii this produces are estimates
along two directions, not enclosures. On the examples they are of the size
of the rigorous ones or smaller, and all twelve acceptance families pass
their independent references.

**What did not help.** Plain interval transport with shorter steps loses
everything (radius 1e114): the growth is cumulative. Lower working precision
does not speed the recurrence up, because the time is in the bookkeeping of
`acb_dot`, not in the limb products. The forward ("factored") method is 1.4
times faster on a double-pentagon sample (164 s) because it moves one column
per leaf integral instead of 69 rows, but it runs toward the singular
endpoints and returned radius 1e41 at 384 bits; it needs 640. A looser
Taylor tail (`--transport-bits 120` instead of 160) removes a fifth of the
charts at an error near 1e-33.

Of a 227 s double-pentagon sample, 192 s is transport and about 150 s of
that the Taylor recurrences of the 69-row and 33-row levels.

## One core

The figures above are on four threads. The following changes are gains on a
single core; numerical seconds on one thread, against the code before any of
the changes in this note (the last column is the transport alone).

| | Before | Now | Transport |
| --- | ---: | ---: | ---: |
| Kite | 12.5 | 5.1 | 9.1 -> 1.9 |
| Planar double box | 40.8 | 18.7 | 34.5 -> 12.3 |
| Four-loop banana, expanded | 240.8 | 49.4 | 218.8 -> 27.6 |
| Four-loop banana, one sample (epsilon = 1/400) | 76.9 | 7.8 | 23.7 -> 4.0 |
| Double pentagon, one sample (four threads) | 1360 | 157 | 413 -> 123 |

**Real-axis paths.** Every leg used to rise to the contour height, run there
and come down, whether or not a singular point lay on its real segment. A leg
of the expanded recursion now stays on the real axis, where chart centres and
hence all shifted coefficients are real, and passes a real singular point of
its level over a triangle on the causal side (neighbouring points share one
rectangle). Passing a pole costs about a dozen charts whatever the size of
the detour, so detours are as wide as the free distance allows. A sampled
recursion keeps the lifted path: its moving poles are real points in nearly
coincident pairs all along each segment.

**Compact recurrence.** One chart of the double box performs 2.5 million
multiply-adds of an exact coefficient, mostly one limb, with a
working-precision value. Through `acb_dot` each cost 120-200 ns, and a
microbenchmark showed why: with operands scattered over a chart-sized
working set (a value is a 96-byte struct pointing at two heap blocks)
`acb_dot` takes 127 ns per term, on cached operands 55 ns, while one
`mpn_addmul_1` over contiguous limbs takes about 6. The recurrence now runs
on values of 2L contiguous limbs with limb-aligned exponents (no bit shifts
to align a product), stores the orders a product sum walks through next to
each other, and keeps one radius per value as a scaled double. Semantics are
those of the ball recurrence, which remains as reference and fallback and is
compared coefficient by coefficient in `tests/test_compact_recurrence.cpp`.
A replacement of `acb_dot` alone, on the old storage, gained 15-30 % in a
microbenchmark and nothing end to end.

**Tail from the leaf accuracy.** A chart's Taylor tail is now admitted at ten
digits below the requested leaf digits (2^-127 by default) instead of 2^-160:
about a fifth fewer charts, with results equal to the old ones to 1e-37.

**What limits a core now.** An instruction profile of the recurrence shows
one third limb arithmetic, a quarter the full-precision multiplications by
the inverse pivots of the quotient streams, and the rest bookkeeping spread
thinly over the product loop; in a sampled recursion the coefficients are
several limbs long and the products are at their arithmetic cost. Outside
the transport, the exact stages of the expanded recursion (endpoint
projection, endpoint series, matching coefficients, all rational-function
arithmetic in epsilon) are untouched and are now 35 to 60 % of a run: 3.2 of
5.1 s for the kite, 21.8 of 49.4 s for the four-loop banana. A fraction-free
rewrite of the projection was tried and was slower, because the stored
coefficients have a different epsilon-denominator at every order.

## The double pentagon, end to end

With a sample at 227 s on four threads (157 s after the single-core work
above) the whole family is within reach of
one machine: twelve samples at epsilon = +-k/400, four at a time on one
thread each, took 48 minutes of wall time (6312 s of CPU; the tool then
still ran its first sample alone, which costs eight of those minutes),
after the five-minute exact preparation. The fit of epsilon^4 times the scalar
integral gives

| Order | Value | Change from the last two samples |
| --- | --- | ---: |
| epsilon^-4 | -0.77777777777777779 + 0 i | 1e-14 |
| epsilon^-3 | -0.71244912607989874 - 17.453292519943297 i | 3e-14 |
| epsilon^-2 | 436.218113787 - 5.202032870 i | 2e-9 |
| epsilon^-1 | -88.950617318 + 4207.353609604 i | 6e-9 |
| epsilon^0 | -29747.6631 - 3507.3533 i | 1e-4 |

There is no independent reference for this integral in the repository (the
published boundary values are for a canonical basis whose top-sector members
carry numerators), so this is not an accepted result. Its leading terms have
the structure they must have: the epsilon^-4 coefficient is -7/9 to 14
digits with no imaginary part, and the epsilon^-3 coefficient matches
-1 - (5/9) log 3 + (14/9) gamma_E - (50/9) pi i to 15 digits, with the
14/9 = -2 * (-7/9) that a two-loop integral in this normalisation carries
and s12 = 3 the only invariant that is not +-1. The sample values are in
[validation/henn-double-pentagon-epsilon-samples.json](validation/henn-double-pentagon-epsilon-samples.json).

## Cost model

Each transport leg now reports its size. For the planar double box:

| Level | Rows x masters | Epsilon orders | Charts (lower / upper) | Seconds (lower / upper) |
| ---: | --- | ---: | --- | --- |
| 5 | 3 x 2 | 11 | 15 / 22 | 0.05 / 0.31 |
| 4 | 7 x 3 | 11 | 16 / 28 | 0.42 / 1.05 |
| 3 | 11 x 7 | 11 | 27 / 26 | 6.7 / 4.9 |
| 2 | 12 x 11 | 9 | 37 / 24 | 6.7 / 12.2 |
| 1 | 1 x 12 | 9 | 34 / 36 | 0.8 / 0.9 |

A leg costs about rows x connection size x epsilon orders x charts x Taylor
order. The row count is the number of integrals the parent level needs, so
the top level (one row) is cheap and the middle levels dominate. For the
kite, transport time is linear in the epsilon orders carried: 9.4 s, 14.0 s
and 19.1 s for 5, 8 and 11 orders.

The double pentagon has 64 rows x 33 masters with epsilon through 12 at level
2 and 33 x 23 through 14 at level 3; the top level needs epsilon^0 only. What
would reduce that, none of it measured:

- fewer epsilon orders at the inner levels, through master bases whose
  connection has fewer epsilon poles (each level currently adds several);
- fewer rows, by transporting the d x d solution map once per chart when a
  level has more rows than masters (at most a factor rows/masters, below 2
  here);
- numeric rather than exact endpoint series inside a sampled recursion,
  which would bring a banana4 sample to about 25 s.
