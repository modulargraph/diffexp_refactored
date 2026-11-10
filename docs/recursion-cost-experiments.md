# What limits the Feynman-trick recursion

Measurements made for 2.2.0 on a 4-core x86-64 Linux container (GCC 13, FLINT
3.4.0), prompted by the two-loop double pentagon, whose numerical recursion
runs for hours. Two candidate remedies were tested and neither is a large
win; the cost model that explains why is at the end. Timings taken while
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

It is not faster. One banana4 sample costs 65 s (shared) against 279 s for
the whole expanded run: the transport drops from 238 s to 20 s, but the exact
endpoint series (21 s) and the matching do not shrink. The series in epsilon
has a small radius of convergence (the first-order coefficient is -439), so
each pair of samples buys about 2.7 digits and 30 digits need about 22
samples: roughly 1200 s of CPU against 279 s. On four cores that is parity.
For the box, ten samples take 5.6 s against 0.1 s.

The reason is the cost model below: the expanded transport is already linear
in the number of epsilon orders, so one sample costs about one order, and
more samples are needed than orders are carried.

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
