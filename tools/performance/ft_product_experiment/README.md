# Individual rational-product recurrence experiment

**Not adopted:** the manufactured chart improved by 1.37×, but the cached sunrise
matrix slowed by 3.35×. Production transport is unchanged. See
[the comparison report](../../../docs/amflow-comparison.md) and
[recorded observations](../../../docs/validation/amflow-product-experiment.json).

This probes a small part of AMFlow's rational off-diagonal forcing strategy,
not its complete SCC algorithm. For each epsilon-expanded entry p(x)/q(x), retain
f=(p/q)y and apply finite polynomial recurrences to f before accumulating y′.
The production baseline already uses row-wise denominator clearing. The candidate
avoids row-LCM growth at the cost of additional per-edge jets and divisions.
It preserves ball arithmetic and incoming uncertainty. It only supports expanded
epsilon entries with real rational x coefficients at center zero; unsupported
inputs delegate to production. Explicit storage/degree/operation caps remain.

The manufactured case has 10 components, 55 couplings, epsilon orders 0–3,
Taylor order 128 and 256 bits. It checks exact Rational Taylor coefficients against
independent dense convolution, enclosure of three perturbed input realizations,
production ball overlap and unsupported-input fallback. Warm times average five
calls. First calls follow small validation and are not cold-process measurements.

The cached sunrise test uses the exact 3×3 closure matrix, epsilon orders 0–2,
Taylor order 80 and 384 bits. Coefficients are translated exactly to center 1/2 and
evaluated over step 1/64. A saved ball vector supplies identical test inputs to all
methods. This is a cached-matrix microbenchmark, **not a replay of an original
physical chart**: the row checkpoint does not encode the prepared stage gauge/path.
It checks both payload hashes, production overlap, an independent dense recurrence
and agreement with the unshifted-coordinate calculation. The two small original
fixtures are retained unchanged, including their certificates.

Neither probe bounds omitted Taylor tails. A speedup on this manufactured system
would not establish a whole-family speedup.

From the repository root on the measured Homebrew setup:

```sh
c++ -std=c++20 -O3 -DNDEBUG -Iinclude -Itools -I/opt/homebrew/include \
  tools/performance/ft_product_experiment/compare.cpp \
  -L/opt/homebrew/lib -lflint -lmpfr -lgmp -lboost_json -lboost_container \
  -o /tmp/diffexp-product-chart
python3 tools/performance/ft_product_experiment/run.py /tmp/diffexp-product-chart

c++ -std=c++20 -O3 -DNDEBUG -Iinclude -Itools -I/opt/homebrew/include \
  tools/performance/ft_product_experiment/cached_chart.cpp \
  -L/opt/homebrew/lib -lflint -lmpfr -lgmp -lboost_json -lboost_container \
  -o /tmp/diffexp-product-cached
python3 tools/performance/ft_product_experiment/run.py /tmp/diffexp-product-cached --cached
```

Run sequentially. Each numerical process has a 60-second timeout. Adapt compiler
include/library directories on other platforms; the runner builds nothing.
