# Mathematica interface

The interface consists of one small Wolfram Language file. It uses
[`RunProcess`](https://reference.wolfram.com/language/ref/RunProcess.html) with
an argument list and imports/exports
[`RawJSON`](https://reference.wolfram.com/language/ref/format/RawJSON.html).
The C++ executable owns all symbolic and numerical computation.

Load `DiffExp.m` from the repository or the installed
`share/diffexp/Mathematica/DiffExp.wl`. There is no paclet installation step.
Executable selection, in order, is the `"Executable"` option,
`$DiffExpExecutable`, the `DIFFEXP_EXECUTABLE` environment variable, the nearby
build/installed binary, then `diffexp` on the process search path.

| Function | Result |
| --- | --- |
| `DiffExpBackendInfo[]` | Backend information as an Association |
| `DiffExpSeries[request]` | Exact regular-series coefficients as a JSON Association |
| `DiffExpFeynmanTrick[family, {arguments}]` | Recursive integral coefficients and tail status |
| `DiffExpRun[{arguments}, input]` | Exit code, stdout and stderr; stdin defaults to empty |

All functions accept `"Executable" -> "/path/to/diffexp"`. A launch failure,
nonzero process exit or malformed JSON response returns `Failure`. Process
failures retain stdout and stderr for inspection. Generic `DiffExpRun` can also
call the original-example commands, `singular-endpoint`, and `prepare`.

`DiffExpSeries` accepts an Association matching `examples/logarithm.json`:
`matrix` is a square list of lists of rational-function strings in `x` and
`eps`; `boundary` is indexed by component then retained epsilon power;
`center` is a rational string. `epsilon_low`, `epsilon_high` and `taylor_order`
are integers. The center must be ordinary. Coefficients are indexed by Taylor
order, component, and epsilon power, each shifted to Mathematica's one-based
indexing. They remain exact rational strings. This function computes a finite
series; it does not certify an omitted infinite tail.

The FT response uses component/epsilon indexing and complex Arb interval
strings, for example `<|"real" -> "[1.234 +/- 0.001]", "imaginary" -> "0"|>`.
Keeping these strings avoids rounding through machine-precision JSON numbers.
Check `omitted_tails_certified` before treating a radius as a bound for the full
integral. Built-in FT family names are listed in [examples.md](examples.md).
Custom families and general transport are available through the C++ API; this
wrapper exposes the executable's commands rather than reproducing the 2.0 API.

Calls are synchronous. For long jobs, use the executable in a terminal to see
progress live. There is no persistent Mathematica session or background retry
loop; expensive work is reused through explicitly selected disk caches.

With a local Wolfram installation, enable the real wrapper regression with:

```sh
cmake -S . -B build -DBUILD_TESTING=ON -DDIFFEXP_WOLFRAMSCRIPT=/path/to/wolframscript
ctest --test-dir build -R mathematica.interface --output-on-failure
```

This is optional and is not required to build or test the native core.

If wolframscript selects a different installation, also set
`-DDIFFEXP_WOLFRAM_KERNEL=/path/to/WolframKernel` when configuring the test.
