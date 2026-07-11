# Feynman-trick examples

These shell examples invoke the implemented DiffExp 2 ladder in
``Scripts/run_ft_stepwise2.m``.  They require FIRE for a cold preparation and
the compiled recurrence backend.  Run them from any directory with ``sh``.

| Example | Command | Status in the documentation snapshot |
| --- | --- | --- |
| bubble | ``sh Examples/FeynmanTrick/Bubble.sh`` | small recommended first run |
| sunrise | ``sh Examples/FeynmanTrick/Sunrise.sh`` | completed parity fixture |
| unequal three-loop banana | ``sh Examples/FeynmanTrick/UnequalBanana.sh`` | completed C++/oracle comparison |
| box-bubble | ``sh Examples/FeynmanTrick/BoxBubble.sh`` | completed Euclidean ladder |
| massive kite | ``sh Examples/FeynmanTrick/Kite.sh`` | completed Euclidean ladder |
| unequal four-loop banana | ``sh Examples/FeynmanTrick/FourLoopUnequalBanana.sh`` | experimental; no source-controlled completed ladder result in this snapshot |

Every script:

- selects ``DE2_RECURRENCE_BACKEND=Cpp`` explicitly;
- accepts an existing ``DE2_CPP_THREADS`` value or uses four workers;
- keeps FIRE preparation under
  ``${DIFFEXP2_CACHE_DIR:-$HOME/.cache/diffexp2}/fire``;
- keeps resumable ladder checkpoints in an example-specific subdirectory;
- records all precision, expansion, epsilon-lookahead, division, and radius
  settings in the invocation.

Set a different cache root without editing a script:

```sh
DIFFEXP2_CACHE_DIR=/fast/local/cache \
DE2_CPP_THREADS=8 \
sh Examples/FeynmanTrick/UnequalBanana.sh
```

To force a new FIRE preparation, add ``FT_REBUILD_PREP=1`` to the environment.
That is normally unnecessary: solver-source changes do not invalidate the
separate prepared-data cache.

The exact propagators and Euclidean points are centralized in
``Scripts/FTExamples.m``, which the runner loads.  Their conventions and all
controls are documented in [Feynman Trick](../../Docs/FeynmanTrick.md).
