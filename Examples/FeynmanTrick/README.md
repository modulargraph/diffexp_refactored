# Feynman-trick examples

These examples use the public `FeynmanTrick` pipeline facade. The shared
`RunExample.m` launcher constructs a validated `PipelinePlan` and executes it
with `RunIntegrationPipeline`; the shell files are only convenient names for
those public calls. They require FIRE for a cold preparation and the compiled
recurrence backend. Run the commands below from the repository root with
`sh`. The scripts may also be invoked by absolute path from another
directory.

| Example | Command | Status in the documentation snapshot |
| --- | --- | --- |
| bubble | ``sh Examples/FeynmanTrick/Bubble.sh`` | small recommended first run |
| sunrise | ``sh Examples/FeynmanTrick/Sunrise.sh`` | completed parity fixture |
| unequal three-loop banana | ``sh Examples/FeynmanTrick/UnequalBanana.sh`` | completed C++/oracle comparison |
| box-bubble | ``sh Examples/FeynmanTrick/BoxBubble.sh`` | completed Euclidean ladder |
| massive kite | ``sh Examples/FeynmanTrick/Kite.sh`` | completed Euclidean ladder |
| unequal four-loop banana | ``sh Examples/FeynmanTrick/FourLoopUnequalBanana.sh`` | experimental; no source-controlled completed ladder result in this snapshot |

The numerical profiles are centralized in `RunExample.m`. Every profile:

- selects the C++ recurrence backend;
- accepts an existing `DE2_CPP_THREADS` value or uses four workers;
- enables regular-chart value transport and native value-hop execution;
- requests paired C++ endpoint-arm prewarming;
- keeps FIRE preparation under
  `${DIFFEXP2_CACHE_DIR:-$HOME/.cache/diffexp2}/fire`;
- keeps resumable ladder checkpoints in an example-specific subdirectory;
- records all precision, expansion, epsilon-lookahead, division, and radius
  settings in the pipeline plan.

Inspect one exact plan without running it:

```sh
sh Examples/FeynmanTrick/Kite.sh --plan
```

Validate every documented profile:

```sh
wolframscript -script Examples/FeynmanTrick/RunExample.m --check
```

The generic launcher is also usable directly:

```sh
wolframscript -script Examples/FeynmanTrick/RunExample.m bubble
```

Set a different cache root without editing a script:

```sh
DIFFEXP2_CACHE_DIR=/fast/local/cache \
DE2_CPP_THREADS=8 \
sh Examples/FeynmanTrick/UnequalBanana.sh
```

Set `FT_FIRE_PATH=/absolute/path/to/fire/FIRE7` when FIRE is not installed at
`Dependencies/fire/FIRE7/FIRE7`; the wrappers preserve that environment
variable.

To force a new FIRE preparation, add ``FT_REBUILD_PREP=1`` to the environment.
That is normally unnecessary: solver-source changes do not invalidate the
separate prepared-data cache.

The exact propagators and Euclidean points are centralized in
`Scripts/FTExamples.m`, which the pipeline implementation loads. Their conventions and all
controls are documented in [Feynman Trick](../../Docs/FeynmanTrick.md).
