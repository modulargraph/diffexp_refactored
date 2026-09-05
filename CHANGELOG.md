# Changelog

## 2.1.1

- Add generic JSON-defined Feynman families, including editable configurations
  for every supplied example and arbitrary requested integral powers.
- Add native rational/canonical/algebraic transport and partial Frobenius
  boundary matching, exposed through the original Mathematica configuration,
  boundary, transport and plotting workflow using RunProcess.
- Preserve boundary precision, report preparation/numerical/wrapper timings,
  and retain explicit error estimates and omitted-tail status.
- Accept explicit algebraic basis prefactors for principal endpoint
  normalization after root continuation, including through the Mathematica wrapper.
- Use guarded grouped-dot rational charts with the scalar fallback, and share
  canonical letter expansions and source convolutions across matrix entries.
- Retain carried boundary uncertainty through local homogeneous maps where
  feedback makes direct interval recurrence poorly conditioned.
- Reuse nearby physical samples for single-coordinate Mathematica plotting.
- Include the original equal-banana reference data needed by a fresh checkout.


## 2.1.0

- Publish the standalone C++20 mathematical backend.
- Retain series recurrences, regular-singular matching, Feynman-trick recursion,
  native FIRE/FIRE7p integration and verified durable checkpoints.
- Preserve shared-boundary correlations across recursion levels; both four-loop
  banana examples pass independent finite-part comparisons.
- Add a small Mathematica `RunProcess` wrapper and JSON CLI output. Remove the
  former Mathematica/LibraryLink application from the current repository tree.
- Install as CMake package `DiffExp`, target `DiffExp::core`, namespace `diffexp`,
  executable `diffexp`, version 2.1.0.
- Keep all original example groups as opt-in data comparisons. Full Henn
  108-component FT reconstruction remains experimental and frozen.

This is an API break from 2.0. Existing mathematical cache formats remain readable. The old application remains available in Git history.
