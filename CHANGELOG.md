# Changelog

## 2.1.0

- Publish the standalone C++20 rewrite previously developed as DiffExp3.
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

This is an API break from 2.0. Existing cache format identifiers retain their
working-name versions. The old application remains available in Git history.
