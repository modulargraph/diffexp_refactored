# IBP Solver core

Unmodified core headers from https://github.com/modulargraph/ibp-solver, commit `8d803554c33ac691877ef72c3dd338929c862abc`. Licensed under GPL-3.0; see LICENSE.

This pinned source snapshot lets DiffExp build offline. The integration uses the public `InputGeometry`, `ParametricProgram`, `ArithmeticTrace`, `Field` and `Solver` interfaces. Update all five headers together from an audited upstream commit. The adapter compiles parameterized equation structure once, learns a guarded CPU arithmetic trace per reconstruction prime, and independently runs full reductions at validation points. GPU batch replay remains in the standalone package.
