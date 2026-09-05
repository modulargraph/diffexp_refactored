Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "DiffExp.m"}]];
(* Automatic finds ../build/diffexp; alternatively set $DiffExpExecutable. *)
Print[DiffExpBackendInfo[]];
result = DiffExpSeries[<|
  "matrix" -> {{"eps/(1-x)"}},
  "boundary" -> {{"1", "0", "0"}},
  "center" -> "0", "taylor_order" -> 8,
  "epsilon_low" -> 0, "epsilon_high" -> 2
|>];
Print[result];
(* coefficients[[n+1, component, epsilon-epsilon_low+1]] are rational strings. *)
(* A recursion example, when a FIRE7 executable is available:
DiffExpFeynmanTrick["sunrise", {"--fire", "/path/to/FIRE7",
  "--cache", "/path/to/cache", "--epsilon-order", "0"}]
*)
