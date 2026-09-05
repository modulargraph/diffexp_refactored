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
family = DiffExpFamilyTemplate["sunrise"];
DiffExpFeynmanTrick[family, {"--fire", "/path/to/FIRE7",
  "--cache", "/path/to/cache", "--epsilon-order", "0"}]
*)

(* Original configuration/boundary/transport workflow for y' = y/(1-t). *)
folder = CreateDirectory[];
Put[{{1/(1-t)}}, FileNameJoin[{folder, "dt_0.m"}]];
LoadConfiguration[{MatrixDirectory -> folder, EpsilonOrder -> 0,
  ExpansionOrder -> 40, WorkingPrecision -> 80, AccuracyGoal -> 20}];
boundary = PrepareBoundaryConditions[{1}, {t -> 0}];
saved = TransportTo[boundary, {t -> x}, 1/2, True];
Print[saved[[1, 2]]]; (* {{2}} *)
Print[DiffExpLastTimings[]];
functions = ToPiecewise[saved];
Print[functions[[1, 1]][1/4]]; (* 4/3 *)
DeleteDirectory[folder, DeleteContents -> True];
