(* Fast reproduction of the failing xx1=1 local solve: transport the level-1
   system from xx1=0.96 (boundary values extracted from a verified regular
   segment) to xx1=1, with the strategy dispatch log enabled.  Saves the
   resulting transport in the same format as the stepwise save for reuse by
   check_transport_ode_residual.m. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["FeynmanTrick/FeynmanTrick.m"], {General::shdw, Symbol::shdw}];
Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];

envOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]
];

bcsData = Get[envOrDefault["BCS_FILE", "/tmp/ft_transport_save/bcs_at_096.m"]];
matrixDir = envOrDefault["MATRIX_DIR",
  Get["/tmp/ft_transport_save/transport_level_1.m"]["MatrixDir"]];
outFile = envOrDefault["OUT_FILE", "/tmp/ft_transport_save/repro_seg33.m"];

bcs = bcsData["BCs"];
pt = bcsData["Point"];
epsOrder = Length[bcs[[1]]] - 1;

Print["repro: transporting ", Length[bcs], " masters from xx1=", N[pt, 6],
  " to 1, epsOrder=", epsOrder];

DiffExp`State`StrategyDispatchLog = {};
DiffExp`State`$LogStrategyDispatch = True;

result = FeynmanTrick`DiffExpIntegration`TransportLevel[
  matrixDir, bcs, epsOrder,
  "FixedParamValue" -> pt,
  "LowerEndpoint" -> None,
  "UpperEndpoint" -> 1,
  "WorkingPrecision" -> 450,
  "ExpansionOrder" -> 50,
  "DivisionOrder" -> 4,
  "Verbosity" -> 1,
  "UseRationalRecurrence" -> True
];

If[!AssociationQ[result],
  Print["TRANSPORT FAILED: ", result];
  Quit[1];
];

(* Enable the dispatch log only for the singular solve trace: the log was
   enabled globally above; print a per-strategy summary now. *)
Print["segments: ", Length[result["SegmentData"]]];
Module[{log = DiffExp`State`StrategyDispatchLog},
  If[ListQ[log] && Length[log] > 0,
    Print["dispatch summary: ", Tally[{#["Label"], #["Strategy"]} & /@ log]];
    ,
    Print["dispatch log empty (enable LogStrategyDispatch)"];
  ];
];

Put[
  <|
    "Level" -> 1,
    "TransportResult" -> result,
    "TransportOrder" -> epsOrder,
    "MatrixDir" -> matrixDir
  |>,
  outFile
];
Print["saved -> ", outFile];
Quit[0];
