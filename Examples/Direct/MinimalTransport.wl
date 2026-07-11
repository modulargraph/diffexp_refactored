(* DiffExp 2: minimal direct transport.

   Run from any directory with:
     wolframscript -file Examples/Direct/MinimalTransport.wl

   The release path uses the compiled recurrence backend. *)

repoRoot = ParentDirectory[
  ParentDirectory[DirectoryName[ExpandFileName[$InputFileName]]]
];

Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

SetAttributes[catchDiffExp2, HoldFirst];
catchDiffExp2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

configuration = catchDiffExp2[
  DiffExp2`Config`LoadConfiguration[{
    "RecurrenceBackend" -> "Cpp",
    "WorkingPrecision" -> 100,
    "ExpansionOrder" -> 30,
    "EpsilonOrder" -> 2,
    "DivisionOrder" -> 3,
    "Verbosity" -> 0
  }]
];

If[FailureQ[configuration],
  Print["Could not configure DiffExp 2: ", configuration];
  Exit[1]
];

x = Global`x;
eps = Global`eps;

(* f'(x) = f(x)/(x-2), f(0)=1/2, hence f(1)=1/4. *)
system = catchDiffExp2[
  DiffExp2`API`LoadSystem[<|
    "Matrix" -> {{1/(x - 2)}},
    "Variable" -> x
  |>]
];

boundary = {{1/2, 0, 0}};

result = If[FailureQ[system], system,
  catchDiffExp2[
    DiffExp2`API`TransportEndpoint[system, boundary, 0, 1]
  ]
];

If[FailureQ[result],
  Print["Transport failed: ", result];
  Exit[1]
];

value = result["Value"];
eps0 = DiffExp2`EpsSeries`ESCoefficient[value, 0];

Print["segments = ", result["SegmentCount"]];
Print["epsilon window = ", DiffExp2`EpsSeries`ESWindow[value]];
Print["value = ", DiffExp2`EpsSeries`ESToExpression[value, eps]];
Print["eps^0 vector = ", eps0];

If[!TrueQ[Abs[N[First[eps0] - 1/4, 50]] < 10^-30],
  Print["Unexpected result; expected {1/4}."];
  Exit[1]
];

Exit[0];
