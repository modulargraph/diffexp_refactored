(* DiffExp 2: inspect chart segments and exact x^(a+b eps) behavior.

   The equation f'(x)=eps f(x)/x with
   f(1/2)=(1/2)^eps has the exact solution x^eps.  The endpoint chart
   therefore exposes a=0, b=1 rather than fitting Log[x] towers.

   Run from the repository root with:
     wolframscript -file Examples/Direct/SingularEndpointAndSegments.wl
*)

repoRoot = ParentDirectory[
  ParentDirectory[DirectoryName[ExpandFileName[$InputFileName]]]
];

Get[FileNameJoin[{repoRoot, "DiffExp2.m"}]];

SetAttributes[catchDiffExp2, HoldFirst];
catchDiffExp2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

configuration = catchDiffExp2[
  DiffExp2`LoadConfiguration[{
    "RecurrenceBackend" -> "Cpp",
    "WorkingPrecision" -> 100,
    "ExpansionOrder" -> 40,
    "EpsilonOrder" -> 3,
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
anchor = 1/2;

system = catchDiffExp2[
  DiffExp2`LoadSystem[<|
    "Matrix" -> {{eps/x}},
    "Variable" -> x
  |>]
];

boundary = DiffExp2`PrepareBoundary[{anchor^eps}];

plan = If[FailureQ[system], system,
  catchDiffExp2[
    DiffExp2`PlanLine[system, {anchor, 0}]
  ]
];

result = If[FailureQ[plan], plan,
  catchDiffExp2[
    DiffExp2`TransportLine[system, boundary, plan]
  ]
];

If[FailureQ[result],
  Print["Transport failed: ", result];
  Exit[1]
];

segments = DiffExp2`LineSegments[result];
chartRows = Map[
  Function[segment, Module[{chart = segment["Chart"]}, <|
    "Name" -> chart["Name"],
    "Center" -> chart["Center"],
    "Scale" -> chart["Scale"],
    "Radius" -> chart["Radius"],
    "IncomingMatchPoint" -> Lookup[chart, "IncomingMatchPoint", None],
    "Singular" -> chart["Singular"]
  |>]],
  segments
];

Print["segment plan:"];
Scan[Print, chartRows];

tags = KeyTake[#, {"a", "b", "p"}] & /@
  DiffExp2`ExactSectors[result];
endpointLimit = catchDiffExp2[DiffExp2`EndpointLimit[result, {1}]];

Print["exact endpoint sector tags = ", tags];
Print["dimreg endpoint limit of f = ", endpointLimit];

If[!MemberQ[tags, <|"a" -> 0, "b" -> 1, "p" -> 0|>],
  Print["Unexpected endpoint sectors; expected {a,b,p}={0,1,0}."];
  Exit[1]
];

If[FailureQ[endpointLimit] ||
    !TrueQ[PossibleZeroQ[DiffExp2`EpsilonCoefficient[endpointLimit, 0]]],
  Print["Unexpected dimensional-regularization endpoint limit; expected 0."];
  Exit[1]
];

(* Plot the eps^1 coefficient.  For x^eps it is Log[x].  The endpoint itself
   is not evaluated: multivalued local solutions use the endpoint-limit API. *)
local = result["Final"];
scale = local["ChartMap", "Scale"];
physicalRadius = N[Abs[scale] local["Radius"], 80];
xmin = N[10^-6, 80];
xmax = Min[N[1/4, 80], 9 physicalRadius/20];
samples = Exp[Subdivide[Log[xmin], Log[xmax], 80]];

eps1At[xx_] := Module[{evaluation, t},
  t = (xx - local["Center"])/scale;
  evaluation = catchDiffExp2[
    DiffExp2`EvaluateLocal[local, t]
  ];
  If[FailureQ[evaluation], Return[Indeterminate, Module]];
  N[First[
    DiffExp2`EpsilonCoefficient[evaluation, 1]
  ], 30]
];

plotData = Select[
  Transpose[{N[samples, 30], eps1At /@ samples}],
  FreeQ[#, Indeterminate] &
];

If[plotData === {},
  Print["No valid local-solution samples were produced."];
  Exit[1]
];

eps1Error = Max[Abs[N[#[[2]] - Log[#[[1]]], 20]] & /@ plotData];
If[!TrueQ[eps1Error < 10^-18],
  Print["Unexpected eps^1 coefficient; max error against Log[x] is ",
    eps1Error];
  Exit[1]
];

plot = ListLinePlot[
  plotData,
  Frame -> True,
  FrameLabel -> {"x", "coefficient of eps"},
  PlotLabel -> "DiffExp 2 local sector: x^eps = 1 + eps Log[x] + ...",
  PlotRange -> All,
  ImageSize -> Large
];

plotFile = FileNameJoin[{
  $TemporaryDirectory,
  "diffexp2-singular-endpoint-eps1.png"
}];
If[FileExistsQ[plotFile], DeleteFile[plotFile]];
exported = Quiet[Check[Export[plotFile, plot], $Failed]];
If[exported === $Failed || !FileExistsQ[plotFile],
  Print["Plot export failed: ", plotFile];
  Exit[1]
];
Print["plot written to ", plotFile];

Exit[0];
