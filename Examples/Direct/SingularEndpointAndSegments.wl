(* DiffExp 2: inspect chart segments and exact x^(a+b eps) behavior.

   The equation f'(x)=eps f(x)/x with
   f(1/2)=(1/2)^eps has the exact solution x^eps.  The endpoint chart
   therefore exposes a=0, b=1 rather than fitting Log[x] towers.

   Run with:
     wolframscript -file Examples/Direct/SingularEndpointAndSegments.wl
*)

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
  DiffExp2`API`LoadSystem[<|
    "Matrix" -> {{eps/x}},
    "Variable" -> x
  |>]
];

boundary = Transpose@Table[
  {SeriesCoefficient[anchor^eps, {eps, 0, k}]},
  {k, 0, 3}
];

plan = If[FailureQ[system], system,
  catchDiffExp2[
    DiffExp2`Transport`SegmentLine[system, {anchor, 0}]
  ]
];

result = If[FailureQ[plan], plan,
  catchDiffExp2[
    DiffExp2`Transport`TransportLine[system, boundary, plan]
  ]
];

If[FailureQ[result],
  Print["Transport failed: ", result];
  Exit[1]
];

chartRows = Map[
  Function[chart, <|
    "Name" -> chart["Name"],
    "Center" -> chart["Center"],
    "Scale" -> chart["Scale"],
    "Radius" -> chart["Radius"],
    "IncomingMatchPoint" -> Lookup[chart, "IncomingMatchPoint", None],
    "Singular" -> chart["Singular"]
  |>],
  plan["Charts"]
];

Print["segment plan:"];
Scan[Print, chartRows];

decomposition = DiffExp2`SectorSeries`SectorDecomposition[result["Final"]];
tags = KeyTake[#, {"a", "b", "p"}] & /@ decomposition["Sectors"];

Print["exact endpoint sector tags = ", tags];
Print["dimreg endpoint limit of f = ",
  catchDiffExp2[DiffExp2`API`EndpointLimitValues[result, {1}]]
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
    DiffExp2`SectorSeries`EvaluateLocalSolution[local, t]
  ];
  If[FailureQ[evaluation], Return[Indeterminate, Module]];
  N[First[
    DiffExp2`EpsSeries`ESCoefficient[evaluation["Value"], 1]
  ], 30]
];

plotData = Select[
  Transpose[{N[samples, 20], eps1At /@ samples}],
  FreeQ[#, Indeterminate] &
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
Export[plotFile, plot];
Print["plot written to ", plotFile];

Exit[0];
