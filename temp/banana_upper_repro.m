(* Fast reproduction of the upper-line eps^-1 impurity: extract clean
   boundary values at x=0.48 from the saved level-1 transport, re-transport
   0.48 -> 1 with the strategy dispatch log and solver residual checks
   enabled, and compare the resulting towers at probe points against the
   saved transport.  Saves the repro transport for iteration. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["FeynmanTrick/FeynmanTrick.m"], {General::shdw, Symbol::shdw}];
Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];

envOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]
];

savedFile = envOrDefault["TRANSPORT_FILE",
  "/tmp/ft_transport_save5/transport_level_1.m"];
startPt = ToExpression[envOrDefault["START_POINT", "12/25"]];
probePts = ToExpression /@ StringSplit[
  envOrDefault["PROBE_POINTS", "7/10,9/10"], ","];
outFile = envOrDefault["OUT_FILE", "/tmp/banana_upper_repro_transport.m"];
wp = ToExpression[envOrDefault["FT_WORKING_PRECISION", "500"]];
expansionOrder = ToExpression[envOrDefault["FT_EXPANSION_ORDER", "50"]];
divisionOrder = ToExpression[envOrDefault["FT_DIVISION_ORDER", "4"]];

xLocal = DiffExp`Symbols`x;

towersAt[transport_Association, pt_] := Module[
  {segs = transport["SegmentData"], numMasters = transport["NumIntegrals"],
   segIdx, seg, loc, series, dir, rules},
  segIdx = SelectFirst[
    Range[Length[segs]],
    Module[{b = DiffExp`RegularizedIntegration`Private`segmentActualBounds[
        segs[[#]]]},
      TrueQ[Min[b] <= pt <= Max[b]]
    ] &
  ];
  seg = segs[[segIdx]];
  loc = DiffExp`RegularizedIntegration`Private`segmentLocalCoordinateForValue[
    seg, pt];
  series = DiffExp`RegularizedIntegration`Private`uncompressSeriesData[
    seg[[5]]];
  dir = If[TrueQ[Re[N[loc]] >= 0], 1, -1];
  rules = DiffExp`RegularizedIntegration`Private`thetaRulesAtPoint[loc, dir];
  Table[
    Module[{expr = If[MatchQ[series[[m, n]], _SeriesData],
        Normal[series[[m, n]]], series[[m, n]]]},
      expr = expr /. rules /. DiffExp`Symbols`Logx -> Log[xLocal];
      SetPrecision[expr /. xLocal -> SetPrecision[loc, wp - 20], wp - 50]
    ],
    {m, numMasters}, {n, Length[series[[1]]]}
  ]
];

Print["Loading saved transport: ", savedFile];
saved = Get[savedFile];
matrixDir = saved["MatrixDir"];
transport0 = saved["TransportResult"];
Print["MatrixDir: ", matrixDir];

bcs = towersAt[transport0, startPt];
Print["BCs at ", N[startPt, 6], " dims=", Dimensions[bcs]];

savedProbes = Table[towersAt[transport0, pt], {pt, probePts}];

Clear[saved, transport0];

DiffExp`State`StrategyDispatchLog = {};
DiffExp`State`$LogStrategyDispatch = True;
DiffExp`State`$DebugFuchsianizedCheck =
  envOrDefault["DEBUG_FUCHS_CHECK", "0"] === "1";
DiffExp`State`$DebugBlockResidualSeries =
  envOrDefault["DEBUG_RESID_SERIES", "0"] === "1";

epsOrder = Length[bcs[[1]]] - 1;
Print["Transporting ", Length[bcs], " masters from ", N[startPt, 6],
  " to 1 at epsOrder ", epsOrder];

result = FeynmanTrick`DiffExpIntegration`TransportLevel[
  matrixDir, bcs, epsOrder,
  "FixedParamValue" -> startPt,
  "LowerEndpoint" -> None,
  "UpperEndpoint" -> 1,
  "WorkingPrecision" -> wp,
  "ExpansionOrder" -> expansionOrder,
  "DivisionOrder" -> divisionOrder,
  "Verbosity" -> 0,
  "UseRationalRecurrence" -> True
];
If[!AssociationQ[result], Print["TRANSPORT FAILED: ", result]; Quit[1]];
Print["segments: ", Length[result["SegmentData"]]];

Module[{log = DiffExp`State`StrategyDispatchLog},
  If[ListQ[log] && Length[log] > 0,
    Print["dispatch tally: ",
      Tally[{#["Label"], #["Strategy"]} & /@ log] // InputForm];
  ];
];

Do[
  Module[{pt = probePts[[i]], reproT, diffs},
    reproT = towersAt[result, pt];
    diffs = Table[
      Max[Abs[N[reproT[[m]] - savedProbes[[i, m]], 30]]],
      {m, Length[reproT]}
    ];
    Print["PROBE ", N[pt, 6], " max tower diff per master: ",
      InputForm[N[diffs, 5]]];
  ],
  {i, Length[probePts]}
];

Put[<|"TransportResult" -> result, "MatrixDir" -> matrixDir,
  "StartPoint" -> startPt|>, outFile];
Print["saved repro -> ", outFile];
Quit[0];
