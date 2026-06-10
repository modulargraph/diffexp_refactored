(* Extract per-master boundary values at a regular point from a saved
   transport result (evaluating the containing segment's local series),
   for use as TransportTo boundary conditions in fast reproductions. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];

envOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]
];

saved = Get[envOrDefault["TRANSPORT_FILE", "/tmp/ft_transport_save/transport_level_1.m"]];
pt = ToExpression[envOrDefault["POINT", "96/100"]];
outFile = envOrDefault["OUT_FILE", "/tmp/ft_transport_save/bcs_at_096.m"];

transport = saved["TransportResult"];
segs = transport["SegmentData"];
numMasters = transport["NumIntegrals"];
xLocal = DiffExp`Symbols`x;

segIdx = SelectFirst[
  Range[Length[segs]],
  Module[{b = DiffExp`RegularizedIntegration`Private`segmentActualBounds[segs[[#]]]},
    TrueQ[Min[b] < pt < Max[b]]
  ] &
];
seg = segs[[segIdx]];
loc = DiffExp`RegularizedIntegration`Private`segmentLocalCoordinateForValue[seg, pt];
Print["point=", N[pt, 10], " segment=", segIdx, " local=", N[loc, 10]];

series = DiffExp`RegularizedIntegration`Private`uncompressSeriesData[seg[[5]]];
dir = If[TrueQ[Re[N[loc]] >= 0], 1, -1];
rules = DiffExp`RegularizedIntegration`Private`thetaRulesAtPoint[loc, dir];

bcs = Table[
  Module[{expr = If[MatchQ[series[[m, n]], _SeriesData],
      Normal[series[[m, n]]], series[[m, n]]], val},
    expr = expr /. rules /. DiffExp`Symbols`Logx -> Log[xLocal];
    val = expr /. xLocal -> SetPrecision[loc, 480];
    SetPrecision[val, 450]
  ],
  {m, numMasters}, {n, Length[series[[1]]]}
];

Print["bcs dims=", Dimensions[bcs], " sample=", N[bcs[[1, 1 ;; 3]], 15]];
Put[<|"Point" -> pt, "BCs" -> bcs|>, outFile];
Print["saved -> ", outFile];
Quit[0];
