(* Evaluate a dumped combined integrand at a chosen main-variable point
   (default: the FT fixed value 11/23) and print the epsilon tower.
   For the banana level-0 dump this should reproduce the pySecDec value
   of the needed integrand at the fixed point: {2,0,1,1} -> 5.4025802965
   at the epsilon order carrying the finite term. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];

envOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]
];

Get[envOrDefault["DUMP_FILE", "/tmp/diffexp_banana_l0_dumps/laurent_integral_0009.m"]];
dump = DiffExp`RegularizedIntegration`Private`laurentIntegralDump;
segs = dump["SavedData"]["SegmentData"];
epsMin = dump["EpsMinPower"];
anchor = ToExpression[envOrDefault["ANCHOR", "11/23"]];
wp = 50;

Print["Call=", dump["CallIndex"], " epsMin=", epsMin, " segs=", Length[segs]];

xLocal = DiffExp`Symbols`x;

(* Find the segment containing the anchor strictly inside *)
segIdx = SelectFirst[
  Range[Length[segs]],
  Module[{b = DiffExp`RegularizedIntegration`Private`segmentActualBounds[segs[[#]]]},
    TrueQ[Min[b] < anchor < Max[b]]
  ] &,
  Missing["None"]
];
If[MissingQ[segIdx],
  segIdx = SelectFirst[
    Range[Length[segs]],
    Module[{b = DiffExp`RegularizedIntegration`Private`segmentActualBounds[segs[[#]]]},
      TrueQ[Min[b] <= anchor <= Max[b]]
    ] &,
    Missing["None"]
  ];
];
If[MissingQ[segIdx], Print["No segment contains the anchor"]; Quit[1]];

seg = segs[[segIdx]];
loc = DiffExp`RegularizedIntegration`Private`segmentLocalCoordinateForValue[seg, anchor];
ser = DiffExp`RegularizedIntegration`Private`uncompressSeriesData[seg[[5]]][[1]];
dir = If[TrueQ[Re[N[loc]] >= 0], 1, -1];
rules = DiffExp`RegularizedIntegration`Private`thetaRulesAtPoint[loc, dir];

Print["segment=", segIdx, " local=", N[loc, 10]];
Do[
  Module[{expr = If[MatchQ[ser[[k]], _SeriesData], Normal[ser[[k]]], ser[[k]]], val},
    expr = expr /. rules /. DiffExp`Symbols`Logx -> Log[xLocal];
    val = N[expr /. xLocal -> SetPrecision[loc, 200], wp];
    Print["ANCHOR eps^", epsMin + k - 1, " = ", InputForm[N[val, 25]]];
  ],
  {k, Length[ser]}
];

Quit[0];
