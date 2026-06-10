(* Profile one epsilon order of a dumped combined integrand at many points.
   Env: DUMP_FILE, POINTS (comma-separated rationals), EPS_ORDER (default -1). *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];

envOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]
];

Get[envOrDefault["DUMP_FILE", "/tmp/diffexp_banana_l0_dumps3/laurent_integral_0001.m"]];
dump = DiffExp`RegularizedIntegration`Private`laurentIntegralDump;
segs = dump["SavedData"]["SegmentData"];
epsMin = dump["EpsMinPower"];
targetOrder = ToExpression[envOrDefault["EPS_ORDER", "-1"]];
points = ToExpression /@ StringSplit[
  envOrDefault["POINTS", "1/2,11/20,3/5,13/20,7/10,3/4,4/5,17/20,9/10,19/20,39/40"], ","];

xLocal = DiffExp`Symbols`x;

Do[
  Module[{segIdx, seg, loc, ser, dir, rules, k, expr, val},
    segIdx = SelectFirst[
      Range[Length[segs]],
      Module[{b = DiffExp`RegularizedIntegration`Private`segmentActualBounds[segs[[#]]]},
        TrueQ[Min[b] <= pt <= Max[b]]
      ] &,
      Missing["None"]
    ];
    If[MissingQ[segIdx], Print["PT ", N[pt], " no segment"]; Continue[]];
    seg = segs[[segIdx]];
    loc = DiffExp`RegularizedIntegration`Private`segmentLocalCoordinateForValue[seg, pt];
    ser = DiffExp`RegularizedIntegration`Private`uncompressSeriesData[seg[[5]]][[1]];
    dir = If[TrueQ[Re[N[loc]] >= 0], 1, -1];
    rules = DiffExp`RegularizedIntegration`Private`thetaRulesAtPoint[loc, dir];
    k = targetOrder - epsMin + 1;
    If[k < 1 || k > Length[ser], Print["PT ", N[pt], " order out of window"]; Continue[]];
    expr = If[MatchQ[ser[[k]], _SeriesData], Normal[ser[[k]]], ser[[k]]];
    expr = expr /. rules /. DiffExp`Symbols`Logx -> Log[xLocal];
    val = N[expr /. xLocal -> SetPrecision[loc, 200], 30];
    Print["PT ", N[pt, 6], " seg=", segIdx, " local=", N[loc, 6],
      " eps^", targetOrder, " = ", InputForm[N[val, 12]]];
  ],
  {pt, points}
];
Quit[0];
