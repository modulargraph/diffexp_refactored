(* Independent check of one segment's contribution: numerically integrate the
   dumped local series per epsilon order (Logx resolved with the saved theta
   branch data) and compare with IntegrateSegmentWithPrefactorLaurent. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];

envOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]
];

Get[envOrDefault["DUMP_FILE", "/tmp/diffexp_banana_l0_dumps/laurent_integral_0009.m"]];
dump = DiffExp`RegularizedIntegration`Private`laurentIntegralDump;
segIndex = ToExpression[envOrDefault["SEG_INDEX", "33"]];
seg = dump["SavedData"]["SegmentData"][[segIndex]];
epsMin = dump["EpsMinPower"];
xLocal = DiffExp`Symbols`x;

bounds = DiffExp`RegularizedIntegration`Private`segmentActualBounds[seg];
overlap = {Max[Min[bounds], dump["Bounds"][[1]]], Min[Max[bounds], dump["Bounds"][[2]]]};
localA = DiffExp`RegularizedIntegration`Private`segmentLocalCoordinateForValue[seg, overlap[[1]]];
localB = DiffExp`RegularizedIntegration`Private`segmentLocalCoordinateForValue[seg, overlap[[2]]];
xMainExpr = DiffExp`Utilities`PChop[Expand[
  DiffExp`RegularizedIntegration`Private`segmentMainExpression[seg]]];
jac = D[xMainExpr, xLocal] /. xLocal -> (localA + localB)/2;

Print["segment=", segIndex, " local: ", N[localA, 10], " -> ", N[localB, 10],
  " jacobian=", N[jac, 10]];

ser = DiffExp`RegularizedIntegration`Private`uncompressSeriesData[seg[[5]]][[1]];

(* code result *)
codeRes = DiffExp`RegularizedIntegration`Private`IntegrateSegmentWithPrefactorLaurent[
  seg, overlap, 1, epsMin,
  dump["PrefactorSpec"]["PowerAtLower"], dump["PrefactorSpec"]["PowerAtUpper"],
  dump["PrefactorSpec"]["RationalFactor"], dump["PrefactorSpec"]["Variable"],
  dump["Bounds"][[1]], dump["Bounds"][[2]]
];

(* numeric reference: jac * NIntegrate of the branch-resolved series order by
   order over the local interval *)
sideDir = If[TrueQ[Re[N[(localA + localB)/2]] < 0], -1, 1];
rules = DiffExp`RegularizedIntegration`Private`thetaRulesAtPoint[(localA + localB)/2, sideDir];

numRes = Table[
  Module[{expr, f, val},
    expr = If[MatchQ[ser[[k]], _SeriesData], Normal[ser[[k]]], ser[[k]]];
    expr = expr /. rules /. DiffExp`Symbols`Logx -> Log[xLocal];
    If[PossibleZeroQ[expr],
      0,
      f = Function @@ {xLocal, expr};
      val = Quiet[NIntegrate[f[t] /. t -> xLocal, {xLocal, localA, localB},
        WorkingPrecision -> 60, MaxRecursion -> 14, AccuracyGoal -> 25,
        PrecisionGoal -> 25]];
      jac * val
    ]
  ],
  {k, Length[ser]}
];

Print["per eps order: code vs numeric"];
Do[
  Module[{cv = If[k - codeRes["MinPower"] + epsMin >= 0 &&
        k - codeRes["MinPower"] + epsMin + 1 <= Length[codeRes["Coefficients"]],
      codeRes["Coefficients"][[epsMin + k - codeRes["MinPower"]]],
      Missing[]
    ]},
    0
  ];
  ,
  {k, 1, 1}
];
(* simpler: print both Laurent lists *)
Print["CODE min=", codeRes["MinPower"]];
Do[
  Print["  CODE eps^", codeRes["MinPower"] + i - 1, " = ",
    InputForm[N[codeRes["Coefficients"][[i]], 22]]],
  {i, Length[codeRes["Coefficients"]]}
];
Print["NUM (termwise, no endpoint resummation) min=", epsMin];
Do[
  Print["  NUM  eps^", epsMin + i - 1, " = ", InputForm[N[numRes[[i]], 22]]],
  {i, Length[numRes]}
];

Quit[0];
