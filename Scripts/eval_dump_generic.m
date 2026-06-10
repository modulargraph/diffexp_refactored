(* Evaluate an isolated Laurent definite-integral dump and print per-segment results. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];

envOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]
];

Get[envOrDefault["DUMP_FILE", "/tmp/diffexp_banana_l1_dumps/laurent_integral_0006.m"]];
dump = DiffExp`RegularizedIntegration`Private`laurentIntegralDump;
withSegments = envOrDefault["PRINT_SEGMENTS", "1"] === "1";

Print["Call=", dump["CallIndex"], " epsMin=", dump["EpsMinPower"],
  " pref=", InputForm[dump["PrefactorSpec"]],
  " segs=", dump["RelevantSegmentCount"]];

res = DiffExp`RegularizedIntegration`DefiniteIntegralWithPrefactorLaurent[
  dump["SavedData"], dump["Bounds"], dump["PrefactorSpec"], dump["EpsMinPower"]
][[1]];
Print["TotalMinPower=", res["MinPower"]];
Do[
  Print["TOTAL eps^", res["MinPower"] + i - 1, " = ",
    InputForm[N[res["Coefficients"][[i]], 30]]],
  {i, Length[res["Coefficients"]]}
];

If[withSegments,
  Do[
    seg = dump["SavedData"]["SegmentData"][[i]];
    actual = DiffExp`RegularizedIntegration`Private`segmentActualBounds[seg];
    overlap = {Max[Min[actual], dump["Bounds"][[1]]], Min[Max[actual], dump["Bounds"][[2]]]};
    If[overlap[[1]] < overlap[[2]],
      sres = DiffExp`RegularizedIntegration`Private`IntegrateSegmentWithPrefactorLaurent[
        seg, overlap, 1, dump["EpsMinPower"],
        dump["PrefactorSpec"]["PowerAtLower"], dump["PrefactorSpec"]["PowerAtUpper"],
        dump["PrefactorSpec"]["RationalFactor"], dump["PrefactorSpec"]["Variable"],
        dump["Bounds"][[1]], dump["Bounds"][[2]]
      ];
      Print["SEG ", i, " bounds=", N[actual, 6], " min=", sres["MinPower"],
        " res=", InputForm[N[sres["Coefficients"], 20]]];
    ];
  , {i, Length[dump["SavedData"]["SegmentData"]]}];
];

Quit[0];
