(* Check continuity of a dumped combined integrand across segment boundaries.
   For each pair of adjacent segments, evaluate both local series (with theta
   branch rules resolved on each side) at the shared main-variable point and
   compare per epsilon order.  Discontinuities localize transport or analytic
   continuation problems. *)

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
wp = 60;

Print["Call=", dump["CallIndex"], " epsMin=", epsMin, " segs=", Length[segs]];

xLocal = DiffExp`Symbols`x;

evalSeriesAt[ser_, localPt_, dir_] := Module[{rules, expr},
  rules = DiffExp`RegularizedIntegration`Private`thetaRulesAtPoint[localPt, dir];
  expr = If[MatchQ[ser, _SeriesData], Normal[ser], ser];
  expr = expr /. rules /. DiffExp`Symbols`Logx -> Log[xLocal];
  N[expr /. xLocal -> SetPrecision[localPt, 200], wp]
];

(* For each adjacent pair, evaluate at the shared main point *)
Do[
  Module[{segA = segs[[i]], segB = segs[[i + 1]], mainPt, locA, locB,
          serA, serB, dirA, dirB, valsA, valsB, diffs, relDiffs, nOrders},
    mainPt = DiffExp`RegularizedIntegration`Private`segmentActualBounds[segA][[2]];
    locA = DiffExp`RegularizedIntegration`Private`segmentLocalCoordinateForValue[segA, mainPt];
    locB = DiffExp`RegularizedIntegration`Private`segmentLocalCoordinateForValue[segB, mainPt];
    If[Min[Abs[N[locA]], Abs[N[locB]]] < 10^-10,
      Print["BOUNDARY ", i, "/", i + 1, " main=", N[mainPt, 8],
        " SKIPPED (singular expansion point on one side)"];
      Continue[];
    ];
    (* direction: approach the shared point from inside each segment *)
    dirA = If[TrueQ[Re[N[locA]] >= 0], 1, -1];
    dirB = If[TrueQ[Re[N[locB]] >= 0], 1, -1];
    serA = DiffExp`RegularizedIntegration`Private`uncompressSeriesData[segA[[5]]][[1]];
    serB = DiffExp`RegularizedIntegration`Private`uncompressSeriesData[segB[[5]]][[1]];
    nOrders = Min[Length[serA], Length[serB]];
    valsA = Table[evalSeriesAt[serA[[k]], locA, dirA], {k, nOrders}];
    valsB = Table[evalSeriesAt[serB[[k]], locB, dirB], {k, nOrders}];
    diffs = valsA - valsB;
    relDiffs = Table[
      Abs[diffs[[k]]] / Max[1, Abs[valsA[[k]]], Abs[valsB[[k]]]],
      {k, nOrders}
    ];
    Print["BOUNDARY ", i, "/", i + 1, " main=", N[mainPt, 8],
      " locA=", N[locA, 6], " locB=", N[locB, 6],
      " relDiffs=", InputForm[N[relDiffs, 3]]];
    If[Max[relDiffs] > 10^-15,
      Do[
        If[relDiffs[[k]] > 10^-15,
          Print["  ORDER eps^", epsMin + k - 1,
            " A=", InputForm[N[valsA[[k]], 25]],
            " B=", InputForm[N[valsB[[k]], 25]]];
        ],
        {k, nOrders}
      ];
    ];
  ],
  {i, Length[segs] - 1}
];

Quit[0];
