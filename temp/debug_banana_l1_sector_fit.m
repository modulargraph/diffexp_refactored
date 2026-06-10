repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];
envOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]
];

Get[envOrDefault[
  "DUMP_FILE",
  "/tmp/diffexp_banana_l1_dumps/laurent_integral_0006.m"
]];

dump = DiffExp`RegularizedIntegration`Private`laurentIntegralDump;
segIndex = ToExpression[envOrDefault["SEG_INDEX", "1"]];
localPowerToInspect = ToExpression[envOrDefault["LOCAL_POWER", "0"]];
seg = dump["SavedData"]["SegmentData"][[segIndex]];
epsMinPower = dump["EpsMinPower"];
xLocal = DiffExp`Symbols`x;
tol = DiffExp`State`FEC[DiffExp`State`RationalizationTolerance];

zeroQ[z_] := TrueQ[PossibleZeroQ[z]] || TrueQ[Abs[N[z, 80]] < tol];
coeffAt[ser_SeriesData, lp_] := Module[{idx = lp * ser[[6]] - ser[[4]] + 1},
  If[IntegerQ[idx] && idx >= 1 && idx <= Length[ser[[3]]], ser[[3, idx]], 0]
];
coeffAt[expr_, lp_] := If[lp === 0, expr, 0];

seriesAtIndex = DiffExp`RegularizedIntegration`Private`uncompressSeriesData[seg[[5]]][[1]];
decomp = DiffExp`SingularityDecomposition`DecomposeSingularity[seriesAtIndex];
xMainExpr = DiffExp`RegularizedIntegration`Private`segmentMainExpression[seg] //
  Expand // DiffExp`Utilities`PChop;
jacobian = D[xMainExpr, xLocal] /. xLocal -> Mean[seg[[4]]];
localA = DiffExp`RegularizedIntegration`Private`segmentLocalCoordinateForValue[
  seg, 0, xMainExpr
];
localB = DiffExp`RegularizedIntegration`Private`segmentLocalCoordinateForValue[
  seg, 0.0480493049491587164576308496538967271824749629414144157378`40,
  xMainExpr
];
branchRules = DiffExp`RegularizedIntegration`Private`thetaRulesAtPoint[localA, -1];

term = decomp[[1]];
g = term["g"];
b = DiffExp`Utilities`PChop[Expand[term["b"] /. branchRules]];
coeffList = Table[coeffAt[g[[idx]], localPowerToInspect], {idx, Length[g]}];
resolved = DiffExp`Utilities`PChop[
  Expand[# /. branchRules]
] & /@ coeffList;
firstIdx = SelectFirst[
  Range[Length[resolved]],
  !DiffExp`RegularizedIntegration`Private`EffectiveZeroExprQ[resolved[[#]], tol] &,
  Missing["None"]
];
maxOffset = Length[resolved] - firstIdx;

logCoeffAt[offset_Integer, logPower_Integer] := If[
  offset < 0 || firstIdx + offset > Length[resolved],
  0,
  DiffExp`SeriesOps`LogxCoeffNS[resolved[[firstIdx + offset]], logPower]
];
moment[k_Integer] := DiffExp`Utilities`PChop[
  Expand[Factorial[k] * logCoeffAt[k, k]]
];

s0 = moment[0]; s1 = moment[1]; s2 = moment[2]; s3 = moment[3]; s4 = moment[4];
det = DiffExp`Utilities`PChop[Expand[s0*s2 - s1^2]];
If[maxOffset < 2 ||
    DiffExp`RegularizedIntegration`Private`EffectiveZeroExprQ[det, tol],
  roots = {If[
      DiffExp`RegularizedIntegration`Private`EffectiveZeroExprQ[s1, tol],
      0,
      DiffExp`Utilities`PChop[Expand[s1/s0]]
    ]};
  weights = {s0};
  e1 = e2 = disc = "one-sector";
  ,
  e1 = DiffExp`Utilities`PChop[Expand[(s0*s3 - s1*s2)/det]];
  e2 = DiffExp`Utilities`PChop[Expand[(s1*s3 - s2^2)/det]];
  disc = DiffExp`Utilities`PChop[Expand[e1^2 - 4*e2]];
  roots = (DiffExp`Utilities`PChop[Expand[#]] &) /@
    {(e1 + Sqrt[disc])/2, (e1 - Sqrt[disc])/2};
  weights = {
    DiffExp`Utilities`PChop[
      Expand[(s1 - s0*roots[[2]])/(roots[[1]] - roots[[2]])]
    ],
    0
  };
  weights[[2]] = DiffExp`Utilities`PChop[Expand[s0 - weights[[1]]]];
];

Print["call=", dump["CallIndex"], " segment=", segIndex,
  " localPower=", localPowerToInspect, " a=", term["a"]];
Print["branchB=", InputForm[N[b, 40]]];
Print["jacobian=", InputForm[N[jacobian, 40]]];
Print["firstIdx=", firstIdx, " epsFirst=", epsMinPower + firstIdx - 1,
  " maxOffset=", maxOffset];
Do[
  Print["resolved[", i, "] freeList=", FreeQ[resolved[[i]], List],
    " logs=", Table[
      {lp, InputForm[N[DiffExp`SeriesOps`LogxCoeffNS[resolved[[i]], lp], 30]]},
      {lp, DiffExp`SeriesOps`LogxPowerRange[resolved[[i]]]}
    ]
  ],
  {i, firstIdx, Min[firstIdx + 5, Length[resolved]]}
];
Print["moments freeList=", FreeQ[{s0, s1, s2, s3, s4}, List]];
Print["s=", InputForm[N[{s0, s1, s2, s3, s4}, 40]]];
Print["e=", InputForm[N[{det, e1, e2, disc}, 40]]];
Print["roots freeList=", FreeQ[roots, List]];
Print["roots=", InputForm[N[roots, 40]]];
Print["weights=", InputForm[N[weights, 40]]];
Print["leading pole candidate=", InputForm[N[jacobian Total[weights/(b + roots)], 40]]];

If[Length[roots] == 2,
  Do[
    Module[{rows = {}, kmax = maxOffset - q},
      Do[
        If[k + 1 <= kmax,
          Module[{m0, m1, c1, c2, denom},
            m0 = Factorial[k] * logCoeffAt[q + k, k];
            m1 = Factorial[k + 1] * logCoeffAt[q + k + 1, k + 1];
            denom = roots[[1]]^k roots[[2]]^(k + 1) -
              roots[[2]]^k roots[[1]]^(k + 1);
            c1 = DiffExp`Utilities`PChop[
              Expand[(m0 roots[[2]]^(k + 1) - m1 roots[[2]]^k) / denom]
            ];
            c2 = DiffExp`Utilities`PChop[
              Expand[(m1 roots[[1]]^k - m0 roots[[1]]^(k + 1)) / denom]
            ];
            AppendTo[rows, {k, N[{c1, c2}, 30]}];
          ]
        ],
        {k, 0, kmax - 1}
      ];
      Print["q ", q, " coeff rows=", InputForm[rows]];
    ],
    {q, 0, Min[3, maxOffset - 1]}
  ];
];

Quit[0];
