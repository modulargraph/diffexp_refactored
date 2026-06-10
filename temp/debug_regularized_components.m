repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];

envOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]
];

Get[envOrDefault[
  "DUMP_FILE",
  "/tmp/diffexp_banana_l1_dumps/laurent_integral_0005.m"
]];

dump = DiffExp`RegularizedIntegration`Private`laurentIntegralDump;
segIndex = ToExpression[envOrDefault["SEG_INDEX", "1"]];
seg = dump["SavedData"]["SegmentData"][[segIndex]];
epsMinPower = dump["EpsMinPower"];
pref = dump["PrefactorSpec"];
xLocal = DiffExp`Symbols`x;
tol = DiffExp`State`FEC[DiffExp`State`RationalizationTolerance];

zeroQ[z_] := TrueQ[PossibleZeroQ[z]] || TrueQ[Abs[N[z, 80]] < tol];
coeffAt[ser_SeriesData, lp_] := Module[{idx = lp * ser[[6]] - ser[[4]] + 1},
  If[IntegerQ[idx] && idx >= 1 && idx <= Length[ser[[3]]], ser[[3, idx]], 0]
];
coeffAt[expr_, lp_] := If[lp === 0, expr, 0];

segmentActualBounds[segment_] :=
  DiffExp`RegularizedIntegration`Private`segmentActualBounds[segment];
segmentLocalCoordinateForValue[segment_, value_, expr_] :=
  DiffExp`RegularizedIntegration`Private`segmentLocalCoordinateForValue[
    segment, value, expr
  ];

actual = segmentActualBounds[seg];
overlap = {
  Max[Min[actual], dump["Bounds"][[1]]],
  Min[Max[actual], dump["Bounds"][[2]]]
};
xMainExpr = DiffExp`Utilities`PChop[
  Expand[DiffExp`RegularizedIntegration`Private`segmentMainExpression[seg]]
];
localA = segmentLocalCoordinateForValue[seg, overlap[[1]], xMainExpr];
localB = segmentLocalCoordinateForValue[seg, overlap[[2]], xMainExpr];
jacobian = D[xMainExpr, xLocal] /. xLocal -> (localA + localB)/2;
seriesAtIndex = DiffExp`RegularizedIntegration`Private`uncompressSeriesData[
  seg[[5]]
][[1]];
decomp = DiffExp`SingularityDecomposition`DecomposeSingularity[seriesAtIndex];

atLowerBound = zeroQ[overlap[[1]] - dump["Bounds"][[1]]];
atUpperBound = zeroQ[overlap[[2]] - dump["Bounds"][[2]]];
singularPowerLower = If[atLowerBound && zeroQ[localA],
  pref["PowerAtLower"],
  0
];
singularPowerUpper = If[atUpperBound && zeroQ[localB],
  pref["PowerAtUpper"],
  0
];
totalSingularPower = singularPowerLower + singularPowerUpper;
prefactorRational = pref["RationalFactor"] /. pref["Variable"] -> xMainExpr;
prefactorLower = If[singularPowerLower =!= 0,
  Abs[jacobian]^pref["PowerAtLower"] *
    DiffExp`RegularizedIntegration`Private`localSidePhase[
      localB - localA, pref["PowerAtLower"]
    ],
  (xMainExpr - dump["Bounds"][[1]])^pref["PowerAtLower"]
];
prefactorUpper = If[singularPowerUpper =!= 0,
  Abs[jacobian]^pref["PowerAtUpper"] *
    DiffExp`RegularizedIntegration`Private`localSidePhase[
      localA - localB, pref["PowerAtUpper"]
    ],
  (dump["Bounds"][[2]] - xMainExpr)^pref["PowerAtUpper"]
];
smoothPrefactor = prefactorLower * prefactorUpper * prefactorRational;
smoothPrefactorIsOne = TrueQ[
  PossibleZeroQ[DiffExp`Utilities`PChop[Expand[smoothPrefactor - 1]]]
];

modified = Table[
  Module[{termA = term["a"], termB = term["b"], termG = term["g"],
      newA, newG, gSeries, smoothSer, expansionOrder = 50},
    newA = termA + totalSingularPower;
    expansionOrder = If[MatchQ[termG[[1]], _SeriesData],
      termG[[1]][[5]] - termG[[1]][[4]],
      50
    ];
    newG = Table[
      If[ord <= Length[termG],
        gSeries = termG[[ord]];
        If[MatchQ[gSeries, _SeriesData],
          If[smoothPrefactorIsOne,
            gSeries,
            smoothSer = Series[
              SetPrecision[
                smoothPrefactor,
                DiffExp`RegularizedIntegration`Private`activeNumericPrecision[]
              ],
              {xLocal, gSeries[[2]], gSeries[[5]] - gSeries[[4]]}
            ];
            DiffExp`SeriesOps`SExpand[gSeries * smoothSer]
          ],
          gSeries * (smoothPrefactor /. xLocal -> 0)
        ],
        0
      ],
      {ord, Length[termG]}
    ];
    Module[{gMinPowers, minPower},
      gMinPowers = Table[
        If[MatchQ[newG[[ord]], _SeriesData],
          newG[[ord]][[4]]/newG[[ord]][[6]],
          0
        ],
        {ord, Length[newG]}
      ];
      minPower = Min[gMinPowers];
      If[TrueQ[minPower < 0],
        newA = newA + minPower;
        newG = Table[
          If[MatchQ[newG[[ord]], _SeriesData],
            Module[{s = newG[[ord]], shiftIndex},
              shiftIndex = minPower * newG[[ord]][[6]];
              SeriesData[s[[1]], s[[2]], s[[3]],
                s[[4]] - shiftIndex, s[[5]] - shiftIndex, s[[6]]]
            ],
            newG[[ord]] * xLocal^(-minPower)
          ],
          {ord, Length[newG]}
        ];
      ];
    ];
    <|"a" -> newA, "b" -> termB, "g" -> newG|>
  ],
  {term, decomp}
];

Print["call=", dump["CallIndex"], " seg=", segIndex, " overlap=", N[overlap, 20],
  " local=", N[{localA, localB}, 20], " jac=", N[jacobian, 30],
  " totalSingularPower=", totalSingularPower,
  " smoothIsOne=", smoothPrefactorIsOne];

Do[
  term = modified[[termIdx]];
  a = term["a"]; b = term["b"]; g = term["g"];
  branchRules = DiffExp`RegularizedIntegration`Private`branchRulesForLocalBounds[
    localA, localB
  ];
  branchB = DiffExp`Utilities`PChop[Expand[b /. branchRules]];
  subtractPowerQ[lp_] :=
    DiffExp`RegularizedIntegration`Private`NumericNegativeQ[a + lp + 1, tol] ||
      DiffExp`RegularizedIntegration`Private`NumericZeroQ[a + lp + 1, tol];
  subtractSeries[ser_SeriesData] := Module[{nmin = ser[[4]], den = ser[[6]], coeffs},
    coeffs = Table[
      If[subtractPowerQ[(nmin + idx - 1)/den], 0, ser[[3, idx]]],
      {idx, Length[ser[[3]]]}
    ];
    SeriesData[ser[[1]], ser[[2]], coeffs, ser[[4]], ser[[5]], ser[[6]]]
  ];
  subtractSeries[expr_] := If[subtractPowerQ[0], 0, expr];
  remainder = subtractSeries /@ g;
  full = DiffExp`RegularizedIntegration`Private`IntegrateAnalyticRegularizedBySubtractionLaurent[
    a, b, epsMinPower, g, {localA, localB},
    If[zeroQ[localA], Sign[N[localB, 50]], False],
    If[zeroQ[localB], Sign[N[localA, 50]], False]
  ];
  base = <|"MinPower" -> epsMinPower, "Coefficients" -> Table[
    DiffExp`RegularizedIntegration`Private`IntegratePowerTimesSeriesAtPower[
      a, branchB, epsMinPower, remainder, targetPower,
      {localA, localB},
      If[zeroQ[localA], Sign[N[localB, 50]], False],
      If[zeroQ[localB], Sign[N[localA, 50]], False]
    ],
    {targetPower, epsMinPower, epsMinPower + Length[g] - 1}
  ]|>;
  powers = DeleteDuplicates[Flatten[Table[
    If[MatchQ[g[[idx]], _SeriesData],
      Table[(g[[idx, 4]] + j - 1)/g[[idx, 6]], {j, Length[g[[idx, 3]]]}],
      {0}
    ],
    {idx, Length[g]}
  ]]];
  subPowers = Select[powers, subtractPowerQ];
  Print["term ", termIdx, " a=", a, " b=", N[b, 30], " branchB=", N[branchB, 30],
    " subPowers=", subPowers];
  Print["  full=", N[DiffExp`RegularizedIntegration`Private`LaurentScale[jacobian, full], 30]];
  Print["  base=", N[DiffExp`RegularizedIntegration`Private`LaurentScale[jacobian, base], 30]];
  Do[
    Print["  lp=", lp, " coeffs=", Table[
      Module[{c = coeffAt[g[[idx]], lp]},
        If[DiffExp`RegularizedIntegration`Private`EffectiveZeroExprQ[c, tol],
          Nothing,
          {epsMinPower + idx - 1, N[c, 20]}
        ]
      ],
      {idx, Length[g]}
    ]];
    onlyG = Table[
      Module[{c = coeffAt[g[[idx]], lp], den = 1, nmin},
        If[MatchQ[g[[idx]], _SeriesData],
          den = g[[idx, 6]];
          nmin = lp * den;
          SeriesData[g[[idx, 1]], g[[idx, 2]], {c}, nmin, nmin + 1, den],
          c
        ]
      ],
      {idx, Length[g]}
    ];
    onlyFull = DiffExp`RegularizedIntegration`Private`IntegrateAnalyticRegularizedBySubtractionLaurent[
      a, b, epsMinPower, onlyG, {localA, localB},
      If[zeroQ[localA], Sign[N[localB, 50]], False],
      If[zeroQ[localB], Sign[N[localA, 50]], False]
    ];
    Print["    onlyFull=", N[
      DiffExp`RegularizedIntegration`Private`LaurentScale[jacobian, onlyFull],
      30
    ]],
    {lp, subPowers}
  ],
  {termIdx, Length[modified]}
];

Quit[0];
