(* Focused tests for singular regularized integration edge cases. *)

SetDirectory[DirectoryName[$InputFileName]];
$Path = Prepend[$Path, FileNameJoin[{ParentDirectory[Directory[]], "DiffExp"}]];

Get["DiffExp.m"];

DiffExp`State`DiffExpConfiguration[RationalizationTolerance] = 10^-30;
DiffExp`State`DiffExpConfiguration[ChopPrecision] = 80;
DiffExp`State`DiffExpConfiguration[WorkingPrecision] = 100;

testsTotal = 0;
testsPassed = 0;

pass[name_] := (testsPassed++; Print["  [PASS] ", name]);
fail[name_, got_] := Print["  [FAIL] ", name, ": ", got];
closeQ[a_, b_, tol_:10^-40] := TrueQ[Abs[N[a - b, 80]] < tol];
laurentCoeff[laur_, power_Integer] := Module[
  {idx = power - laur["MinPower"] + 1},
  If[idx >= 1 && idx <= Length[laur["Coefficients"]],
    laur["Coefficients"][[idx]],
    0
  ]
];
laurentMaxPower[laur_] := laur["MinPower"] + Length[laur["Coefficients"]] - 1;
expectedLaurentCoeffs[expr_, var_, min_Integer, max_Integer] :=
  Table[SeriesCoefficient[expr, {var, 0, power}], {power, min, max}];
laurentCloseQ[laur_, min_Integer, coeffs_List, tol_:10^-35] :=
  And @@ Table[
    closeQ[laurentCoeff[laur, min + idx - 1], coeffs[[idx]], tol],
    {idx, Length[coeffs]}
  ];
formulaAgreementQ[a_, b_, epsMin_Integer, gList_List, maxPower_:Automatic] :=
  Module[{old, new, min, max},
    old = DiffExp`RegularizedIntegration`Private`LaurentTrim[
      DiffExp`RegularizedIntegration`Private`IntegrateAnalyticRegularizedByIBPLaurent[
        a, b, epsMin, gList, {0, 1}, 1, False
      ]
    ];
    new = DiffExp`RegularizedIntegration`Private`LaurentTrim[
      DiffExp`RegularizedIntegration`IntegrateSingularTermLaurent[
        a, b, epsMin, gList, {0, 1}
      ]
    ];
    min = Min[old["MinPower"], new["MinPower"]];
    max = If[maxPower === Automatic,
      Max[laurentMaxPower[old], laurentMaxPower[new]],
      maxPower
    ];
    And @@ Table[
      closeQ[laurentCoeff[old, power], laurentCoeff[new, power]],
      {power, min, max}
    ]
  ];

x = DiffExp`Symbols`x;
logx = DiffExp`Symbols`Logx;
ee = Unique["eps"];

Print["==========================================="];
Print["Regularized Integration Edge Case Tests"];
Print["==========================================="];

testsTotal++;
pureLogSeries = SeriesData[x, 0, {logx}, 0, 1, 1];
pureLogDecomp = DiffExp`SingularityDecomposition`DecomposeSingularity[{pureLogSeries}];
If[Length[pureLogDecomp] == 1,
  pass["pure Logx series is not dropped by decomposition"],
  fail["pure Logx series is not dropped by decomposition", pureLogDecomp]
];

testsTotal++;
directNonzeroInterval =
  DiffExp`RegularizedIntegration`IntegrateSingularTermLaurent[
    -1, 1, 0, {1, 0}, {1, 2}
  ];
If[directNonzeroInterval["MinPower"] == 0 &&
    Length[directNonzeroInterval["Coefficients"]] >= 1 &&
    closeQ[directNonzeroInterval["Coefficients"][[1]], Log[2]],
  pass["no analytic-regularization pole away from x=0"],
  fail["no analytic-regularization pole away from x=0", directNonzeroInterval]
];

testsTotal++;
topOrderPole =
  DiffExp`RegularizedIntegration`IntegrateSingularTermLaurent[
    -1, 1, 0, {1}, {0, -1}
  ];
If[topOrderPole["MinPower"] == -1 &&
    Length[topOrderPole["Coefficients"]] >= 2 &&
    closeQ[topOrderPole["Coefficients"][[1]], 1] &&
    closeQ[topOrderPole["Coefficients"][[2]], I Pi],
  pass["a=-1 regularization preserves finite part with one input order"],
  fail["a=-1 regularization preserves finite part with one input order", topOrderPole]
];

testsTotal++;
explicitLogIntegral =
  DiffExp`RegularizedIntegration`IntegrateSingularTermLaurent[
    0, 0, 0, {pureLogSeries}, {0, 1}
  ];
If[explicitLogIntegral["MinPower"] == 0 &&
    Length[explicitLogIntegral["Coefficients"]] >= 1 &&
    closeQ[explicitLogIntegral["Coefficients"][[1]], -1],
  pass["explicit Logx in Taylor coefficient integrates by parts correctly"],
  fail["explicit Logx in Taylor coefficient integrates by parts correctly", explicitLogIntegral]
];

testsTotal++;
subtractionAminus1 =
  DiffExp`RegularizedIntegration`IntegrateSingularTermLaurent[
    -1, 2, 0, {Series[3 + 5*x + 7*x^2, {x, 0, 5}], 0, 0}, {0, 1}
  ];
subtractionAminus1Expected = expectedLaurentCoeffs[
  3/(2*ee) + 5/(1 + 2*ee) + 7/(2 + 2*ee),
  ee, -1, 2
];
If[laurentCloseQ[subtractionAminus1, -1, subtractionAminus1Expected],
  pass["Taylor subtraction a=-1 uses g(0) pole and finite remainder"],
  fail["Taylor subtraction a=-1 uses g(0) pole and finite remainder",
    {subtractionAminus1, subtractionAminus1Expected}]
];

testsTotal++;
subtractionAminus2 =
  DiffExp`RegularizedIntegration`IntegrateSingularTermLaurent[
    -2, 2, 0, {Series[3 + 5*x + 7*x^2, {x, 0, 5}], 0, 0}, {0, 1}
  ];
subtractionAminus2Expected = expectedLaurentCoeffs[
  3/(-1 + 2*ee) + 5/(2*ee) + 7/(1 + 2*ee),
  ee, -1, 2
];
If[laurentCloseQ[subtractionAminus2, -1, subtractionAminus2Expected],
  pass["Taylor subtraction a=-2 uses the linear Taylor coefficient pole"],
  fail["Taylor subtraction a=-2 uses the linear Taylor coefficient pole",
    {subtractionAminus2, subtractionAminus2Expected}]
];

testsTotal++;
splitPoleCancellation =
  DiffExp`RegularizedIntegration`IntegrateDecompositionLaurent[
    {
      <|"a" -> -1, "b" -> 1, "g" -> {1, 7}|>,
      <|"a" -> -1, "b" -> 1, "g" -> {-1, -2}|>
    },
    {0, 1},
    0
  ];
wholePoleCancellation =
  DiffExp`RegularizedIntegration`IntegrateDecompositionLaurent[
    {<|"a" -> -1, "b" -> 1, "g" -> {0, 5}|>},
    {0, 1},
    0
  ];
If[laurentCloseQ[splitPoleCancellation, 0, {5}] &&
    laurentCloseQ[wholePoleCancellation, 0, {5}],
  pass["Taylor subtraction poles cancel consistently across split terms"],
  fail["Taylor subtraction poles cancel consistently across split terms",
    {splitPoleCancellation, wholePoleCancellation}]
];

testsTotal++;
splitLogPoleCancellation =
  DiffExp`RegularizedIntegration`IntegrateDecompositionLaurent[
    {
      <|"a" -> -1, "b" -> 1,
        "g" -> {Series[logx + 2*x, {x, 0, 6}], 0, 0}|>,
      <|"a" -> -1, "b" -> 1,
        "g" -> {Series[-logx + 3*x, {x, 0, 6}], 0, 0}|>
    },
    {0, 1},
    0
  ];
wholeLogPoleCancellation =
  DiffExp`RegularizedIntegration`IntegrateDecompositionLaurent[
    {<|"a" -> -1, "b" -> 1,
      "g" -> {Series[5*x, {x, 0, 6}], 0, 0}|>},
    {0, 1},
    0
  ];
If[laurentCloseQ[splitLogPoleCancellation, 0, {5, -5, 5}] &&
    laurentCloseQ[wholeLogPoleCancellation, 0, {5, -5, 5}],
  pass["endpoint Logx poles cancel consistently across split terms"],
  fail["endpoint Logx poles cancel consistently across split terms",
    {splitLogPoleCancellation, wholeLogPoleCancellation}]
];

testsTotal++;
splitMixedLogPoleCancellation =
  DiffExp`RegularizedIntegration`IntegrateDecompositionLaurent[
    {
      <|"a" -> -2, "b" -> 1,
        "g" -> {Series[3 + 2*x + 5*x*logx, {x, 0, 6}], 0, 0}|>,
      <|"a" -> -2, "b" -> 1,
        "g" -> {Series[-3 - 2*x - 5*x*logx + 7*x^2, {x, 0, 6}], 0, 0}|>
    },
    {0, 1},
    0
  ];
wholeMixedLogPoleCancellation =
  DiffExp`RegularizedIntegration`IntegrateDecompositionLaurent[
    {<|"a" -> -2, "b" -> 1,
      "g" -> {Series[7*x^2, {x, 0, 6}], 0, 0}|>},
    {0, 1},
    0
  ];
If[laurentCloseQ[splitMixedLogPoleCancellation, 0, {7, -7, 7}] &&
    laurentCloseQ[wholeMixedLogPoleCancellation, 0, {7, -7, 7}],
  pass["mixed endpoint Logx poles cancel after Laurent summation"],
  fail["mixed endpoint Logx poles cancel after Laurent summation",
    {splitMixedLogPoleCancellation, wholeMixedLogPoleCancellation}]
];

testsTotal++;
endpointLogSubtraction =
  DiffExp`RegularizedIntegration`IntegrateSingularTermLaurent[
    -1, 2, 0, {Series[logx, {x, 0, 1}], 0, 0}, {0, 1}
  ];
endpointLogSubtractionExpected = expectedLaurentCoeffs[
  -1/(2*ee)^2,
  ee, -2, 0
];
If[laurentCloseQ[endpointLogSubtraction, -2, endpointLogSubtractionExpected],
  pass["endpoint Logx jet is integrated by analytic monomial subtraction"],
  fail["endpoint Logx jet is integrated by analytic monomial subtraction",
    {endpointLogSubtraction, endpointLogSubtractionExpected}]
];

testsTotal++;
endpointLogSquaredSubtraction =
  DiffExp`RegularizedIntegration`IntegrateSingularTermLaurent[
    -1, 2, 0, {Series[logx^2, {x, 0, 1}], 0, 0}, {0, 1}
  ];
endpointLogSquaredSubtractionExpected = expectedLaurentCoeffs[
  2/(2*ee)^3,
  ee, -3, 0
];
If[laurentCloseQ[
    endpointLogSquaredSubtraction, -3, endpointLogSquaredSubtractionExpected],
  pass["endpoint Logx^2 jet produces the expected higher Laurent pole"],
  fail["endpoint Logx^2 jet produces the expected higher Laurent pole",
    {endpointLogSquaredSubtraction, endpointLogSquaredSubtractionExpected}]
];

testsTotal++;
mixedEndpointLogSubtraction =
  DiffExp`RegularizedIntegration`IntegrateSingularTermLaurent[
    -2, 2, 0, {Series[3 + 5*x*logx + 7*x^2, {x, 0, 5}], 0, 0}, {0, 1}
  ];
mixedEndpointLogSubtractionExpected = expectedLaurentCoeffs[
  3/(-1 + 2*ee) - 5/(2*ee)^2 + 7/(1 + 2*ee),
  ee, -2, 2
];
If[laurentCloseQ[
    mixedEndpointLogSubtraction, -2, mixedEndpointLogSubtractionExpected],
  pass["nonresonant and resonant endpoint logs subtract together"],
  fail["nonresonant and resonant endpoint logs subtract together",
    {mixedEndpointLogSubtraction, mixedEndpointLogSubtractionExpected}]
];

testsTotal++;
topOrderLogGuard =
  Quiet[
    DiffExp`RegularizedIntegration`IntegrateSingularTermLaurent[
      -1, 2, 0,
      {Series[3, {x, 0, 3}], Series[5, {x, 0, 3}],
       Series[7 + 11*logx + 13*logx^2, {x, 0, 3}]},
      {0, 1}
    ]
  ];
topOrderLogGuardExpected = expectedLaurentCoeffs[
  3/(2*ee) + 5*ee/(2*ee) + 7*ee^2/(2*ee),
  ee, -1, 2
];
If[laurentCloseQ[topOrderLogGuard, -1, topOrderLogGuardExpected],
  pass["top-order endpoint Logx terms do not back-propagate Laurent poles"],
  fail["top-order endpoint Logx terms do not back-propagate Laurent poles",
    {topOrderLogGuard, topOrderLogGuardExpected}]
];

testsTotal++;
If[formulaAgreementQ[
    -1, 2, 0,
    {Series[3 + 5*x + 7*x^2, {x, 0, 6}], 0, 0}
  ],
  pass["monomial subtraction agrees with recursive IBP for a=-1 Taylor jets"],
  fail["monomial subtraction agrees with recursive IBP for a=-1 Taylor jets",
    "mismatch"]
];

testsTotal++;
If[formulaAgreementQ[
    -2, 2, 0,
    {Series[3 + 5*x + 7*x^2, {x, 0, 6}], 0, 0}
  ],
  pass["monomial subtraction agrees with recursive IBP for a=-2 Taylor jets"],
  fail["monomial subtraction agrees with recursive IBP for a=-2 Taylor jets",
    "mismatch"]
];

testsTotal++;
If[Quiet[
    formulaAgreementQ[
      -1, 2, -1,
      {Series[1 + x, {x, 0, 6}], Series[2 - x, {x, 0, 6}],
       Series[3 + 5*logx + 7*logx^2, {x, 0, 6}]},
      0
    ]
  ],
  pass["top-order endpoint Logx guard preserves lower IBP-equivalent orders"],
  fail["top-order endpoint Logx guard preserves lower IBP-equivalent orders",
    "mismatch"]
];

testsTotal++;
integrableLogRemainder =
  DiffExp`RegularizedIntegration`IntegrateSingularTermLaurent[
    -1, 1, 0, {Series[x*logx, {x, 0, 3}], 0, 0}, {0, 1}
  ];
integrableLogRemainderExpected = expectedLaurentCoeffs[
  -1/(1 + ee)^2,
  ee, 0, 2
];
If[laurentCloseQ[integrableLogRemainder, 0, integrableLogRemainderExpected],
  pass["integrable Logx remainder is preserved by subtraction formula"],
  fail["integrable Logx remainder is preserved by subtraction formula",
    {integrableLogRemainder, integrableLogRemainderExpected}]
];

testsTotal++;
negativeSegmentSeries = Table[
  SeriesData[x, 0, {logx^k/k!}, -1, 0, 1],
  {k, 0, 2}
];
negativeSegmentData = {
  {<|x -> x|>, x -> x, {-1, 0}, {-1, 0}, {negativeSegmentSeries}}
};
prefactorPhaseResult =
  DiffExp`RegularizedIntegration`DefiniteIntegralWithPrefactorLaurent[
    negativeSegmentData,
    {-1, 0},
    <|"PowerAtLower" -> 0, "PowerAtUpper" -> 1,
      "RationalFactor" -> 1, "Variable" -> x|>,
    0
  ][[1]];
If[prefactorPhaseResult["MinPower"] == 0 &&
    Length[prefactorPhaseResult["Coefficients"]] >= 1 &&
    closeQ[prefactorPhaseResult["Coefficients"][[1]], -1],
  pass["upper endpoint prefactor keeps negative-side phase"],
  fail["upper endpoint prefactor keeps negative-side phase", prefactorPhaseResult]
];

testsTotal++;
artifactG = {
  0,
  Series[3 + 6*x, {x, 0, 2}],
  0,
  0,
  Series[(8/3)*x*logx^3, {x, 0, 2}]
};
artifactSeries = Table[
  Series[
    x^-2 * Sum[
      artifactG[[k]] * (2*logx)^(n - k)/Factorial[n - k],
      {k, 1, n}
    ],
    {x, 0, 2}
  ],
  {n, Length[artifactG]}
];
artifactSegmentData = {
  {<|x -> x|>, x -> x, {0, 1}, {0, 1}, {artifactSeries}}
};
artifactResult =
  DiffExp`RegularizedIntegration`DefiniteIntegralWithPrefactorLaurent[
    artifactSegmentData,
    {0, 1},
    <|"PowerAtLower" -> 0, "PowerAtUpper" -> 0,
      "RationalFactor" -> 1, "Variable" -> x|>,
    -2
  ][[1]];
If[artifactResult["MinPower"] == -2 &&
    closeQ[laurentCoeff[artifactResult, -2], 3],
  pass["reconstructed endpoint Logx artifacts do not generate extra poles"],
  fail["reconstructed endpoint Logx artifacts do not generate extra poles",
    artifactResult]
];

(* ----------------------------------------------------------------- *)
(* Residual endpoint sector recovery (multi-sector x^(a + b eps))     *)
(* ----------------------------------------------------------------- *)

(* Build the eps-order series list of sum_i w_i(eps) x^a x^(b_i eps),
   i.e. the towers DecomposeSingularity sees when several sectors share
   the same integer power a.  ws[[i]] is the eps-coefficient list of the
   weight w_i(eps). *)
sectorPow[base_, 0] := 1;
sectorPow[base_, k_Integer] := base^k;
sectorTowerSeries[a_, bs_List, ws_List, epsOrders_Integer] := Table[
  SeriesData[x, 0,
    {Sum[
       Sum[
         If[q + 1 <= Length[ws[[i]]],
           ws[[i, q + 1]] *
             sectorPow[bs[[i]] * logx, n - 1 - q] / Factorial[n - 1 - q],
           0
         ],
         {q, 0, n - 1}
       ],
       {i, Length[bs]}
     ]},
    a, a + 1, 1
  ],
  {n, 1, epsOrders}
];

splitSectorIntegral[a_, bs_List, ws_List, bounds_List] := Module[
  {parts},
  parts = Table[
    DiffExp`RegularizedIntegration`IntegrateSingularTermLaurent[
      a, bs[[i]], 0, ws[[i]], bounds
    ],
    {i, Length[bs]}
  ];
  Fold[
    DiffExp`RegularizedIntegration`Private`LaurentAdd,
    First[parts],
    Rest[parts]
  ]
];

combinedSectorIntegral[a_, bs_List, ws_List, epsOrders_Integer, bounds_List] :=
  DiffExp`RegularizedIntegration`IntegrateDecompositionLaurent[
    DiffExp`SingularityDecomposition`DecomposeSingularity[
      sectorTowerSeries[a, bs, ws, epsOrders]
    ],
    bounds,
    0
  ];

laurentAgreeQ[laurA_, laurB_, min_Integer, max_Integer, tol_:10^-35] :=
  And @@ Table[
    closeQ[laurentCoeff[laurA, power], laurentCoeff[laurB, power], tol],
    {power, min, max}
  ];

testsTotal++;
twoSectorConstW = Quiet[combinedSectorIntegral[
  -1, {2, 1}, {{3, 0, 0, 0, 0, 0, 0}, {5, 0, 0, 0, 0, 0, 0}}, 7, {0, 1/2}
]];
twoSectorConstWExpected = expectedLaurentCoeffs[
  3*(1/2)^(2*ee)/(2*ee) + 5*(1/2)^ee/ee,
  ee, -1, 3
];
If[laurentCloseQ[twoSectorConstW, -1, twoSectorConstWExpected],
  pass["two sectors with equal a and different b are recovered exactly"],
  fail["two sectors with equal a and different b are recovered exactly",
    {twoSectorConstW, twoSectorConstWExpected}]
];

testsTotal++;
twoSectorEpsW = Quiet[combinedSectorIntegral[
  -1, {2, 1},
  {{3, 2, -1, 0, 0, 0, 0, 0}, {5, -1, 4, 0, 0, 0, 0, 0}}, 8, {0, 1/2}
]];
twoSectorEpsWExpected = expectedLaurentCoeffs[
  (3 + 2*ee - ee^2)*(1/2)^(2*ee)/(2*ee) +
    (5 - ee + 4*ee^2)*(1/2)^ee/ee,
  ee, -1, 3
];
If[laurentCloseQ[twoSectorEpsW, -1, twoSectorEpsWExpected],
  pass["epsilon-dependent sector weights are resummed exactly"],
  fail["epsilon-dependent sector weights are resummed exactly",
    {twoSectorEpsW, twoSectorEpsWExpected}]
];

testsTotal++;
(* Leading sector weights cancel: the first nonzero epsilon order is pure
   Logx, so the moment s0 vanishes.  This is the banana {1,0,2,1} failure
   mode; the two-sector solve must proceed through the s0 = 0 Hankel
   system instead of bailing out to the explicit-log fallback. *)
cancellingW = Quiet[combinedSectorIntegral[
  -1, {2, 1}, {{4, 0, 0, 0, 0, 0, 0}, {-4, 0, 0, 0, 0, 0, 0}}, 7, {0, 1/2}
]];
cancellingWExpected = expectedLaurentCoeffs[
  4*(1/2)^(2*ee)/(2*ee) - 4*(1/2)^ee/ee,
  ee, -1, 3
];
If[laurentCloseQ[cancellingW, -1, cancellingWExpected],
  pass["cancelling leading weights (s0 = 0) still resolve two sectors"],
  fail["cancelling leading weights (s0 = 0) still resolve two sectors",
    {cancellingW, cancellingWExpected}]
];

testsTotal++;
(* Same two-sector structure on a negative local interval, where the
   integration endpoint carries analytic-continuation phases.  Compare the
   tower-recovered result against integrating each known sector directly,
   so the test is independent of the branch convention. *)
negativeCombined = Quiet[combinedSectorIntegral[
  -1, {2, 1}, {{3, 1, 0, 0, 0, 0, 0}, {5, -2, 0, 0, 0, 0, 0}}, 7, {-1, 0}
]];
negativeSplit = Quiet[splitSectorIntegral[
  -1, {2, 1}, {{3, 1, 0, 0, 0, 0, 0}, {5, -2, 0, 0, 0, 0, 0}}, {-1, 0}
]];
If[laurentAgreeQ[negativeCombined, negativeSplit, -1, 3],
  pass["negative-axis two-sector recovery matches direct sector integration"],
  fail["negative-axis two-sector recovery matches direct sector integration",
    {negativeCombined, negativeSplit}]
];

testsTotal++;
(* Single sector with an epsilon-dependent weight: the residual fit sees a
   zero root, which must go through the zero-power-safe convention instead
   of generating 0^0 and falling back to the explicit-log path. *)
zeroRootSector = Quiet[combinedSectorIntegral[
  -1, {3}, {{1, 5, 0, 0, 0, 0}}, 6, {0, 1}
]];
zeroRootSectorExpected = expectedLaurentCoeffs[
  (1 + 5*ee)/(3*ee),
  ee, -1, 3
];
If[laurentCloseQ[zeroRootSector, -1, zeroRootSectorExpected],
  pass["zero residual root (single sector) validates without 0^0"],
  fail["zero residual root (single sector) validates without 0^0",
    {zeroRootSector, zeroRootSectorExpected}]
];

testsTotal++;
(* Truncated data: zero out the top epsilon order, as happens when an
   incompletely combined order slips through.  The validated-prefix logic
   must keep all lower Laurent orders exact and confine the damage to the
   top of the window. *)
truncatedSeries = sectorTowerSeries[
  -1, {2, 1}, {{3, 0, 0, 0, 0, 0, 0, 0}, {5, 0, 0, 0, 0, 0, 0, 0}}, 8
];
truncatedSeries[[8]] = SeriesData[x, 0, {0}, -1, 0, 1];
truncatedResult = Quiet[
  DiffExp`RegularizedIntegration`IntegrateDecompositionLaurent[
    DiffExp`SingularityDecomposition`DecomposeSingularity[truncatedSeries],
    {0, 1/2},
    0
  ]
];
truncatedExpected = expectedLaurentCoeffs[
  3*(1/2)^(2*ee)/(2*ee) + 5*(1/2)^ee/ee,
  ee, -1, 2
];
If[laurentCloseQ[truncatedResult, -1, truncatedExpected],
  pass["a truncated top epsilon order does not corrupt lower Laurent orders"],
  fail["a truncated top epsilon order does not corrupt lower Laurent orders",
    {truncatedResult, truncatedExpected}]
];

(* ----------------------------------------------------------------- *)
(* Sector-aware endpoint limits (EvaluateEndpointLimitSectors)        *)
(* ----------------------------------------------------------------- *)

(* Like sectorTowerSeries, but with sectors at several integer powers:
   specs is a list of <|"a" -> integer, "bs" -> {..}, "ws" -> {..}|>. *)
multiPowerSectorSeries[specs_List, epsOrders_Integer] := Module[
  {amin, amax},
  amin = Min[#["a"] & /@ specs];
  amax = Max[#["a"] & /@ specs];
  Table[
    SeriesData[x, 0,
      Table[
        Total[Table[
          If[spec["a"] === p,
            Sum[
              Sum[
                If[q + 1 <= Length[spec["ws"][[i]]],
                  spec["ws"][[i, q + 1]] *
                    sectorPow[spec["bs"][[i]] * logx, n - 1 - q] /
                    Factorial[n - 1 - q],
                  0
                ],
                {q, 0, n - 1}
              ],
              {i, Length[spec["bs"]]}
            ],
            0
          ],
          {spec, specs}
        ]],
        {p, amin, amax}
      ],
      amin, amax + 1, 1
    ],
    {n, 1, epsOrders}
  ]
];

testsTotal++;
(* The banana {1,0,0,1} limitUpper mechanism: endpoint data mixing
   x^(-1 + eps), x^(2 eps), and x^0 sectors.  DecomposeSingularity
   extracts the most negative power's exponent (a = -1, b = 1), burying
   the genuine b = 0 sector inside g; the naive "drop the term unless
   a >= 0 and b == 0" rule then returns 0.  The sector-aware limit must
   recover exactly the b = 0 sector's weight tower. *)
limitTowerW0 = {3, 1, -4, 0, 0, 0, 0, 0};
limitMultiSector = Quiet[
  DiffExp`RegularizedIntegration`EvaluateEndpointLimitSectors[
    multiPowerSectorSeries[{
      <|"a" -> -1, "bs" -> {1}, "ws" -> {{7, -2, 5, 0, 0, 0, 0, 0}}|>,
      <|"a" -> 0, "bs" -> {0, 2},
        "ws" -> {limitTowerW0, {-6, 2, 9, 0, 0, 0, 0, 0}}|>
    }, 8],
    1
  ]
];
If[ListQ[limitMultiSector] &&
    And @@ Table[
      closeQ[limitMultiSector[[k]], limitTowerW0[[k]]],
      {k, 5}
    ],
  pass["multi-sector endpoint limit keeps exactly the b = 0 sector"],
  fail["multi-sector endpoint limit keeps exactly the b = 0 sector",
    limitMultiSector]
];

testsTotal++;
(* Single sector x^(2 eps) at a = 0 (clean b extraction): regulated to
   zero in the limit. *)
limitPureBSector = Quiet[
  DiffExp`RegularizedIntegration`EvaluateEndpointLimitSectors[
    sectorTowerSeries[0, {2}, {{5, 3, 0, 0, 0, 0}}, 6],
    1
  ]
];
If[ListQ[limitPureBSector] &&
    And @@ Table[closeQ[limitPureBSector[[k]], 0], {k, 4}],
  pass["pure b != 0 sector limit is regulated to zero"],
  fail["pure b != 0 sector limit is regulated to zero", limitPureBSector]
];

testsTotal++;
(* Plain Taylor data must pass through unchanged (the single-sector
   b = 0 case every other example exercises). *)
limitTaylor = Quiet[
  DiffExp`RegularizedIntegration`EvaluateEndpointLimitSectors[
    {
      SeriesData[x, 0, {4, 11}, 0, 2, 1],
      SeriesData[x, 0, {-9, 2}, 0, 2, 1],
      SeriesData[x, 0, {13/7, 1}, 0, 2, 1]
    },
    1
  ]
];
If[ListQ[limitTaylor] &&
    closeQ[limitTaylor[[1]], 4] &&
    closeQ[limitTaylor[[2]], -9] &&
    closeQ[limitTaylor[[3]], 13/7],
  pass["plain Taylor endpoint limit passes through unchanged"],
  fail["plain Taylor endpoint limit passes through unchanged", limitTaylor]
];

Print["==========================================="];
Print["Results: ", testsPassed, " / ", testsTotal, " tests passed"];
If[testsPassed === testsTotal,
  Print["All tests PASSED!"],
  Print["Some tests FAILED"];
  Exit[1]
];
Print["==========================================="];
