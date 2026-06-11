(* Tests for the N-root / confluent residual endpoint sector fitter
   (FitResidualEndpointSectors with N >= 3 plain sectors and repeated
   roots), plus the over-fit displacement regression: truncation-doctored
   towers must keep the true root count, keep the leading Laurent orders
   exact, shrink the trustworthy window, and warn about omitted content. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];
Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];

DiffExp`State`DiffExpConfiguration[RationalizationTolerance] = 10^-30;
DiffExp`State`DiffExpConfiguration[ChopPrecision] = 80;
DiffExp`State`DiffExpConfiguration[WorkingPrecision] = 100;

passed = 0; failed = 0;
test[name_String, ok_] := If[TrueQ[ok],
  passed++; Print["  PASS: ", name],
  failed++; Print["  FAIL: ", name]
];

x = DiffExp`Symbols`x;
logx = DiffExp`Symbols`Logx;
ee = Unique["eps"];

closeQ[a_, b_, tol_:10^-20] := TrueQ[Abs[N[a - b, 80]] < tol];
laurentCoeff[laur_, power_Integer] := Module[
  {idx = power - laur["MinPower"] + 1},
  If[idx >= 1 && idx <= Length[laur["Coefficients"]],
    laur["Coefficients"][[idx]],
    0
  ]
];
expectedLaurentCoeffs[expr_, var_, min_Integer, max_Integer] :=
  Table[SeriesCoefficient[expr, {var, 0, power}], {power, min, max}];
laurentCloseQ[laur_, min_Integer, coeffs_List, tol_:10^-20] :=
  And @@ Table[
    closeQ[laurentCoeff[laur, min + idx - 1], coeffs[[idx]], tol],
    {idx, Length[coeffs]}
  ];

(* Local replica of the sector-tower builders from
   test_regularized_integration_edge_cases.m, generalized to confluent
   specs {b_i, p_i}: each sector contributes
   w_i(eps) x^(b_i eps) (eps Logx)^(p_i), so the eps^(n-1) coefficient is
   sum_q w_{i,q} (b_i Logx)^k / k! * Logx^(p_i) with k = n-1-q-p_i. *)
sectorPow[base_, 0] := 1;
sectorPow[base_, k_Integer] := base^k;
specTowerCoeffs[specs_List, ws_List, epsOrders_Integer] := Table[
  Sum[
    Module[{bi = specs[[i, 1]], pi = specs[[i, 2]]},
      Sum[
        Module[{k = n - 1 - q - pi},
          If[k >= 0 && q + 1 <= Length[ws[[i]]],
            ws[[i, q + 1]] * sectorPow[bi * logx, k] / Factorial[k] *
              sectorPow[logx, pi],
            0
          ]
        ],
        {q, 0, n - 1}
      ]
    ],
    {i, Length[specs]}
  ],
  {n, 1, epsOrders}
];
specTowerSeries[a_, coeffs_List] :=
  Map[SeriesData[x, 0, {#}, a, a + 1, 1] &, coeffs];
towerIntegral[a_, coeffs_List, bounds_List] :=
  DiffExp`RegularizedIntegration`IntegrateDecompositionLaurent[
    DiffExp`SingularityDecomposition`DecomposeSingularity[
      specTowerSeries[a, coeffs]
    ],
    bounds,
    0
  ];

(* PrintWarning is a plain Print; capture honesty warnings by shadowing
   it inside a Block so the omitted-content reporting can be asserted. *)
$capturedWarnings = {};
SetAttributes[captureWarnings, HoldFirst];
captureWarnings[body_] := Module[{res},
  $capturedWarnings = {};
  Block[{DiffExp`Utilities`PrintWarning},
    DiffExp`Utilities`PrintWarning[args___] := AppendTo[
      $capturedWarnings, StringJoin[ToString /@ {args}]
    ];
    res = body;
  ];
  res
];

Print["==========================================="];
Print["Multi-Sector / Confluent Fit Tests"];
Print["==========================================="];

(* ----------------------------------------------------------------- *)
(* (a) Three plain sectors: x^-1 (w1 x^eps + w2 x^(2 eps) + w3 x^(3 eps)) *)
(* ----------------------------------------------------------------- *)

bsA = {1, 2, 3};
wA = {
  {3, 1, -4, 2, -1, 5, 1, -3, 2, 1, -2, 4},
  {7, -2, 5, -3, 1, -1, 2, 4, -5, 3, 1, -1},
  {-6, 2, 9, 1, -2, 3, -1, 1, 2, -3, 4, 1}
};
epsOrdersA = 12;
specsA = {{1, 0}, {2, 0}, {3, 0}};
towerA = specTowerCoeffs[specsA, wA, epsOrdersA];

fitA = DiffExp`RegularizedIntegration`FitResidualEndpointSectors[
  towerA, {}
];
fitASectors = If[AssociationQ[fitA], fitA["Sectors"], {}];

test["(a) fit returns exactly 3 sectors",
  Length[fitASectors] === 3
];
test["(a) residual exponents are {1, 2, 3} with LogPower 0",
  Length[fitASectors] === 3 &&
    (#["ResidualB"] & /@ fitASectors) === {1, 2, 3} &&
    (#["LogPower"] & /@ fitASectors) === {0, 0, 0}
];
(* usable weight orders for 12 offsets / 3 sectors: q = 0..9 *)
test["(a) sector weights match the input Laurents to 20+ digits",
  Length[fitASectors] === 3 &&
    And @@ Flatten[Table[
      closeQ[fitASectors[[i]]["Coefficients"][[q + 1]], wA[[i, q + 1]]],
      {i, 3}, {q, 0, 9}
    ]]
];

wPolyA[i_] := Sum[wA[[i, q + 1]] * ee^q, {q, 0, epsOrdersA - 1}];
closedFormA = Sum[
  wPolyA[i] * (1/2)^(bsA[[i]] * ee) / (bsA[[i]] * ee),
  {i, 3}
];

resA = captureWarnings[towerIntegral[-1, towerA, {0, 1/2}]];
test["(a) integral matches the closed form at orders ee^-1..ee^3",
  AssociationQ[resA] && resA["MinPower"] == -1 &&
    laurentCloseQ[resA, -1, expectedLaurentCoeffs[closedFormA, ee, -1, 3]]
];

(* ----------------------------------------------------------------- *)
(* (b) Confluent double root: w1 x^eps + w2 x^eps (eps Logx)          *)
(* ----------------------------------------------------------------- *)

wB = {
  {5, -1, 3, 2, -2, 1, 4, -3, 1, 2, -1, 1},
  {3, 2, -4, 1, 5, -2, 1, 3, -1, -2, 2, 1}
};
epsOrdersB = 12;
specsB = {{1, 0}, {1, 1}};
towerB = specTowerCoeffs[specsB, wB, epsOrdersB];

fitB = DiffExp`RegularizedIntegration`FitResidualEndpointSectors[
  towerB, {}
];
fitBSectors = If[AssociationQ[fitB], fitB["Sectors"], {}];

test["(b) confluent fit returns exactly 2 sectors",
  Length[fitBSectors] === 2
];
test["(b) double root r = 1 recovered with LogPower {0, 1}",
  Length[fitBSectors] === 2 &&
    (#["ResidualB"] & /@ fitBSectors) === {1, 1} &&
    (#["LogPower"] & /@ fitBSectors) === {0, 1}
];
(* usable weight orders for 12 offsets / 2 sectors: q = 0..10 *)
test["(b) confluent sector weights match the input Laurents to 20+ digits",
  Length[fitBSectors] === 2 &&
    And @@ Flatten[Table[
      closeQ[fitBSectors[[i]]["Coefficients"][[q + 1]], wB[[i, q + 1]]],
      {i, 2}, {q, 0, 10}
    ]]
];

(* closed form: Int_0^c x^(-1 + b ee) dx = c^(b ee)/(b ee) =: F(b); the
   (eps Logx) sector integrates to the b-derivative at b = 1,
   d/db F(b) = c^(b ee) (Log[c]/b - 1/(b^2 ee)). *)
wPolyB[i_] := Sum[wB[[i, q + 1]] * ee^q, {q, 0, epsOrdersB - 1}];
closedFormB = wPolyB[1] * (1/2)^ee / ee +
  wPolyB[2] * (1/2)^ee * (Log[1/2] - 1/ee);

resB = captureWarnings[towerIntegral[-1, towerB, {0, 1/2}]];
test["(b) integral matches the b-derivative closed form at ee^-1..ee^3",
  AssociationQ[resB] && resB["MinPower"] == -1 &&
    laurentCloseQ[resB, -1, expectedLaurentCoeffs[closedFormB, ee, -1, 3]]
];

(* ----------------------------------------------------------------- *)
(* (c) Adversarial truncation: zero the highest-log-power slices of   *)
(* the top two eps orders of (a)'s tower.  Regression for the         *)
(* over-fit displacement mechanism: the fitter must keep N = 3 (the   *)
(* extra-root candidates are unfalsifiable / non-dominant), keep the  *)
(* leading Laurent orders exact, shrink the trustworthy window, and   *)
(* report the omitted content.                                        *)
(* ----------------------------------------------------------------- *)

dropTopLog[expr_, k_Integer] := Expand[
  expr - Coefficient[expr, logx, k] * logx^k
];
towerC = towerA;
towerC[[epsOrdersA]] = dropTopLog[towerC[[epsOrdersA]], epsOrdersA - 1];
towerC[[epsOrdersA - 1]] = dropTopLog[towerC[[epsOrdersA - 1]], epsOrdersA - 2];

fitC = DiffExp`RegularizedIntegration`FitResidualEndpointSectors[
  towerC, {}
];
fitCSectors = If[AssociationQ[fitC], fitC["Sectors"], {}];

test["(c) doctored tower still selects exactly 3 sectors (not 4+)",
  Length[fitCSectors] === 3 &&
    (#["ResidualB"] & /@ fitCSectors) === {1, 2, 3}
];
test["(c) trustworthy window shrank (salvage starts at offset 7, was 10)",
  AssociationQ[fitA] && AssociationQ[fitC] &&
    Min[fitA["SalvageOffsets"]] === 10 &&
    Min[fitC["SalvageOffsets"]] === 7
];

resC = captureWarnings[towerIntegral[-1, towerC, {0, 1/2}]];
warningsC = $capturedWarnings;
test["(c) the three leading Laurent orders remain exact",
  AssociationQ[resC] && resC["MinPower"] == -1 &&
    laurentCloseQ[resC, -1, expectedLaurentCoeffs[closedFormA, ee, -1, 1]]
];
test["(c) omitted-content trust warning fires",
  AnyTrue[warningsC, StringContainsQ[#, "NOT trustworthy"] &]
];

Print["\nResults: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
