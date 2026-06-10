(* ::Package:: *)

(* DiffExp SeriesOps Subpackage *)
(* This package provides series manipulation operations *)

BeginPackage["DiffExp`SeriesOps`", {"DiffExp`Symbols`", "DiffExp`State`", "DiffExp`Utilities`"}];

(* Series operations *)
SApply::usage = "SApply[f_,a_] applies function f to series coefficients.";
SExpand::usage = "SExpand[a_] expands a series with chopping.";
SN::usage = "SN applies N to series.";
SSN::usage = "SSN applies N to series.";
SMultiply::usage = "SMultiply[a_,b_] multiplies two series with expansion.";
SeriesCoefficientMinus::usage = "SeriesCoefficientMinus[a_,k_] gets coefficient from end of series.";
ApplyAnalyticContinuation::usage = "ApplyAnalyticContinuation[s_] applies analytic continuation to series.";
SafeReplaceSeries11::usage = "SafeReplaceSeries11[a_,b_] safely replaces in series.";

(* Log operations *)
MaxLogxPower::usage = "MaxLogxPower[ex_] returns maximum power of Logx in expression.";
LogxPowerRange::usage = "LogxPowerRange[ex_] returns the non-negative integer Logx powers to scan.";
LogxCoeff::usage = "LogxCoeff[Ser_,Which_] extracts coefficient of Logx^Which from series.";
LogxCoeffNS::usage = "LogxCoeffNS[Ser_,Which_] extracts coefficient of Logx^Which (non-series).";
LogxCoeffList::usage = "LogxCoeffList[Ser_] returns list of Logx coefficients.";

(* Matrix operations *)
MatrixMultiplySExpand::usage = "MatrixMultiplySExpand[MatA_,MatB_] multiplies matrices with series expansion.";
MatrixPowerSExpand::usage = "MatrixPowerSExpand[a_,n_] raises matrix to power with expansion.";

(* DiffExp series functions *)
DiffExpSeries::usage = "DiffExpSeries[Ser_,ord_] creates a DiffExp series.";
SeriesAlways::usage = "SeriesAlways[term_,{a_,b_,c_},ex_] creates series even for constant terms.";
LeadingCoefficientSeries::usage = "LeadingCoefficientSeries[Ser_] extracts leading coefficient as series.";
SeriesMinPower::usage = "SeriesMinPower[Ser_] returns minimum power in series.";
SeriesMaxPower::usage = "SeriesMaxPower[Ser_] returns maximum power in series.";
DecreaseSeriesOrderBy::usage = "DecreaseSeriesOrderBy[a_,k_] decreases series order.";

(* Series derivative *)
SD::usage = "SD[a_,b_] takes derivative avoiding Log[x] terms.";

Begin["`Private`"];

activeNumericPrecision[] := Module[{precision},
  precision = Quiet[Check[DiffExp`State`FEWorkingPrecision, 500]];
  If[IntegerQ[precision] && precision > 0, precision, 500]
];

numericAtActivePrecision[expr_, precision_:Automatic] := Module[
  {p = If[precision === Automatic, activeNumericPrecision[], precision], val},
  val = Quiet[Check[N[expr, p], $Failed]];
  If[val === $Failed, val, SetPrecision[val, p]]
];

(* Set attributes as in original code *)
(* Note: SEval1, SEval2 are defined in Pade.m and should have their attributes set there *)
SetAttributes[{SApply, SN, SExpand, SN, SSN, SMultiply, SEval, DecreaseSeriesOrderBy, SeriesCoefficientMinus, SplitTimes, ApplyAnalyticContinuation}, Listable];
SetAttributes[DiffExpSeries, Listable];

(* Series application *)
SApply[f_, 0] := 0;
SApply[f_, a_SeriesData] := MapAt[f, a, 3];
SApply[f_, a_] := f[a];

(* Series expand with chopping *)
SExpand[0] := 0;
SExpand[a_SeriesData] := SApply[DiffExp`Utilities`PChop@*Expand, a];
SExpand[a_] := Expand[a];

(* Numeric conversions *)
SN := SApply[N, #] &;
SSN := SN;

(* Series multiplication *)
SMultiply[a_, b_] := SExpand[a] SExpand[b] // SExpand;

(* Get coefficient from end of series *)
SeriesCoefficientMinus[a_SeriesData, k_: 1] := a[[3]][[Max[-k, -Length[a[[3]]]]]] DiffExp`Symbols`x^((a[[5]] - k)/a[[6]]);
SeriesCoefficientMinus[a_, k_: 1] := 0;

(* Apply analytic continuation *)
ApplyAnalyticContinuation[s_SeriesData] := Series[
  Normal[s] /. DiffExp`State`AnalyticContinuationReplacements,
  {DiffExp`Symbols`x, 0, Floor[s[[5]]/s[[6]]]}
];

(* Safe replacement in series *)
SafeReplaceSeries11[a_SeriesData, b_] := Block[{ex = Normal[a], ex2},
  ex2 = Series[(ex /. b), {DiffExp`Symbols`x, 0, Ceiling[a[[5]]/a[[6]]]}];
  ex2 + O[DiffExp`Symbols`x]^(a[[5]]/a[[6]])
];
SafeReplaceSeries11[a_List, b_] := SafeReplaceSeries11[#, b] & /@ a;
SafeReplaceSeries11[a_, b_] := a /. b;

(* Log operations *)
NormalizeLogPower[p_] := Module[
  {tol = DiffExp`State`FEC[RationalizationTolerance], np, rounded, exact},
  If[IntegerQ[p], Return[p, Module]];
  np = numericAtActivePrecision[p];
  If[np =!= $Failed && NumericQ[np] &&
      TrueQ[Abs[np - Round[np]] < tol],
    rounded = Round[np];
    exact = Rationalize[rounded, 0];
    If[IntegerQ[exact], exact, rounded],
    p
  ]
];

MaxLogxPower[ex_] := Module[{rawPowers, powers},
  rawPowers = Append[
    DiffExp`Utilities`GetCases[
      ex // SExpand, DiffExp`Symbols`Logx^(k_: 1) :> k
    ],
    0
  ];
  powers = Cases[
    NormalizeLogPower /@ Flatten[{rawPowers}],
    _Integer?NonNegative
  ];
  If[powers === {}, 0, Max[powers]]
];

LogxPowerRange[ex_] := Module[{maxpow = NormalizeLogPower[MaxLogxPower[ex]]},
  If[IntegerQ[maxpow] && maxpow >= 0, Range[0, maxpow], {0}]
];

LogxCoeff[Ser_, Which_] := If[Which === 0,
  SApply[(# /. DiffExp`Symbols`Logx -> 0) &, Ser],
  SApply[(Coefficient[#, DiffExp`Symbols`Logx^Which]) &, Ser]
];

LogxCoeffNS[Ser_, Which_] := If[Which === 0,
  Ser /. DiffExp`Symbols`Logx -> 0,
  Coefficient[Ser, DiffExp`Symbols`Logx^Which]
];

LogxCoeffList[Ser_] := LogxCoeff[Ser, #] & /@ LogxPowerRange[Ser];

(* Matrix multiplication with series expansion *)
MatrixMultiplySExpand[MatA_, MatB_] := Module[{Dim1 = Dimensions[MatA], Dim2 = Dimensions[MatB]},
  If[!(Dim1[[2]] === Dim2[[1]]),
    DiffExp`State`LastErrorContext = {MatA, MatB};
    DiffExp`Utilities`ReportError["Matrix dimensions don't match."];
  ];

  (* Non-parallel version - parallel code commented out in original *)
  Table[
    SExpand@Sum[SMultiply[MatA[[iind, kind]], MatB[[kind, jind]]], {kind, Dim1[[2]]}],
    {iind, Dim1[[1]]}, {jind, Dim2[[2]]}
  ]
];

(* Matrix power with series expansion *)
MatrixPowerSExpand[a_, 0] := Block[{dim = a // Dimensions},
  If[!(dim[[1]] === dim[[2]]),
    DiffExp`Utilities`ReportError["Matrix is not square!"];
  ];
  IdentityMatrix[dim[[1]]]
];
MatrixPowerSExpand[a_, n_] := Nest[MatrixMultiplySExpand[a, #] &, a, n - 1];

(* DiffExp series creation *)
DiffExpSeries[Ser_, ord_] := SeriesAlways[
  Ser // Quiet[N[#, DiffExp`State`FEWorkingPrecision]] &,
  {DiffExp`Symbols`x, 0, ord}
];
DiffExpSeries[Ser_] := DiffExpSeries[Ser, DiffExp`State`ExpansionOrderVal];

SeriesAlways[term_, {a_, b_, c_}, ex_: 1] := If[DiffExp`Utilities`DependsQ[term, a],
  Series[term, {a, b, c}, Assumptions -> DiffExp`Symbols`x > 0],
  SeriesData[a, b, {term}, 0, c ex + 1, ex]
];

(* Leading coefficient series *)
LeadingCoefficientSeries[Ser_, AddTo2_: 1] := Module[{leadingSeries, truncatedSeries, combinedDenominator, Den1, Den2, fractionalPower},
  If[!DiffExp`Utilities`DependsQ[Ser // Normal, DiffExp`Symbols`x],
    Ser + O[DiffExp`Symbols`x]^(1/AddTo2),
    Assuming[DiffExp`Symbols`x > 0,
      leadingSeries = Ser /. DiffExp`Symbols`x -> (DiffExp`Symbols`x + O[DiffExp`Symbols`x]^2);
      truncatedSeries = leadingSeries /. a_List :> Take[a, 1];

      Den1 = Association[#[[1]] -> #[[2]] & /@ (leadingSeries[[6]] // FactorInteger)];
      Den2 = Association[#[[1]] -> #[[2]] & /@ (AddTo2 // FactorInteger)];
      combinedDenominator = Times @@ KeyValueMap[#1^#2 &, Merge[{Den1, Den2}, Max]];

      fractionalPower = leadingSeries[[4]]/leadingSeries[[6]] + 1/combinedDenominator;
      truncatedSeries + SeriesData[DiffExp`Symbols`x, 0, List[], fractionalPower // Numerator, fractionalPower // Numerator, fractionalPower // Denominator]
    ]
  ]
];

(* Series min/max power *)
SeriesMinPower[Ser_] := Ser[[4]] Ser[[6]];
SeriesMaxPower[Ser_] := Ser[[5]] Ser[[6]];

(* Decrease series order *)
DecreaseSeriesOrderBy[a_, k_: 1] := Block[{tmp},
  If[a[[0]] === SeriesData,
    tmp = a;
    tmp[[5]] = tmp[[5]] - k tmp[[6]];
    tmp,
    a
  ]
];

(* Series derivative avoiding Log[x] terms *)
SD[a_, b_] := Block[{Tmp},
  Sum[
    Tmp = LogxCoeff[a, logxord];
    D[Tmp, b] DiffExp`Symbols`Logx^logxord + (D[Log[DiffExp`Symbols`x]^logxord, b] /. Log[DiffExp`Symbols`x] -> DiffExp`Symbols`Logx) Tmp,
    {logxord, LogxPowerRange[a]}
  ] // SExpand
];
SD[a_, b__] := SD[SD[a, {b} // First], Sequence @@ Delete[{b}, 1]];
SD[a_] := a;

End[];

EndPackage[];
