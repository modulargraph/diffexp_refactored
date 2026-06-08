(* ::Package:: *)

(* DiffExp Pade Subpackage *)
(* This package provides Pade approximant functions *)

BeginPackage["DiffExp`Pade`", {"DiffExp`Symbols`", "DiffExp`State`", "DiffExp`Utilities`", "DiffExp`SeriesOps`"}];

(* Pade functions *)
GetPade::usage = "GetPade[a_] returns the Pade approximant of a series.";
SEval1::usage = "SEval1[a_] evaluates series, optionally using Pade approximant.";
SEval2::usage = "SEval2[a_,at_] evaluates expression at a point with analytic continuation.";
SEval::usage = "SEval[a_,at_] evaluates series at a point.";

Begin["`Private`"];

SafeSeriesCoefficient[s_SeriesData, ord_] := Module[
  {order = Round[ord], minOrder, maxOrder, step},

  minOrder = s[[4]];
  maxOrder = s[[5]];
  step = s[[6]];

  If[order < minOrder || order >= maxOrder || !IntegerQ[(order - minOrder)/step],
    0,
    SeriesCoefficient[s, order]
  ]
];

(* GetPade for various inputs *)
GetPade[0] = 0;
GetPade[a_?NumericQ] := a;
GetPade[a_SeriesData] := Block[
  {
    $MinPrecision = DiffExp`State`FEWorkingPrecision,
    maxPadeOrder = {Floor[((a[[5]] - a[[4]])/a[[6]] + 1)/2], Floor[((a[[5]] - a[[4]])/a[[6]] + 1)/2]},
    MaxLogOrder = DiffExp`SeriesOps`MaxLogxPower[a],
    tempPadeApprox
  },

  Sum[
    DiffExp`Symbols`Logx^ind (
      (Quiet[
        tempPadeApprox = PadeApproximant[Chop[#, 10^-DiffExp`State`ChopPrecisionVal], {DiffExp`Symbols`x, 0, maxPadeOrder}];
        If[Length[DiffExp`Utilities`GetCases[tempPadeApprox, PadeApproximant[__]]] > 0,
          DiffExp`Utilities`PrintWarning["Error in Pade approximant of ", # // N, ". Evaluating normally."][1];
          DiffExp`Utilities`PChop@#,
          tempPadeApprox
        ]
      ]) &@(Normal@DiffExp`SeriesOps`LogxCoeff[a, ind])
    ),
    {ind, 0, MaxLogOrder}
  ]
];

(* SEval1 - get evaluable form of series *)
SetAttributes[SEval1, Listable];
SEval1[a_SeriesData] := If[DiffExp`State`FEC[UsePade] === True, GetPade[a], Normal[a]];
SEval1[0] := 0;
SEval1[a_ /; NumericQ[a]] := a;
NormalizeEmbeddedSeries[a_] := a /. {
      HoldPattern[SeriesCoefficient[s_SeriesData, {_, _, ord_}]] :>
        SafeSeriesCoefficient[s, ord],
      HoldPattern[SeriesCoefficient[s_SeriesData, ord_]] :>
        SafeSeriesCoefficient[s, ord]
    } /. s_SeriesData :> SEval1[s];
SEval1[a_] := NormalizeEmbeddedSeries[a];

(* SEval2 - evaluate at a point with analytic continuation *)
SetAttributes[SEval2, Listable];
SEval2[a_, at_] := Block[{$MinPrecision = DiffExp`State`FEWorkingPrecision},
  a /. (DiffExp`State`AnalyticContinuationReplacements /.
    If[at >= 0,
      {DiffExp`Symbols`\[Theta]p -> 1, DiffExp`Symbols`\[Theta]m -> 0},
      {DiffExp`Symbols`\[Theta]p -> 0, DiffExp`Symbols`\[Theta]m -> 1}
    ]) /. DiffExp`Symbols`Logx -> Log[DiffExp`Symbols`x] /. DiffExp`Symbols`x -> at //
    SetPrecision[#, DiffExp`State`FEWorkingPrecision] & // Expand
];

(* SEval - evaluate series at a point *)
SEval[a_SeriesData, at_] := Block[{$MinPrecision = DiffExp`State`FEWorkingPrecision},
  If[DiffExp`State`FEC[UsePade] === True, GetPade[a], Normal[a]] /.
    (DiffExp`State`AnalyticContinuationReplacements /.
      If[at >= 0,
        {DiffExp`Symbols`\[Theta]p -> 1, DiffExp`Symbols`\[Theta]m -> 0},
        {DiffExp`Symbols`\[Theta]p -> 0, DiffExp`Symbols`\[Theta]m -> 1}
      ]) /. DiffExp`Symbols`Logx -> Log[DiffExp`Symbols`x] /. DiffExp`Symbols`x -> at //
    SetPrecision[#, DiffExp`State`FEWorkingPrecision] & // Expand
];
SEval[0, at_] := 0;
SEval[a_ /; NumericQ[a], at_] := a;
SEval[a_, at_] := SEval2[NormalizeEmbeddedSeries[a], at];

End[];

EndPackage[];
