(* ::Package:: *)

(* DiffExp Integration Subpackage *)
(* This package provides series integration functions *)

BeginPackage["DiffExp`Integration`", {"DiffExp`Symbols`", "DiffExp`State`", "DiffExp`Utilities`", "DiffExp`SeriesOps`"}];

(* Integration functions *)
DiffExpIntegrate::usage = "DiffExpIntegrate[a__] integrates series expressions.";
DiffExpIntegrate1::usage = "DiffExpIntegrate1[a_] internal integration function.";
UpdateIntReps::usage = "UpdateIntReps[MaxOrd_] updates integration replacement rules.";
IntReps::usage = "Integration replacement rules.";

Begin["`Private`"];

(* Set attributes as in original code *)
SetAttributes[DiffExpIntegrate1, Listable];

(* Integration replacement rules - starts empty, gets populated by UpdateIntReps *)
IntReps = {};

(* Main integration function *)
DiffExpIntegrate[a__] := Block[{LogOrd, rawLogPowers, logPowers},
  (* If maximum power of Logx has increased, update the replacement relations. *)
  rawLogPowers = Append[
    DiffExp`Utilities`GetCases[{a} // DiffExp`SeriesOps`SExpand, DiffExp`Symbols`Logx^(k_: 1) | Log[DiffExp`Symbols`x]^(k_: 1) :> k],
    1
  ];
  logPowers = Cases[
    DiffExp`SeriesOps`NormalizeLogPower /@ Flatten[{rawLogPowers}],
    _Integer?NonNegative
  ];
  LogOrd = If[logPowers === {}, 1, Max[logPowers]];

  If[LogOrd > DiffExp`State`IMaxLogOrder,
    DiffExp`Utilities`PrintDebug["Encountered ", Log[DiffExp`Symbols`x]^LogOrd, ". Updating IntReps."][2];
    UpdateIntReps[LogOrd];
  ];

  DiffExpIntegrate1[a]
];

(* Internal integration function for single argument *)
DiffExpIntegrate1[a_] := DiffExpIntegrate1[a, DiffExp`Symbols`x];

(* Integration for SeriesData *)
DiffExpIntegrate1[exp0_SeriesData, var_] := Block[
  {exp = (Expand@Normal@exp0) /. DiffExp`Symbols`Logx -> Log[DiffExp`Symbols`x], Out, MaxPow, Const},
  Out = exp /. var -> b /. IntReps;
  Const = (Out /. a -> 0);
  MaxPow = Floor[(exp0[[5]]/exp0[[6]])];
  DiffExp`SeriesOps`SeriesAlways[Out - Const + b Const /. a -> 1 /. b -> var /. Log[DiffExp`Symbols`x] -> DiffExp`Symbols`Logx, {DiffExp`Symbols`x, 0, MaxPow}] // DiffExp`SeriesOps`SExpand
];

(* Integration for numeric values *)
DiffExpIntegrate1[exp0_ /; NumericQ[exp0], var_] := DiffExp`SeriesOps`SeriesAlways[var exp0, {DiffExp`Symbols`x, 0, DiffExp`State`ExpansionOrderVal}];

(* Error case *)
DiffExpIntegrate1[exp0_, var_] := DiffExp`Utilities`ReportError["DiffExpIntegrate called with unsupported arguments. ", exp0, " ", var];

(* Update integration replacement rules for given max log order *)
UpdateIntReps[MaxOrd_] := Block[{},
  DiffExp`State`IMaxLogOrder = MaxOrd;

  IntReps = Join[
    Table[Log[DiffExp`Symbols`x]^n -> a Integrate[Log[DiffExp`Symbols`x]^n, DiffExp`Symbols`x] /. DiffExp`Symbols`x -> b, {n, DiffExp`State`IMaxLogOrder}],
    Table[Log[DiffExp`Symbols`x]^n DiffExp`Symbols`x -> a Integrate[Log[DiffExp`Symbols`x]^n DiffExp`Symbols`x, DiffExp`Symbols`x] /. DiffExp`Symbols`x -> b, {n, DiffExp`State`IMaxLogOrder}],
    Table[Log[DiffExp`Symbols`x]^n DiffExp`Symbols`x^m_ /; m != -1 -> a Integrate[Log[DiffExp`Symbols`x]^n DiffExp`Symbols`x^m, DiffExp`Symbols`x] /. DiffExp`Symbols`x -> b, {n, DiffExp`State`IMaxLogOrder}],
    Table[Log[DiffExp`Symbols`x]^n/DiffExp`Symbols`x -> a Integrate[Log[DiffExp`Symbols`x]^n DiffExp`Symbols`x^-1, DiffExp`Symbols`x] /. DiffExp`Symbols`x -> b, {n, DiffExp`State`IMaxLogOrder}],
    {DiffExp`Symbols`x^m_ /; m != -1 -> a Integrate[DiffExp`Symbols`x^m, DiffExp`Symbols`x] /. DiffExp`Symbols`x -> b},
    {DiffExp`Symbols`x -> a Integrate[DiffExp`Symbols`x, DiffExp`Symbols`x] /. DiffExp`Symbols`x -> b},
    {1/DiffExp`Symbols`x -> a Integrate[1/DiffExp`Symbols`x, DiffExp`Symbols`x] /. DiffExp`Symbols`x -> b}
  ] // Reverse // Expand;
];

End[];

EndPackage[];
