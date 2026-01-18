(* ::Package:: *)

(* DiffExp Utilities Subpackage *)
(* This package provides small helper functions used throughout DiffExp *)

BeginPackage["DiffExp`Utilities`", {"DiffExp`Symbols`", "DiffExp`State`"}];

(* Printing functions *)
PrintDebug::usage = "PrintDebug[args__][lev_] prints debug messages at specified verbosity level.";
PrintInfo::usage = "PrintInfo[args__][lev_] prints info messages at specified verbosity level.";
PrintWarning::usage = "PrintWarning[args__] prints warning messages.";
ReportError::usage = "ReportError[mes__] reports an error and aborts.";

(* General utilities *)
AllSameQ::usage = "AllSameQ[l_,b_] checks if all elements in l are equal to b.";
CA::usage = "CA is an alias for ConstantArray.";
GetCases::usage = "GetCases[expr_,case_] extracts matching cases from expression.";
DependsQ::usage = "DependsQ[a_,b_] checks if a depends on b.";
ZeroQ::usage = "ZeroQ[a_] checks if a is exactly zero.";
R::usage = "R is an alias for ReplaceAll.";
FirstOrNull::usage = "FirstOrNull[l_] returns first element or Null if empty.";
FindPivots::usage = "FindPivots[Matrix_] finds pivot positions in a matrix.";
SplitTimes::usage = "SplitTimes[Expr_] splits a Times expression into a list.";
SplitSum::usage = "SplitSum[Expr_] splits a Plus expression into a list.";

(* Chopping functions *)
PChop::usage = "PChop chops numbers below ChopPrecision.";
LSPChop::usage = "LSPChop chops numbers below LinearSolveChopPrecision.";
CPChop::usage = "CPChop chops numbers below CrosscheckChopPrecision.";

(* Line/Point detection *)
IsPoint::usage = "IsPoint[line_] checks if line is actually a point (no x dependence).";
IsLine::usage = "IsLine[line_] checks if line depends on x.";

(* Interval functions *)
IntervalOverlapQ::usage = "IntervalOverlapQ[intv1_,intv2_] checks if intervals overlap.";
IntervalIntersec::usage = "IntervalIntersec[intv1_,intv2_] returns intersection of intervals.";
IntervalContainsQ::usage = "IntervalContainsQ[intv_,point_] checks if interval contains point.";

(* Line type checking *)
ExactLineQ::usage = "ExactLineQ[line_Association] checks if line has exact coefficients.";
FactorOrTogether::usage = "FactorOrTogether[line_Association] returns Factor or Together based on line type.";

Begin["`Private`"];

(* Import state variables *)
(* These are accessed via DiffExp`State` context *)

(* Printing functions *)
PrintDebug[args__][lev_] := If[lev <= DiffExp`State`FEC["VerbosityDebug"],
  Print[Sequence @@ ConstantArray["   ", Max[lev - 1, 0]], "DiffExp Debug: ", args]
];

PrintInfo[args__][lev_] := If[lev <= DiffExp`State`FEC[Verbosity],
  Print[Sequence @@ ConstantArray["   ", Max[lev - 1, 0]], "DiffExp: ", args]
];

PrintWarning[args__] := Print["DiffExp Warning: ", args];

ReportError[mes__] := (Print["DiffExp Error: ", mes]; Abort[];);
ReportError[mes__, False] := Print["DiffExp Error: ", mes];

(* General utilities *)
AllSameQ[l_, b_] := If[Length[l] > 0, SameQ[l] && l[[1]] === b, True];

CA = ConstantArray;

GetCases[expr_, case_] := expr // Cases[{#}, case, Infinity] & // DeleteDuplicates // Sort;

DependsQ[a_, b_] := Length[GetCases[a, b]] > 0;

ZeroQ[a_] := a === 0;

R = ReplaceAll;

FirstOrNull[l_] := If[l === {}, Null, First[l]];

FindPivots[Matrix_] := DeleteCases[
  Flatten[FirstOrNull[Position[#, a_ /; a != 0]] & /@ Matrix],
  Null
];

SetAttributes[SplitTimes, Listable];
SplitTimes[Expr_] := If[Expr[[0]] === Times, List @@ Expr, {Expr}];

SplitSum[Expr_] := If[Expr[[0]] === Plus, List @@ Expr, {Expr}];

(* Chopping functions - using state variables *)
PChop := Chop[#, 10^-DiffExp`State`FEC[ChopPrecision]] &;
LSPChop := Chop[#, 10^-DiffExp`State`FEC["LinearSolveChopPrecision"]] &;
CPChop := Chop[#, 10^-DiffExp`State`CrosscheckChopPrecision] &;

(* Line/Point detection - uses x from Symbols *)
IsPoint[line_] := !DependsQ[line // Values, DiffExp`Symbols`x];
IsLine[line_] := !IsPoint[line];

(* Interval functions *)
IntervalOverlapQ[intv1_, intv2_] := !(IntervalIntersection[Interval[intv1], Interval[intv2]] === Interval[]);
IntervalIntersec[intv1_, intv2_] := IntervalIntersection[Interval[intv1], Interval[intv2]][[1]];
IntervalContainsQ[intv_, point_] := intv[[1]] <= point <= intv[[2]];

(* Line type checking *)
ExactLineQ[line_Association] := And @@ (ExactNumberQ /@ Flatten[({
  Numerator[#] /. DiffExp`Symbols`x -> 0,
  Coefficient[Numerator[#], DiffExp`Symbols`x],
  Denominator[#] /. DiffExp`Symbols`x -> 0,
  Coefficient[Denominator[#], DiffExp`Symbols`x]
} & /@ (Together /@ Values[line]))]);

FactorOrTogether[line_Association] := If[ExactLineQ[line], Factor, Together];
FactorOrTogether[line1_Association, line2_Association] := If[ExactLineQ[line1] && ExactLineQ[line2], Factor, Together];

End[];

EndPackage[];
