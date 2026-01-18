(* ::Package:: *)

(* DiffExp Symbols Subpackage *)
(* This package defines the core symbols used throughout DiffExp *)
(* All other subpackages should load this package to ensure consistent symbol usage *)

BeginPackage["DiffExp`Symbols`"];

(* Export symbols *)
\[Epsilon]::usage = "Symbol representing the dimensional regulator.";
eps::usage = "Alternative symbol for the dimensional regulator, which is converted into \[Epsilon].";
Logx::usage = "Symbol representing Log[x].";
\[Theta]p::usage = "Symbol representing HeavisideTheta[x].";
\[Theta]m::usage = "Symbol representing HeavisideTheta[-x].";

(* Internal x symbol - will be set to the line parameter *)
x::usage = "Internal line parameter symbol.";

Begin["`Private`"];

(* eps is an alias for \[Epsilon] - this is a delayed rule so it always returns the current value *)
eps := \[Epsilon];

End[];

EndPackage[];
