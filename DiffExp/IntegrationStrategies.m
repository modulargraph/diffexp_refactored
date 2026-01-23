(* ::Package:: *)

(* DiffExp Integration Strategies Subpackage *)
(* This package provides different integration strategy implementations *)
(* for solving systems of differential equations *)

BeginPackage["DiffExp`IntegrationStrategies`", {
  "DiffExp`Symbols`",
  "DiffExp`State`",
  "DiffExp`Utilities`",
  "DiffExp`SeriesOps`",
  "DiffExp`Integration`",
  "DiffExp`Frobenius`",
  "DiffExp`Wronskian`"
}];

(* Strategy functions *)
SolveSimple::usage = "SolveSimple[intind, bVec, line, epsord] solves a simple integration case where the integral has no homogeneous components.";
SolveDefault::usage = "SolveDefault[intind, bVec, line, epsord, BufferedData] solves using the default Frobenius/Wronskian strategy.";
SolveVOP::usage = "SolveVOP[intind, bVec, line, epsord, BufferedData] solves using the variation of parameters strategy.";
SolveVOPAlt::usage = "SolveVOPAlt[intind, bVec, line, epsord, BufferedData] solves using the alternative variation of parameters strategy.";
SolveRationalRecurrence::usage = "SolveRationalRecurrence[intind, bVec, line, epsord, BufferedData] solves using the rational recurrence method for non-singular points with rational matrices.";
RationalRecurrenceApplicableQ::usage = "RationalRecurrenceApplicableQ[intind, line] checks if the rational recurrence method is applicable.";
SolveSingularRecurrence::usage = "SolveSingularRecurrence[intind, bVec, line, epsord, BufferedData] solves using the Frobenius recurrence method for regular singular points (simple pole) with non-resonant, diagonalizable residue.";
SingularRecurrenceApplicableQ::usage = "SingularRecurrenceApplicableQ[intind, line] checks if the singular recurrence method is applicable (simple pole, diagonalizable, non-resonant).";
DispatchStrategy::usage = "DispatchStrategy[intind, bVec, line, epsord, BufferedData] dispatches to the appropriate integration strategy based on configuration.";

Begin["`Private`"];

Get[FileNameJoin[{DirectoryName[$InputFileName], "IntegrationStrategies", "Helpers.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "IntegrationStrategies", "Default.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "IntegrationStrategies", "VOP.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "IntegrationStrategies", "Recurrence.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "IntegrationStrategies", "Dispatch.m"}]];

End[];

EndPackage[];
