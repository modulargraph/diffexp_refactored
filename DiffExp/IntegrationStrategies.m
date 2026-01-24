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
SolveSimple::usage = "SolveSimple[ctx, bVec, epsord] solves a simple integration case where the integral has no homogeneous components. ctx is a SegmentContext association.";
SolveDefault::usage = "SolveDefault[ctx, bVec, epsord, cache] solves using the default Frobenius/Wronskian strategy. Returns {cIndices, fGeneral, cache}.";
SolveVOP::usage = "SolveVOP[ctx, bVec, epsord, cache] solves using the variation of parameters strategy. Returns {cIndices, fGeneral, cache}.";
SolveVOPAlt::usage = "SolveVOPAlt[ctx, bVec, epsord, cache] solves using the alternative variation of parameters strategy. Returns {cIndices, fGeneral, cache}.";
SolveRationalRecurrence::usage = "SolveRationalRecurrence[ctx, bVec, epsord, cache] solves using the rational recurrence method for non-singular points with rational matrices. Returns {cIndices, fGeneral, cache}.";
RationalRecurrenceApplicableQ::usage = "RationalRecurrenceApplicableQ[ctx] checks if the rational recurrence method is applicable for the given SegmentContext.";
SolveSingularRecurrence::usage = "SolveSingularRecurrence[ctx, bVec, epsord, cache] solves using the Frobenius recurrence method for regular singular points (simple pole) with non-resonant, diagonalizable residue. Returns {cIndices, fGeneral, cache}.";
SingularRecurrenceApplicableQ::usage = "SingularRecurrenceApplicableQ[ctx] checks if the singular recurrence method is applicable (simple pole, diagonalizable, non-resonant) for the given SegmentContext.";
DispatchStrategy::usage = "DispatchStrategy[ctx, bVec, epsord, cache] dispatches to the appropriate integration strategy based on the SegmentContext configuration. Returns {cIndices, fGeneral, cache}.";
SolveGeneralSingularRecurrence::usage = "SolveGeneralSingularRecurrence[ctx, bVec, epsord, cache] solves using the general singular recurrence method for regular singular points with resonant or non-diagonalizable residue. Returns {cIndices, fGeneral, cache}.";
GeneralSingularRecurrenceApplicableQ::usage = "GeneralSingularRecurrenceApplicableQ[ctx] checks if the general singular recurrence method is applicable (simple pole, any eigenvalue structure).";

Begin["`Private`"];

Get[FileNameJoin[{DirectoryName[$InputFileName], "IntegrationStrategies", "Helpers.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "IntegrationStrategies", "Default.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "IntegrationStrategies", "VOP.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "IntegrationStrategies", "Recurrence.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "IntegrationStrategies", "ResonantRecurrence.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "IntegrationStrategies", "Dispatch.m"}]];

End[];

EndPackage[];
