(* ::Package:: *)
(* FeynmanTrick - Feynman Parameter Integration Package *)
(* Implements iterative Feynman trick with FIRE6 IBP reductions *)

BeginPackage["FeynmanTrick`"];

(* --- Public symbols --- *)
xx::usage = "Feynman parameter variable used in FeynmanTrick (avoids conflict with DiffExp's x).";
FTeps::usage = "Dimensional regulator symbol: d = 4 - 2*FTeps.";

FTConfiguration::usage = "FTConfiguration[] returns the current FeynmanTrick configuration.";
SetFTOption::usage = "SetFTOption[key, value] sets a configuration option.";

(* Load subpackages *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "PropagatorAlgebra.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "FIREInterface.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "MatrixExport.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "EpsPrefactors.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "FeynmanTrickIteration.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "BoundaryConditions.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "DiffExpIntegration.m"}]];

Begin["`Private`"];

(* --- Configuration --- *)
$FTConfig = <|
  "FIREPath" -> FileNameJoin[{
    ParentDirectory[DirectoryName[$InputFileName]],
    "Dependencies", "fire", "FIRE6"
  }],
  "WorkDirectory" -> FileNameJoin[{$TemporaryDirectory, "FeynmanTrick"}],
  "FeynmanParameter" -> FeynmanTrick`xx,
  "FixedParameterValue" -> 11/23,
  "DimensionVariable" -> Global`d,
  "EpsilonSymbol" -> FeynmanTrick`FTeps,
  "Threads" -> 4,
  "FThreads" -> 4,
  "WorkingPrecision" -> 500,
  "Verbosity" -> 1
|>;

FTConfiguration[] := $FTConfig;

SetFTOption[key_String, value_] := ($FTConfig[key] = value);

End[];
EndPackage[];

Print["FeynmanTrick package loaded."];
