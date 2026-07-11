(* ::Package:: *)
(* FeynmanTrick - Feynman Parameter Integration Package *)
(* Implements iterative Feynman trick with FIRE6 IBP reductions *)

BeginPackage["FeynmanTrick`"];

(* --- Public symbols --- *)
xx::usage = "Feynman parameter variable used in FeynmanTrick (avoids conflict with DiffExp's x).";
FTeps::usage = "Dimensional regulator symbol. By default d = 4 - 2*FTeps.";

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
  "DimensionExpression" -> Automatic,
  "EpsilonSymbol" -> FeynmanTrick`FTeps,
  "Threads" -> 1,        (* FIRE6 multi-threading can be unstable; default to 1 *)
  "FThreads" -> 1,
  "FIRETimeoutSeconds" -> 600,  (* watchdog timeout for a single FIRE6 run *)
  "AutoDetectRestrictions" -> False,
  "ReductionCache" -> True,
  (* Experimental, opt-in persistence of FIRE's internal sector databases
     across successive reduction requests for one exact topology.  The
     implementation uses content-addressed immutable generations; Automatic
     keeps them below the package temporary root. *)
  "PersistentFIREStorage" -> False,
  "FIREStorageDirectory" -> Automatic,
  "WorkingPrecision" -> 500,
  (* Extra epsilon orders requested from each level transport beyond the
     orders needed by the IBP coefficients alone.  Regularized endpoint
     integration consumes lookahead: each x^(a + b eps) endpoint pole
     deepens the Laurent window by one order, and recovering several
     x^(a + b_i eps) sectors that DecomposeSingularity collapsed into one
     average exponent needs sectorCount extra orders for the residual
     moment solve plus one guard order for the truncation boundary. *)
  "IntegrationPoleAllowance" -> 4,
  "Verbosity" -> 1
|>;

FTConfiguration[] := $FTConfig;

SetFTOption[key_String, value_] := ($FTConfig[key] = value);

DimensionExpression[] := Module[
  {expr = Lookup[$FTConfig, "DimensionExpression", Automatic],
   eps = $FTConfig["EpsilonSymbol"]},
  If[expr === Automatic, 4 - 2*eps, expr]
];

End[];
EndPackage[];

Print["FeynmanTrick package loaded."];
