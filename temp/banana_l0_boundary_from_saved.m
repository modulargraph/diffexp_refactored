(* Recompute the banana L0 boundary from a SAVED level-1 transport
   (FT_SAVED_TRANSPORT, default /tmp/ft_transport_save5/transport_level_1.m).
   Rebuilds ftData (FIRE reductions) but skips the expensive transport.
   Set DIFFEXP_DUMP_LAURENT_DIR to capture the combined Laurent integrand
   dump for per-segment dissection. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["FeynmanTrick/FeynmanTrick.m"], {General::shdw, Symbol::shdw}];
Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];

envOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]
];

wp = ToExpression[envOrDefault["FT_WORKING_PRECISION", "500"]];
epsOrder = ToExpression[envOrDefault["FT_EPS_ORDER", "0"]];
savedFile = envOrDefault["FT_SAVED_TRANSPORT",
  "/tmp/ft_transport_save5/transport_level_1.m"];

FeynmanTrick`SetFTOption["Threads", 1];
FeynmanTrick`SetFTOption["FThreads", 1];
FeynmanTrick`SetFTOption["Verbosity", 0];
FeynmanTrick`SetFTOption["WorkingPrecision", wp];
FeynmanTrick`SetFTOption["DimensionExpression", 2 - 2*FeynmanTrick`FTeps];
FeynmanTrick`SetFTOption["ReductionCache", False];
FeynmanTrick`SetFTOption["FIRETimeoutSeconds", 1800];

topology = FeynmanTrick`FIREInterface`DefineTopology[
  "step_banana",
  {Global`l1, Global`l2, Global`l3},
  {Global`p},
  {
    1 - Global`l1^2,
    1 - Global`l2^2,
    1 - Global`l3^2,
    1 - (-Global`l1 - Global`l2 - Global`l3 + Global`p)^2
  },
  {Global`p^2 -> -1}
];

ftData = FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
  topology, {{1, 2}, {1, 3}, {1, 4}}, {}
];
outputDir = FileNameJoin[{$TemporaryDirectory,
  "FT_l0_from_saved_" <> ToString[$ProcessID]}];
If[DirectoryQ[outputDir], DeleteDirectory[outputDir, DeleteContents -> True]];
CreateDirectory[outputDir, CreateIntermediateDirectories -> True];
ftData = FeynmanTrick`FeynmanTrickIteration`RunFullIteration[ftData, outputDir];
If[ftData === $Failed, Print["RunFullIteration failed"]; Exit[10]];

Print["Loading saved transport: ", savedFile];
saved = Get[savedFile];
transportResult = saved["TransportResult"];
Print["Saved level: ", saved["Level"], "  TransportOrder: ",
  saved["TransportOrder"]];

levelBoundary = FeynmanTrick`DiffExpIntegration`ComputeLevelBoundary[
  ftData, 0, transportResult, epsOrder
];
If[!AssociationQ[levelBoundary], Print["ComputeLevelBoundary failed"]; Exit[11]];

raw = levelBoundary["RawBoundaryValues"][[1]];
Print["L0 raw MinPower: ", raw["MinPower"]];
Do[
  Print["L0 eps^", raw["MinPower"] + k - 1, " = ",
    InputForm[N[raw["Coefficients"][[k]], 30]]],
  {k, Length[raw["Coefficients"]]}
];
Print["DONE"];
