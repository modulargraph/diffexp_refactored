Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "DiffExp.m"}]];
family = DiffExpFamilyTemplate["sunrise"];
family["name"] = "unregistered-native-ibp-sunrise";
cache = CreateDirectory[];
result = DiffExpFeynmanTrick[family, {"--ibp-provider", "ibp-solver", "--cache", cache, "--no-numerical-cache"}];
passed = AssociationQ[result] && result["status"] === "completed" &&
 result["ibp_provider"] === "ibp-solver" && result["ibp_statistics"]["fresh_probes"] > 0 &&
 StringContainsQ[Last[First[result["coefficients"]]]["real"], "2.2367927002126465"];
If[passed, Print["Mathematica native IBP sunrise passed; timings: ", result["timings"]], Print[result]];
DeleteDirectory[cache, DeleteContents -> True];
Exit[If[TrueQ[passed], 0, 1]];
