Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "DiffExp.m"}]];
failures = {};
check[condition_, name_] := If[!TrueQ[condition], AppendTo[failures, name]];
request = <|"matrix" -> {{"eps/(1-x)"}}, "boundary" -> {{"1", "0", "0"}},
  "center" -> "0", "taylor_order" -> 8, "epsilon_low" -> 0, "epsilon_high" -> 2|>;
result = DiffExpSeries[request];
check[AssociationQ[result], "series response"];
If[AssociationQ[result],
  check[result["coefficients"][[4, 1, 2]] === "1/3", "exact logarithm coefficient"];
  check[result["coefficients"][[3, 1, 3]] === "1/2", "epsilon convolution"];
  check[result["infinite_tail_certified"] === False, "tail status"]];
check[DiffExpBackendInfo[]["standalone_cpp"] === True, "standalone backend"];
check[StringContainsQ[DiffExpRun[{"--version"}]["StandardOutput"], "2.1.1"], "version"];
testProgram = Environment["DIFFEXP_EXECUTABLE"];
If[!StringQ[testProgram] || testProgram === "", testProgram = ExpandFileName[
  FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "build", "diffexp"}]]];
relativeProgram = FileNameJoin[{DirectoryName[testProgram], "..",
  FileNameTake[DirectoryName[testProgram]], FileNameTake[testProgram]}];
check[AssociationQ[DiffExpRun[{"--version"}, "Executable" -> relativeProgram]], "normalized executable path"];
check[And @@ Table[AssociationQ[DiffExpRun[{"--version"}]], {5}], "fast process without stdin"];
check[FailureQ[DiffExpRun[{"unknown-command"}]], "nonzero process exit"];
check[FailureQ[DiffExpRun[{"--version"}, "Executable" -> "/nonexistent/diffexp"]], "missing executable"];
check[FailureQ[DiffExpSeries[Join[request, <|"center" -> "1"|>]]], "singular ordinary center"];
check[FailureQ[DiffExpFeynmanTrick["not-a-family"]], "FT failure propagation"];
(* Spaces and shell metacharacters must remain a single literal argument. *)
path = FileNameJoin[{$TemporaryDirectory, "diffexp test ; literal " <> CreateUUID[] <> ".json"}];
Export[path, request, "RawJSON"];
process = DiffExpRun[{"series", path}];
DeleteFile[path];
check[AssociationQ[process] && process["ExitCode"] === 0, "literal argument path"];
If[failures === {}, Print["Mathematica RunProcess interface: all checks passed"]; Exit[0], Print[failures]; Exit[1]];
