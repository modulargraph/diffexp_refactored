(* DiffExp 2.1: a process wrapper; all mathematics runs in the C++ executable. *)
BeginPackage["DiffExp`"];
$DiffExpExecutable::usage = "$DiffExpExecutable is Automatic or the path to the diffexp executable.";
DiffExpRun::usage = "DiffExpRun[{arguments}, input] runs diffexp and returns its exit code, stdout and stderr. input defaults to an empty string. A failed process returns Failure.";
DiffExpSeries::usage = "DiffExpSeries[request] sends an Association describing a regular Taylor-series problem to the C++ backend and returns an Association of exact coefficient strings.";
DiffExpFeynmanTrick::usage = "DiffExpFeynmanTrick[family, {arguments}] evaluates a built-in family through the C++ Feynman-trick recursion and returns JSON data as an Association.";
DiffExpBackendInfo::usage = "DiffExpBackendInfo[] returns information about the C++ backend.";
Begin["`Private`"];
$wrapperDirectory = DirectoryName[$InputFileName];
If[!ValueQ[$DiffExpExecutable], $DiffExpExecutable = Automatic];
Options[DiffExpRun] = {"Executable" -> Automatic};
Options[DiffExpSeries] = Options[DiffExpRun];
Options[DiffExpFeynmanTrick] = Options[DiffExpRun];
Options[DiffExpBackendInfo] = Options[DiffExpRun];

executable[option_] := Module[{selected = option, environment, candidates},
  If[selected === Automatic, selected = $DiffExpExecutable];
  If[selected =!= Automatic, Return[selected]];
  environment = Environment["DIFFEXP_EXECUTABLE"];
  If[StringQ[environment] && environment =!= "", Return[environment]];
  candidates = {
    FileNameJoin[{$wrapperDirectory, "..", "build", "diffexp"}],
    FileNameJoin[{$wrapperDirectory, "..", "..", "..", "bin", "diffexp"}]
  };
  SelectFirst[candidates, FileExistsQ, "diffexp"]
];

DiffExpRun[arguments : {___String}, input_String : "", OptionsPattern[]] :=
 Module[{program = executable[OptionValue["Executable"]], result},
  If[!StringQ[program] || program === "",
    Return[Failure["Executable", <|"MessageTemplate" -> "Set a nonempty executable path."|>]]];
  (* RunProcess requires normalized executable paths on some platforms. *)
  If[FileExistsQ[program], program = ExpandFileName[program]];
  (* Omit stdin for empty input: a fast process may exit before a pipe write. *)
  result = Quiet[If[input === "",
    RunProcess[Prepend[arguments, program], All],
    RunProcess[Prepend[arguments, program], All, input]]];
  If[!AssociationQ[result], Return[Failure["ProcessLaunch", <|
    "MessageTemplate" -> "Could not launch the DiffExp executable.", "Executable" -> program|>]]];
  If[Lookup[result, "ExitCode", -1] =!= 0,
    Return[Failure["ProcessExit", Join[<|"MessageTemplate" -> "DiffExp reported an error."|>, result]]]];
  result
];

jsonRun[arguments_, input_, program_] := Module[{process, data},
  process = DiffExpRun[arguments, input, "Executable" -> program];
  If[FailureQ[process], Return[process]];
  data = Quiet[Check[ImportString[process["StandardOutput"], "RawJSON"], $Failed]];
  If[!AssociationQ[data], Return[Failure["JSONResponse", Join[<|
    "MessageTemplate" -> "DiffExp did not return a JSON object."|>, process]]]];
  data
];

DiffExpSeries[request_Association, OptionsPattern[]] := Module[{input},
  input = Quiet[Check[ExportString[request, "RawJSON", "Compact" -> True], $Failed]];
  If[!StringQ[input], Return[Failure["JSONRequest", <|
    "MessageTemplate" -> "The request must contain JSON-compatible values."|>]]];
  jsonRun[{"series", "-"}, input, OptionValue["Executable"]]
];
DiffExpFeynmanTrick[family_String, arguments : {___String} : {}, OptionsPattern[]] :=
  jsonRun[Join[{"ft", family, "--json"}, arguments], "", OptionValue["Executable"]];
DiffExpBackendInfo[OptionsPattern[]] :=
  jsonRun[{"backend-info"}, "", OptionValue["Executable"]];
End[];
EndPackage[];
