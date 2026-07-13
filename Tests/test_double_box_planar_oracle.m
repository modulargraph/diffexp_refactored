(* Structural and synthetic-output checks for the planar double-box oracle.
   No FIRE preparation or differential-equation solve occurs. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
SetEnvironment["DOUBLE_BOX_ORACLE_DEFINITIONS_ONLY" -> "1"];
Get[FileNameJoin[{repoRoot, "Scripts", "verify_double_box_planar.m"}]];

passed = 0; failed = 0;
assert[label_, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

reference = dbDoubleBoxReference[];
exactRow = <|
  "Example" -> "double_box_planar", "Level" -> 0,
  "Master" -> dbDoubleBoxMaster[], "EpsPrefactor" -> 0,
  "RawMinPower" -> -4,
  "Coefficients" -> ({#, reference[#]} & /@ Keys[reference])|>;
exactOutput = StringRiffle[{
  "unrelated diagnostic",
  "FINAL " <> ExportString[
    <|"Example" -> "double_box_planar", "Finite" -> 999|>,
    "RawJSON", "Compact" -> True],
  "STEPWISE " <> ExportString[exactRow, "RawJSON", "Compact" -> True]},
  "\n"];

exactResult = dbVerifyDoubleBoxOutput[exactOutput, 40];
finiteOnlyResult = dbVerifyDoubleBoxOutput[
  "FINAL " <> ExportString[
    <|"Example" -> "double_box_planar", "Finite" -> reference[0]|>,
    "RawJSON", "Compact" -> True], 20];
truncatedResult = dbVerifyDoubleBoxRow[
  Join[exactRow, <|"RawMinPower" -> -3,
    "Coefficients" -> Rest[exactRow["Coefficients"]]|>], 20];
wrongValueRow = exactRow;
wrongCoefficients = exactRow["Coefficients"];
wrongCoefficients[[3, 2]] += 1/1000;
wrongValueRow = Join[exactRow, <|"Coefficients" -> wrongCoefficients|>];
wrongValueResult = dbVerifyDoubleBoxRow[wrongValueRow, 20];
fileArguments = dbScriptArguments[{
  "wolframscript", "-file", "Scripts/verify_double_box_planar.m",
  "positional.log", "37"}];
directArguments = dbScriptArguments[{
  "wolframscript", "Scripts/verify_double_box_planar.m",
  "direct.log", "29"}];
positionalInputs = dbResolveVerifierInputs[
  fileArguments, "environment.log", "17"];
environmentInputs = dbResolveVerifierInputs[{}, "environment.log", "23"];
defaultDigitsInputs = dbResolveVerifierInputs[{"one.log"}, "", ""];
badDigitInputs = dbResolveVerifierInputs[
  {"bad.log", "20;Quit[]"}, "", ""];
tooManyInputs = dbResolveVerifierInputs[
  {"bad.log", "20", "unexpected"}, "", ""];

assert["double_box_oracle_definitions_only",
  $dbDoubleBoxVerifierRan === False];
assert["double_box_oracle_reference_window_and_normalization_pin",
  Keys[reference] === Range[-4, 0] && reference[-4] == 12 &&
  Round[10^40 reference[0]] ===
    -2470036782205925670094377082844956237592781];
assert["double_box_oracle_accepts_full_terminal_stepwise_row",
  exactResult["Status"] === "PASS" &&
  exactResult["Powers"] === Range[-4, 0] &&
  PossibleZeroQ[exactResult["MaximumScaledError"]]];
assert["double_box_oracle_does_not_treat_final_as_full_result",
  finiteOnlyResult["Status"] === "FAIL" &&
  finiteOnlyResult["MatchingRows"] === 0];
assert["double_box_oracle_rejects_truncated_laurent_window",
  truncatedResult["Status"] === "FAIL"];
assert["double_box_oracle_rejects_wrong_laurent_coefficient",
  wrongValueResult["Status"] === "FAIL" &&
  !wrongValueResult["Comparisons"][[3, "Pass"]]];
assert["double_box_oracle_extracts_file_and_direct_positional_arguments",
  fileArguments === {"positional.log", "37"} &&
  directArguments === {"direct.log", "29"}];
assert["double_box_oracle_positional_path_and_digits_override_environment",
  positionalInputs === <|"Path" -> "positional.log", "Digits" -> 37|>];
assert["double_box_oracle_environment_and_default_inputs_remain_supported",
  environmentInputs === <|"Path" -> "environment.log", "Digits" -> 23|> &&
  defaultDigitsInputs === <|"Path" -> "one.log", "Digits" -> 20|>];
assert["double_box_oracle_rejects_nonnumeric_or_extra_arguments_without_evaluation",
  FailureQ[badDigitInputs] && FailureQ[tooManyInputs]];

SetEnvironment["DOUBLE_BOX_ORACLE_DEFINITIONS_ONLY" -> None];
Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Quit[1], Quit[0]];
