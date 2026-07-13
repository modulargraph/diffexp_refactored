(* Numerical oracle for the terminal double_box_planar STEPWISE row.

   The runner's FINAL row is intentionally only a presentation of the finite
   coefficient.  This verifier instead checks the complete level-0 Laurent
   list, eps^-4 through eps^0.

   Usage:
     wolframscript -file Scripts/verify_double_box_planar.m run.log [digits]

   DOUBLE_BOX_RESULT_LOG may be used instead of the positional path, and
   DOUBLE_BOX_VERIFY_DIGITS may be used instead of positional digits.  Explicit
   positional values take precedence over their environment counterparts.

   Normalization: with x=t/s=1/3 at s=-1, t=-1/3, the conventions of the
   FeynmanTrick fixture and Smirnov's planar-double-box K-function are related
   by

       I_FT(eps) = -3 Exp[-2 EulerGamma eps] K(1/3, eps).

   Expanding that expression gives the constants pinned below. *)

dbEnvOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]];

dbDoubleBoxReference[] := <|
  -4 -> 12.000000000000000000000000000000000000000000000000,
  -3 -> 2.6260083723848567163723883918602272247005343047831,
  -2 -> -77.808155752086180303874072147668656652688133264584,
  -1 -> -201.73675884525049357623701386093884968969770043649,
   0 -> -247.00367822059256700943770828449562375927814454946|>;

dbDoubleBoxMaster[] := ConstantArray[1, 7];

dbJsonRows[output_String, prefix_String] := DeleteCases[
  Map[Function[line, If[StringStartsQ[line, prefix],
      Quiet[Check[ImportString[StringDrop[line, StringLength[prefix]],
        "RawJSON"], Missing["InvalidJSON", line]]], Nothing]],
    StringSplit[output, {"\r\n", "\n", "\r"}]],
  _Missing];

dbNumericValue[value_] := Which[
  NumericQ[value], value,
  AssociationQ[value] &&
      And @@ (KeyExistsQ[value, #] & /@ {"Re", "Im"}) &&
      NumericQ[value["Re"]] && NumericQ[value["Im"]],
    value["Re"] + I value["Im"],
  True, $Failed];

dbFailure[detail_String, extra_:<||>] := Join[
  <|"Example" -> "double_box_planar", "Status" -> "FAIL",
    "Detail" -> detail|>, extra];

dbVerifyDoubleBoxRow[row_Association, digits_Integer:20] := Module[
  {reference, expectedPowers, coefficients, powers, values, comparisons,
    scaledErrors},
  If[digits < 1, Return[dbFailure["target digits must be positive"], Module]];
  If[Lookup[row, "Example", None] =!= "double_box_planar" ||
      Lookup[row, "Level", None] =!= 0,
    Return[dbFailure[
      "the selected STEPWISE row is not double_box_planar level 0"], Module]];
  If[Lookup[row, "Master", None] =!= dbDoubleBoxMaster[],
    Return[dbFailure["the level-0 row is not the seven-line scalar master",
      <|"ObservedMaster" -> Lookup[row, "Master", Missing["Absent"]]|>],
      Module]];
  reference = dbDoubleBoxReference[];
  expectedPowers = Keys[reference];
  coefficients = Lookup[row, "Coefficients", Missing["Absent"]];
  If[!ListQ[coefficients] ||
      !AllTrue[coefficients, MatchQ[#, {_, _}] &],
    Return[dbFailure["Coefficients is not a list of {power,value} pairs"],
      Module]];
  powers = coefficients[[All, 1]];
  If[powers =!= expectedPowers ||
      Lookup[row, "RawMinPower", None] =!= First[expectedPowers],
    Return[dbFailure[
      "the terminal Laurent list must contain exactly eps^-4 through eps^0",
      <|"ObservedPowers" -> powers,
        "ObservedRawMinPower" -> Lookup[row, "RawMinPower", Missing["Absent"]]|>],
      Module]];
  values = dbNumericValue /@ coefficients[[All, 2]];
  If[MemberQ[values, $Failed],
    Return[dbFailure["at least one Laurent coefficient is nonnumeric"],
      Module]];
  comparisons = MapThread[Function[{power, observed}, Module[
      {expected = reference[power], difference, scaledError},
      difference = Abs[observed - expected];
      scaledError = difference/Max[1, Abs[expected]];
      <|"Power" -> power, "Expected" -> expected,
        "Observed" -> If[PossibleZeroQ[Im[observed]], Re[observed],
          <|"Re" -> Re[observed], "Im" -> Im[observed]|>],
        "AbsoluteDifference" -> difference,
        "ScaledError" -> scaledError,
        "Pass" -> TrueQ[scaledError < 10^-digits]|>]],
    {expectedPowers, values}];
  scaledErrors = comparisons[[All, "ScaledError"]];
  <|"Example" -> "double_box_planar",
    "Status" -> If[And @@ comparisons[[All, "Pass"]], "PASS", "FAIL"],
    "TargetDigits" -> digits, "Powers" -> expectedPowers,
    "MaximumScaledError" -> Max[scaledErrors],
    "Comparisons" -> comparisons|>];

dbVerifyDoubleBoxOutput[output_String, digits_Integer:20] := Module[
  {rows, candidates},
  rows = dbJsonRows[output, "STEPWISE "];
  candidates = Select[rows,
    AssociationQ[#] && Lookup[#, "Example", None] === "double_box_planar" &&
      Lookup[#, "Level", None] === 0 &];
  If[Length[candidates] =!= 1,
    Return[dbFailure[
      "expected exactly one double_box_planar level-0 STEPWISE row; FINAL is not a Laurent-series substitute",
      <|"MatchingRows" -> Length[candidates]|>], Module]];
  dbVerifyDoubleBoxRow[First[candidates], digits]];

dbScriptArguments[commandLine_List] := Module[
  {filePosition, scriptPosition},
  filePosition = FirstPosition[commandLine, "-file", Missing[]];
  If[!MissingQ[filePosition],
    (* Drop the executable, -file, and the script path itself. *)
    Return[Drop[commandLine, First[filePosition] + 1], Module]];
  (* Also support direct executable/script command lines, whose runtime does
     not retain the -file marker in $ScriptCommandLine. *)
  scriptPosition = FirstPosition[commandLine,
    value_String /; FileNameTake[value] === "verify_double_box_planar.m",
    Missing[]];
  If[MissingQ[scriptPosition], {}, Drop[commandLine, First[scriptPosition]]]];

dbParsePositiveInteger[text_] := Module[{trim},
  If[!StringQ[text], Return[$Failed, Module]];
  trim = StringTrim[text];
  If[!StringMatchQ[trim, DigitCharacter ..], Return[$Failed, Module]];
  FromDigits[trim]];

dbResolveVerifierInputs[args_List, environmentPath_String,
    environmentDigits_String] := Module[{path, digitsText, digits},
  If[Length[args] > 2 || !AllTrue[args, StringQ],
    Return[Failure["DoubleBoxVerifierArguments", <|
      "Detail" -> "expected result-log path and optional positive digits",
      "Arguments" -> args|>], Module]];
  path = If[Length[args] >= 1 && StringLength[StringTrim[args[[1]]]] > 0,
    StringTrim[args[[1]]], StringTrim[environmentPath]];
  digitsText = If[Length[args] >= 2, args[[2]],
    If[StringLength[StringTrim[environmentDigits]] > 0,
      environmentDigits, "20"]];
  digits = dbParsePositiveInteger[digitsText];
  If[digits === $Failed || digits < 1,
    Return[Failure["DoubleBoxVerifierArguments", <|
      "Detail" -> "verification digits must be a positive integer",
      "Value" -> digitsText|>], Module]];
  <|"Path" -> path, "Digits" -> digits|>];

dbVerifierInputs[] := dbResolveVerifierInputs[
  dbScriptArguments[$ScriptCommandLine],
  dbEnvOrDefault["DOUBLE_BOX_RESULT_LOG", ""],
  dbEnvOrDefault["DOUBLE_BOX_VERIFY_DIGITS", "20"]];

runDoubleBoxVerifier[] := Module[{inputs, path, digits, output, result},
  inputs = dbVerifierInputs[];
  If[FailureQ[inputs],
    Print[inputs[[2, "Detail"]]];
    Print["Usage: wolframscript -file Scripts/verify_double_box_planar.m run.log [digits]"];
    Return[2, Module]];
  path = inputs["Path"];
  digits = inputs["Digits"];
  If[path === "" || !FileExistsQ[path],
    Print["Usage: wolframscript -file Scripts/verify_double_box_planar.m run.log [digits]"];
    Return[2, Module]];
  output = Quiet[Check[Import[path, "Text"], $Failed]];
  If[!StringQ[output], Print["Could not read result log: ", path];
    Return[2, Module]];
  result = dbVerifyDoubleBoxOutput[output, digits];
  Print["DOUBLE_BOX_ORACLE ",
    ExportString[result, "RawJSON", "Compact" -> True]];
  If[result["Status"] === "PASS", 0, 1]];

$dbDoubleBoxVerifierRan = False;
If[dbEnvOrDefault["DOUBLE_BOX_ORACLE_DEFINITIONS_ONLY", "0"] =!= "1",
  $dbDoubleBoxVerifierRan = True;
  Exit[runDoubleBoxVerifier[]]];
