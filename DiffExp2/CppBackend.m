(* DiffExp2/CppBackend.m -- coarse-grained LibraryLink bridge to the
   FLINT/Arb framed recurrence kernel.  The bridge intentionally accepts
   only already-prepared numeric, rational, or rational-function coefficient
   tensors.  All
   indicial, resonance, epsilon-window, and branch decisions stay in
   Solve.m and are serialized explicitly. *)

BeginPackage["DiffExp2`CppBackend`",
  {"DiffExp2`Tolerances`", "DiffExp2`Config`"}];

BackendAvailableQ::usage = "BackendAvailableQ[] reports whether the compiled DiffExp2 LibraryLink recurrence backend can be loaded.";
BackendInformation::usage = "BackendInformation[] returns the compiled backend version record, or a Failure if it cannot be loaded.";
EncodeScalar::usage = "EncodeScalar[z, digits] converts a numeric scalar to the canonical C++ JSON pair of Arb scalar strings, retaining inexact input uncertainty. Non-numeric analytic regulators are rejected.";
EncodeSymbolicScalar::usage = "EncodeSymbolicScalar[z, vars] converts an exact rational function of named analytic regulators to the FLINT symbolic coefficient-field syntax.";
RunRequest::usage = "RunRequest[jsonReadyAssociation] executes one coarse-grained compiled recurrence request and returns its decoded JSON response.";
DecodeScalar::usage = "DecodeScalar[encoded, precision] reconstructs a Wolfram scalar from a C++ rational or Acb midpoint/radius record with honest Accuracy.";
ResetBackend::usage = "ResetBackend[] unloads the cached LibraryFunction handles so a rebuilt library can be loaded.";

Begin["`Private`"];

$moduleDirectory = DirectoryName[$InputFileName];
$repositoryRoot = ParentDirectory[$moduleDirectory];
$backendLibrary = None;
$runFunction = None;
$infoFunction = None;

libraryCandidates[] := DeleteDuplicates[Select[{
    Quiet[Environment["DE2_CPP_LIBRARY"]],
    FileNameJoin[{$repositoryRoot, "build", "cpp",
      "diffexp2_librarylink." <> Switch[$OperatingSystem,
        "MacOSX", "dylib", "Windows", "dll", _, "so"]}]
  }, StringQ[#] && StringLength[StringTrim[#]] > 0 &]];

loadBackend[] := Module[{lib},
  If[Head[$runFunction] === LibraryFunction &&
      Head[$infoFunction] === LibraryFunction, Return[True, Module]];
  lib = SelectFirst[libraryCandidates[], FileExistsQ, None];
  If[lib === None, Return[False, Module]];
  Quiet[Check[
    $runFunction = LibraryFunctionLoad[lib, "de2RunRecurrence",
      {"UTF8String"}, LibraryDataType[NumericArray, "UnsignedInteger8"]];
    $infoFunction = LibraryFunctionLoad[lib, "de2BackendInfo", {},
      LibraryDataType[NumericArray, "UnsignedInteger8"]];
    $backendLibrary = lib;
    True,
    $runFunction = None; $infoFunction = None; $backendLibrary = None;
    False]]];

ResetBackend[] := Module[{},
  If[Head[$runFunction] === LibraryFunction, Quiet[LibraryFunctionUnload[$runFunction]]];
  If[Head[$infoFunction] === LibraryFunction, Quiet[LibraryFunctionUnload[$infoFunction]]];
  $runFunction = None; $infoFunction = None; $backendLibrary = None; Null];

BackendAvailableQ[] := TrueQ[loadBackend[]];

bytesToJSON[bytes_NumericArray] := ImportString[
  FromCharacterCode[Normal[bytes], "UTF-8"], "RawJSON"];

BackendInformation[] := If[loadBackend[],
  Quiet[Check[bytesToJSON[$infoFunction[]],
    Failure["CppBackend", <|"Detail" -> "compiled backend information call failed"|>]]],
  Failure["CppBackend", <|"Detail" -> "compiled backend library was not found or could not be loaded",
    "Candidates" -> libraryCandidates[]|>]];

decimalString[x_?InexactNumberQ, digits_Integer] := Module[
  {available, used, y, sign, ds, exponent},
  available = Precision[x];
  used = If[NumericQ[available], Max[1, Min[digits, Floor[available]]], digits];
  y = N[x, used];
  If[TrueQ[y == 0], Return["0", Module]];
  sign = If[TrueQ[y < 0], "-", ""];
  {ds, exponent} = RealDigits[Abs[y], 10, used];
  If[!VectorQ[ds, IntegerQ], Return[Failure["UnsupportedScalar", <|
    "Scalar" -> x, "Detail" ->
      "could not obtain a finite fixed-precision decimal expansion"|>], Module]];
  sign <> ToString[First[ds]] <> "." <>
    StringJoin[ToString /@ Rest[ds]] <> "e" <> ToString[exponent - 1]];

arbInexactString[x_?InexactNumberQ, digits_Integer] := Module[
  {midpoint, accuracy, radiusExponent, radius},
  midpoint = decimalString[x, digits];
  If[FailureQ[midpoint], Return[midpoint, Module]];
  accuracy = Accuracy[x];
  If[accuracy === Infinity, Return[midpoint, Module]];
  If[!NumericQ[accuracy], Return[Failure["UnsupportedScalar", <|
    "Scalar" -> x, "Detail" ->
      "inexact coefficient has no finite accuracy estimate"|>], Module]];
  (* Arb interval syntax preserves the uncertainty already tracked by the
     Wolfram handoff.  The factor two covers midpoint decimalization as well
     as the conventional 10^-Accuracy source-radius estimate. *)
  radiusExponent = -Floor[accuracy];
  radius = "2e" <> If[radiusExponent >= 0, "+", ""] <>
    ToString[radiusExponent];
  "[" <> midpoint <> " +/- " <> radius <> "]"];

realScalarString[x_, digits_Integer] := Which[
  IntegerQ[x] || Head[x] === Rational, ToString[x, InputForm],
  InexactNumberQ[x] && TrueQ[Im[x] == 0], arbInexactString[Re[x], digits],
  NumericQ[x] && TrueQ[Im[N[x, digits]] == 0],
    arbInexactString[Re[N[x, digits]], digits],
  True, Failure["UnsupportedScalar", <|"Scalar" -> x,
    "Detail" -> "coefficient is not a real numeric scalar"|>]];

EncodeScalar[z_, digits_Integer] := Module[{re, im},
  If[!NumericQ[z], Return[Failure["UnsupportedScalar", <|"Scalar" -> z,
    "Detail" -> "symbolic analytic-regulator coefficient is not yet supported by the Acb backend"|>], Module]];
  re = realScalarString[Re[z], digits];
  im = realScalarString[Im[z], digits];
  If[FailureQ[re], Return[re, Module]];
  If[FailureQ[im], Return[im, Module]];
  {re, im}];

EncodeSymbolicScalar[z_, vars_List] := Module[{value, num, den, names, extra},
  value = Cancel[Together[z]];
  names = SymbolName /@ vars;
  If[!FreeQ[value, _?InexactNumberQ] || !FreeQ[value, Complex[0, _]] ||
      !FreeQ[value, I] || !FreeQ[value, _Root] ||
      !FreeQ[value, Power[_, r_Rational /; Denominator[r] =!= 1]],
    Return[Failure["UnsupportedScalar", <|"Scalar" -> z,
      "Detail" -> "symbolic C++ coefficients must be exact rational functions over Q"|>],
      Module]];
  num = Numerator[value]; den = Denominator[value];
  extra = Complement[DeleteDuplicates[Cases[{num, den},
    s_Symbol /; Context[s] =!= "System`", Infinity]], vars];
  If[!PolynomialQ[num, vars] || !PolynomialQ[den, vars] ||
      extra =!= {},
    Return[Failure["UnsupportedScalar", <|"Scalar" -> z,
      "Variables" -> names,
      "Detail" -> "coefficient lies outside the declared rational-function field"|>],
      Module]];
  ToString[value, InputForm]];

parseDecimal[s_String, precision_Integer] := Module[{held, value},
  held = StringReplace[s, {"e+" -> "*^", "e-" -> "*^-", "e" -> "*^"}];
  value = Quiet[Check[ToExpression[held], $Failed]];
  If[value === $Failed || !NumericQ[value],
    Failure["CppBackend", <|"Detail" -> "invalid scalar returned by compiled backend",
      "Scalar" -> s|>],
    SetPrecision[value, precision]]];

parseRadiusExponent["zero"] := -Infinity;
parseRadiusExponent[s_String] := Module[
  {value = Quiet[Check[ToExpression[s], $Failed]]},
  If[value === $Failed || !IntegerQ[value],
    Failure["CppBackend", <|"Detail" ->
      "invalid Arb radius exponent returned by compiled backend",
      "RadiusExponent" -> s|>], value]];

applyBallAccuracy[mid_, radiusExponent_, requested_Integer] := Module[
  {numeric, printedAccuracy, ballAccuracy, targetAccuracy, adjusted},
  If[radiusExponent === -Infinity,
    Return[If[TrueQ[mid === 0], 0, SetPrecision[mid, requested]], Module]];
  numeric = N[mid, requested];
  printedAccuracy = Accuracy[numeric];
  (* C++ supplies the exact integer e for the bound radius < 2^e. *)
  ballAccuracy = N[-radiusExponent*Log[10, 2], 30];
  targetAccuracy = Min[printedAccuracy, ballAccuracy];
  adjusted = SetAccuracy[numeric, targetAccuracy];
  If[NumericQ[Precision[adjusted]] && Precision[adjusted] > requested,
    SetPrecision[adjusted, requested], adjusted]];

DecodeScalar[s_String, _Integer] := Quiet[Check[ToExpression[s],
  Failure["CppBackend", <|"Detail" -> "invalid exact scalar returned by compiled backend",
    "Scalar" -> s|>]]];
DecodeScalar[data_List, precision_Integer] /; Length[data] === 4 := Module[
  {re = parseDecimal[data[[1]], precision], im = parseDecimal[data[[2]], precision],
   reRadius, imRadius},
  If[FailureQ[re], Return[re, Module]];
  If[FailureQ[im], Return[im, Module]];
  reRadius = parseRadiusExponent[data[[3]]];
  imRadius = parseRadiusExponent[data[[4]]];
  If[FailureQ[reRadius], Return[reRadius, Module]];
  If[FailureQ[imRadius], Return[imRadius, Module]];
  applyBallAccuracy[re, reRadius, precision] +
    I applyBallAccuracy[im, imRadius, precision]];
DecodeScalar[x_, _Integer] := Failure["CppBackend", <|
  "Detail" -> "malformed scalar record returned by compiled backend", "Scalar" -> x|>];

RunRequest[request_Association] := Module[{json, bytes, result},
  If[!loadBackend[], Return[BackendInformation[], Module]];
  json = Quiet[Check[ExportString[request, "RawJSON", "Compact" -> True], $Failed]];
  If[json === $Failed, Return[Failure["CppBackend", <|
    "Detail" -> "could not serialize recurrence request"|>], Module]];
  bytes = Quiet[Check[$runFunction[json], $Failed]];
  If[bytes === $Failed, Return[Failure["CppBackend", <|
    "Detail" -> "compiled recurrence LibraryLink call failed"|>], Module]];
  result = Quiet[Check[bytesToJSON[bytes], $Failed]];
  If[!AssociationQ[result], Failure["CppBackend", <|
    "Detail" -> "compiled recurrence returned malformed JSON"|>], result]];

End[];
EndPackage[];
