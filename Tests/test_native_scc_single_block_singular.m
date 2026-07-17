(* One strongly connected singular Jordan block must use the same retained
   composite/basis lifecycle as a multi-block SCC DAG. *)

repo = DirectoryName[DirectoryName[$InputFileName]];
Get[FileNameJoin[{repo, "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["SKIP: compiled DiffExp2 backend is not available"];
  Exit[If[Environment["DE2_REQUIRE_CPP"] === "1", 1, 0]]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

decodeEvaluation[evaluated_] := Module[
  {decoded, min, max, dimension},
  If[!AssociationQ[evaluated] ||
      Lookup[evaluated, "status", "error"] =!= "ok",
    Return[<||>, Module]];
  min = Lookup[evaluated["value"], "min", 1];
  max = Lookup[evaluated["value"], "max", 0];
  dimension = Lookup[evaluated["value"], "dimension", 0];
  decoded = DiffExp2`CppBackend`DecodeScalars[
    evaluated["value", "coefficients"], 60];
  If[!IntegerQ[min] || !IntegerQ[max] || min > max ||
      !IntegerQ[dimension] || dimension < 1 ||
      !ListQ[decoded] || !AllTrue[decoded, NumericQ] ||
      Length[decoded] =!= (max - min + 1) dimension,
    Return[<||>, Module]];
  <|"Min" -> min, "Max" -> max, "Dimension" -> dimension,
    "Table" -> ArrayReshape[
      decoded, {max - min + 1, dimension}]|>];

frameDifference[left_Association, right_Association] := Module[
  {min = Max[left["Min"], right["Min"]],
   max = Min[left["Max"], right["Max"]], l, r},
  If[left["Dimension"] =!= right["Dimension"] || min > max,
    Return[Infinity, Module]];
  l = Flatten[left["Table"][[
    min - left["Min"] + 1 ;; max - left["Min"] + 1]]];
  r = Flatten[right["Table"][[
    min - right["Min"] + 1 ;; max - right["Min"] + 1]]];
  If[l === {} || Length[l] =!= Length[r], Infinity,
    Max[Abs[N[l - r, 50]]]]];

x = Global`x; t = Global`t; eps = Global`eps;
lambda = 1/2 + eps/3;
system = <|"Matrix" -> {{lambda/x, 1/x}, {1, lambda/x}},
  "Variable" -> x|>;
chart = <|"ChartVar" -> t, "Center" -> 0, "Scale" -> 1,
  "Radius" -> Sqrt[2], "LocalRadius" -> Sqrt[2],
  "Name" -> "native-single-scc-singular-jordan",
  "Prescriptions" -> {}, "UseSCCSkeleton" -> True|>;
request = <|"EpsWindow" -> <|"Min" -> -3, "CompleteMax" -> 0|>,
  "TOrder" -> 2|>;

runDomain[exactDomain_] := Block[
  {DiffExp2`Solve`Private`$cppExactDomain = exactDomain},
  Module[{cs, prepared, stats, basis, evaluations, values, result},
    catchDE2[DiffExp2`Config`LoadConfiguration[{
      "WorkingPrecision" -> 80, "ExpansionOrder" -> 10,
      "EpsilonOrder" -> 0, "RecurrenceBackend" -> "Cpp",
      "Variables" -> {}, "Verbosity" -> 0}]];
    DiffExp2`Solve`ClearSolveCaches[];
    cs = catchDE2[DiffExp2`Solve`PrepareChart[system, chart]];
    prepared = If[FailureQ[cs], cs,
      catchDE2[DiffExp2`Solve`PrepareNativeSCCComposite[cs, request]]];
    stats = If[FailureQ[prepared], prepared,
      DiffExp2`CppBackend`PersistentSCCStatistics[prepared]];
    basis = If[FailureQ[prepared], prepared,
      catchDE2[DiffExp2`Solve`SolveNativeSCCBasis[
        cs, request, 2]]];
    evaluations = If[FailureQ[basis], {basis},
      DiffExp2`CppBackend`EvaluatePersistentLocal[#,
          <|"exact" -> "1/4"|>, <|"tail_estimate" -> False|>, 60] & /@
        basis["Columns"]];
    values = decodeEvaluation /@ evaluations;
    result = <|"Chart" -> cs, "Prepared" -> prepared,
      "Stats" -> stats, "Basis" -> basis,
      "Evaluations" -> evaluations, "Values" -> values|>;
    If[AssociationQ[basis] && ListQ[Lookup[basis, "Columns", None]],
      Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &,
        basis["Columns"]]];
    If[AssociationQ[prepared] && StringQ[Lookup[prepared, "SCC", None]],
      Quiet[DiffExp2`CppBackend`ReleasePersistentSCC[prepared]]];
    DiffExp2`Solve`ClearSolveCaches[];
    result]];

rational = runDomain[True];
acb = runDomain[False];

domainOK[result_, capability_] := Module[
  {cs = result["Chart"], stats = result["Stats"],
   basis = result["Basis"], evaluations = result["Evaluations"],
   values = result["Values"]},
  !AnyTrue[{cs, stats, basis}, FailureQ] &&
  AssociationQ[cs] && !TrueQ[Lookup[cs, "SCCSkeleton", False]] &&
  Lookup[cs["IntegrationSequence"], "Components", None] === {{1, 2}} &&
  !TrueQ[Lookup[cs["IndicialData"], "Regular", True]] &&
  Lookup[stats, "status", "error"] === "ok" &&
  Lookup[stats, "blocks", None] === 1 &&
  Lookup[stats, "coupling_groups", None] === 0 &&
  TrueQ[Lookup[stats, "execution_implemented", False]] &&
  Lookup[stats, "execution_scope", None] === capability &&
  Lookup[basis, "Type", None] === "DiffExp2NativeSCCBasis" &&
  Lookup[basis, "Dimension", None] === 2 &&
  Lookup[Lookup[basis, "Columns", {}], "BasisIndex", {}] === {1, 2} &&
  TrueQ[Lookup[basis["NativeSummary"], "atomic_retention", False]] &&
  AllTrue[Lookup[basis, "Columns", {}],
    Lookup[# ["NativeSummary"], "execution_capability", None] ===
      capability &] &&
  Length[evaluations] === 2 &&
  AllTrue[evaluations,
    AssociationQ[#] && Lookup[#, "status", "error"] === "ok" &] &&
  Length[values] === 2 && AllTrue[values, AssociationQ[#] &&
    # =!= <||> &]];

rationalCapability =
  "exact-rational-regular-singular-jordan-block-dag-column-v2";
acbCapability = "acb-regular-singular-jordan-block-dag-column-v1";
parity = If[AllTrue[Join[rational["Values"], acb["Values"]],
    AssociationQ[#] && # =!= <||> &],
  Max[MapThread[frameDifference,
    {acb["Values"], rational["Values"]}]] < 10^-40,
  False];

(* A Jordan/log column can honestly retain eps^-1 even when the public
   request starts at eps^0.  The Rational shadow and Acb owner still share
   the same compact work rectangle; specialization must preserve that useful
   lower Laurent coefficient instead of demanding an identical public
   minimum. *)
shadowRequest = <|
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 0|>,
  "TOrder" -> 2|>;
shadowCs = catchDE2[DiffExp2`Solve`PrepareChart[system, chart]];
targetPrepared = Block[
  {DiffExp2`Solve`Private`$cppExactDomain = False},
  catchDE2[DiffExp2`Solve`PrepareNativeSCCComposite[
    shadowCs, shadowRequest]]];
targetStats = If[AssociationQ[targetPrepared],
  DiffExp2`CppBackend`PersistentSCCStatistics[targetPrepared],
  targetPrepared];
sourceBasis = Block[
  {DiffExp2`Solve`Private`$cppExactDomain = True},
  catchDE2[DiffExp2`Solve`SolveNativeSCCBasis[
    shadowCs, shadowRequest, 2]]];
shadowIdentity = If[AssociationQ[targetStats],
  Lookup[targetStats, "rational_shadow_identity", None], None];
shadowImports = If[AssociationQ[sourceBasis] &&
    StringQ[shadowIdentity],
  MapIndexed[
    DiffExp2`CppBackend`SpecializePersistentRationalSCCColumn[
      #1, targetPrepared, shadowIdentity,
      "single-scc-extra-lower-shadow:" <>
        ToString[First[#2]]] &,
    sourceBasis["Columns"]], {}];
shadowSourceHasExtraLower = AssociationQ[sourceBasis] &&
  AnyTrue[Lookup[sourceBasis, "Columns", {}],
    Lookup[Lookup[#, "EpsWindow", <||>], "Min", 0] <
      shadowRequest["EpsWindow", "Min"] &];
shadowImportsOK = Length[shadowImports] === 2 &&
  AllTrue[shadowImports, AssociationQ[#] &&
    Lookup[#, "status", "error"] === "ok" &&
    Lookup[#, "specialization_capability", None] ===
      "exact-rational-shadow-to-acb-local-v1" &];
If[AssociationQ[sourceBasis],
  Quiet[DiffExp2`CppBackend`ClosePersistentSession[sourceBasis]]];
Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &,
  Select[shadowImports, AssociationQ]];
If[AssociationQ[targetPrepared],
  Quiet[DiffExp2`CppBackend`ReleasePersistentSCC[targetPrepared]]];

ok = domainOK[rational, rationalCapability] &&
  domainOK[acb, acbCapability] && TrueQ[parity] &&
  TrueQ[shadowSourceHasExtraLower] && TrueQ[shadowImportsOK];

DiffExp2`CppBackend`ClearPersistentSessions[];

If[TrueQ[ok],
  Print["PASS: Rational/Acb one-block singular Jordan native basis"],
  Print["FAIL: ", InputForm[<|
    "Rational" -> KeyTake[rational, {"Chart", "Stats", "Basis"}],
    "Acb" -> KeyTake[acb, {"Chart", "Stats", "Basis"}],
    "Parity" -> parity,
    "ShadowSourceHasExtraLower" -> shadowSourceHasExtraLower,
    "ShadowImports" -> shadowImports|>]];
  Exit[1]];
