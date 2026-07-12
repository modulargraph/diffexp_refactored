(* One strongly connected singular Jordan block must use the same retained
   composite/basis lifecycle as a multi-block SCC DAG. *)

repo = DirectoryName[DirectoryName[$InputFileName]];
Get[FileNameJoin[{repo, "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["SKIP: compiled DiffExp2 backend is not available"];
  Exit[If[Environment["DE2_REQUIRE_CPP"] === "1", 1, 0]]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

decodeEvaluation[evaluated_] := Module[{decoded},
  If[!AssociationQ[evaluated] ||
      Lookup[evaluated, "status", "error"] =!= "ok",
    Return[{}, Module]];
  decoded = DiffExp2`CppBackend`DecodeScalars[
    evaluated["value", "coefficients"], 60];
  If[ListQ[decoded] && AllTrue[decoded, NumericQ], decoded, {}]];

x = Global`x; t = Global`t; eps = Global`eps;
lambda = 1/2 + eps/3;
system = <|"Matrix" -> {{lambda/x, 1/x}, {1, lambda/x}},
  "Variable" -> x|>;
chart = <|"ChartVar" -> t, "Center" -> 0, "Scale" -> 1,
  "Radius" -> 2, "LocalRadius" -> 2,
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
  Length[values] === 2 && AllTrue[values, Length[#] > 0 &]];

rationalCapability =
  "exact-rational-regular-singular-jordan-block-dag-column-v2";
acbCapability = "acb-regular-singular-jordan-block-dag-column-v1";
parity = If[AllTrue[Join[rational["Values"], acb["Values"]],
    ListQ[#] && # =!= {} &],
  Max[Abs[N[Flatten[acb["Values"] - rational["Values"]], 50]]] < 10^-40,
  False];
ok = domainOK[rational, rationalCapability] &&
  domainOK[acb, acbCapability] && TrueQ[parity];

DiffExp2`CppBackend`ClearPersistentSessions[];

If[TrueQ[ok],
  Print["PASS: Rational/Acb one-block singular Jordan native basis"],
  Print["FAIL: ", InputForm[<|
    "Rational" -> KeyTake[rational, {"Chart", "Stats", "Basis"}],
    "Acb" -> KeyTake[acb, {"Chart", "Stats", "Basis"}],
    "Parity" -> parity|>]];
  Exit[1]];
