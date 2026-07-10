(* Exact cross-lag denominator-fusion benchmark on the banana L1 endpoint
   charts.  Compares the production fused path with the private unfused
   parity seam; both use the same adaptive lower frames and solver code.

   Environment:
     RDF_ORDERS=10,20  RDF_WP=100  RDF_EPS_ORDER=2  RDF_REPS=1 *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

envOr[name_, def_] := Module[{v = Environment[name]},
  If[StringQ[v] && StringLength[StringTrim[v]] > 0, v, def]];
orders = ToExpression /@ StringSplit[envOr["RDF_ORDERS", "10,20"], ","];
wp = ToExpression[envOr["RDF_WP", "100"]];
epsOrder = ToExpression[envOr["RDF_EPS_ORDER", "2"]];
reps = ToExpression[envOr["RDF_REPS", "1"]];
If[!AllTrue[orders, IntegerQ[#] && # >= 10 &] ||
    !IntegerQ[reps] || reps < 1,
  Print["invalid RDF_ORDERS/RDF_REPS"]; Quit[2]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[e_] := Quiet[Catch[e, "DiffExp2Error"]];
catchDE2[DiffExp2`Config`LoadConfiguration[{
  "WorkingPrecision" -> wp, "ExpansionOrder" -> Max[orders],
  "EpsilonOrder" -> epsOrder, "DivisionOrder" -> 3,
  "StepDivisionOrder" -> 3, "Variables" -> {}}]];

fix = Get[FileNameJoin[{repoRoot, "Tests", "refs", "bench", "banana_L1.m"}]];
sys = catchDE2[DiffExp2`API`LoadSystem[
  <|"Matrix" -> fix["Matrix"], "Variable" -> fix["Variable"]|>]];
If[FailureQ[sys], Print[sys]; Quit[1]];

req[n_] := <|"EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> epsOrder|>,
  "TOrder" -> n|>;
emit[a_] := Print["RDF ", ExportString[a, "RawJSON", "Compact" -> True]];

timeMode[cs_, n_, disable_] := Module[{times, result},
  times = Table[
    DiffExp2`Solve`ClearSolveCaches[];
    First@AbsoluteTiming[result = Block[{
      DiffExp2`Solve`Private`$disableRationalDenominatorFusion = disable},
      catchDE2[DiffExp2`Solve`SolveChart[cs, req[n]]]]],
    {reps}];
  {Min[times], result}];

divisionCounts[prep_, n_] := Module[
  {groups = prep["NhatRationalGroups"], dN = prep["dN"], unfused, fused},
  unfused = Sum[(n - j + 1) Length[groups[[j + 1]]],
    {j, 1, Min[n, dN]}];
  fused = Sum[Module[{active =
        Flatten[groups[[2 ;; Min[n0, dN] + 1]], 1]},
      If[active === {}, 0,
        Length[DeleteDuplicates[Lookup[active, "DenominatorIndex"]]]]],
    {n0, 1, n}];
  {unfused, fused}];

failed = False;
Do[
  plan = catchDE2[DiffExp2`Transport`SegmentLine[sys, {11/23, endpoint}]];
  chart = SelectFirst[Reverse[plan["Charts"]], TrueQ[# ["Singular"]] &];
  cs = catchDE2[DiffExp2`Solve`PrepareChart[sys, chart]];
  symbolic = DiffExp2`Solve`Private`clearedSymbolic[cs];
  (* Denominator identities and lag occupancy are frame-independent here;
     this compact frame keeps benchmark setup negligible. *)
  prep = catchDE2[DiffExp2`Solve`Private`prepareCleared[
    cs, -20, 45, symbolic]];
  Do[
    fusedRun = timeMode[cs, n, False];
    unfusedRun = timeMode[cs, n, True];
    fused = fusedRun[[2]]; unfused = unfusedRun[[2]];
    parity = !FailureQ[fused] && !FailureQ[unfused] &&
      fused["Basis", "Columns"] === unfused["Basis", "Columns"] &&
      fused["Basis", "Specs"] === unfused["Basis", "Specs"];
    If[!TrueQ[parity], failed = True];
    counts = divisionCounts[prep, n];
    emit[<|"Endpoint" -> endpoint, "TaylorOrder" -> n,
      "WorkingPrecision" -> wp, "EpsilonOrder" -> epsOrder,
      "ExactParity" -> parity,
      "FusedSeconds" -> N[fusedRun[[1]], 8],
      "UnfusedSeconds" -> N[unfusedRun[[1]], 8],
      "Speedup" -> N[unfusedRun[[1]]/fusedRun[[1]], 8],
      "DistinctDenominators" -> Length[prep["NhatRationalDenominators"]],
      "TheoreticalUnfusedDivisionsPerColumnLog" -> counts[[1]],
      "TheoreticalFusedDivisionsPerColumnLog" -> counts[[2]]|>],
    {n, orders}],
  {endpoint, {0, 1}}];

If[failed, Quit[1], Quit[0]];
