(* Focused parity tests for the private, default-off singular match-frame
   normalization seam. *)
repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];
catchDE2[DiffExp2`Config`LoadConfiguration[{
  "WorkingPrecision" -> 100, "ExpansionOrder" -> 20,
  "EpsilonOrder" -> 2}]];

esNew = DiffExp2`EpsSeries`ESNew;
esAdd = DiffExp2`EpsSeries`ESAdd;
esScale = DiffExp2`EpsSeries`ESScale;
esTimes = DiffExp2`EpsSeries`ESTimes;
esCoeff = DiffExp2`EpsSeries`ESCoefficient;
esMin = DiffExp2`EpsSeries`ESMinPower;
esCM = DiffExp2`EpsSeries`ESCompleteMax;

(* One t^a Exp[b eps Log[t]] Log[t]^p sector.  A coefficient vector is
   inserted at each listed epsilon order; all Taylor orders above zero are
   absent, which is sufficient for exact match-point frame tests. *)
mkSector[a_, b_, p_, cm_Integer, byOrder_Association, ncomp_Integer] := <|
  "a" -> a, "b" -> b, "p" -> p,
  "Coeffs" -> Table[{Lookup[byOrder, k, ConstantArray[0, ncomp]]},
    {k, 0, cm}]|>;

sharedPrescriptions = {<|"Factor" -> Global`x, "Sign" -> 1,
  "Multiplicity" -> 1, "LeadingCoeffSign" -> 1|>};
mkTaggedColumn[cm_Integer, main_List, extras_List] := <|
  "Center" -> 0, "ChartMap" -> <|"Center" -> 0, "Scale" -> 1|>,
  "Radius" -> 1,
  "Sectors" -> Join[{
    mkSector[0, 0, 0, cm, <|0 -> main|>, Length[main]]}, extras],
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> cm|>,
  "TWindow" -> <|"CompleteMax" -> 0|>,
  "ErrorEstimate" -> ConstantArray[0, cm + 1],
  "Prescriptions" -> sharedPrescriptions|>;

matchPoint = 1/4;
basisValues[bb_List] := Module[{ev},
  ev = Map[DiffExp2`SectorSeries`EvaluateLocalSolution[#,
      matchPoint, "UsePade" -> False, "ImSign" -> 1,
      "ComputeTailEstimates" -> False]["Value"] &, bb];
  Table[esNew[esMin[ev[[i]]],
      Table[esCoeff[ev[[i]], k][[c]], {k, esMin[ev[[i]]], esCM[ev[[i]]]}]],
    {c, Length[esCoeff[ev[[1]], esMin[ev[[1]]]]]},
    {i, Length[bb]}]];

matVec[M_List, w_List] := Table[
  Fold[esAdd, Table[esTimes[M[[r, j]], w[[j]]], {j, Length[w]}]],
  {r, Length[M]}];

mapRightWeights[P_List, w_List] := Table[
  Fold[esAdd, Table[esScale[P[[r, j]], w[[j]]], {j, Length[w]}]],
  {r, Length[P]}];

esCommonSameQ[a_, b_] := Module[{lo = Min[esMin[a], esMin[b]],
    hi = Min[esCM[a], esCM[b]]},
  hi >= lo && AllTrue[Range[lo, hi],
    AllTrue[Flatten[{esCoeff[a, #] - esCoeff[b, #]}],
      Together[Expand[#]] === 0 &] &]];

vectorCommonSameQ[a_List, b_List] := Length[a] === Length[b] &&
  And @@ MapThread[esCommonSameQ, {a, b}];

(* The seam must remain opt-in while its banana replay is under review. *)
assert["singular_match_precondition_default_off",
  DiffExp2`Transport`Private`$enableSingularMatchPrecondition === False];
singularCS = <|"IndicialData" -> <|"Regular" -> False|>|>;
regularCS = <|"IndicialData" -> <|"Regular" -> True|>|>;
assert["singular_match_precondition_toggle_routes_only_when_enabled",
  DiffExp2`Transport`Private`mwUseConstantMatchPreconditionQ[regularCS] &&
  !DiffExp2`Transport`Private`mwUseConstantMatchPreconditionQ[singularCS] &&
  Block[{DiffExp2`Transport`Private`$enableSingularMatchPrecondition = True},
    DiffExp2`Transport`Private`mwUseConstantMatchPreconditionQ[singularCS]]];

(* Exact two-column singular frame.  Its leading matrix
     {{K,K+1},{1,1}}
   has a dense exact inverse.  The columns deliberately have unequal honest
   tops and different fractional/resonant/log tags.  Dense constant mixing
   may impose the common top 2, but it must neither shrink below that top nor
   discard any tag. *)
K = 10^8;
col1 = mkTaggedColumn[4, {K, 1}, {
  mkSector[1/2, 0, 0, 4, <|1 -> {2, 3}|>, 2]}];
col2 = mkTaggedColumn[2, {K + 1, 1}, {
  mkSector[0, 0, 1, 2, <|1 -> {5, 7}|>, 2],
  mkSector[0, 1, 0, 2, <|1 -> {11, 13}|>, 2]}];
basis = {col1, col2};
F = basisValues[basis];
pre = catchDE2[DiffExp2`Transport`Private`mwConstantMatchPrecondition[
  basis, F, basisValues, 2, "singular-tagged-exact"]];

allTags = Sort[DeleteDuplicates[Flatten[
  ({#["a"], #["b"], #["p"]} & /@ # ["Sectors"]) & /@ basis, 1]]];
preTags = If[FailureQ[pre], {},
  (Sort[DeleteDuplicates[{#["a"], #["b"], #["p"]} & /@
      # ["Sectors"]]] &) /@ pre["Basis"]];
assert["singular_constant_frame_preserves_unequal_window_top",
  !FailureQ[pre] && TrueQ[pre["Applied"]] &&
  Min[# ["EpsWindow", "CompleteMax"] & /@ pre["Basis"]] === 2];
assert["singular_constant_frame_preserves_fractional_resonant_log_tags",
  !FailureQ[pre] && AllTrue[preTags, # === allTags &] &&
  AllTrue[pre["Basis"], # ["Prescriptions"] === sharedPrescriptions &]];
assert["singular_constant_frame_rechecks_identity_and_saturation",
  !FailureQ[pre] &&
  DiffExp2`Transport`Private`mwIdentityMatrixQ[
    DiffExp2`Transport`Private`mwLeadingValueMatrix[
      pre["Matrix"], "singular-tagged-check"], "singular-tagged-check"] &&
  With[{verify = catchDE2[DiffExp2`Transport`Private`mwSaturationPlan[
      pre["Matrix"], "singular-tagged-saturation"]]},
    !FailureQ[verify] && verify["Steps"] === 0 &&
      verify["InitialShifts"] === {0, 0}]];
topFailure = catchDE2[DiffExp2`Transport`Private`mwConstantMatchPrecondition[
  basis, F, basisValues, 3, "singular-tagged-short-top"]];
assert["singular_constant_frame_insufficient_common_top_loud",
  FailureQ[topFailure] && topFailure["ID"] === "E5" &&
  topFailure["RequiredTop"] === 3];

(* Literal d=1 exact Laurent parity: the normalized coordinate is seven
   times the original one, and mapping it back by P=1/7 must recover every
   complete Laurent coefficient exactly. *)
scalarBasis = {mkTaggedColumn[3, {7}, {}]};
scalarF = basisValues[scalarBasis];
scalarPre = catchDE2[
  DiffExp2`Transport`Private`mwConstantMatchPrecondition[
    scalarBasis, scalarF, basisValues, 2, "singular-scalar-exact"]];
scalarTruth = {esNew[-1, {2, 3, 5, 0, 0}]};
scalarV = matVec[scalarF, scalarTruth];
scalarNormalizedW = If[FailureQ[scalarPre], scalarPre,
  catchDE2[DiffExp2`Transport`MatchWeights[
    scalarPre["Matrix"], scalarV, "singular-scalar-exact"]]];
scalarMappedW = If[FailureQ[scalarNormalizedW], scalarNormalizedW,
  mapRightWeights[scalarPre["LeadingInverse"], scalarNormalizedW]];
assert["singular_right_frame_exact_scalar_laurent_parity",
  !FailureQ[scalarMappedW] &&
  vectorCommonSameQ[scalarMappedW, scalarTruth]];

(* Laurent weights make the coordinate mapping observable.  Solving in the
   normalized frame and mapping by P must reproduce both the raw exact solve
   and the known original coordinates on every common complete order.  The
   downstream LocalSolution made from normalized columns must evaluate to
   the same vector as the original tagged basis. *)
wTruth = {esNew[-1, {2, 3, 5, 0, 0}], esNew[0, {7, 11, 0, 0}]};
v = matVec[F, wTruth];

(* Expensive-run replay is a no-op unless a path is explicitly supplied.
   When enabled it writes the saturated basis seam atomically and publishes
   a versioned payload that can be loaded without rerunning the recurrence. *)
fixtureDir = CreateDirectory[FileNameJoin[{$TemporaryDirectory,
  "de2-match-fixture-" <> StringReplace[CreateUUID[], "-" -> ""]}]];
fixtureFile = FileNameJoin[{fixtureDir, "match.mx"}];
oldFixtureEnv = Quiet[Environment["DE2_MATCH_FIXTURE_FILE"]];
Clear[Global`$DE2MatchFixture];
Internal`WithLocalSettings[
  SetEnvironment["DE2_MATCH_FIXTURE_FILE" -> None],
  fixtureOffResult = catchDE2[
    DiffExp2`Transport`Private`mwMaybeDumpMatchFixture[
      basis, F, v, matchPoint, 2,
      <|"Name" -> "fixture-off", "Center" -> 0,
        "Scale" -> 1, "Singular" -> True|>, matchPoint]],
  SetEnvironment["DE2_MATCH_FIXTURE_FILE" ->
    If[StringQ[oldFixtureEnv], oldFixtureEnv, None]]];
assert["match_fixture_default_off_has_no_side_effect",
  fixtureOffResult === Null && !FileExistsQ[fixtureFile] &&
  !ValueQ[Global`$DE2MatchFixture]];

Internal`WithLocalSettings[
  SetEnvironment["DE2_MATCH_FIXTURE_FILE" -> fixtureFile],
  fixtureWriteResult = catchDE2[
    DiffExp2`Transport`Private`mwMaybeDumpMatchFixture[
      basis, F, v, matchPoint, 2,
      <|"Name" -> "fixture-on", "Center" -> 0,
        "Scale" -> 1, "Singular" -> True|>, matchPoint]],
  SetEnvironment["DE2_MATCH_FIXTURE_FILE" ->
    If[StringQ[oldFixtureEnv], oldFixtureEnv, None]]];
fixturePayload = If[FileExistsQ[fixtureFile],
  Clear[Global`$DE2MatchFixture]; Get[fixtureFile];
  Global`$DE2MatchFixture, $Failed];
Clear[Global`$DE2MatchFixture];
assert["match_fixture_atomic_write_and_schema",
  fixtureWriteResult === ExpandFileName[fixtureFile] &&
  AssociationQ[fixturePayload] &&
  fixturePayload["Schema"] ===
    DiffExp2`Transport`Private`$matchFixtureSchema &&
  fixturePayload["Label"] === "fixture-on" &&
  fixturePayload["RequiredTop"] === 2 &&
  And @@ (fixturePayload["Config", #] === DiffExp2`Config`CFG[#] & /@
    {"WorkingPrecision", "ExpansionOrder", "EpsilonOrder",
      "DivisionOrder"}) &&
  fixturePayload["Basis"] === basis && fixturePayload["F"] === F &&
  fixturePayload["V"] === v &&
  FileNames["*.tmp-*.mx", fixtureDir] === {}];
DeleteDirectory[fixtureDir, DeleteContents -> True];

wRaw = catchDE2[DiffExp2`Transport`MatchWeights[F, v,
  "singular-tagged-raw"]];
wPre = If[FailureQ[pre], pre,
  catchDE2[DiffExp2`Transport`MatchWeights[pre["Matrix"], v,
    "singular-tagged-preconditioned"]]];
wMapped = If[FailureQ[wPre], wPre,
  mapRightWeights[pre["LeadingInverse"], wPre]];
assert["singular_right_frame_exact_laurent_weight_mapping",
  !FailureQ[wRaw] && !FailureQ[wPre] && !FailureQ[wMapped] &&
  vectorCommonSameQ[wRaw, wTruth] &&
  vectorCommonSameQ[wMapped, wRaw]];

rawLS = If[FailureQ[wRaw], wRaw,
  DiffExp2`SectorSeries`CombineLocalSolutions[wRaw, basis]];
preLS = If[FailureQ[wPre], wPre,
  DiffExp2`SectorSeries`CombineLocalSolutions[wPre, pre["Basis"]]];
rawValue = If[FailureQ[rawLS], rawLS,
  DiffExp2`SectorSeries`EvaluateLocalSolution[rawLS, matchPoint,
    "UsePade" -> False, "ImSign" -> 1]["Value"]];
preValue = If[FailureQ[preLS], preLS,
  DiffExp2`SectorSeries`EvaluateLocalSolution[preLS, matchPoint,
    "UsePade" -> False, "ImSign" -> 1]["Value"]];
vVector = esNew[Min[esMin /@ v], Table[
  Table[esCoeff[v[[c]], k], {c, Length[v]}],
  {k, Min[esMin /@ v], Min[esCM /@ v]}]];
assert["singular_right_frame_tagged_local_solution_parity",
  !FailureQ[rawValue] && !FailureQ[preValue] &&
  esCommonSameQ[rawValue, preValue] &&
  esCommonSameQ[preValue, vVector]];

(* Numeric Schur-cancellation witness.  Every input coefficient has 200-digit
   precision and the incoming vector still has about 98 digits of absolute
   accuracy, but the raw Laurent Gaussian back-substitution manufactures an
   underresolved centered zero from 10^100-sized columns.  Constant right
   normalization has a certified identity product and solves the unchanged
   input without relaxing the zero or residual contracts. *)
mkSeriesColumn[rows_List] := Module[{cm = Length[rows] - 1}, <|
  "Center" -> 0, "ChartMap" -> <|"Center" -> 0, "Scale" -> 1|>,
  "Radius" -> 1,
  "Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
    "Coeffs" -> ({#} & /@ rows)|>},
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> cm|>,
  "TWindow" -> <|"CompleteMax" -> 0|>,
  "ErrorEstimate" -> ConstantArray[0, cm + 1],
  "Prescriptions" -> {}|>];

Knum = N[10^100, 200];
numericBasis = {
  mkSeriesColumn[{{Knum, 1}, {3 Knum, 2}, {-2 Knum, 5}}],
  mkSeriesColumn[{{Knum + 1, 1}, {3 Knum + 1, 2},
    {-2 Knum + 4, 5}}]};
numericF = basisValues[numericBasis];
numericTruth = {esNew[0, {1, 2, 3}], esNew[0, {-1, 4, 7}]};
numericV = matVec[numericF, numericTruth];
numericRaw = catchDE2[DiffExp2`Transport`MatchWeights[
  numericF, numericV, "singular-schur-raw"]];
numericPre = catchDE2[
  DiffExp2`Transport`Private`mwConstantMatchPrecondition[
    numericBasis, numericF, basisValues, 2, "singular-schur-precondition"]];
numericW = If[FailureQ[numericPre], numericPre,
  catchDE2[DiffExp2`Transport`MatchWeights[
    numericPre["Matrix"], numericV, "singular-schur-normalized"]]];
assert["singular_raw_match_reproduces_underresolved_schur_cancellation",
  FailureQ[numericRaw] && numericRaw["ID"] === "E5" &&
  StringContainsQ[numericRaw["Context"], "back substitution"]];
assert["singular_right_normalization_keeps_underresolved_identity_loud",
  FailureQ[numericPre] && numericPre["ID"] === "E5" &&
  StringContainsQ[numericPre["Detail"],
    "re-evaluated preconditioned basis is not identity"]];

(* Actual TransportLine route through a resonant Jordan residue.  The
   second Frobenius solution contains x Log[x], so this exercises the
   private singular toggle beyond the mock basis and compares the physical
   value with both the default path and the closed-form solution. *)
x = Global`x;
logSys = <|"Matrix" -> {{1/x, 1/x}, {0, 1/x}}, "Variable" -> x,
  "SingularFactors" -> {x}|>;
logPlan = DiffExp2`Transport`SegmentLine[logSys, {1/2, 0}];
logBoundary = {{1, 0, 0}, {1, 0, 0}};
logOff = catchDE2[DiffExp2`Transport`TransportLine[
  logSys, logBoundary, logPlan]];
logOn = Block[{
    DiffExp2`Transport`Private`$enableSingularMatchPrecondition = True},
  catchDE2[DiffExp2`Transport`TransportLine[
    logSys, logBoundary, logPlan]]];
logOffValue = If[FailureQ[logOff], logOff,
  DiffExp2`SectorSeries`EvaluateLocalSolution[logOff["Final"], 1/4,
    "UsePade" -> False, "ImSign" -> 1]["Value"]];
logOnValue = If[FailureQ[logOn], logOn,
  DiffExp2`SectorSeries`EvaluateLocalSolution[logOn["Final"], 1/4,
    "UsePade" -> False, "ImSign" -> 1]["Value"]];
numericMax[z_] := Max[0, Sequence @@
  (Abs[N[#, 80]] & /@ Flatten[{z}])];
assert["singular_toggle_transportline_resonant_log_parity",
  !FailureQ[logOffValue] && !FailureQ[logOnValue] &&
  AllTrue[Range[Min[esMin[logOffValue], esMin[logOnValue]],
      Min[esCM[logOffValue], esCM[logOnValue]]],
    numericMax[esCoeff[logOffValue, #] -
      esCoeff[logOnValue, #]] < 10^-70 &] &&
  numericMax[esCoeff[logOnValue, 0] -
    {1/2 - Log[2]/2, 1/2}] < 10^-12];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
