(* Public IntegrateLine selection of the persistent concurrent-arm solver. *)

repo = DirectoryName[DirectoryName[$InputFileName]];
Get[FileNameJoin[{repo, "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["SKIP: compiled backend unavailable"];
  Exit[0]];

DiffExp2`LoadConfiguration[{
  "RecurrenceBackend" -> "Cpp",
  "WorkingPrecision" -> 100,
  "ChopPrecision" -> 50,
  "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 2,
  "DivisionOrder" -> 3,
  "Verbosity" -> 0}];

x = Global`x;
sys = DiffExp2`LoadSystem[<|"Matrix" -> {{0}}, "Variable" -> x|>];
boundary = {{1, 0, 0}};
coefficients = {1};
coefficient = DiffExp2`EpsilonCoefficient;
nearQ[a_, b_] := TrueQ[Abs[N[a - b, 40]] < 10^-30];

DiffExp2`CppBackend`ClearPersistentSessions[];
nativeValue = Catch[
  DiffExp2`IntegrateLine[sys, boundary, 0, {-1, 1}, coefficients,
    "ExtraSingularFactors" -> {x - 2}],
  "DiffExp2Error"];
nativeSessionInfo = DiffExp2`CppBackend`PersistentSessionInformation[];
nativeStats = Values[nativeSessionInfo];

nativeOK = !FailureQ[nativeValue] &&
  nearQ[coefficient[nativeValue, 0], 2] &&
  nearQ[coefficient[nativeValue, 1], 0] &&
  nearQ[coefficient[nativeValue, 2], 0] &&
  Length[nativeStats] === 1 &&
  Lookup[First[nativeStats], "tile_plans_created", 0] === 1 &&
  Lookup[First[nativeStats], "line_integrations", 0] > 0 &&
  Lookup[First[nativeStats], "line_exports", 0] === 1 &&
  Lookup[First[nativeStats], "tile_plans", -1] === 0 &&
  Lookup[First[nativeStats], "line_results", -1] === 0 &&
  Lookup[First[nativeStats], "locals", -1] === 0 &&
  Lookup[First[nativeStats], "matches", -1] === 0;

(* PrecomputedCharts has an existing object-level meaning: integrate that
   exact supplied chain.  It is therefore an explicit non-native dispatch
   reason, not permission to silently discard and rebuild the chain. *)
DiffExp2`CppBackend`ClearPersistentSessions[];
DiffExp2`UpdateConfiguration[{"RecurrenceBackend" -> "Wolfram"}];
lower = Catch[DiffExp2`TransportLine[sys, boundary, {0, -1}],
  "DiffExp2Error"];
upper = Catch[DiffExp2`TransportLine[sys, boundary, {0, 1}],
  "DiffExp2Error"];
precomputed = If[FailureQ[lower] || FailureQ[upper], {},
  Join[lower["Charts"], upper["Charts"]]];
DiffExp2`UpdateConfiguration[{"RecurrenceBackend" -> "Cpp"}];
DiffExp2`CppBackend`ClearPersistentSessions[];

decision = DiffExp2`API`Private`lineIntegralDispatchDecision[
  sys, coefficients, 0, {-1, 1}, precomputed];
precomputedValue = Catch[
  DiffExp2`IntegrateLine[sys, boundary, 0, {-1, 1}, coefficients,
    "PrecomputedCharts" -> precomputed],
  "DiffExp2Error"];
precomputedSessionInfo =
  DiffExp2`CppBackend`PersistentSessionInformation[];

precomputedOK = precomputed =!= {} &&
  decision === <|"Mode" -> "EstablishedOrchestration",
    "Reason" -> "PrecomputedChartsOwnTheIntegrationChain"|> &&
  !FailureQ[precomputedValue] &&
  nearQ[coefficient[precomputedValue, 0], 2] &&
  nearQ[coefficient[precomputedValue, 1], 0] &&
  nearQ[coefficient[precomputedValue, 2], 0] &&
  precomputedSessionInfo === <||>;

DiffExp2`CppBackend`ClearPersistentSessions[];

If[TrueQ[nativeOK && precomputedOK],
  Print["PASS: public persistent-native line dispatch"],
  Print["FAIL: nativeValue=", InputForm[nativeValue],
    "; nativeStats=", InputForm[nativeStats],
    "; decision=", InputForm[decision],
    "; precomputedValue=", InputForm[precomputedValue],
    "; precomputedSessions=", InputForm[precomputedSessionInfo]];
  Exit[1]];
