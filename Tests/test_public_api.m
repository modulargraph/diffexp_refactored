(* Contract tests for the DiffExp2` release umbrella.  These tests exercise
   only object accessors and local evaluation; no recurrence backend, FIRE,
   or transport solve is launched. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
Get[FileNameJoin[{repoRoot, "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

exports = {
  DiffExp2`LoadConfiguration, DiffExp2`UpdateConfiguration,
  DiffExp2`CurrentConfiguration, DiffExp2`LoadSystem,
  DiffExp2`LoadCanonicalSystem,
  DiffExp2`CanonicalLineChartGeometry,
  DiffExp2`TransportCanonicalLine,
  DiffExp2`PrepareBoundary, DiffExp2`PrepareLaurentBoundary,
  DiffExp2`PlanLine, DiffExp2`TransportEndpoint, DiffExp2`TransportLine,
  DiffExp2`LineSegments, DiffExp2`LineSegment, DiffExp2`EvaluateLine,
  DiffExp2`PiecewiseSolution, DiffExp2`EvaluateLocal,
  DiffExp2`LocalBehavior, DiffExp2`ExactSectors, DiffExp2`EndpointLimit,
  DiffExp2`IntegrateLine, DiffExp2`EpsilonWindow,
  DiffExp2`EpsilonCoefficient, DiffExp2`EpsilonCoefficientList,
  DiffExp2`EpsilonExpression
};
assert["umbrella exports have usage text",
  AllTrue[exports, StringQ[MessageName[#, "usage"]] &]];
assert["compiled backend is the release default",
  DiffExp2`CurrentConfiguration[]["RecurrenceBackend"] === "Cpp"];
loadSystemDefinitions = DownValues[DiffExp2`LoadSystem];
Get[FileNameJoin[{repoRoot, "DiffExp2.m"}]];
assert["root loader is idempotent and keeps only the umbrella context",
  DownValues[DiffExp2`LoadSystem] === loadSystemDefinitions &&
  First[$ContextPath] === "DiffExp2`" &&
  !MemberQ[$ContextPath, "DiffExp2`API`"]];

DiffExp2`LoadConfiguration["RecurrenceBackend" -> "Wolfram"];
assert["explicit Wolfram diagnostic selection is respected",
  DiffExp2`CurrentConfiguration[]["RecurrenceBackend"] === "Wolfram"];
DiffExp2`LoadConfiguration["WorkingPrecision" -> 100];
assert["release default is reinjected on configuration reset",
  DiffExp2`CurrentConfiguration[]["RecurrenceBackend"] === "Cpp"];

epsValue = DiffExp2`EpsSeries`ESNew[-1, {2, 3, 5}];
assert["honest epsilon window",
  DiffExp2`EpsilonWindow[epsValue] ===
    <|"Min" -> -1, "CompleteMax" -> 1|>];
assert["epsilon coefficient",
  DiffExp2`EpsilonCoefficient[epsValue, 0] === 3];
assert["epsilon coefficient list",
  DiffExp2`EpsilonCoefficientList[epsValue, -1, 1] === {2, 3, 5}];
assert["epsilon presentation expression",
  Expand[DiffExp2`EpsilonExpression[epsValue, Global`z]] ===
    2/Global`z + 3 + 5 Global`z];

preparedBoundary = DiffExp2`PrepareBoundary[
  {Exp[Global`eps], (1/2)^Global`\[Epsilon]},
  "EpsilonOrder" -> 2];
assert["closed-form regular boundary preparation",
  preparedBoundary === {
    {1, 1, 1/2}, {1, -Log[2], Log[2]^2/2}}];
poleBoundary = Catch[
  DiffExp2`PrepareBoundary[{1/Global`eps}, "EpsilonOrder" -> 2],
  "DiffExp2Error"];
assert["closed-form boundary preparation preserves the Laurent guard",
  FailureQ[poleBoundary]];
emptyBoundary = Catch[DiffExp2`PrepareBoundary[{}], "DiffExp2Error"];
assert["empty closed-form boundaries fail at the public seam",
  FailureQ[emptyBoundary]];

laurentBoundary = DiffExp2`PrepareLaurentBoundary[
  {1/Global`eps + 2 + 3 Global`eps, 1/Global`eps^2 + 5},
  7, "EpsilonOrder" -> 2];
assert["Laurent boundary constructor preserves the common honest window",
  laurentBoundary["EpsWindow"] ===
    <|"Min" -> -2, "CompleteMax" -> 2|> &&
  laurentBoundary["PointDatum"] === True &&
  laurentBoundary["Center"] === 7 &&
  laurentBoundary["Radius"] === Infinity];
assert["Laurent boundary constructor aligns component coefficients",
  laurentBoundary["Sectors"][[1, "Coeffs"]] === {
    {{0, 1}}, {{1, 0}}, {{2, 5}}, {{3, 0}}, {{0, 0}}
  } &&
  DiffExp2`SectorSeries`ValidateLocalSolution[laurentBoundary] ===
    laurentBoundary];
badLaurentRadius = Catch[
  DiffExp2`PrepareLaurentBoundary[{1/Global`eps}, 0,
    "Radius" -> 0],
  "DiffExp2Error"];
assert["Laurent boundary constructor rejects a nonpositive radius",
  FailureQ[badLaurentRadius]];

aliasSystem = DiffExp2`LoadSystem[<|
  "Matrix" -> {{Global`\[Epsilon]/(1 - Global`x)}},
  "Variable" -> Global`x|>];
assert["public ingestion canonicalizes epsilon aliases",
  FreeQ[aliasSystem["Matrix"], Global`\[Epsilon]] &&
  !FreeQ[aliasSystem["Matrix"], Global`eps]];

ls = <|
  "Center" -> 0,
  "ChartMap" -> <|"Center" -> 0, "Scale" -> 1|>,
  "Radius" -> 2,
  "Sectors" -> {<|
    "a" -> -1/2, "b" -> 3, "p" -> 1,
    "Coeffs" -> {{{7}}, {{11}}}|>},
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 1|>,
  "TWindow" -> <|"CompleteMax" -> 0|>,
  "ErrorEstimate" -> {0, 0},
  "Prescriptions" -> {}
|>;

behavior = DiffExp2`LocalBehavior[ls];
sectors = DiffExp2`ExactSectors[ls];
assert["local behavior keeps exact tags",
  KeyTake[behavior["Sectors"][[1]], {"a", "b", "p"}] ===
    <|"a" -> -1/2, "b" -> 3, "p" -> 1|>];
assert["exact x^(a+b eps) exponent is exposed",
  Together[sectors[[1, "Exponent"]] - (-1/2 + 3 Global`eps)] === 0 &&
  sectors[[1, "LogPower"]] === 1];

mockResult = <|
  "Schema" -> "DiffExp2.TransportResult/v1", "From" -> 0, "To" -> 1,
  "Plan" -> <|"From" -> 0, "To" -> 1|>,
  "Charts" -> {<|"Chart" -> <|"Center" -> 0, "Radius" -> 2|>,
    "LocalSolution" -> ls|>},
  "Final" -> ls, "EndpointIsSingular" -> False
|>;
segments = DiffExp2`LineSegments[mockResult];
piecewise = DiffExp2`PiecewiseSolution[mockResult];
assert["named line segment schema",
  Length[segments] === 1 &&
  segments[[1, "Schema"]] === "DiffExp2.LineSegment/v1" &&
  segments[[1, "Domain"]] === {0, 1}];
assert["piecewise solution is inspectable",
  piecewise["Schema"] === "DiffExp2.PiecewiseSolution/v1" &&
  Length[piecewise["Segments"]] === 1 &&
  Head[piecewise["Function"]] === Function];
assert["transport result exposes its honest epsilon window",
  DiffExp2`EpsilonWindow[mockResult] === ls["EpsWindow"]];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Exit[1], Exit[0]];
