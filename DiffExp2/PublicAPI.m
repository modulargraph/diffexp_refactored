(* DiffExp2/PublicAPI.m -- the small, stable DiffExp2` user surface.

   The solver modules deliberately live in implementation subcontexts.  This
   file turns their typed objects into a coherent release API without copying
   any recurrence, matching, or regularisation mathematics.  In particular,
   exact local sectors remain the source of truth: no exponent fitting is
   reintroduced at this layer. *)

BeginPackage["DiffExp2`"];

$DiffExp2Version::usage = "$DiffExp2Version is the public API version string.";

LoadConfiguration::usage = "LoadConfiguration[rules] resets DiffExp2 to its validated defaults and applies rules. The compiled C++ recurrence backend is the default.";
UpdateConfiguration::usage = "UpdateConfiguration[rules] atomically validates and updates the current DiffExp2 configuration.";
CurrentConfiguration::usage = "CurrentConfiguration[] returns the resolved, String-keyed DiffExp2 configuration.";

LoadSystem::usage = "LoadSystem[spec] loads an exact epsilon-rational differential system. spec may be an Association, a d<var>_full.m file, or a directory containing exactly one such file.";
BackendInformation::usage = "BackendInformation[] reports the compiled recurrence backend version and availability.";

PlanLine::usage = "PlanLine[sys, {from, to}, opts] constructs and validates a transport plan without solving it. Option: \"ExtraSingularFactors\".";
TransportEndpoint::usage = "TransportEndpoint[sys, boundary, from, to, opts] transports boundary data to one endpoint and returns a named TransportResult.";
TransportLine::usage = "TransportLine[sys, boundary, {from, to}, opts] plans and transports a complete line. TransportLine[sys, boundary, plan] reuses an existing plan.";
LineSegments::usage = "LineSegments[result] returns named, individually inspectable local segments with domains and exact LocalSolution objects.";
LineSegment::usage = "LineSegment[result, i] returns segment i from LineSegments[result].";
EvaluateLine::usage = "EvaluateLine[result, point, opts] evaluates a transported line at a covered regular point and returns the local evaluation record.";
PiecewiseSolution::usage = "PiecewiseSolution[result] returns a PiecewiseSolution record with inspectable segments and an evaluation function suitable for sampling and plotting.";

EvaluateLocal::usage = "EvaluateLocal[local, t, opts] evaluates an exact LocalSolution at local chart coordinate t.";
LocalBehavior::usage = "LocalBehavior[localOrResult] returns the exact canonical local-sector decomposition, including all x^(a+b eps) and log sectors.";
ExactSectors::usage = "ExactSectors[localOrResult] returns user-facing exact sector records with Exponent = a+b eps, LogPower, and coefficient tensors.";
EndpointLimit::usage = "EndpointLimit[result] returns the dimensionally regularized endpoint limit of every component. EndpointLimit[result, weights] combines components before taking the limit.";
IntegrateLine::usage = "IntegrateLine[sys, boundary, from, {lo, hi}, coefficients, opts] integrates a rational linear combination over a transported line.";

EpsilonWindow::usage = "EpsilonWindow[value] returns its honest Laurent window <|\"Min\", \"CompleteMax\"|>. It accepts EpsSeries values, evaluation records, transport results, and lists thereof.";
EpsilonCoefficient::usage = "EpsilonCoefficient[value, k] returns the coefficient of eps^k, loudly rejecting reads above the complete window.";
EpsilonCoefficientList::usage = "EpsilonCoefficientList[value, k1, k2] returns coefficients eps^k1 through eps^k2 with the honest-window drop guard.";

Begin["`Private`"];

$DiffExp2Version = "2.0.0-prototype";

de2Error[id_String, detail_, payload_:<||>] :=
  DiffExp2`Tolerances`DE2Error[id, Join[
    <|"Module" -> "PublicAPI", "Detail" -> detail|>, payload]];

(* Configuration is stateful by design, but only Config.m owns the state.
   The low-level implementation loader retains its Wolfram-reference default
   for its focused tests.  The RELEASE umbrella injects Cpp unless the user
   explicitly selected a backend; this keeps the default at the right API
   boundary without changing developer-module semantics. *)
withReleaseBackend[rules_List] := Module[{names},
  names = DiffExp2`Config`CanonicalKey[First[#]] & /@ rules;
  If[MemberQ[names, "RecurrenceBackend"], rules,
    Append[rules, "RecurrenceBackend" -> "Cpp"]]];

DiffExp2`LoadConfiguration[] :=
  DiffExp2`Config`LoadConfiguration[
    {"RecurrenceBackend" -> "Cpp"}];
DiffExp2`LoadConfiguration[rules : {___Rule}] :=
  DiffExp2`Config`LoadConfiguration[withReleaseBackend[rules]];
DiffExp2`LoadConfiguration[assoc_Association] :=
  DiffExp2`LoadConfiguration[Normal[assoc]];
DiffExp2`LoadConfiguration[rules__Rule] :=
  DiffExp2`LoadConfiguration[{rules}];
DiffExp2`LoadConfiguration[args___] :=
  DiffExp2`Config`LoadConfiguration[args];
DiffExp2`UpdateConfiguration[args___] :=
  DiffExp2`Config`UpdateConfiguration[args];
DiffExp2`CurrentConfiguration[] :=
  DiffExp2`Config`CurrentConfiguration[];

DiffExp2`BackendInformation[] :=
  DiffExp2`CppBackend`BackendInformation[];

(* ---- systems ---- *)

Options[DiffExp2`LoadSystem] = {"Variable" -> Automatic};

systemRecord[sys_Association, source_] := Join[sys, <|
  "Schema" -> "DiffExp2.System/v1",
  "Source" -> source,
  "Dimension" -> Length[sys["Matrix"]]
|>];

variableFromFullFile[file_String] := Module[{base, names},
  base = FileBaseName[file];
  names = StringCases[base,
    StartOfString ~~ "d" ~~ name__ ~~ "_full" ~~ EndOfString :> name];
  If[Length[names] =!= 1 || StringLength[First[names]] === 0,
    de2Error["E6", "could not infer the differential variable from the exact matrix filename",
      <|"File" -> file, "Expected" -> "d<variable>_full.m"|>]];
  Symbol["Global`" <> First[names]]];

loadSystemFile[file_String, variable_] := Module[{var},
  If[!FileExistsQ[file],
    de2Error["E6", "exact matrix file does not exist", <|"File" -> file|>]];
  var = Replace[variable, Automatic :> variableFromFullFile[file]];
  If[!MatchQ[var, _Symbol],
    de2Error["E6", "\"Variable\" must be a Symbol", <|"Variable" -> var|>]];
  systemRecord[DiffExp2`API`LoadSystem[
    <|"FullMatrixFile" -> ExpandFileName[file], "Variable" -> var|>],
    ExpandFileName[file]]];

DiffExp2`LoadSystem[spec_Association, OptionsPattern[]] := Module[{sys},
  sys = DiffExp2`API`LoadSystem[spec];
  systemRecord[sys, Lookup[spec, "FullMatrixFile", "InMemory"]]];

DiffExp2`LoadSystem[path_String, OptionsPattern[]] := Module[{fullFiles},
  If[DirectoryQ[path],
    fullFiles = Sort[FileNames["d*_full.m", ExpandFileName[path], 1]];
    If[Length[fullFiles] =!= 1,
      de2Error["E6",
        "a release-v1 directory must contain exactly one exact d<variable>_full.m file",
        <|"Directory" -> ExpandFileName[path], "Files" -> fullFiles|>]];
    loadSystemFile[First[fullFiles], OptionValue["Variable"]],
    loadSystemFile[path, OptionValue["Variable"]]]];

(* ---- line planning and transport ---- *)

Options[DiffExp2`PlanLine] = {"ExtraSingularFactors" -> {}};

systemWithExtraFactors[sys_Association, extra_] := Module[{var, clean},
  var = sys["Variable"];
  If[!ListQ[extra],
    de2Error["E8", "\"ExtraSingularFactors\" must be a list",
      <|"Value" -> extra|>]];
  clean = Select[extra, !FreeQ[#, var] &];
  Join[sys, <|"ExtraSingularFactors" -> clean|>]];

DiffExp2`PlanLine[sys_Association, {from_, to_}, OptionsPattern[]] := Module[
  {transportSystem, plan},
  transportSystem = systemWithExtraFactors[sys,
    OptionValue["ExtraSingularFactors"]];
  plan = DiffExp2`Transport`SegmentLine[transportSystem, {from, to}];
  Join[plan, <|"Schema" -> "DiffExp2.LinePlan/v1",
    "ExtraSingularFactors" -> transportSystem["ExtraSingularFactors"]|>]];

Options[DiffExp2`TransportLine] = Options[DiffExp2`PlanLine];

transportResult[sys_, boundary_, plan_] := Module[{transportSystem, raw},
  transportSystem = systemWithExtraFactors[sys,
    Lookup[plan, "ExtraSingularFactors", {}]];
  raw = DiffExp2`Transport`TransportLine[transportSystem, boundary, plan];
  Join[raw, <|
    "Schema" -> "DiffExp2.TransportResult/v1",
    "From" -> plan["From"], "To" -> plan["To"], "Plan" -> plan
  |>]];

DiffExp2`TransportLine[sys_Association, boundary_, {from_, to_},
    OptionsPattern[]] := Module[{plan},
  plan = DiffExp2`PlanLine[sys, {from, to},
    "ExtraSingularFactors" -> OptionValue["ExtraSingularFactors"]];
  transportResult[sys, boundary, plan]];

DiffExp2`TransportLine[sys_Association, boundary_, plan_Association] /;
    KeyExistsQ[plan, "Charts"] := transportResult[sys, boundary, plan];

Options[DiffExp2`TransportEndpoint] = Options[DiffExp2`TransportLine];
DiffExp2`TransportEndpoint[sys_Association, boundary_, from_, to_,
    OptionsPattern[]] := DiffExp2`TransportLine[sys, boundary, {from, to},
  "ExtraSingularFactors" -> OptionValue["ExtraSingularFactors"]];

localSolutionOf[obj_Association] := Which[
  KeyExistsQ[obj, "Sectors"], obj,
  KeyExistsQ[obj, "LocalSolution"], obj["LocalSolution"],
  KeyExistsQ[obj, "Final"], obj["Final"],
  True, de2Error["E9", "object does not contain a LocalSolution",
    <|"Keys" -> Keys[obj]|>]];

segmentRecords[result_Association] := Module[
  {entries, plan, from, to, lo, hi},
  entries = Lookup[result, "Charts", Missing["NotAvailable"]];
  If[!ListQ[entries],
    de2Error["E9", "transport result has no chart list",
      <|"Keys" -> Keys[result]|>]];
  plan = Lookup[result, "Plan", <||>];
  {from, to} = Lookup[plan, {"From", "To"}, {-Infinity, Infinity}];
  {lo, hi} = If[NumericQ[from] && NumericQ[to],
    {Min[from, to], Max[from, to]}, {-Infinity, Infinity}];
  MapIndexed[Function[{entry, index}, Module[
    {ls = entry["LocalSolution"], center, scale, halfWidth, domain},
    center = ls["Center"];
    scale = ls["ChartMap", "Scale"];
    halfWidth = Abs[scale]*ls["Radius"]/2;
    domain = {Max[lo, center - halfWidth], Min[hi, center + halfWidth]};
    <|"Schema" -> "DiffExp2.LineSegment/v1", "Index" -> First[index],
      "Center" -> center, "Scale" -> scale, "Radius" -> ls["Radius"],
      "Domain" -> domain, "Chart" -> entry["Chart"],
      "LocalSolution" -> ls|>]], entries]];

DiffExp2`LineSegments[result_Association] := segmentRecords[result];

DiffExp2`LineSegment[result_Association, i_Integer] := Module[{segments},
  segments = DiffExp2`LineSegments[result];
  If[!(1 <= i <= Length[segments]),
    de2Error["E9", "line-segment index is out of range",
      <|"Index" -> i, "SegmentCount" -> Length[segments]|>]];
  segments[[i]]];

Options[DiffExp2`EvaluateLocal] = Options[
  DiffExp2`SectorSeries`EvaluateLocalSolution];
DiffExp2`EvaluateLocal[obj_Association, t_, opts:OptionsPattern[]] :=
  DiffExp2`SectorSeries`EvaluateLocalSolution[localSolutionOf[obj], t, opts];

insideDomainQ[x_, {lo_, hi_}] := TrueQ[
  N[lo - x, 50] <= 0 && N[x - hi, 50] <= 0];

Options[DiffExp2`EvaluateLine] = Options[DiffExp2`EvaluateLocal];
DiffExp2`EvaluateLine[result_Association, point_, opts:OptionsPattern[]] :=
 Module[{segments, candidates, chosen, ls, t},
  If[TrueQ[Lookup[result, "EndpointIsSingular", False]] &&
      SameQ[point, Lookup[result, "To", Missing["Unknown"]]],
    de2Error["E9",
      "a singular endpoint cannot be sampled directly; use EndpointLimit or LocalBehavior",
      <|"Point" -> point|>]];
  segments = DiffExp2`LineSegments[result];
  candidates = Select[segments, insideDomainQ[point, #["Domain"]] &];
  If[candidates === {},
    de2Error["E9", "point is outside the certified piecewise line coverage",
      <|"Point" -> point, "Domains" -> (# ["Domain"] & /@ segments)|>]];
  chosen = First[SortBy[candidates,
    Abs[N[(point - #["Center"])/#["Scale"], 50]] &]];
  ls = chosen["LocalSolution"];
  t = Together[(point - chosen["Center"])/chosen["Scale"]];
  DiffExp2`EvaluateLocal[ls, t, opts]];

DiffExp2`PiecewiseSolution[result_Association] := Module[{segments},
  segments = DiffExp2`LineSegments[result];
  <|"Schema" -> "DiffExp2.PiecewiseSolution/v1",
    "Domain" -> Lookup[result, {"From", "To"}, Missing["NotAvailable"]],
    "Segments" -> segments,
    "Function" -> Function[point, DiffExp2`EvaluateLine[result, point]]|>];

(* ---- local sectors, limits, and integration ---- *)

DiffExp2`LocalBehavior[obj_Association] :=
  DiffExp2`SectorSeries`SectorDecomposition[localSolutionOf[obj]];

DiffExp2`ExactSectors[obj_Association] := Module[{decomp, eps},
  decomp = DiffExp2`LocalBehavior[obj];
  eps = DiffExp2`Config`CanonicalEps[];
  Map[Function[sector, <|
    "a" -> sector["a"], "b" -> sector["b"], "p" -> sector["p"],
    "Exponent" -> Together[sector["a"] + sector["b"]*eps],
    "LogPower" -> sector["p"], "Coefficients" -> sector["Coeffs"]|>],
    decomp["Sectors"]]];

DiffExp2`EndpointLimit[result_Association, weights_:Automatic] := Module[{ls},
  ls = localSolutionOf[result];
  If[!TrueQ[Lookup[result, "EndpointIsSingular", False]],
    de2Error["E9", "EndpointLimit requires a singular-endpoint transport result"]];
  If[weights === Automatic,
    DiffExp2`Integrate`EndpointSectorLimit[ls],
    DiffExp2`API`EndpointLimitValues[result, weights]]];

Options[DiffExp2`IntegrateLine] = Options[DiffExp2`API`LineIntegral];
DiffExp2`IntegrateLine[sys_Association, boundary_, from_, {lo_, hi_},
    coefficients_List, opts:OptionsPattern[]] :=
  DiffExp2`API`LineIntegral[sys, boundary, from, {lo, hi}, coefficients, opts];

(* ---- honest epsilon-window accessors ---- *)

epsValue[obj_] := Which[
  DiffExp2`EpsSeries`ESQ[obj], obj,
  AssociationQ[obj] && KeyExistsQ[obj, "Value"] &&
      DiffExp2`EpsSeries`ESQ[obj["Value"]], obj["Value"],
  True, obj];

DiffExp2`EpsilonWindow[obj_List] := DiffExp2`EpsilonWindow /@ obj;
DiffExp2`EpsilonWindow[obj_] := Module[{value = epsValue[obj]},
  If[!DiffExp2`EpsSeries`ESQ[value],
    de2Error["ERR-WINDOW-READ", "object is not an EpsSeries or evaluation record"]];
  DiffExp2`EpsSeries`ESWindow[value]];

DiffExp2`EpsilonCoefficient[obj_List, k_Integer] :=
  DiffExp2`EpsilonCoefficient[#, k] & /@ obj;
DiffExp2`EpsilonCoefficient[obj_, k_Integer] := Module[{value = epsValue[obj]},
  If[!DiffExp2`EpsSeries`ESQ[value],
    de2Error["ERR-WINDOW-READ", "object is not an EpsSeries or evaluation record"]];
  DiffExp2`EpsSeries`ESCoefficient[value, k]];

DiffExp2`EpsilonCoefficientList[obj_List, k1_Integer, k2_Integer] :=
  DiffExp2`EpsilonCoefficientList[#, k1, k2] & /@ obj;
DiffExp2`EpsilonCoefficientList[obj_, k1_Integer, k2_Integer] :=
 Module[{value = epsValue[obj]},
  If[!DiffExp2`EpsSeries`ESQ[value],
    de2Error["ERR-WINDOW-READ", "object is not an EpsSeries or evaluation record"]];
  DiffExp2`EpsSeries`ESCoefficientList[value, k1, k2]];

End[];
EndPackage[];
