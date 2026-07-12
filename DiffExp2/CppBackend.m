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
RunPersistentRequest::usage = "RunPersistentRequest[schema1Request, metadata] executes a recurrence through the persistent schema-2 session, preparing immutable chart/operator and SCC data once and sending only run-dependent frames on later calls.";
RunPersistentRequests::usage = "RunPersistentRequests[schema1Requests, metadata, threads] executes several runs sharing one retained operator through the persistent native worker pool and returns ordered per-run responses.";
RunPersistentRequestGroups::usage = "RunPersistentRequestGroups[groups, threads] prepares several chart groups in one solver session and executes all of their dynamic runs through one ordered session.solve_many worker pool. Each group contains Requests and Metadata.";
PreparePersistentSCC::usage = "PreparePersistentSCC[groups, manifest] prepares one retained chart for each Requests/Metadata group, binds those charts and one full-parent physical q/C equation owner into a typed schema-2 SCC manifest, and returns an opaque session-owned SCC handle.";
PersistentSCCStatistics::usage = "PersistentSCCStatistics[handle] returns statistics and exact metadata for a retained native SCC chart.";
RunPersistentSCCColumn::usage = "RunPersistentSCCColumn[handle, seed, targets, checkpointIdentity] executes one strict schema-2 native SCC basis column and returns an opaque retained local-solution summary. seed and every target contain exactly block, run, and metadata; no coefficient slab is returned.";
RunPersistentSCCColumns::usage = "RunPersistentSCCColumns[handle, columns, threads] executes an ordered all-or-nothing batch of retained SCC basis columns through one native worker pool. Each column contains Seed, Targets, and CheckpointIdentity; successful locals are retained atomically and no coefficient slab is returned.";
SpecializePersistentRationalSCCColumn::usage = "SpecializePersistentRationalSCCColumn[source,targetSCC,rationalShadowIdentity,checkpointIdentity] imports one completed exact Rational SCC column into the paired target Acb SCC. Native code proves the shared domain-independent owner identity, physical q/C containment, exact geometry/prescriptions, and public epsilon/Taylor contract before coefficientwise specialization.";
ReleasePersistentSCC::usage = "ReleasePersistentSCC[handle] releases one retained native SCC chart and removes its Wolfram collision certificate.";
RunPersistentLocalSolve::usage = "RunPersistentLocalSolve[schema1Request, metadata, localMetadata] executes recurrence plus retained native assembly and returns an opaque session-owned local-solution handle without returning its coefficient slab.";
EvaluatePersistentLocal::usage = "EvaluatePersistentLocal[handle, point, options, outputDigits] evaluates a retained native local solution at the JSON-ready exact rational point record. The handle is the response returned by RunPersistentLocalSolve or an association containing session/local keys.";
CertifyPersistentLocalResidual::usage = "CertifyPersistentLocalResidual[handle, request, outputDigits] evaluates a retained local and certifies the native Acb residual of its immutable physical q(t,eps), C(t,eps), and homogeneous-source payload. request carries point, relative_tolerance, and the operator, owner-signature, physical-payload, source, local-checkpoint, analytic-metadata, and provenance identities advertised by residual_binding; caller-supplied operators or sources are rejected.";
RunPersistentLocalMatch::usage = "RunPersistentLocalMatch[basis, incoming, request] matches retained exact-rational regular locals entirely in their common native session and retains the exact lattice transformation and Laurent weights. request binds chart/checkpoint identities, local match points, the work epsilon window, and the required residual CompleteMax.";
RunPersistentAcbLocalMatch::usage = "RunPersistentAcbLocalMatch[basis, incoming, request] evaluates retained Acb locals at one exact physical point and retains an exact-Rational-lattice-guided Laurent match with bounded refinement and residual diagnostics. request supplies exact_lattice and refinement records in addition to chart, branch, checkpoint, point, and epsilon identities.";
RunPersistentPlannedMatch::usage = "RunPersistentPlannedMatch[plan, arm, index, basis, incoming, policy] performs one plan-derived retained Rational or Acb match. The one-based match index selects exact physical/local coordinates and branch data from the plan; policy supplies epsilon/checkpoint fields and, for Acb, exact_lattice/refinement.";
MaterializePersistentLocalMatch::usage = "MaterializePersistentLocalMatch[match, checkpointIdentity] combines a retained plan-driven match's Laurent weights with its strongly owned receiving basis entirely in C++ and returns the next opaque retained local without coefficient JSON.";
ApplyPersistentRationalRow::usage = "ApplyPersistentRationalRow[local, row, checkpointIdentity] applies one JSON-ready exact rational coefficient row to a retained Rational or Acb vector local entirely in C++. The prepared row carries one multiplier per active zero-based component, exact epsilon/Taylor coverage, and structural identities; the result is an opaque retained scalar local with unchanged chart prescriptions and no coefficient slab.";
PersistentLocalMatchStatistics::usage = "PersistentLocalMatchStatistics[handle] returns the opaque summary and exact provenance of one retained native local match.";
ReleasePersistentLocalMatch::usage = "ReleasePersistentLocalMatch[handle] releases one retained native local match state. A second release is a loud native error.";
RunPersistentEndpointLimit::usage = "RunPersistentEndpointLimit[local, request] is the unplanned low-level endpoint seam. It applies the native sector endpoint gate to a retained local using caller-supplied approach_direction and optional rim, and returns an opaque session-owned endpoint-result handle. Prefer RunPersistentPlannedEndpointLimit for retained transport arms.";
RunPersistentPlannedEndpointLimit::usage = "RunPersistentPlannedEndpointLimit[plan, arm, local, policy] applies the native endpoint gate to the final retained local of arm \"lower\" or \"upper\". The retained tile plan alone supplies the endpoint, final chart, local approach side, and nullable exact prescription rim; policy contains exactly checkpoint_identity and cancellation. No caller direction/rim or coefficient slab is accepted.";
RunPersistentWeightedPlannedEndpointLimit::usage = "RunPersistentWeightedPlannedEndpointLimit[plan, arm, local, row, checkpointIdentity, cancellation] first applies one exact prepared rational row to the retained vector local in C++, then applies the plan-bound endpoint gate to that retained scalar projection. Distinct deterministic row/endpoint checkpoint identities are derived from checkpointIdentity. The projected local is never serialized, is released publicly only after the endpoint strongly owns it, and is cleaned up on endpoint failure.";
PersistentEndpointStatistics::usage = "PersistentEndpointStatistics[handle] returns the opaque endpoint summary, analytic-regularization/branch provenance, and export counters.";
ExportPersistentEndpoint::usage = "ExportPersistentEndpoint[handle, checkpointIdentity, outputDigits] explicitly exports the retained specialized Acb epsilon vector for final compatibility.";
ReleasePersistentEndpoint::usage = "ReleasePersistentEndpoint[handle] releases one retained native endpoint result. A second release is a loud native error.";
CreatePersistentTilePlan::usage = "CreatePersistentTilePlan[owner, lower, upper, checkpointIdentity, divisionOrder] asks the persistent C++ session to build and retain exact lower/upper arm plans from retained chart handles. Exact match coordinates and physical/local tile intervals come from the native planner; DivisionOrder defaults to 3.";
CreatePersistentArmTilePlan::usage = "CreatePersistentArmTilePlan[owner, arm, checkpointIdentity, divisionOrder] asks the persistent C++ session to build and retain one exact lower or upper arm plan without fabricating an opposite arm. The native planner derives the arm name from its exact direction; DivisionOrder defaults to 3.";
PersistentTilePlanStatistics::usage = "PersistentTilePlanStatistics[handle] returns the retained immutable lower/upper native plan, including exact match and tile intervals plus branch/prescription provenance.";
PersistentTileMatchInterval::usage = "PersistentTileMatchInterval[handle, arm, index] returns one exact native match handoff for arm \"lower\" or \"upper\". Wolfram indices are one-based.";
PersistentTileIntegrationInterval::usage = "PersistentTileIntegrationInterval[handle, arm, index] returns one exact physical/local tile interval selected by the retained native plan. Wolfram indices are one-based.";
ReleasePersistentTilePlan::usage = "ReleasePersistentTilePlan[handle] releases one public native tile-plan token. Already-retained line results keep strong ownership of their immutable plan snapshot.";
RunPersistentTileIntegral::usage = "RunPersistentTileIntegral[plan, arm, tile, local, epsilon, checkpointIdentity, certifyTail] integrates the retained local over the exact tile selected by plan and applies the exact affine Jacobian. With certifyTail True (the seventh Boolean argument), an attached regular-tail model may promote the result to FullLocalWithCertifiedTail; unsupported or inconclusive requests remain StoredTruncation. Wolfram tile indices are one-based.";
RunPersistentNativeArms::usage = "RunPersistentNativeArms[plan, anchor, arms, epsilon, checkpointRoot, refinement, certifyTail] marches both retained tile-plan arms concurrently inside one persistent C++ request. arms has exactly lower/upper records, each with receiving_basis (one retained-local basis list per match) and one precomputed integrand_rows entry per tile. epsilon contains min, max, required_complete_max for projected/public lines, and match_required_complete_max for the source/match halo. C++ derives every live match window and, for ordinary Acb bases, certifies identity saturation from the actual evaluated basis (exact-zero negative powers and a certified full-rank epsilon-zero frame); it applies each rational row to a hidden scalar local, matches the unprojected vector local, integrates and aggregates natively, and returns only opaque final-local/line handles plus timing.";
RunPersistentTransportArms::usage = "RunPersistentTransportArms[plan,anchor,arms,epsilon,checkpointRoot,refinement] atomically marches both retained lower/upper arms in exactly two native workers and returns exactly two public transport-state tokens. arms has exactly lower/upper records, each containing only receiving_basis (one retained-local basis list per match). Each state strongly owns its hidden match/materialized-local/final-local closure. A returned final_local record is dependency-only: it is not a public local token and must not be passed to ReleasePersistentLocal; release the corresponding transport state instead.";
RunPersistentConsumingTransportArms::usage = "RunPersistentConsumingTransportArms[plan,anchor,arms,epsilon,checkpointRoot,refinement] marches lower then upper sequentially and consumes each receiving-basis handle after its last successful exact match. Returned compact v4 transport states retain immutable match/basis provenance and materialized tile locals, but not receiving-basis coefficient slabs.";
ConsumePersistentTransportHop::usage = "ConsumePersistentTransportHop[plan,arm,index,basis,incoming,epsilon,checkpointRoot,refinement] performs one positive one-based plan match, materializes and seals its next local, then consumes the complete receiving basis. epsilon has exactly min,max,required_complete_max. The response publishes next_local plus compact basis_reference and match_reference records for final state publication.";
PublishPersistentConsumedTransportStates::usage = "PublishPersistentConsumedTransportStates[plan,anchor,arms,epsilon,checkpointRoot,refinement] atomically publishes lower/upper compact certificate-only states after streaming hops. Each arm contains only ordered tile_sources (including the common anchor); C++ validates the already-sealed per-hop lineage directly, without echoing basis/operator/match provenance through Wolfram. Non-anchor tile-source public tokens are consumed only after both states validate.";
RunPersistentTransportArm::usage = "RunPersistentTransportArm[plan,arm,anchor,receivingBasis,epsilon,checkpointRoot,refinement] marches one retained lower or upper arm entirely in C++ without projecting or integrating observables. It returns an opaque transport-state handle that strongly owns its plan, anchor, receiving bases, hidden planned matches, one unprojected source local per tile, and final local; no coefficient slab is serialized.";
ContractPersistentTransportObservables::usage = "ContractPersistentTransportObservables[state,observables,checkpointRoot] contracts an ordered list of zero, one, or many scalar observables against one retained native transport-arm state without rematching. Each observable has exactly Identity, CheckpointIdentity, IntegrandRows, Epsilon, and TailPolicy; Epsilon has exactly Min, Max, and RequiredCompleteMax. IntegrandRows contains one prepared rational row per retained tile. TailPolicy is \"stored\", \"attempt\", or \"require\": stored never requests tail certification, attempt may remain stored-truncation, and require fails atomically unless every tile aggregates with a certified full-local tail. The result retains input order and returns directly usable opaque line handles; an empty observable list succeeds without publishing lines.";
ContractPersistentTransportPairObservables::usage = "ContractPersistentTransportPairObservables[lowerState,upperState,observables,checkpointRoot] contracts an ordered list of zero, one, or many scalar observables against exact retained lower/upper states in one native paired request. Each observable has exactly Identity, CheckpointIdentity, LowerIntegrandRows, UpperIntegrandRows, Epsilon, and optionally TailPolicy; omitted TailPolicy means \"stored\". Epsilon has exactly Min, Max, and RequiredCompleteMax. The wrapper accepts no caller signs, arms, points, rims, or cancellation data: native transport.contract_pair always combines -lower+upper. Results retain request order and are opaque paired line handles.";
RunPersistentTransportEndpointBatch::usage = "RunPersistentTransportEndpointBatch[state,observables,checkpointRoot] atomically contracts an ordered list of zero, one, or many prepared scalar rows against the final retained local of one native transport-arm state and returns opaque endpoint handles. Each observable has exactly Identity, CheckpointIdentity, IntegrandRow, and Epsilon; Epsilon has exactly Min, Max, and RequiredCompleteMax. The retained state and plan derive the arm, endpoint, local coordinate, approach direction, and analytic prescription; callers cannot override them.";
PersistentTransportArmStatistics::usage = "PersistentTransportArmStatistics[state] returns the opaque retained arm-state topology, exact provenance, ownership counts, epsilon/refinement contract, final-local handle, and statistics.";
ReleasePersistentTransportArm::usage = "ReleasePersistentTransportArm[state] releases one public transport-state token. A second release is a loud native error; independently published final locals remain governed by their own tokens.";
PersistentLineIntegralStatistics::usage = "PersistentLineIntegralStatistics[handle] returns one retained physical-tile integral summary, exact provenance, stored-or-certified-tail scope diagnostics, and export counters.";
ExportPersistentLineIntegral::usage = "ExportPersistentLineIntegral[handle, checkpointIdentity, outputDigits] explicitly exports one retained physical-tile epsilon vector for compatibility.";
ReleasePersistentLineIntegral::usage = "ReleasePersistentLineIntegral[handle] releases one retained line-integral result. A second release is a loud native error.";
PersistentLocalStatistics::usage = "PersistentLocalStatistics[handle] returns statistics and exact metadata for a retained native local solution.";
ReleasePersistentLocal::usage = "ReleasePersistentLocal[handle] releases one retained native local solution. A second release is a loud native error.";
SavePersistentCheckpoint::usage = "SavePersistentCheckpoint[owner, path, identity] atomically writes an opaque versioned native checkpoint for prepared charts/SCCs, retained Rational or Acb locals, ordinary and plan-driven exact/Acb matches, match-materialized locals, raw and plan-bound endpoint results, exact tile plans, retained transport-arm states, and completed line results. The complete strong-owner closure is serialized separately from public registry visibility; symbolic locals remain the sole deferred scalar kind.";
RestorePersistentCheckpoint::usage = "RestorePersistentCheckpoint[path, expectedIdentity] validates and restores an opaque native checkpoint into a new persistent C++ session without replaying Wolfram preprocessing. It returns stable public chart, SCC, local, ordinary/planned-match, endpoint, tile-plan, transport-state, and line-result handle maps while restoring dependency-only owners without making them public.";
ClosePersistentSession::usage = "ClosePersistentSession[owner] closes exactly one persistent native session and removes its process-local restored-session registration. owner may be a session token or any opaque handle carrying that session.";
ReleasePersistentPreparedToken::usage = "ReleasePersistentPreparedToken[token] releases retained native charts certified by one prepared-operator token and removes its collision certificate.";
ClearPersistentSessions::usage = "ClearPersistentSessions[] closes every process-local native solver session owned by this Wolfram kernel and clears its chart and SCC handle registries.";
PersistentSessionInformation::usage = "PersistentSessionInformation[] returns native statistics for the live persistent solver sessions owned by this Wolfram kernel.";
DecodeScalar::usage = "DecodeScalar[encoded, precision] reconstructs a Wolfram scalar from a C++ rational or Acb midpoint/radius record with honest Accuracy.";
DecodeScalars::usage = "DecodeScalars[encodedList, precision] reconstructs a list of C++ rational or Acb scalar records in bounded parser batches, preserving DecodeScalar semantics and falling back elementwise for mixed or malformed input.";
ResetBackend::usage = "ResetBackend[] unloads the cached LibraryFunction handles so a rebuilt library can be loaded.";

Begin["`Private`"];

$moduleDirectory = DirectoryName[$InputFileName];
$repositoryRoot = ParentDirectory[$moduleDirectory];
$backendLibrary = None;
$runFunction = None;
$infoFunction = None;
$persistentSessionCache = <||>;
$persistentChartCache = <||>;
$persistentSCCCache = <||>;
$persistentPreparedTokenCache = <||>;
$persistentRestoredSessionHandles = <||>;
$persistentChartCacheMax = 1024;
$persistentSCCCacheMax = 128;

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
  Quiet[ClearPersistentSessions[]];
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

decimalRecord[x_?InexactNumberQ, digits_Integer] := Module[
  {available, used, y, sign, ds, exponent, rd, tail},
  available = Precision[x];
  used = If[NumericQ[available], Max[1, Min[digits, Floor[available]]], digits];
  (* A BigReal may report p reliable binary-derived decimal digits while
     RealDigits at exactly Floor[p] exposes one or more uncertain tail cells
     as Indeterminate.  Retry monotonically at a shorter honest midpoint;
     never replace those cells or stamp the requested precision onto them. *)
  While[used >= 1,
    y = N[x, used];
    If[TrueQ[y == 0], Return[<|"String" -> "0", "Digits" -> used,
      "Exponent" -> 0, "DecimalErrorExponent" -> -Infinity|>, Module]];
    sign = If[TrueQ[y < 0], "-", ""];
    rd = Quiet[Check[RealDigits[Abs[y], 10, used], $Failed]];
    If[MatchQ[rd, {_List, _Integer}],
      {ds, exponent} = rd;
      If[VectorQ[ds, IntegerQ],
        tail = If[Rest[ds] === {}, "0",
          StringJoin[ToString /@ Rest[ds]]];
        Return[<|"String" -> sign <> ToString[First[ds]] <> "." <>
            tail <> "e" <> ToString[exponent - 1],
          "Digits" -> used, "Exponent" -> exponent,
          (* The decimal midpoint differs from x by at most one unit at
             the first omitted base-10 place. *)
          "DecimalErrorExponent" -> exponent - used|>, Module]]];
    used--];
  Failure["UnsupportedScalar", <|"Scalar" -> x, "Detail" ->
    "could not obtain a finite fixed-precision decimal expansion"|>]];

decimalString[x_?InexactNumberQ, digits_Integer] := Module[{record},
  record = decimalRecord[x, digits];
  If[FailureQ[record], record, record["String"]]];

arbInexactString[x_?InexactNumberQ, digits_Integer] := Module[
  {record, midpoint, accuracy, sourceRadiusExponent,
   radiusExponent, radius},
  record = decimalRecord[x, digits];
  If[FailureQ[record], Return[record, Module]];
  midpoint = record["String"];
  accuracy = Accuracy[x];
  If[accuracy =!= Infinity && !NumericQ[accuracy],
    Return[Failure["UnsupportedScalar", <|
    "Scalar" -> x, "Detail" ->
      "inexact coefficient has no finite accuracy estimate"|>], Module]];
  (* Arb interval syntax preserves the uncertainty already tracked by the
     Wolfram handoff.  A shortened retry midpoint may have a larger decimal
     rounding error than the source radius; include the larger exponent.
     The factor two bounds the sum of both contributions. *)
  sourceRadiusExponent = If[accuracy === Infinity,
    -Infinity, -Floor[accuracy]];
  radiusExponent = Max[sourceRadiusExponent,
    record["DecimalErrorExponent"]];
  If[radiusExponent === -Infinity, Return[midpoint, Module]];
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

(* Parsing every returned coefficient with a separate ToExpression call is
   disproportionately expensive for the large framed tensors produced by the
   native recurrence.  Keep batches bounded: this amortizes parser setup
   without constructing a single potentially enormous Wolfram expression at
   high working precision.  A failed/shape-changing batch is deliberately
   retried through DecodeScalar so malformed public input retains the scalar
   decoder's diagnostics instead of poisoning neighbouring records. *)
$decodeScalarBatchSize = 2048;

parseExpressionBatch[strings_List] := Module[{values},
  values = Quiet[Check[ToExpression[
      "{" <> StringRiffle[strings, ","] <> "}"], $Failed]];
  If[ListQ[values] && Length[values] === Length[strings], values, $Failed]];

parseDecimalBatch[strings_List, precision_Integer] := Module[
  {held, values},
  held = StringReplace[strings,
    {"e+" -> "*^", "e-" -> "*^-", "e" -> "*^"}];
  values = parseExpressionBatch[held];
  If[values === $Failed || !AllTrue[values, NumericQ], Return[$Failed, Module]];
  SetPrecision[#, precision] & /@ values];

parseRadiusExponentBatch[strings_List] := Module[
  {values = ConstantArray[-Infinity, Length[strings]], positions, parsed},
  positions = Flatten[Position[strings, s_String /; s =!= "zero", {1}]];
  If[positions === {}, Return[values, Module]];
  parsed = parseExpressionBatch[strings[[positions]]];
  If[parsed === $Failed || !AllTrue[parsed, IntegerQ], Return[$Failed, Module]];
  values[[positions]] = parsed;
  values];

decodeExactBatch[strings_List, precision_Integer] := Module[{values},
  values = parseExpressionBatch[strings];
  If[values === $Failed,
    DecodeScalar[#, precision] & /@ strings,
    values]];

decodeAcbBatch[records_List, precision_Integer] := Module[
  {n = Length[records], midpoints, radii, re, im, reRadius, imRadius,
   decodedRe, decodedIm},
  midpoints = parseDecimalBatch[
    Join[records[[All, 1]], records[[All, 2]]], precision];
  radii = parseRadiusExponentBatch[
    Join[records[[All, 3]], records[[All, 4]]]];
  If[midpoints === $Failed || radii === $Failed,
    Return[DecodeScalar[#, precision] & /@ records, Module]];
  re = Take[midpoints, n];
  im = Take[midpoints, -n];
  reRadius = Take[radii, n];
  imRadius = Take[radii, -n];
  decodedRe = MapThread[applyBallAccuracy[#1, #2, precision] &,
    {re, reRadius}];
  decodedIm = MapThread[applyBallAccuracy[#1, #2, precision] &,
    {im, imRadius}];
  decodedRe + I decodedIm];

DecodeScalars[encoded_List, precision_Integer] := Module[{decoder},
  If[encoded === {}, Return[{}, Module]];
  decoder = Which[
    AllTrue[encoded, StringQ], decodeExactBatch,
    AllTrue[encoded, MatchQ[#, {_String, _String, _String, _String}] &],
      decodeAcbBatch,
    True, Function[{batch, ignoredPrecision},
      DecodeScalar[#, ignoredPrecision] & /@ batch]];
  Join @@ (decoder[#, precision] & /@
    Partition[encoded, UpTo[$decodeScalarBatchSize]])];

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

(* ---- persistent schema-2 recurrence sessions -------------------------
   A recurrence request has a large immutable prepared operator and a much
   smaller per-column/per-source state.  Retain the former as typed FLINT/Arb
   data in C++ and send only the latter after the first solve.  Hashes are
   cache indices only: every hit is checked against its complete signature. *)

$persistentStaticKeys = {"domain", "symbols", "precision_bits", "d", "fb",
  "w", "d_lags", "denominators", "nhat_lags", "d0_inverse", "blocks",
  "assembly", "physical_ode", "chop_digits"};
$persistentRequiredStaticKeys = DeleteCases[
  $persistentStaticKeys, "physical_ode"];
$persistentRunKeys = {"nmax", "p", "has_initial", "adaptive_probe",
  "a_target", "b_target", "a_shift_min", "a_shifts", "schedule",
  "initial", "initial_validity", "source", "return_u"};

persistentCacheLookup[cache_Association, key_, signature_, label_String] :=
  Module[{entry = Lookup[cache, key, None]},
    Which[
      entry === None, None,
      AssociationQ[entry] && SameQ[entry["Signature"], signature], entry,
      True, Failure["CppBackend", <|"Detail" ->
        ("persistent " <> label <>
          " cache-key collision with unequal full identity"),
        "Key" -> key|>]]];

persistentCommandOKQ[response_] := AssociationQ[response] &&
  Lookup[response, "status", "error"] === "ok";

persistentIdentityString[value_] := ToString[value, InputForm];

persistentPreparedTokenCertificate[token_String, static_Association] := Module[
  {certified = Lookup[$persistentPreparedTokenCache, token, None]},
  Which[
    certified === None,
      AssociateTo[$persistentPreparedTokenCache, token -> static]; token,
    SameQ[certified, static], token,
    True, Failure["CppBackend", <|"Detail" ->
      "persistent prepared token was reused with unequal static operator data",
      "PreparedToken" -> token|>]]];

persistentReleaseChartKey[key_] := Module[
  {entry = Lookup[$persistentChartCache, key, None], token},
  If[!AssociationQ[entry], Return[Null, Module]];
  Quiet[Check[RunRequest[<|"schema" -> 2, "op" -> "chart.release",
      "session" -> entry["Session"], "chart" -> entry["Handle"]|>], Null]];
  token = Lookup[entry, "PreparedToken", None];
  KeyDropFrom[$persistentChartCache, key];
  If[StringQ[token] && !AnyTrue[Values[$persistentChartCache],
      Lookup[#, "PreparedToken", None] === token &],
    KeyDropFrom[$persistentPreparedTokenCache, token]];
  Null];

persistentTouchChart[key_, entry_Association] := (
  KeyDropFrom[$persistentChartCache, key];
  AssociateTo[$persistentChartCache, key -> entry];
  entry);

ReleasePersistentPreparedToken[token_String] := Module[{keys},
  keys = Keys@Select[$persistentChartCache,
    Lookup[#, "PreparedToken", None] === token &];
  Scan[persistentReleaseChartKey, keys];
  KeyDropFrom[$persistentPreparedTokenCache, token];
  Null];

persistentForgetSessionHandle[handle_String] := Module[
  {sessionKeys, chartKeys, sccKeys, activeTokens},
  sessionKeys = Keys@Select[$persistentSessionCache,
    Lookup[#, "Handle", None] === handle &];
  chartKeys = Keys@Select[$persistentChartCache,
    Lookup[#, "Session", None] === handle &];
  sccKeys = Keys@Select[$persistentSCCCache,
    Lookup[#, "Session", None] === handle &];
  KeyDropFrom[$persistentSessionCache, sessionKeys];
  KeyDropFrom[$persistentRestoredSessionHandles, handle];
  KeyDropFrom[$persistentChartCache, chartKeys];
  KeyDropFrom[$persistentSCCCache, sccKeys];
  activeTokens = DeleteDuplicates@Select[
    Map[Lookup[#, "PreparedToken", None] &,
      Values[$persistentChartCache]], StringQ];
  $persistentPreparedTokenCache = KeyTake[
    $persistentPreparedTokenCache, activeTokens];
  Null];

persistentCloseSessionHandle[handle_String] := (
  Quiet[Check[RunRequest[<|"schema" -> 2, "op" -> "session.close",
      "session" -> handle|>], Null]];
  persistentForgetSessionHandle[handle]);

persistentCloseIncompatibleSymbolicSessions[symbols_List] := Module[
  {entries, handles},
  entries = Select[Values[$persistentSessionCache],
    Module[{sig = Lookup[#, "Signature", {}]},
      Length[sig] >= 3 && sig[[1]] === "symbolic" &&
        sig[[3]] =!= symbols] &];
  handles = DeleteDuplicates[Lookup[entries, "Handle", {}]];
  Scan[persistentCloseSessionHandle, handles];
  Null];

preparePersistentRequest[request_Association, metadata_Association] := Module[
  {missingStatic, missingRun, domain, symbols, precisionBits, outputDigits,
   systemIdentity, sessionAnalytic, chartAnalytic, scc, static, run,
   sessionSignature, sessionKey, sessionEntry, sessionResponse, session,
   chartIdentity, chartIdentityString, chartSignature, chartKey, chartEntry,
   chartResponse, chart,
   preparedToken, createRequest, prepareRequest},
  missingStatic = Complement[$persistentRequiredStaticKeys, Keys[request]];
  missingRun = Complement[$persistentRunKeys, Keys[request]];
  If[missingStatic =!= {} || missingRun =!= {},
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent recurrence request is missing required fields",
      "MissingStatic" -> missingStatic, "MissingRun" -> missingRun|>], Module]];
  If[!KeyExistsQ[metadata, "SystemIdentity"] ||
      !KeyExistsQ[metadata, "ChartIdentity"] ||
      !AssociationQ[Lookup[metadata, "SCC", None]],
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent recurrence metadata requires SystemIdentity, ChartIdentity, and SCC"|>], Module]];
  domain = request["domain"];
  symbols = Lookup[request, "symbols", {}];
  precisionBits = Lookup[request, "precision_bits", 256];
  outputDigits = Lookup[request, "output_digits", 50];
  systemIdentity = metadata["SystemIdentity"];
  sessionAnalytic = Lookup[metadata, "SessionAnalytic", <||>];
  chartAnalytic = Lookup[metadata, "ChartAnalytic", <||>];
  scc = metadata["SCC"];

  sessionSignature = {domain, precisionBits, symbols, systemIdentity,
    sessionAnalytic};
  (* SymbolicRational owns one FLINT multivariate context.  Values prepared
     in a different variable field must be destroyed before that context can
     be reconfigured; keep same-field sessions, close only incompatible ones. *)
  If[domain === "symbolic",
    persistentCloseIncompatibleSymbolicSessions[symbols]];
  sessionKey = Hash[sessionSignature, "SHA256"];
  sessionEntry = persistentCacheLookup[$persistentSessionCache, sessionKey,
    sessionSignature, "session"];
  If[FailureQ[sessionEntry], Return[sessionEntry, Module]];
  If[sessionEntry === None,
    createRequest = Join[<|"schema" -> 2, "op" -> "session.create",
        "domain" -> domain, "output_digits" -> outputDigits,
        "chart_capacity" -> $persistentChartCacheMax,
        "scc_capacity" -> $persistentSCCCacheMax,
        "analytic" -> sessionAnalytic|>,
      If[domain === "acb", <|"precision_bits" -> precisionBits|>, <||>],
      If[domain === "symbolic", <|"symbols" -> symbols|>, <||>]];
    sessionResponse = RunRequest[createRequest];
    If[!persistentCommandOKQ[sessionResponse], Return[sessionResponse, Module]];
    session = sessionResponse["session"];
    AssociateTo[$persistentSessionCache, sessionKey -> <|
      "Signature" -> sessionSignature, "Handle" -> session|>],
    session = sessionEntry["Handle"]];

  static = KeyTake[request, $persistentStaticKeys];
  run = KeyTake[request, $persistentRunKeys];
  chartIdentity = metadata["ChartIdentity"];
  chartIdentityString = persistentIdentityString[chartIdentity];
  preparedToken = Lookup[metadata, "PreparedToken", None];
  (* Solve's static-payload cache supplies a stable content token only after
     proving equality against its complete prepared tensor.  This module also
     binds every token to the full encoded static record, so even a direct
     caller cannot reuse a token for unequal operator data.  Callers without
     such a token retain the full-signature path. *)
  If[StringQ[preparedToken],
    preparedToken = persistentPreparedTokenCertificate[
      preparedToken, static];
    If[FailureQ[preparedToken], Return[preparedToken, Module]]];
  (* The native chart's public exact identity is the same collision-certified
     immutable-operator token embedded in physical_ode.  The richer Wolfram
     ChartIdentity remains in chartSignature as an independently compared
     structural cache certificate. *)
  If[StringQ[preparedToken], chartIdentityString = preparedToken];
  chartSignature = If[StringQ[preparedToken],
    {session, preparedToken, chartIdentity, chartAnalytic, scc},
    {session, static, chartIdentity, chartAnalytic, scc}];
  chartKey = Hash[chartSignature, "SHA256"];
  chartEntry = persistentCacheLookup[$persistentChartCache, chartKey,
    chartSignature, "chart"];
  If[FailureQ[chartEntry], Return[chartEntry, Module]];
  If[chartEntry === None,
    If[Length[$persistentChartCache] >= $persistentChartCacheMax,
      persistentReleaseChartKey[First[Keys[$persistentChartCache]]]];
    prepareRequest = <|"schema" -> 2, "op" -> "chart.prepare",
      "session" -> session,
      "key" -> ("chart:" <> IntegerString[chartKey, 16, 64]),
      "identity" -> chartIdentityString,
      "analytic" -> chartAnalytic, "scc" -> scc,
      "problem" -> static|>;
    chartResponse = RunRequest[prepareRequest];
    If[!persistentCommandOKQ[chartResponse], Return[chartResponse, Module]];
    chart = chartResponse["chart"];
    AssociateTo[$persistentChartCache, chartKey -> <|
      "Signature" -> chartSignature, "Session" -> session,
      "Handle" -> chart, "Identity" -> chartIdentityString,
      "PreparedToken" -> preparedToken|>],
    chartEntry = persistentTouchChart[chartKey, chartEntry];
    chart = chartEntry["Handle"];
    chartIdentityString = Lookup[chartEntry, "Identity",
      chartIdentityString]];

  <|"Session" -> session, "Chart" -> chart, "Run" -> run,
    "ChartIdentity" -> chartIdentityString,
    "Static" -> static, "OutputDigits" -> outputDigits|>];

RunPersistentRequest[request_Association, metadata_Association] := Module[
  {prepared = preparePersistentRequest[request, metadata]},
  If[FailureQ[prepared] || !AssociationQ[prepared] ||
      !KeyExistsQ[prepared, "Session"], Return[prepared, Module]];
  RunRequest[<|"schema" -> 2, "op" -> "chart.solve",
    "session" -> prepared["Session"], "chart" -> prepared["Chart"],
    "output_digits" -> prepared["OutputDigits"],
    "run" -> prepared["Run"]|>]];

RunPersistentRequests[requests_List, metadata_Association,
    threads_Integer] := Module[
  {first, prepared, static, outputDigits, incompatible, runs, response},
  If[requests === {}, Return[<|"status" -> "ok", "results" -> {}|>, Module]];
  If[threads < 1 || !AllTrue[requests, AssociationQ],
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent recurrence batch requires associations and a positive thread count"|>], Module]];
  first = First[requests];
  static = KeyTake[first, $persistentStaticKeys];
  outputDigits = Lookup[first, "output_digits", 50];
  incompatible = Select[Rest[requests],
    !SameQ[KeyTake[#, $persistentStaticKeys], static] ||
      Lookup[#, "output_digits", 50] =!= outputDigits &];
  If[incompatible =!= {},
    Return[Failure["CppBackend", <|"Detail" ->
      "one persistent recurrence batch cannot mix prepared operators or output precision"|>], Module]];
  prepared = preparePersistentRequest[first, metadata];
  If[FailureQ[prepared] || !AssociationQ[prepared] ||
      !KeyExistsQ[prepared, "Session"], Return[prepared, Module]];
  runs = KeyTake[#, $persistentRunKeys] & /@ requests;
  response = RunRequest[<|"schema" -> 2, "op" -> "chart.solve_batch",
    "session" -> prepared["Session"], "chart" -> prepared["Chart"],
    "output_digits" -> outputDigits, "threads" -> threads,
    "runs" -> runs|>];
  response];

preparePersistentRequestGroup[group_Association] := Module[
  {requests = Lookup[group, "Requests", None],
   metadata = Lookup[group, "Metadata", None], first, static,
   outputDigits, incompatible, prepared, jobs},
  If[!ListQ[requests] || requests === {} ||
      !AllTrue[requests, AssociationQ] || !AssociationQ[metadata],
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent request group requires nonempty Requests and Metadata"|>],
      Module]];
  first = First[requests];
  static = KeyTake[first, $persistentStaticKeys];
  outputDigits = Lookup[first, "output_digits", 50];
  incompatible = Select[Rest[requests],
    !SameQ[KeyTake[#, $persistentStaticKeys], static] ||
      Lookup[#, "output_digits", 50] =!= outputDigits &];
  If[incompatible =!= {},
    Return[Failure["CppBackend", <|"Detail" ->
      "one persistent request group cannot mix prepared operators or output precision"|>],
      Module]];
  prepared = preparePersistentRequest[first, metadata];
  If[FailureQ[prepared] || !AssociationQ[prepared] ||
      !KeyExistsQ[prepared, "Session"], Return[prepared, Module]];
  jobs = Map[<|"chart" -> prepared["Chart"],
      "run" -> KeyTake[#, $persistentRunKeys],
      "output_digits" -> outputDigits|> &, requests];
  <|"Session" -> prepared["Session"], "Chart" -> prepared["Chart"],
    "ChartIdentity" -> prepared["ChartIdentity"], "Jobs" -> jobs|>];

RunPersistentRequestGroups[groups_List, threads_Integer] := Module[
  {prepared, failures, sessions, jobs},
  If[groups === {}, Return[<|"status" -> "ok", "results" -> {}|>, Module]];
  If[threads < 1 || !AllTrue[groups, AssociationQ],
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent request groups require associations and a positive thread count"|>],
      Module]];
  prepared = preparePersistentRequestGroup /@ groups;
  failures = Select[prepared, FailureQ[#] || !AssociationQ[#] ||
      !KeyExistsQ[#, "Session"] &];
  If[failures =!= {}, Return[First[failures], Module]];
  sessions = DeleteDuplicates[prepared[[All, "Session"]]];
  If[Length[sessions] =!= 1,
    Return[Failure["CppBackend", <|"Detail" ->
      "one cross-chart persistent batch must belong to exactly one solver session",
      "Sessions" -> sessions|>], Module]];
  jobs = Flatten[prepared[[All, "Jobs"]], 1];
  RunRequest[<|"schema" -> 2, "op" -> "session.solve_many",
    "session" -> First[sessions], "threads" -> threads,
    "jobs" -> jobs|>]];

(* Native SCC keys are hashes of the complete JSON-bound manifest, including
   the actual retained chart handles and exact chart identity strings.  Sort
   object keys recursively so caller Association insertion order cannot
   perturb that key. *)
persistentCanonicalJSONValue[value_Association] := Association@Map[
  Function[key, key -> persistentCanonicalJSONValue[value[key]]],
  Sort[Keys[value]]];
persistentCanonicalJSONValue[value_List] :=
  persistentCanonicalJSONValue /@ value;
persistentCanonicalJSONValue[value_] := value;

persistentSCCOpaqueHandle[entry_Association] := <|
  "Session" -> entry["Session"], "SCC" -> entry["Handle"],
  "Key" -> entry["NativeKey"]|>;

persistentSCCHandles[handle_Association] := Module[{session, scc},
  session = Lookup[handle, "session", Lookup[handle, "Session", None]];
  scc = Lookup[handle, "scc", Lookup[handle, "SCC", None]];
  If[!StringQ[session] || !StringQ[scc],
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent SCC handle requires exact session and SCC tokens"|>],
      Module]];
  <|"Session" -> session, "SCC" -> scc|>];

PreparePersistentSCC[groups_List, manifest_Association] := Module[
  {requiredKeys = {"identity", "parent", "blocks", "couplings",
      "physical_ode"}, allowedKeys,
   reservedBlockKeys = {"chart", "principal_identity"}, blocks,
   badBlockPositions, preparedGroups = {}, prepared, sessions, session,
   filledBlocks, filledManifest, canonicalManifest, manifestJSON,
   nativeKey, cacheKey, cacheSignature, cacheEntry, response, scc},
  If[groups === {} || !AllTrue[groups, AssociationQ],
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent SCC preparation requires a nonempty list of request groups"|>],
      Module]];
  allowedKeys = Append[requiredKeys, "rational_shadow_identity"];
  If[!MemberQ[Sort /@ {requiredKeys, allowedKeys}, Sort[Keys[manifest]]],
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent SCC manifest must contain exactly the lowercase native fields",
      "ExpectedKeys" -> {requiredKeys, allowedKeys},
      "ActualKeys" -> Keys[manifest]|>],
      Module]];
  blocks = manifest["blocks"];
  If[!StringQ[manifest["identity"]] ||
      StringLength[manifest["identity"]] == 0 ||
      !StringQ[Lookup[manifest, "rational_shadow_identity",
        manifest["identity"]]] ||
      StringLength[Lookup[manifest, "rational_shadow_identity",
        manifest["identity"]]] == 0 ||
      !AssociationQ[manifest["parent"]] || !ListQ[blocks] ||
      !AllTrue[blocks, AssociationQ] || !ListQ[manifest["couplings"]] ||
      !AssociationQ[manifest["physical_ode"]],
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent SCC identity, parent, blocks, couplings, or physical ODE has the wrong native JSON type"|>],
      Module]];
  If[Length[blocks] =!= Length[groups],
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent SCC requires exactly one request group per manifest block",
      "Groups" -> Length[groups], "Blocks" -> Length[blocks]|>], Module]];
  badBlockPositions = Select[Range[Length[blocks]], Function[index,
    AnyTrue[reservedBlockKeys,
      Function[key, KeyExistsQ[blocks[[index]], key]]]]];
  If[badBlockPositions =!= {},
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent SCC input blocks must omit chart and principal_identity",
      "BlockPositions" -> badBlockPositions|>], Module]];

  Do[
    prepared = preparePersistentRequestGroup[group];
    If[FailureQ[prepared] || !AssociationQ[prepared] ||
        !StringQ[Lookup[prepared, "Session", None]] ||
        !StringQ[Lookup[prepared, "Chart", None]] ||
        !StringQ[Lookup[prepared, "ChartIdentity", None]],
      Return[prepared, Module]];
    preparedGroups = Append[preparedGroups, prepared],
    {group, groups}];
  sessions = DeleteDuplicates[preparedGroups[[All, "Session"]]];
  If[Length[sessions] =!= 1,
    Return[Failure["CppBackend", <|"Detail" ->
      "all persistent SCC diagonal charts must belong to one solver session",
      "Sessions" -> sessions|>], Module]];
  session = First[sessions];
  filledBlocks = MapThread[Join[#1, <|"chart" -> #2["Chart"],
        "principal_identity" -> #2["ChartIdentity"]|>] &,
    {blocks, preparedGroups}];
  filledManifest = <|"identity" -> manifest["identity"],
    "rational_shadow_identity" -> Lookup[manifest,
      "rational_shadow_identity", manifest["identity"]],
    "parent" -> manifest["parent"], "blocks" -> filledBlocks,
    "couplings" -> manifest["couplings"],
    "physical_ode" -> manifest["physical_ode"]|>;
  canonicalManifest = persistentCanonicalJSONValue[filledManifest];
  manifestJSON = Quiet[Check[ExportString[canonicalManifest, "RawJSON",
      "Compact" -> True], $Failed]];
  If[manifestJSON === $Failed,
    Return[Failure["CppBackend", <|"Detail" ->
      "could not serialize the complete persistent SCC manifest"|>], Module]];
  nativeKey = "scc-manifest:" <>
    IntegerString[Hash[manifestJSON, "SHA256"], 16, 64];
  cacheKey = Hash[{session, nativeKey}, "SHA256"];
  cacheSignature = {session, nativeKey, manifestJSON};
  cacheEntry = persistentCacheLookup[$persistentSCCCache, cacheKey,
    cacheSignature, "SCC"];
  If[FailureQ[cacheEntry], Return[cacheEntry, Module]];
  If[AssociationQ[cacheEntry],
    Return[persistentSCCOpaqueHandle[cacheEntry], Module]];
  (* SCC handles are public lifetime tokens.  Do not invalidate an older
     opaque handle by silently evicting it as the internal chart cache does. *)
  If[Length[$persistentSCCCache] >= $persistentSCCCacheMax,
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent SCC cache capacity is exhausted; release an SCC handle before preparing another",
      "Capacity" -> $persistentSCCCacheMax|>], Module]];
  response = RunRequest[Join[<|"schema" -> 2, "op" -> "scc.prepare",
      "session" -> session, "key" -> nativeKey|>, canonicalManifest]];
  If[!persistentCommandOKQ[response], Return[response, Module]];
  scc = Lookup[response, "scc", None];
  If[!StringQ[scc],
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent SCC preparation returned no native SCC handle"|>], Module]];
  cacheEntry = <|"Signature" -> cacheSignature, "Session" -> session,
    "Handle" -> scc, "NativeKey" -> nativeKey|>;
  AssociateTo[$persistentSCCCache, cacheKey -> cacheEntry];
  persistentSCCOpaqueHandle[cacheEntry]];

PersistentSCCStatistics[handle_Association] := Module[
  {tokens = persistentSCCHandles[handle]},
  If[FailureQ[tokens], Return[tokens, Module]];
  RunRequest[<|"schema" -> 2, "op" -> "scc.stats",
    "session" -> tokens["Session"], "scc" -> tokens["SCC"]|>]];

RunPersistentSCCColumn[handle_Association, seed_Association,
    targets_List, checkpointIdentity_String] := Module[
  {tokens = persistentSCCHandles[handle], entries, badEntries},
  If[FailureQ[tokens], Return[tokens, Module]];
  If[StringLength[checkpointIdentity] == 0,
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent SCC column checkpoint identity must be nonempty"|>],
      Module]];
  entries = Prepend[targets, seed];
  badEntries = Select[Range[Length[entries]], Function[index,
    !AssociationQ[entries[[index]]] ||
      Sort[Keys[entries[[index]]]] =!= Sort[{"block", "run", "metadata"}] ||
      !IntegerQ[Lookup[entries[[index]], "block", None]] ||
      Lookup[entries[[index]], "block", -1] < 0 ||
      !AssociationQ[Lookup[entries[[index]], "run", None]] ||
      Sort[Keys[Lookup[entries[[index]], "run", <||>]]] =!=
        Sort[$persistentRunKeys] ||
      !AssociationQ[Lookup[entries[[index]], "metadata", None]]]];
  If[badEntries =!= {},
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent SCC column entries require exactly a zero-based block, canonical recurrence run, and local metadata",
      "EntryPositions" -> badEntries|>], Module]];
  RunRequest[<|"schema" -> 2, "op" -> "scc.solve_column",
    "session" -> tokens["Session"], "scc" -> tokens["SCC"],
    "seed" -> seed, "targets" -> targets,
    "checkpoint_identity" -> checkpointIdentity|>]];

RunPersistentSCCColumns[handle_Association, columns_List,
    threads_Integer] := Module[
  {tokens = persistentSCCHandles[handle], native, entries, badColumns,
   badEntries},
  If[FailureQ[tokens], Return[tokens, Module]];
  If[columns === {} || threads < 1,
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent SCC column batch requires at least one column and a positive thread count"|>], Module]];
  badColumns = Select[Range[Length[columns]], Function[index,
    !AssociationQ[columns[[index]]] ||
      Sort[Keys[columns[[index]]]] =!=
        Sort[{"Seed", "Targets", "CheckpointIdentity"}] ||
      !AssociationQ[Lookup[columns[[index]], "Seed", None]] ||
      !ListQ[Lookup[columns[[index]], "Targets", None]] ||
      !StringQ[Lookup[columns[[index]], "CheckpointIdentity", None]] ||
      StringLength[Lookup[columns[[index]], "CheckpointIdentity", ""]] == 0]];
  If[badColumns =!= {},
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent SCC batch columns require exactly Seed, Targets, and a nonempty CheckpointIdentity",
      "ColumnPositions" -> badColumns|>], Module]];
  entries = Flatten[Map[
      Prepend[# ["Targets"], # ["Seed"]] &, columns], 1];
  badEntries = Select[Range[Length[entries]], Function[index,
    !AssociationQ[entries[[index]]] ||
      Sort[Keys[entries[[index]]]] =!= Sort[{"block", "run", "metadata"}] ||
      !IntegerQ[Lookup[entries[[index]], "block", None]] ||
      Lookup[entries[[index]], "block", -1] < 0 ||
      !AssociationQ[Lookup[entries[[index]], "run", None]] ||
      Sort[Keys[Lookup[entries[[index]], "run", <||>]]] =!=
        Sort[$persistentRunKeys] ||
      !AssociationQ[Lookup[entries[[index]], "metadata", None]]]];
  If[badEntries =!= {},
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent SCC batch entries require exactly a zero-based block, canonical recurrence run, and local metadata",
      "EntryPositions" -> badEntries|>], Module]];
  native = Map[<|"seed" -> # ["Seed"],
      "targets" -> # ["Targets"],
      "checkpoint_identity" -> # ["CheckpointIdentity"]|> &, columns];
  RunRequest[<|"schema" -> 2, "op" -> "scc.solve_columns",
    "session" -> tokens["Session"], "scc" -> tokens["SCC"],
    "columns" -> native, "threads" -> threads|>]];

ReleasePersistentSCC[handle_Association] := Module[
  {tokens = persistentSCCHandles[handle], response, keys},
  If[FailureQ[tokens], Return[tokens, Module]];
  response = RunRequest[<|"schema" -> 2, "op" -> "scc.release",
    "session" -> tokens["Session"], "scc" -> tokens["SCC"]|>];
  If[persistentCommandOKQ[response],
    keys = Keys@Select[$persistentSCCCache,
      Lookup[#, "Session", None] === tokens["Session"] &&
        Lookup[#, "Handle", None] === tokens["SCC"] &];
    KeyDropFrom[$persistentSCCCache, keys]];
  response];

(* A native local handle is deliberately represented only by its owning
   session and process-local local token.  Accept both the lower-case native
   response and the capitalized record used by Solve's opaque seam, but never
   guess either token from chart/cache identity. *)
persistentLocalHandles[handle_Association] := Module[{session, local},
  session = Lookup[handle, "session", Lookup[handle, "Session", None]];
  local = Lookup[handle, "local", Lookup[handle, "Local", None]];
  If[!StringQ[session] || !StringQ[local],
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent local handle requires exact session and local tokens"|>],
      Module]];
  <|"Session" -> session, "Local" -> local|>];

SpecializePersistentRationalSCCColumn[source_Association,
    targetSCC_Association, rationalShadowIdentity_String,
    checkpointIdentity_String] := Module[
  {sourceTokens = persistentLocalHandles[source],
   targetTokens = persistentSCCHandles[targetSCC], sourceCheckpoint},
  If[FailureQ[sourceTokens], Return[sourceTokens, Module]];
  If[FailureQ[targetTokens], Return[targetTokens, Module]];
  sourceCheckpoint = Lookup[source, "checkpoint_identity",
    Lookup[source, "CheckpointIdentity", None]];
  If[!StringQ[sourceCheckpoint] || StringLength[sourceCheckpoint] == 0 ||
      StringLength[rationalShadowIdentity] == 0 ||
      StringLength[checkpointIdentity] == 0,
    Return[Failure["CppBackend", <|"Detail" ->
      "Rational-shadow specialization requires nonempty source/new checkpoint and shadow identities"|>], Module]];
  RunRequest[<|"schema" -> 2,
    "op" -> "local.specialize_rational_shadow",
    "session" -> targetTokens["Session"],
    "source_session" -> sourceTokens["Session"],
    "source_local" -> sourceTokens["Local"],
    "source_checkpoint_identity" -> sourceCheckpoint,
    "target_scc" -> targetTokens["SCC"],
    "rational_shadow_identity" -> rationalShadowIdentity,
    "checkpoint_identity" -> checkpointIdentity|>]];

RunPersistentLocalSolve[request_Association, metadata_Association,
    localMetadata_Association] := Module[{prepared, response},
  prepared = preparePersistentRequest[request, metadata];
  If[FailureQ[prepared] || !AssociationQ[prepared] ||
      !KeyExistsQ[prepared, "Session"], Return[prepared, Module]];
  response = RunRequest[<|"schema" -> 2, "op" -> "local.solve",
    "session" -> prepared["Session"], "chart" -> prepared["Chart"],
    "run" -> prepared["Run"], "metadata" -> localMetadata|>];
  response];

EvaluatePersistentLocal[handle_Association, point_Association,
    options_Association:<||>, outputDigits_:Automatic] := Module[
  {tokens = persistentLocalHandles[handle], request},
  If[FailureQ[tokens], Return[tokens, Module]];
  If[outputDigits =!= Automatic &&
      (!IntegerQ[outputDigits] || outputDigits < 1),
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent local output digits must be a positive integer"|>], Module]];
  request = <|"schema" -> 2, "op" -> "local.evaluate",
    "session" -> tokens["Session"], "local" -> tokens["Local"],
    "point" -> point, "options" -> options|>;
  If[IntegerQ[outputDigits],
    request = Append[request, "output_digits" -> outputDigits]];
  RunRequest[request]];

CertifyPersistentLocalResidual[handle_Association,
    certificateRequest_Association, outputDigits_:Automatic] := Module[
  {tokens = persistentLocalHandles[handle], reserved, missing, payload,
   request},
  If[FailureQ[tokens], Return[tokens, Module]];
  If[outputDigits =!= Automatic &&
      (!IntegerQ[outputDigits] || outputDigits < 1),
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent residual output digits must be a positive integer"|>],
      Module]];
  reserved = Intersection[Keys[certificateRequest],
    {"schema", "op", "session", "local", "output_digits"}];
  missing = Select[{"point", "relative_tolerance", "operator_identity",
      "source_identity", "checkpoint_identity", "analytic_metadata",
      "owner_signature_identity", "physical_payload_identity",
      "provenance_identity"},
    !KeyExistsQ[certificateRequest, #] &];
  If[reserved =!= {} || missing =!= {},
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent residual request is missing required fields or contains reserved protocol keys",
      "Missing" -> missing, "Reserved" -> reserved|>], Module]];
  payload = Join[<|"scope" -> "stored_truncation", "include_residual" -> False,
      "options" -> <|"tail_estimate" -> False|>|>,
    certificateRequest];
  request = Join[payload, <|"schema" -> 2,
    "op" -> "local.certify_residual", "session" -> tokens["Session"],
    "local" -> tokens["Local"]|>];
  If[IntegerQ[outputDigits],
    request = Append[request, "output_digits" -> outputDigits]];
  RunRequest[request]];

RunPersistentLocalMatch[basis_List, incoming_Association,
    matchRequest_Association] := Module[
  {basisTokens, incomingTokens = persistentLocalHandles[incoming], bad,
   sessions, reserved, required, missing, request},
  If[basis === {} || !AllTrue[basis, AssociationQ],
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent local matching requires a nonempty list of local handles"|>],
      Module]];
  basisTokens = persistentLocalHandles /@ basis;
  bad = Select[Range[Length[basisTokens]],
    FailureQ[basisTokens[[#]]] &];
  If[bad =!= {},
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent local match basis contains malformed handles",
      "Positions" -> bad|>], Module]];
  If[FailureQ[incomingTokens], Return[incomingTokens, Module]];
  sessions = DeleteDuplicates[Join[Lookup[basisTokens, "Session"],
    {incomingTokens["Session"]}]];
  If[Length[sessions] =!= 1,
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent local matching requires every local in one native session",
      "Sessions" -> sessions|>], Module]];
  reserved = Intersection[Keys[matchRequest],
    {"schema", "op", "session", "basis", "incoming"}];
  required = {"basis_chart", "incoming_chart", "basis_point",
    "incoming_point", "epsilon", "basis_checkpoint_identities",
    "incoming_checkpoint_identity", "checkpoint_identity"};
  missing = Select[required, !KeyExistsQ[matchRequest, #] &];
  If[reserved =!= {} || missing =!= {},
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent local match request is missing required fields or contains reserved protocol keys",
      "Missing" -> missing, "Reserved" -> reserved|>], Module]];
  request = Join[matchRequest, <|"schema" -> 2, "op" -> "local.match",
    "session" -> First[sessions],
    "basis" -> Lookup[basisTokens, "Local"],
    "incoming" -> incomingTokens["Local"]|>];
  RunRequest[request]];

RunPersistentAcbLocalMatch[basis_List, incoming_Association,
    matchRequest_Association] := Module[
  {basisTokens, incomingTokens = persistentLocalHandles[incoming], bad,
   sessions, reserved, required, missing, request},
  If[basis === {} || !AllTrue[basis, AssociationQ],
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent Acb local matching requires a nonempty list of local handles"|>],
      Module]];
  basisTokens = persistentLocalHandles /@ basis;
  bad = Select[Range[Length[basisTokens]],
    FailureQ[basisTokens[[#]]] &];
  If[bad =!= {},
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent Acb local match basis contains malformed handles",
      "Positions" -> bad|>], Module]];
  If[FailureQ[incomingTokens], Return[incomingTokens, Module]];
  sessions = DeleteDuplicates[Join[Lookup[basisTokens, "Session"],
    {incomingTokens["Session"]}]];
  If[Length[sessions] =!= 1,
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent Acb local matching requires every local in one native session",
      "Sessions" -> sessions|>], Module]];
  reserved = Intersection[Keys[matchRequest],
    {"schema", "op", "session", "basis", "incoming"}];
  required = {"basis_chart", "incoming_chart", "basis_point",
    "incoming_point", "epsilon", "basis_checkpoint_identities",
    "incoming_checkpoint_identity", "checkpoint_identity",
    "exact_lattice", "refinement"};
  missing = Select[required, !KeyExistsQ[matchRequest, #] &];
  If[reserved =!= {} || missing =!= {},
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent Acb local match request is missing required fields or contains reserved protocol keys",
      "Missing" -> missing, "Reserved" -> reserved|>], Module]];
  request = Join[matchRequest, <|"schema" -> 2,
    "op" -> "local.match_acb", "session" -> First[sessions],
    "basis" -> Lookup[basisTokens, "Local"],
    "incoming" -> incomingTokens["Local"]|>];
  RunRequest[request]];

persistentMatchHandles[handle_Association] := Module[{session, match},
  session = Lookup[handle, "session", Lookup[handle, "Session", None]];
  match = Lookup[handle, "match", Lookup[handle, "Match", None]];
  If[!StringQ[session] || !StringQ[match],
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent match handle requires exact session and match tokens"|>],
      Module]];
  <|"Session" -> session, "Match" -> match|>];

PersistentLocalMatchStatistics[handle_Association] := Module[
  {tokens = persistentMatchHandles[handle]},
  If[FailureQ[tokens], Return[tokens, Module]];
  RunRequest[<|"schema" -> 2, "op" -> "match.stats",
    "session" -> tokens["Session"], "match" -> tokens["Match"]|>]];

ReleasePersistentLocalMatch[handle_Association] := Module[
  {tokens = persistentMatchHandles[handle]},
  If[FailureQ[tokens], Return[tokens, Module]];
  RunRequest[<|"schema" -> 2, "op" -> "match.release",
    "session" -> tokens["Session"], "match" -> tokens["Match"]|>]];

RunPersistentEndpointLimit[handle_Association,
    endpointRequest_Association] := Module[
  {tokens = persistentLocalHandles[handle], reserved, required, missing,
   request},
  If[FailureQ[tokens], Return[tokens, Module]];
  reserved = Intersection[Keys[endpointRequest],
    {"schema", "op", "session", "local", "output_digits",
      "include_coefficients"}];
  required = {"source_checkpoint_identity", "checkpoint_identity",
    "approach_direction", "cancellation"};
  missing = Select[required, !KeyExistsQ[endpointRequest, #] &];
  If[reserved =!= {} || missing =!= {},
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent endpoint request is missing required fields or contains reserved protocol keys",
      "Missing" -> missing, "Reserved" -> reserved|>], Module]];
  request = Join[endpointRequest, <|"schema" -> 2,
    "op" -> "local.endpoint_limit", "session" -> tokens["Session"],
    "local" -> tokens["Local"]|>];
  RunRequest[request]];

persistentEndpointHandles[handle_Association] := Module[{session, endpoint},
  session = Lookup[handle, "session", Lookup[handle, "Session", None]];
  endpoint = Lookup[handle, "endpoint", Lookup[handle, "Endpoint", None]];
  If[!StringQ[session] || !StringQ[endpoint],
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent endpoint handle requires exact session and endpoint tokens"|>],
      Module]];
  <|"Session" -> session, "Endpoint" -> endpoint|>];

PersistentEndpointStatistics[handle_Association] := Module[
  {tokens = persistentEndpointHandles[handle]},
  If[FailureQ[tokens], Return[tokens, Module]];
  RunRequest[<|"schema" -> 2, "op" -> "endpoint.stats",
    "session" -> tokens["Session"],
    "endpoint" -> tokens["Endpoint"]|>]];

ExportPersistentEndpoint[handle_Association, checkpointIdentity_String,
    outputDigits_:Automatic] := Module[
  {tokens = persistentEndpointHandles[handle], request},
  If[FailureQ[tokens], Return[tokens, Module]];
  If[StringLength[checkpointIdentity] == 0,
    Return[Failure["CppBackend", <|"Detail" ->
      "endpoint export checkpoint identity must be nonempty"|>], Module]];
  If[outputDigits =!= Automatic &&
      (!IntegerQ[outputDigits] || outputDigits < 1),
    Return[Failure["CppBackend", <|"Detail" ->
      "endpoint export digits must be a positive integer"|>], Module]];
  request = <|"schema" -> 2, "op" -> "endpoint.export",
    "session" -> tokens["Session"], "endpoint" -> tokens["Endpoint"],
    "checkpoint_identity" -> checkpointIdentity|>;
  If[IntegerQ[outputDigits],
    request = Append[request, "output_digits" -> outputDigits]];
  RunRequest[request]];

ReleasePersistentEndpoint[handle_Association] := Module[
  {tokens = persistentEndpointHandles[handle]},
  If[FailureQ[tokens], Return[tokens, Module]];
  RunRequest[<|"schema" -> 2, "op" -> "endpoint.release",
    "session" -> tokens["Session"],
    "endpoint" -> tokens["Endpoint"]|>]];

persistentTilePlanHandles[handle_Association] := Module[
  {session, plan, checkpoint},
  session = Lookup[handle, "session", Lookup[handle, "Session", None]];
  plan = Lookup[handle, "tile_plan", Lookup[handle, "TilePlan", None]];
  checkpoint = Lookup[handle, "checkpoint_identity",
    Lookup[handle, "CheckpointIdentity", None]];
  If[!StringQ[session] || !StringQ[plan] || !StringQ[checkpoint] ||
      StringLength[checkpoint] == 0,
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent tile-plan handle requires exact session, tile-plan, and checkpoint tokens"|>],
      Module]];
  <|"Session" -> session, "TilePlan" -> plan,
    "CheckpointIdentity" -> checkpoint|>];

persistentTransportArmHandles[handle_Association] := Module[
  {session, state, checkpoint},
  session = Lookup[handle, "session", Lookup[handle, "Session", None]];
  state = Lookup[handle, "transport_state",
    Lookup[handle, "TransportState", None]];
  checkpoint = Lookup[handle, "checkpoint_identity",
    Lookup[handle, "CheckpointIdentity", None]];
  If[!StringQ[session] || !StringQ[state] || !StringQ[checkpoint] ||
      StringLength[checkpoint] == 0,
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent transport-arm handle requires exact session, state, and checkpoint tokens"|>],
      Module]];
  <|"Session" -> session, "TransportState" -> state,
    "CheckpointIdentity" -> checkpoint|>];

RunPersistentPlannedEndpointLimit[plan_Association, arm_String,
    local_Association, policy_Association] := Module[
  {planTokens = persistentTilePlanHandles[plan],
   localTokens = persistentLocalHandles[local], sourceCheckpoint},
  If[FailureQ[planTokens], Return[planTokens, Module]];
  If[FailureQ[localTokens], Return[localTokens, Module]];
  sourceCheckpoint = Lookup[local, "checkpoint_identity",
    Lookup[local, "CheckpointIdentity", None]];
  If[planTokens["Session"] =!= localTokens["Session"] ||
      !MemberQ[{"lower", "upper"}, arm] ||
      !StringQ[sourceCheckpoint] || StringLength[sourceCheckpoint] == 0 ||
      Sort[Keys[policy]] =!= Sort[{"checkpoint_identity", "cancellation"}] ||
      !StringQ[Lookup[policy, "checkpoint_identity", None]] ||
      StringLength[Lookup[policy, "checkpoint_identity", ""]] == 0 ||
      !AssociationQ[Lookup[policy, "cancellation", None]],
    Return[Failure["CppBackend", <|"Detail" ->
      "plan-bound endpoint evaluation requires one session, arm lower/upper, a retained local checkpoint, and exactly checkpoint_identity plus cancellation policy fields"|>], Module]];
  RunRequest[<|"schema" -> 2, "op" -> "tile.endpoint_limit",
    "session" -> planTokens["Session"],
    "tile_plan" -> planTokens["TilePlan"],
    "tile_plan_checkpoint_identity" -> planTokens["CheckpointIdentity"],
    "arm" -> arm, "local" -> localTokens["Local"],
    "source_checkpoint_identity" -> sourceCheckpoint,
    "checkpoint_identity" -> policy["checkpoint_identity"],
    "cancellation" -> policy["cancellation"]|>]];

RunPersistentWeightedPlannedEndpointLimit[plan_Association, arm_String,
    local_Association, row_Association, checkpointIdentity_String,
    cancellation_Association] := Module[
  {planTokens = persistentTilePlanHandles[plan],
   localTokens = persistentLocalHandles[local], rowKeys,
   rowCheckpoint, endpointCheckpoint, projected = None,
   endpoint = None, cleanup, released},
  If[FailureQ[planTokens], Return[planTokens, Module]];
  If[FailureQ[localTokens], Return[localTokens, Module]];
  rowKeys = {"schema", "columns", "exact_identity", "entries"};
  If[planTokens["Session"] =!= localTokens["Session"] ||
      !MemberQ[{"lower", "upper"}, arm] ||
      StringLength[checkpointIdentity] == 0 ||
      Sort[Keys[row]] =!= Sort[rowKeys] ||
      Sort[Keys[cancellation]] =!= {"mode"} ||
      !MemberQ[{"exact-coefficient-field", "exact-or-acb-singleton"},
        Lookup[cancellation, "mode", None]],
    Return[Failure["CppBackend", <|"Detail" ->
      "weighted plan-bound endpoint evaluation requires one session, arm lower/upper, a nonempty checkpoint root, an exact prepared rational row, and one supported exact cancellation mode"|>], Module]];
  rowCheckpoint = checkpointIdentity <> ":weighted-row";
  endpointCheckpoint = checkpointIdentity <> ":weighted-endpoint";
  cleanup[] := Module[{result},
    If[!AssociationQ[projected], Return[Null, Module]];
    result = Quiet[ReleasePersistentLocal[projected]];
    projected = None;
    result];
  projected = ApplyPersistentRationalRow[
    local, row, rowCheckpoint];
  If[!persistentCommandOKQ[projected], Return[projected, Module]];
  If[Lookup[projected, "json_coefficients", None] =!= 0 ||
      FailureQ[persistentLocalHandles[projected]],
    cleanup[];
    Return[Failure["CppBackend", <|"Detail" ->
      "weighted endpoint projection violated its opaque retained-local contract"|>], Module]];
  endpoint = RunPersistentPlannedEndpointLimit[
    plan, arm, projected, <|
      "checkpoint_identity" -> endpointCheckpoint,
      "cancellation" -> cancellation|>];
  If[!persistentCommandOKQ[endpoint],
    released = cleanup[];
    If[AssociationQ[released] && !persistentCommandOKQ[released],
      Return[Failure["CppBackend", <|"Detail" ->
        "weighted endpoint failure also failed to release its projected scalar local",
        "BackendFailure" -> endpoint,
        "CleanupFailure" -> released|>], Module]];
    Return[endpoint, Module]];
  If[Lookup[endpoint, "json_coefficients", None] =!= 0 ||
      FailureQ[persistentEndpointHandles[endpoint]],
    Quiet[ReleasePersistentEndpoint[endpoint]];
    cleanup[];
    Return[Failure["CppBackend", <|"Detail" ->
      "weighted endpoint gate violated its opaque retained-result contract"|>], Module]];
  released = cleanup[];
  If[!AssociationQ[released] || !persistentCommandOKQ[released],
    Quiet[ReleasePersistentEndpoint[endpoint]];
    Return[Failure["CppBackend", <|"Detail" ->
      "weighted endpoint was created but its projected scalar public token could not be released",
      "CleanupFailure" -> released|>], Module]];
  Join[endpoint, <|"weighted_composition" -> <|
    "capability" ->
      "retained-native-weighted-plan-bound-endpoint-v1",
    "order" -> "rational-row-before-endpoint-gate",
    "coefficient_transport" -> "native-retained-only",
    "row_checkpoint_identity" -> rowCheckpoint,
    "endpoint_checkpoint_identity" -> endpointCheckpoint,
    "projected_local_public_token_released" -> True|>|>]];

CreatePersistentTilePlan[owner_, lower_Association, upper_Association,
    checkpointIdentity_String, divisionOrder_:3] := Module[
  {session = persistentCheckpointSession[owner], required,
   malformed, response},
  If[FailureQ[session], Return[session, Module]];
  If[StringLength[checkpointIdentity] == 0 || !IntegerQ[divisionOrder] ||
      divisionOrder < 2,
    Return[Failure["CppBackend", <|"Detail" ->
      "native tile planning requires a nonempty checkpoint identity and integer DivisionOrder >= 2"|>],
      Module]];
  required = Sort[{"from_exact", "to_exact", "charts", "topology"}];
  malformed = Select[{lower, upper}, Sort[Keys[#]] =!= required &];
  If[malformed =!= {},
    Return[Failure["CppBackend", <|"Detail" ->
      "each native tile arm must contain exactly from_exact, to_exact, charts, and topology"|>],
      Module]];
  response = RunRequest[<|"schema" -> 2, "op" -> "tile.plan",
    "session" -> session, "checkpoint_identity" -> checkpointIdentity,
    "division_order" -> divisionOrder, "lower" -> lower,
    "upper" -> upper|>];
  Which[
    FailureQ[response], response,
    persistentCommandOKQ[response], response,
    True, Failure["CppBackend", <|
      "Operation" -> "tile.plan", "BackendResponse" -> response,
      "Detail" -> "persistent native tile planning returned a non-ok backend status"|>]]];

CreatePersistentArmTilePlan[owner_, arm_Association,
    checkpointIdentity_String, divisionOrder_:3] := Module[
  {session = persistentCheckpointSession[owner], required},
  If[FailureQ[session], Return[session, Module]];
  If[StringLength[checkpointIdentity] == 0 || !IntegerQ[divisionOrder] ||
      divisionOrder < 2,
    Return[Failure["CppBackend", <|"Detail" ->
      "native single-arm tile planning requires a nonempty checkpoint identity and integer DivisionOrder >= 2"|>],
      Module]];
  required = Sort[{"from_exact", "to_exact", "charts", "topology"}];
  If[Sort[Keys[arm]] =!= required,
    Return[Failure["CppBackend", <|"Detail" ->
      "the native tile arm must contain exactly from_exact, to_exact, charts, and topology"|>],
      Module]];
  RunRequest[<|"schema" -> 2, "op" -> "tile.plan_arm",
    "session" -> session, "checkpoint_identity" -> checkpointIdentity,
    "division_order" -> divisionOrder, "arm" -> arm|>]];

PersistentTilePlanStatistics[handle_Association] := Module[
  {tokens = persistentTilePlanHandles[handle]},
  If[FailureQ[tokens], Return[tokens, Module]];
  RunRequest[<|"schema" -> 2, "op" -> "tile.stats",
    "session" -> tokens["Session"],
    "tile_plan" -> tokens["TilePlan"]|>]];

PersistentTileMatchInterval[handle_Association, arm_String,
    index_Integer] := Module[{tokens = persistentTilePlanHandles[handle]},
  If[FailureQ[tokens], Return[tokens, Module]];
  If[!MemberQ[{"lower", "upper"}, arm] || index < 1,
    Return[Failure["CppBackend", <|"Detail" ->
      "native tile match lookup requires arm lower/upper and a positive one-based index"|>],
      Module]];
  RunRequest[<|"schema" -> 2, "op" -> "tile.match_interval",
    "session" -> tokens["Session"],
    "tile_plan" -> tokens["TilePlan"], "arm" -> arm,
    "match" -> index - 1|>]];

PersistentTileIntegrationInterval[handle_Association, arm_String,
    index_Integer] := Module[{tokens = persistentTilePlanHandles[handle]},
  If[FailureQ[tokens], Return[tokens, Module]];
  If[!MemberQ[{"lower", "upper"}, arm] || index < 1,
    Return[Failure["CppBackend", <|"Detail" ->
      "native tile interval lookup requires arm lower/upper and a positive one-based index"|>],
      Module]];
  RunRequest[<|"schema" -> 2, "op" -> "tile.integration_interval",
    "session" -> tokens["Session"],
    "tile_plan" -> tokens["TilePlan"], "arm" -> arm,
    "tile" -> index - 1|>]];

RunPersistentPlannedMatch[plan_Association, arm_String,
    index_Integer, basis_List, incoming_Association,
    policy_Association] := Module[
  {planTokens = persistentTilePlanHandles[plan], basisTokens,
   incomingTokens = persistentLocalHandles[incoming], sessions, bad,
   reserved, required, missing, payload},
  If[FailureQ[planTokens], Return[planTokens, Module]];
  If[!MemberQ[{"lower", "upper"}, arm] || index < 1 ||
      basis === {} || !AllTrue[basis, AssociationQ],
    Return[Failure["CppBackend", <|"Detail" ->
      "planned matching requires lower/upper, a one-based positive match index, and a nonempty basis"|>], Module]];
  basisTokens = persistentLocalHandles /@ basis;
  bad = Select[Range[Length[basisTokens]], FailureQ[basisTokens[[#]]] &];
  If[bad =!= {},
    Return[Failure["CppBackend", <|"Detail" ->
      "planned match basis contains malformed retained locals",
      "Positions" -> bad|>], Module]];
  If[FailureQ[incomingTokens], Return[incomingTokens, Module]];
  sessions = DeleteDuplicates[Join[{planTokens["Session"],
      incomingTokens["Session"]}, Lookup[basisTokens, "Session"]]];
  If[Length[sessions] =!= 1,
    Return[Failure["CppBackend", <|"Detail" ->
      "planned matching requires the plan, basis, and incoming local in one native session",
      "Sessions" -> sessions|>], Module]];
  reserved = Intersection[Keys[policy],
    {"schema", "op", "session", "tile_plan", "arm", "match",
     "basis", "incoming"}];
  required = {"epsilon", "checkpoint_identity"};
  missing = Select[required, !KeyExistsQ[policy, #] &];
  If[reserved =!= {} || missing =!= {},
    Return[Failure["CppBackend", <|"Detail" ->
      "planned match policy is missing required fields or contains reserved protocol keys",
      "Missing" -> missing, "Reserved" -> reserved|>], Module]];
  payload = Join[policy, <|"schema" -> 2,
    "op" -> "tile.match_advance", "session" -> First[sessions],
    "tile_plan" -> planTokens["TilePlan"], "arm" -> arm,
    "match" -> index - 1, "basis" -> Lookup[basisTokens, "Local"],
    "incoming" -> incomingTokens["Local"]|>];
  RunRequest[payload]];

MaterializePersistentLocalMatch[match_Association,
    checkpointIdentity_String] := Module[
  {tokens = persistentMatchHandles[match]},
  If[FailureQ[tokens], Return[tokens, Module]];
  If[StringLength[checkpointIdentity] == 0,
    Return[Failure["CppBackend", <|"Detail" ->
      "match materialization checkpoint identity must be nonempty"|>],
      Module]];
  RunRequest[<|"schema" -> 2, "op" -> "match.materialize_local",
    "session" -> tokens["Session"], "match" -> tokens["Match"],
    "checkpoint_identity" -> checkpointIdentity|>]];

ApplyPersistentRationalRow[local_Association, row_Association,
    checkpointIdentity_String] := Module[
  {tokens = persistentLocalHandles[local], sourceCheckpoint, rowKeys},
  If[FailureQ[tokens], Return[tokens, Module]];
  sourceCheckpoint = Lookup[local, "checkpoint_identity",
    Lookup[local, "CheckpointIdentity", None]];
  rowKeys = {"schema", "columns", "exact_identity", "entries"};
  If[!StringQ[sourceCheckpoint] || StringLength[sourceCheckpoint] == 0 ||
      StringLength[checkpointIdentity] == 0 ||
      Sort[Keys[row]] =!= Sort[rowKeys],
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent rational-row application requires nonempty source/result checkpoint identities and exactly schema, columns, exact_identity, entries in its prepared row"|>],
      Module]];
  RunRequest[<|"schema" -> 2, "op" -> "local.apply_rational_row",
    "session" -> tokens["Session"], "local" -> tokens["Local"],
    "row" -> row, "source_checkpoint_identity" -> sourceCheckpoint,
    "checkpoint_identity" -> checkpointIdentity|>]];

ReleasePersistentTilePlan[handle_Association] := Module[
  {tokens = persistentTilePlanHandles[handle]},
  If[FailureQ[tokens], Return[tokens, Module]];
  RunRequest[<|"schema" -> 2, "op" -> "tile.release",
    "session" -> tokens["Session"],
    "tile_plan" -> tokens["TilePlan"]|>]];

RunPersistentTileIntegral[plan_Association, arm_String, tile_Integer,
    local_Association, epsilon_Association,
    checkpointIdentity_String, certifyTail_:False] := Module[
  {planTokens = persistentTilePlanHandles[plan],
   localTokens = persistentLocalHandles[local], sourceCheckpoint, request},
  If[FailureQ[planTokens], Return[planTokens, Module]];
  If[FailureQ[localTokens], Return[localTokens, Module]];
  sourceCheckpoint = Lookup[local, "checkpoint_identity",
    Lookup[local, "CheckpointIdentity", None]];
  If[planTokens["Session"] =!= localTokens["Session"] ||
      !MemberQ[{"lower", "upper"}, arm] || tile < 1 ||
      Sort[Keys[epsilon]] =!= Sort[{"min", "max"}] ||
      !StringQ[sourceCheckpoint] || StringLength[sourceCheckpoint] == 0 ||
      StringLength[checkpointIdentity] == 0 || !BooleanQ[certifyTail],
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent tile integration requires one session, arm lower/upper, a positive one-based tile, exact min/max epsilon window, nonempty source/result checkpoint identities, and a Boolean certifyTail request"|>],
      Module]];
  request = <|"schema" -> 2, "op" -> "integration.line",
    "session" -> planTokens["Session"],
    "tile_plan" -> planTokens["TilePlan"],
    "tile_plan_checkpoint_identity" -> planTokens["CheckpointIdentity"],
    "local" -> localTokens["Local"], "arm" -> arm,
    "tile" -> tile - 1, "epsilon" -> epsilon,
    "source_checkpoint_identity" -> sourceCheckpoint,
    "checkpoint_identity" -> checkpointIdentity|>;
  If[TrueQ[certifyTail], request = Append[request, "certify_tail" -> True]];
  RunRequest[request]];

RunPersistentNativeArms[plan_Association, anchor_Association,
    arms_Association, epsilon_Association, checkpointRoot_String,
    refinement_Association, certifyTail_:False] := Module[
  {planTokens = persistentTilePlanHandles[plan],
   anchorTokens = persistentLocalHandles[anchor], normalizeArm,
   normalized, sessions, request},
  If[FailureQ[planTokens], Return[planTokens, Module]];
  If[FailureQ[anchorTokens], Return[anchorTokens, Module]];
  If[Sort[Keys[arms]] =!= Sort[{"lower", "upper"}] ||
      Sort[Keys[epsilon]] =!=
        Sort[{"min", "max", "required_complete_max",
          "match_required_complete_max"}] ||
      Sort[Keys[refinement]] =!=
        Sort[{"relative_tolerance", "max_steps"}] ||
      StringLength[checkpointRoot] == 0 || !BooleanQ[certifyTail],
    Return[Failure["CppBackend", <|"Detail" ->
      "native arm marching requires exact lower/upper arm records, epsilon and bounded-refinement contracts, a nonempty checkpoint root, and a Boolean certifyTail request"|>], Module]];
  normalizeArm[raw_Association] := Module[
    {sets, tokens, bad},
    If[Sort[Keys[raw]] =!= Sort[{"receiving_basis", "integrand_rows"}] ||
        !ListQ[raw["receiving_basis"]] ||
        !ListQ[raw["integrand_rows"]] ||
        Length[raw["integrand_rows"]] =!=
          Length[raw["receiving_basis"]] + 1 ||
        !AllTrue[raw["integrand_rows"], AssociationQ],
      Return[Failure["CppBackend", <|"Detail" ->
        "each native arm requires one receiving basis per match and one prepared integrand row per tile"|>], Module]];
    sets = raw["receiving_basis"];
    If[!AllTrue[sets, ListQ[#] && # =!= {} &&
          AllTrue[#, AssociationQ] &],
      Return[Failure["CppBackend", <|"Detail" ->
        "every native receiving basis must be a nonempty list of retained-local associations"|>], Module]];
    tokens = Map[persistentLocalHandles, sets, {2}];
    bad = Cases[tokens, _Failure, Infinity];
    If[bad =!= {}, Return[First[bad], Module]];
    <|"receiving_basis" -> Map[Lookup[#, "Local"] &, tokens],
      "integrand_rows" -> raw["integrand_rows"],
      "sessions" -> Flatten[Map[Lookup[#, "Session"] &, tokens]]|>];
  normalized = Map[normalizeArm, arms];
  If[AnyTrue[Values[normalized], FailureQ],
    Return[First[Select[Values[normalized], FailureQ]], Module]];
  sessions = DeleteDuplicates[Join[
    {planTokens["Session"], anchorTokens["Session"]},
    normalized["lower", "sessions"],
    normalized["upper", "sessions"]]];
  If[Length[sessions] =!= 1,
    Return[Failure["CppBackend", <|"Detail" ->
      "the retained plan, anchor, and every receiving basis local must belong to one persistent session",
      "Sessions" -> sessions|>], Module]];
  request = <|"schema" -> 2, "op" -> "integration.run_arms",
    "session" -> First[sessions],
    "tile_plan" -> planTokens["TilePlan"],
    "tile_plan_checkpoint_identity" ->
      planTokens["CheckpointIdentity"],
    "anchor" -> anchorTokens["Local"],
    "anchor_checkpoint_identity" ->
      Lookup[anchor, "checkpoint_identity",
        Lookup[anchor, "CheckpointIdentity", ""]],
    "epsilon" -> epsilon,
    "refinement" -> refinement,
    "checkpoint_policy" -> <|
      "schema" -> "diffexp2-deterministic-arm-checkpoints-v1",
      "root" -> checkpointRoot|>,
    "lower" -> KeyDrop[normalized["lower"], "sessions"],
    "upper" -> KeyDrop[normalized["upper"], "sessions"]|>;
  If[TrueQ[certifyTail], request = Append[request, "certify_tail" -> True]];
  RunRequest[request]];

runPersistentTransportArmsOperation[operation_String,
    plan_Association, anchor_Association,
    arms_Association, epsilon_Association, checkpointRoot_String,
    refinement_Association] := Module[
  {planTokens = persistentTilePlanHandles[plan],
   anchorTokens = persistentLocalHandles[anchor], normalizeArm,
   normalized, bad, sessions, anchorCheckpoint, epsilonKeys,
   refinementKeys},
  If[FailureQ[planTokens], Return[planTokens, Module]];
  If[FailureQ[anchorTokens], Return[anchorTokens, Module]];
  epsilonKeys = {"min", "max", "required_complete_max",
    "match_required_complete_max"};
  refinementKeys = {"relative_tolerance", "max_steps"};
  anchorCheckpoint = Lookup[anchor, "checkpoint_identity",
    Lookup[anchor, "CheckpointIdentity", None]];
  If[Sort[Keys[arms]] =!= Sort[{"lower", "upper"}] ||
      Sort[Keys[epsilon]] =!= Sort[epsilonKeys] ||
      !AllTrue[Lookup[epsilon, epsilonKeys], IntegerQ] ||
      !TrueQ[epsilon["min"] <= epsilon["required_complete_max"] <=
        epsilon["match_required_complete_max"] <= epsilon["max"]] ||
      Sort[Keys[refinement]] =!= Sort[refinementKeys] ||
      !StringQ[Lookup[refinement, "relative_tolerance", None]] ||
      StringLength[Lookup[refinement, "relative_tolerance", ""]] == 0 ||
      !IntegerQ[Lookup[refinement, "max_steps", None]] ||
      !TrueQ[0 <= refinement["max_steps"] <= 32] ||
      !StringQ[anchorCheckpoint] || StringLength[anchorCheckpoint] == 0 ||
      StringLength[StringTrim[checkpointRoot]] == 0,
    Return[Failure["CppBackend", <|"Detail" ->
      "two-arm retained transport marching requires exact lower/upper arm records, ordered integer epsilon bounds, bounded refinement, and nonempty source/run checkpoint identities"|>], Module]];
  normalizeArm[raw_] := Module[{sets, tokens},
    If[!AssociationQ[raw] ||
        Sort[Keys[raw]] =!= {"receiving_basis"} ||
        !ListQ[raw["receiving_basis"]],
      Return[Failure["CppBackend", <|"Detail" ->
        "each two-arm transport record must contain exactly receiving_basis"|>], Module]];
    sets = raw["receiving_basis"];
    If[!AllTrue[sets,
        ListQ[#] && # =!= {} && AllTrue[#, AssociationQ] &],
      Return[Failure["CppBackend", <|"Detail" ->
        "every two-arm transport receiving basis must be a nonempty list of retained-local associations"|>], Module]];
    tokens = Map[persistentLocalHandles, sets, {2}];
    If[Cases[tokens, _Failure, Infinity] =!= {},
      Return[First[Cases[tokens, _Failure, Infinity]], Module]];
    <|"receiving_basis" -> Map[Lookup[#, "Local"] &, tokens, {2}],
      "sessions" -> Flatten[Map[Lookup[#, "Session"] &, tokens, {2}]]|>];
  normalized = Map[normalizeArm, arms];
  bad = Select[Values[normalized], FailureQ];
  If[bad =!= {}, Return[First[bad], Module]];
  sessions = DeleteDuplicates[Join[
    {planTokens["Session"], anchorTokens["Session"]},
    normalized["lower", "sessions"],
    normalized["upper", "sessions"]]];
  If[Length[sessions] =!= 1,
    Return[Failure["CppBackend", <|"Detail" ->
      "the retained plan, anchor, and every two-arm receiving-basis local must belong to one persistent session",
      "Sessions" -> sessions|>], Module]];
  RunRequest[<|"schema" -> 2, "op" -> operation,
    "session" -> First[sessions],
    "tile_plan" -> planTokens["TilePlan"],
    "tile_plan_checkpoint_identity" -> planTokens["CheckpointIdentity"],
    "anchor" -> anchorTokens["Local"],
    "anchor_checkpoint_identity" -> anchorCheckpoint,
    "epsilon" -> epsilon, "refinement" -> refinement,
    "checkpoint_policy" -> <|
      "schema" -> "diffexp2-deterministic-arm-checkpoints-v1",
      "root" -> checkpointRoot|>,
    "lower" -> KeyDrop[normalized["lower"], "sessions"],
    "upper" -> KeyDrop[normalized["upper"], "sessions"]|>]];

RunPersistentTransportArms[plan_Association, anchor_Association,
    arms_Association, epsilon_Association, checkpointRoot_String,
    refinement_Association] := runPersistentTransportArmsOperation[
  "transport.run_arms", plan, anchor, arms, epsilon, checkpointRoot,
  refinement];

RunPersistentConsumingTransportArms[plan_Association,
    anchor_Association, arms_Association, epsilon_Association,
    checkpointRoot_String, refinement_Association] :=
  runPersistentTransportArmsOperation[
    "transport.run_arms_consuming", plan, anchor, arms, epsilon,
    checkpointRoot, refinement];

ConsumePersistentTransportHop[plan_Association, arm_String,
    index_Integer, basis_List, incoming_Association,
    epsilon_Association, checkpointRoot_String,
    refinement_Association] := Module[
  {planTokens = persistentTilePlanHandles[plan], basisTokens,
   incomingTokens = persistentLocalHandles[incoming], sessions,
   incomingCheckpoint, epsilonKeys, refinementKeys},
  If[FailureQ[planTokens], Return[planTokens, Module]];
  If[!MemberQ[{"lower", "upper"}, arm] || index < 1 || basis === {} ||
      !AllTrue[basis, AssociationQ],
    Return[Failure["CppBackend", <|"Detail" ->
      "consuming transport hop requires lower/upper, a positive one-based match index, and one nonempty retained basis"|>], Module]];
  basisTokens = persistentLocalHandles /@ basis;
  If[Cases[basisTokens, _Failure, Infinity] =!= {},
    Return[First[Cases[basisTokens, _Failure, Infinity]], Module]];
  If[FailureQ[incomingTokens], Return[incomingTokens, Module]];
  incomingCheckpoint = Lookup[incoming, "checkpoint_identity",
    Lookup[incoming, "CheckpointIdentity", None]];
  epsilonKeys = {"min", "max", "required_complete_max"};
  refinementKeys = {"relative_tolerance", "max_steps"};
  sessions = DeleteDuplicates[Join[{planTokens["Session"],
      incomingTokens["Session"]}, Lookup[basisTokens, "Session"]]];
  If[Length[sessions] =!= 1 ||
      Sort[Keys[epsilon]] =!= Sort[epsilonKeys] ||
      !AllTrue[Lookup[epsilon, epsilonKeys], IntegerQ] ||
      !TrueQ[epsilon["min"] <= epsilon["required_complete_max"] <=
        epsilon["max"]] ||
      Sort[Keys[refinement]] =!= Sort[refinementKeys] ||
      !StringQ[Lookup[refinement, "relative_tolerance", None]] ||
      !IntegerQ[Lookup[refinement, "max_steps", None]] ||
      !TrueQ[0 <= refinement["max_steps"] <= 32] ||
      !StringQ[incomingCheckpoint] || StringLength[checkpointRoot] == 0,
    Return[Failure["CppBackend", <|"Detail" ->
      "consuming transport hop received inconsistent sessions, epsilon/refinement bounds, or checkpoint identities"|>], Module]];
  RunRequest[<|"schema" -> 2, "op" -> "transport.consume_hop",
    "session" -> First[sessions],
    "tile_plan" -> planTokens["TilePlan"],
    "tile_plan_checkpoint_identity" -> planTokens["CheckpointIdentity"],
    "arm" -> arm, "match" -> index - 1,
    "receiving_basis" -> Lookup[basisTokens, "Local"],
    "incoming" -> incomingTokens["Local"],
    "incoming_checkpoint_identity" -> incomingCheckpoint,
    "epsilon" -> epsilon, "refinement" -> refinement,
    "checkpoint_policy" -> <|
      "schema" -> "diffexp2-deterministic-arm-checkpoints-v1",
      "root" -> checkpointRoot|>|>]];

PublishPersistentConsumedTransportStates[plan_Association,
    anchor_Association, arms_Association, epsilon_Association,
    checkpointRoot_String, refinement_Association] := Module[
  {planTokens = persistentTilePlanHandles[plan],
   anchorTokens = persistentLocalHandles[anchor], anchorCheckpoint,
   normalizeArm, normalized, failures, sessions, epsilonKeys,
   refinementKeys},
  If[FailureQ[planTokens], Return[planTokens, Module]];
  If[FailureQ[anchorTokens], Return[anchorTokens, Module]];
  anchorCheckpoint = Lookup[anchor, "checkpoint_identity",
    Lookup[anchor, "CheckpointIdentity", None]];
  normalizeArm[record_] := Module[{keys, tileTokens},
    keys = {"tile_sources"};
    If[!AssociationQ[record] || Sort[Keys[record]] =!= Sort[keys] ||
        !ListQ[record["tile_sources"]] || record["tile_sources"] === {} ||
        !AllTrue[record["tile_sources"], AssociationQ],
      Return[Failure["CppBackend", <|"Detail" ->
        "each consumed-state arm requires exactly ordered tile_sources"|>], Module]];
    tileTokens = persistentLocalHandles /@ record["tile_sources"];
    If[Cases[tileTokens, _Failure, Infinity] =!= {},
      Return[First[Cases[tileTokens, _Failure, Infinity]], Module]];
    <|"tile_sources" -> Lookup[tileTokens, "Local"],
      "sessions" -> Lookup[tileTokens, "Session"]|>];
  If[Sort[Keys[arms]] =!= {"lower", "upper"},
    Return[Failure["CppBackend", <|"Detail" ->
      "consumed-state publication requires exactly lower and upper arms"|>], Module]];
  normalized = Map[normalizeArm, arms];
  failures = Cases[normalized, _Failure, Infinity];
  If[failures =!= {}, Return[First[failures], Module]];
  sessions = DeleteDuplicates[Join[{planTokens["Session"],
      anchorTokens["Session"]}, Flatten[Lookup[Values[normalized],
        "sessions"]]]];
  epsilonKeys = {"min", "max", "required_complete_max",
    "match_required_complete_max"};
  refinementKeys = {"relative_tolerance", "max_steps"};
  If[Length[sessions] =!= 1 || !StringQ[anchorCheckpoint] ||
      Sort[Keys[epsilon]] =!= Sort[epsilonKeys] ||
      !AllTrue[Lookup[epsilon, epsilonKeys], IntegerQ] ||
      !TrueQ[epsilon["min"] <= epsilon["required_complete_max"] <=
        epsilon["match_required_complete_max"] <= epsilon["max"]] ||
      Sort[Keys[refinement]] =!= Sort[refinementKeys] ||
      !StringQ[Lookup[refinement, "relative_tolerance", None]] ||
      !IntegerQ[Lookup[refinement, "max_steps", None]] ||
      !TrueQ[0 <= refinement["max_steps"] <= 32] ||
      StringLength[checkpointRoot] == 0,
    Return[Failure["CppBackend", <|"Detail" ->
      "consumed-state publication received inconsistent sessions, epsilon/refinement bounds, or checkpoint identities"|>], Module]];
  RunRequest[<|"schema" -> 2,
    "op" -> "transport.publish_consumed_states",
    "session" -> First[sessions],
    "tile_plan" -> planTokens["TilePlan"],
    "tile_plan_checkpoint_identity" -> planTokens["CheckpointIdentity"],
    "anchor" -> anchorTokens["Local"],
    "anchor_checkpoint_identity" -> anchorCheckpoint,
    "epsilon" -> epsilon, "refinement" -> refinement,
    "checkpoint_policy" -> <|
      "schema" -> "diffexp2-deterministic-arm-checkpoints-v1",
      "root" -> checkpointRoot|>,
    "lower" -> KeyDrop[normalized["lower"], "sessions"],
    "upper" -> KeyDrop[normalized["upper"], "sessions"]|>]];

RunPersistentTransportArm[plan_Association, arm_String,
    anchor_Association, receivingBasis_List, epsilon_Association,
    checkpointRoot_String, refinement_Association] := Module[
  {planTokens = persistentTilePlanHandles[plan],
   anchorTokens = persistentLocalHandles[anchor], basisTokens, bad,
   sessions, anchorCheckpoint, epsilonKeys, refinementKeys},
  If[FailureQ[planTokens], Return[planTokens, Module]];
  If[FailureQ[anchorTokens], Return[anchorTokens, Module]];
  epsilonKeys = {"min", "max", "required_complete_max",
    "match_required_complete_max"};
  refinementKeys = {"relative_tolerance", "max_steps"};
  anchorCheckpoint = Lookup[anchor, "checkpoint_identity",
    Lookup[anchor, "CheckpointIdentity", None]];
  If[!MemberQ[{"lower", "upper"}, arm] ||
      !AllTrue[receivingBasis,
        ListQ[#] && # =!= {} && AllTrue[#, AssociationQ] &] ||
      Sort[Keys[epsilon]] =!= Sort[epsilonKeys] ||
      !AllTrue[Lookup[epsilon, epsilonKeys], IntegerQ] ||
      !TrueQ[epsilon["min"] <= epsilon["required_complete_max"] <=
        epsilon["match_required_complete_max"] <= epsilon["max"]] ||
      Sort[Keys[refinement]] =!= Sort[refinementKeys] ||
      !StringQ[Lookup[refinement, "relative_tolerance", None]] ||
      StringLength[Lookup[refinement, "relative_tolerance", ""]] == 0 ||
      !IntegerQ[Lookup[refinement, "max_steps", None]] ||
      !TrueQ[0 <= refinement["max_steps"] <= 32] ||
      !StringQ[anchorCheckpoint] || StringLength[anchorCheckpoint] == 0 ||
      StringLength[checkpointRoot] == 0,
    Return[Failure["CppBackend", <|"Detail" ->
      "retained transport marching requires lower/upper, zero or more complete retained-local bases, ordered integer epsilon bounds, bounded refinement, and nonempty source/run checkpoint identities"|>], Module]];
  basisTokens = Map[persistentLocalHandles, receivingBasis, {2}];
  bad = Cases[basisTokens, _Failure, Infinity];
  If[bad =!= {}, Return[First[bad], Module]];
  sessions = DeleteDuplicates[Join[
    {planTokens["Session"], anchorTokens["Session"]},
    Flatten[Map[Lookup[#, "Session"] &, basisTokens, {2}]]]];
  If[Length[sessions] =!= 1,
    Return[Failure["CppBackend", <|"Detail" ->
      "the retained plan, anchor, and every receiving-basis local must belong to one persistent session",
      "Sessions" -> sessions|>], Module]];
  RunRequest[<|"schema" -> 2, "op" -> "transport.run_arm",
    "session" -> First[sessions],
    "tile_plan" -> planTokens["TilePlan"],
    "tile_plan_checkpoint_identity" -> planTokens["CheckpointIdentity"],
    "anchor" -> anchorTokens["Local"],
    "anchor_checkpoint_identity" -> anchorCheckpoint,
    "arm" -> arm,
    "receiving_basis" -> Map[Lookup[#, "Local"] &, basisTokens, {2}],
    "epsilon" -> epsilon, "refinement" -> refinement,
    "checkpoint_policy" -> <|
      "schema" -> "diffexp2-deterministic-arm-checkpoints-v1",
      "root" -> checkpointRoot|>|>]];

persistentTransportContractStateHandles[state_Association] := Module[
  {tokens = persistentTransportArmHandles[state], provenance, tiles,
   epsilon, requiredCompleteMax},
  If[FailureQ[tokens], Return[tokens, Module]];
  provenance = Lookup[state, "provenance_identity",
    Lookup[state, "ProvenanceIdentity", None]];
  tiles = Lookup[state, "tiles", Lookup[state, "Tiles", Automatic]];
  epsilon = Lookup[state, "epsilon", Lookup[state, "Epsilon", <||>]];
  requiredCompleteMax = If[AssociationQ[epsilon],
    Lookup[epsilon, "required_complete_max",
      Lookup[epsilon, "RequiredCompleteMax", Automatic]], Automatic];
  If[StringLength[tokens["Session"]] == 0 ||
      StringLength[tokens["TransportState"]] == 0 ||
      !StringQ[provenance] || StringLength[provenance] == 0 ||
      (tiles =!= Automatic && (!IntegerQ[tiles] || tiles < 1)) ||
      (requiredCompleteMax =!= Automatic &&
        !IntegerQ[requiredCompleteMax]),
    Return[Failure["CppBackend", <|"Detail" ->
      "transport contraction requires nonempty session, state, checkpoint, and provenance tokens plus valid optional state topology metadata"|>],
      Module]];
  Join[tokens, <|"ProvenanceIdentity" -> provenance,
    "Tiles" -> tiles,
    "RequiredCompleteMax" -> requiredCompleteMax|>]];

persistentNonemptyStringQ[value_] := StringQ[value] &&
  StringLength[value] > 0;

persistentPreparedRationalMultiplierQ[multiplier_] := Module[
  {kernels, keys, baseKeys = {"epsilon_shift", "center_pole_order",
    "kernels", "exact_identity", "proven_zero"}, analytic},
  keys = If[AssociationQ[multiplier], Sort[Keys[multiplier]], {}];
  If[!AssociationQ[multiplier] ||
      !MemberQ[{Sort[baseKeys], Sort[Append[baseKeys,
        "analytic_coefficients"]]}, keys], Return[False, Module]];
  kernels = multiplier["kernels"];
  analytic = Lookup[multiplier, "analytic_coefficients", None];
  IntegerQ[multiplier["epsilon_shift"]] &&
    IntegerQ[multiplier["center_pole_order"]] &&
    multiplier["center_pole_order"] >= 0 &&
    ListQ[kernels] && kernels =!= {} &&
    AllTrue[kernels, ListQ[#] && # =!= {} &] &&
    (analytic === None || (ListQ[analytic] &&
      Length[analytic] === Length[kernels] &&
      AllTrue[analytic, AssociationQ[#] &&
        Sort[Keys[#]] === Sort[{"numerator", "denominator"}] &&
        ListQ[# ["numerator"]] && # ["numerator"] =!= {} &&
        ListQ[# ["denominator"]] && # ["denominator"] =!= {} &])) &&
    persistentNonemptyStringQ[multiplier["exact_identity"]] &&
    multiplier["proven_zero"] === False];

persistentPreparedRationalRowQ[row_] := Module[
  {columns, entries, indices},
  If[!AssociationQ[row] ||
      Sort[Keys[row]] =!= Sort[{
        "schema", "columns", "exact_identity", "entries"}] ||
      Lookup[row, "schema", None] =!=
        "diffexp2-prepared-rational-local-row-v1",
    Return[False, Module]];
  columns = row["columns"];
  entries = row["entries"];
  If[!IntegerQ[columns] || columns < 1 ||
      !persistentNonemptyStringQ[row["exact_identity"]] ||
      !ListQ[entries] ||
      !AllTrue[entries, AssociationQ[#] &&
          Sort[Keys[#]] === Sort[{"column", "multiplier"}] &&
          IntegerQ[#["column"]] && 0 <= #["column"] < columns &&
          persistentPreparedRationalMultiplierQ[#["multiplier"]] &],
    Return[False, Module]];
  indices = If[entries === {}, {}, Lookup[entries, "column"]];
  SameQ[indices, Sort[DeleteDuplicates[indices]]]];

normalizePersistentTransportObservable[observable_, tiles_] := Module[
  {keys = {"Identity", "CheckpointIdentity", "IntegrandRows",
     "Epsilon", "TailPolicy"}, epsilonKeys = {
     "Min", "Max", "RequiredCompleteMax"}, rows, epsilon, columns},
  If[!AssociationQ[observable] ||
      Sort[Keys[observable]] =!= Sort[keys],
    Return[Failure["CppBackend", <|"Detail" ->
      "each transport observable requires exactly Identity, CheckpointIdentity, IntegrandRows, Epsilon, and TailPolicy"|>], Module]];
  rows = observable["IntegrandRows"];
  epsilon = observable["Epsilon"];
  If[!persistentNonemptyStringQ[observable["Identity"]] ||
      !persistentNonemptyStringQ[observable["CheckpointIdentity"]] ||
      !ListQ[rows] || rows === {} ||
      !AllTrue[rows, persistentPreparedRationalRowQ] ||
      (IntegerQ[tiles] && Length[rows] =!= tiles),
    Return[Failure["CppBackend", <|"Detail" ->
      "transport observable identities must be nonempty and IntegrandRows must contain one valid prepared rational row per retained tile",
      "Identity" -> Lookup[observable, "Identity", None]|>], Module]];
  columns = Lookup[rows, "columns"];
  If[Length[DeleteDuplicates[columns]] =!= 1,
    Return[Failure["CppBackend", <|"Detail" ->
      "all prepared rows of one transport observable must have one physical dimension",
      "Identity" -> observable["Identity"]|>], Module]];
  If[!AssociationQ[epsilon] ||
      Sort[Keys[epsilon]] =!= Sort[epsilonKeys] ||
      !AllTrue[Lookup[epsilon, epsilonKeys], IntegerQ] ||
      !TrueQ[epsilon["Min"] <= epsilon["RequiredCompleteMax"] <=
        epsilon["Max"]],
    Return[Failure["CppBackend", <|"Detail" ->
      "transport observable Epsilon requires exact integer Min, Max, and RequiredCompleteMax with Min <= RequiredCompleteMax <= Max",
      "Identity" -> observable["Identity"]|>], Module]];
  If[!MemberQ[{"stored", "attempt", "require"},
      observable["TailPolicy"]],
    Return[Failure["CppBackend", <|"Detail" ->
      "transport observable TailPolicy must be stored, attempt, or require",
      "Identity" -> observable["Identity"]|>], Module]];
  <|"identity" -> observable["Identity"],
    "checkpoint_identity" -> observable["CheckpointIdentity"],
    "integrand_rows" -> rows,
    "epsilon" -> <|"min" -> epsilon["Min"],
      "max" -> epsilon["Max"],
      "required_complete_max" -> epsilon["RequiredCompleteMax"]|>,
    "tail_policy" -> observable["TailPolicy"]|>];

ContractPersistentTransportObservables[state_Association,
    observables_List, checkpointRoot_String] := Module[
  {tokens = persistentTransportContractStateHandles[state], normalized,
   bad, identities, checkpoints, rowCounts, dimensions,
   requiredCompleteMax, request},
  If[FailureQ[tokens], Return[tokens, Module]];
  If[StringLength[StringTrim[checkpointRoot]] == 0,
    Return[Failure["CppBackend", <|"Detail" ->
      "transport contraction checkpoint root must be nonempty"|>],
      Module]];
  normalized = normalizePersistentTransportObservable[
      #, tokens["Tiles"]] & /@ observables;
  bad = Select[normalized, FailureQ];
  If[bad =!= {}, Return[First[bad], Module]];
  identities = Lookup[normalized, "identity"];
  checkpoints = Lookup[normalized, "checkpoint_identity"];
  If[Length[DeleteDuplicates[identities]] =!= Length[identities] ||
      Length[DeleteDuplicates[checkpoints]] =!= Length[checkpoints],
    Return[Failure["CppBackend", <|"Detail" ->
      "transport observable identities and checkpoint identities must each be pairwise unique and retain request order"|>], Module]];
  If[normalized =!= {},
    rowCounts = Length /@ Lookup[normalized, "integrand_rows"];
    dimensions = Flatten[
      Lookup[#, "columns"] & /@ Lookup[normalized, "integrand_rows"]];
    If[Length[DeleteDuplicates[rowCounts]] =!= 1 ||
        Length[DeleteDuplicates[dimensions]] =!= 1,
      Return[Failure["CppBackend", <|"Detail" ->
        "all transport observables must reproduce one retained tile count and physical dimension"|>], Module]];
    requiredCompleteMax = Lookup[
      Lookup[normalized, "epsilon"], "required_complete_max"];
    If[IntegerQ[tokens["RequiredCompleteMax"]] &&
        Max[requiredCompleteMax] > tokens["RequiredCompleteMax"],
      Return[Failure["CppBackend", <|"Detail" ->
        "a transport observable required epsilon maximum exceeds the retained state's public target",
        "StateRequiredCompleteMax" -> tokens["RequiredCompleteMax"],
        "ObservableRequiredCompleteMax" -> Max[requiredCompleteMax]|>],
        Module]]];
  request = <|"schema" -> 2, "op" -> "transport.contract",
    "session" -> tokens["Session"],
    "transport_state" -> tokens["TransportState"],
    "transport_state_checkpoint_identity" ->
      tokens["CheckpointIdentity"],
    "transport_state_provenance_identity" ->
      tokens["ProvenanceIdentity"],
    "checkpoint_policy" -> <|
      "schema" ->
        "diffexp2-deterministic-transport-contraction-checkpoints-v1",
      "root" -> checkpointRoot|>,
    "observables" -> normalized|>;
  RunRequest[request]];

normalizePersistentTransportPairObservable[observable_, lowerTiles_,
    upperTiles_] := Module[
  {baseKeys = {"Identity", "CheckpointIdentity", "LowerIntegrandRows",
      "UpperIntegrandRows", "Epsilon"}, keys, tailPolicy, lower,
   upper, dimensions},
  If[!AssociationQ[observable],
    Return[Failure["CppBackend", <|"Detail" ->
      "each transport-pair observable must be an association"|>], Module]];
  keys = Sort[Keys[observable]];
  If[keys =!= Sort[baseKeys] &&
      keys =!= Sort[Append[baseKeys, "TailPolicy"]],
    Return[Failure["CppBackend", <|"Detail" ->
      "each transport-pair observable requires exactly Identity, CheckpointIdentity, LowerIntegrandRows, UpperIntegrandRows, Epsilon, and optional TailPolicy; signs, arms, points, rims, and cancellation are not accepted"|>], Module]];
  tailPolicy = Lookup[observable, "TailPolicy", "stored"];
  lower = normalizePersistentTransportObservable[<|
    "Identity" -> observable["Identity"],
    "CheckpointIdentity" -> observable["CheckpointIdentity"],
    "IntegrandRows" -> observable["LowerIntegrandRows"],
    "Epsilon" -> observable["Epsilon"],
    "TailPolicy" -> tailPolicy|>, lowerTiles];
  If[FailureQ[lower], Return[lower, Module]];
  upper = normalizePersistentTransportObservable[<|
    "Identity" -> observable["Identity"],
    "CheckpointIdentity" -> observable["CheckpointIdentity"],
    "IntegrandRows" -> observable["UpperIntegrandRows"],
    "Epsilon" -> observable["Epsilon"],
    "TailPolicy" -> tailPolicy|>, upperTiles];
  If[FailureQ[upper], Return[upper, Module]];
  dimensions = Lookup[Join[lower["integrand_rows"],
    upper["integrand_rows"]], "columns"];
  If[Length[DeleteDuplicates[dimensions]] =!= 1,
    Return[Failure["CppBackend", <|"Detail" ->
      "lower and upper prepared rows of one transport-pair observable must have one common physical dimension",
      "Identity" -> observable["Identity"]|>], Module]];
  <|"identity" -> lower["identity"],
    "checkpoint_identity" -> lower["checkpoint_identity"],
    "lower_integrand_rows" -> lower["integrand_rows"],
    "upper_integrand_rows" -> upper["integrand_rows"],
    "epsilon" -> lower["epsilon"],
    "tail_policy" -> tailPolicy|>];

ContractPersistentTransportPairObservables[lowerState_Association,
    upperState_Association, observables_List,
    checkpointRoot_String] := Module[
  {lower = persistentTransportContractStateHandles[lowerState],
   upper = persistentTransportContractStateHandles[upperState],
   lowerArm, upperArm, normalized, bad, identities, checkpoints,
   lowerRowCounts, upperRowCounts, dimensions, requiredCompleteMax,
   stateLimits},
  If[FailureQ[lower], Return[lower, Module]];
  If[FailureQ[upper], Return[upper, Module]];
  lowerArm = Lookup[lowerState, "arm",
    Lookup[lowerState, "Arm", Automatic]];
  upperArm = Lookup[upperState, "arm",
    Lookup[upperState, "Arm", Automatic]];
  If[lower["Session"] =!= upper["Session"] ||
      lower["TransportState"] === upper["TransportState"] ||
      (lowerArm =!= Automatic && lowerArm =!= "lower") ||
      (upperArm =!= Automatic && upperArm =!= "upper"),
    Return[Failure["CppBackend", <|"Detail" ->
      "paired transport contraction requires distinct exact lower/upper retained states from one persistent session"|>], Module]];
  If[StringLength[StringTrim[checkpointRoot]] == 0,
    Return[Failure["CppBackend", <|"Detail" ->
      "transport-pair contraction checkpoint root must be nonempty"|>],
      Module]];
  normalized = normalizePersistentTransportPairObservable[
      #, lower["Tiles"], upper["Tiles"]] & /@ observables;
  bad = Select[normalized, FailureQ];
  If[bad =!= {}, Return[First[bad], Module]];
  identities = Lookup[normalized, "identity"];
  checkpoints = Lookup[normalized, "checkpoint_identity"];
  If[Length[DeleteDuplicates[identities]] =!= Length[identities] ||
      Length[DeleteDuplicates[checkpoints]] =!= Length[checkpoints],
    Return[Failure["CppBackend", <|"Detail" ->
      "transport-pair observable identities and checkpoint identities must each be pairwise unique and retain request order"|>], Module]];
  If[normalized =!= {},
    lowerRowCounts = Length /@
      Lookup[normalized, "lower_integrand_rows"];
    upperRowCounts = Length /@
      Lookup[normalized, "upper_integrand_rows"];
    dimensions = Flatten[
      Lookup[Join[#1, #2], "columns"] & @@@
        Transpose[{Lookup[normalized, "lower_integrand_rows"],
          Lookup[normalized, "upper_integrand_rows"]}]];
    If[Length[DeleteDuplicates[lowerRowCounts]] =!= 1 ||
        Length[DeleteDuplicates[upperRowCounts]] =!= 1 ||
        Length[DeleteDuplicates[dimensions]] =!= 1,
      Return[Failure["CppBackend", <|"Detail" ->
        "all transport-pair observables must reproduce the retained lower/upper tile counts and one common physical dimension"|>], Module]];
    requiredCompleteMax = Lookup[
      Lookup[normalized, "epsilon"], "required_complete_max"];
    stateLimits = {lower["RequiredCompleteMax"],
      upper["RequiredCompleteMax"]};
    If[AnyTrue[Select[stateLimits, IntegerQ],
        Max[requiredCompleteMax] > # &],
      Return[Failure["CppBackend", <|"Detail" ->
        "a transport-pair observable required epsilon maximum exceeds a retained arm state's public target",
        "StateRequiredCompleteMax" -> stateLimits,
        "ObservableRequiredCompleteMax" -> Max[requiredCompleteMax]|>],
        Module]]];
  RunRequest[<|"schema" -> 2, "op" -> "transport.contract_pair",
    "session" -> lower["Session"],
    "lower" -> <|
      "transport_state" -> lower["TransportState"],
      "checkpoint_identity" -> lower["CheckpointIdentity"],
      "provenance_identity" -> lower["ProvenanceIdentity"]|>,
    "upper" -> <|
      "transport_state" -> upper["TransportState"],
      "checkpoint_identity" -> upper["CheckpointIdentity"],
      "provenance_identity" -> upper["ProvenanceIdentity"]|>,
    "checkpoint_policy" -> <|
      "schema" ->
        "diffexp2-deterministic-transport-pair-contraction-checkpoints-v1",
      "root" -> checkpointRoot|>,
    "observables" -> normalized|>]];

normalizePersistentTransportEndpointObservable[observable_] := Module[
  {keys = {"Identity", "CheckpointIdentity", "IntegrandRow",
     "Epsilon"}, epsilonKeys = {
     "Min", "Max", "RequiredCompleteMax"}, epsilon},
  If[!AssociationQ[observable] ||
      Sort[Keys[observable]] =!= Sort[keys],
    Return[Failure["CppBackend", <|"Detail" ->
      "each transport endpoint observable requires exactly Identity, CheckpointIdentity, IntegrandRow, and Epsilon"|>], Module]];
  epsilon = observable["Epsilon"];
  If[!persistentNonemptyStringQ[observable["Identity"]] ||
      !persistentNonemptyStringQ[observable["CheckpointIdentity"]] ||
      !persistentPreparedRationalRowQ[observable["IntegrandRow"]],
    Return[Failure["CppBackend", <|"Detail" ->
      "transport endpoint identities must be nonempty and IntegrandRow must be one valid prepared rational row",
      "Identity" -> Lookup[observable, "Identity", None]|>], Module]];
  If[!AssociationQ[epsilon] ||
      Sort[Keys[epsilon]] =!= Sort[epsilonKeys] ||
      !AllTrue[Lookup[epsilon, epsilonKeys], IntegerQ] ||
      !TrueQ[epsilon["Min"] <= epsilon["RequiredCompleteMax"] <=
        epsilon["Max"]],
    Return[Failure["CppBackend", <|"Detail" ->
      "transport endpoint Epsilon requires exact integer Min, Max, and RequiredCompleteMax with Min <= RequiredCompleteMax <= Max",
      "Identity" -> observable["Identity"]|>], Module]];
  <|"identity" -> observable["Identity"],
    "checkpoint_identity" -> observable["CheckpointIdentity"],
    "integrand_row" -> observable["IntegrandRow"],
    "epsilon" -> <|"min" -> epsilon["Min"],
      "max" -> epsilon["Max"],
      "required_complete_max" -> epsilon["RequiredCompleteMax"]|>|>];

RunPersistentTransportEndpointBatch[state_Association,
    observables_List, checkpointRoot_String] := Module[
  {tokens = persistentTransportContractStateHandles[state], normalized,
   bad, identities, checkpoints, requiredCompleteMax},
  If[FailureQ[tokens], Return[tokens, Module]];
  If[StringLength[StringTrim[checkpointRoot]] == 0,
    Return[Failure["CppBackend", <|"Detail" ->
      "transport endpoint checkpoint root must be nonempty"|>],
      Module]];
  normalized = normalizePersistentTransportEndpointObservable /@
    observables;
  bad = Select[normalized, FailureQ];
  If[bad =!= {}, Return[First[bad], Module]];
  identities = Lookup[normalized, "identity"];
  checkpoints = Lookup[normalized, "checkpoint_identity"];
  If[Length[DeleteDuplicates[identities]] =!= Length[identities] ||
      Length[DeleteDuplicates[checkpoints]] =!= Length[checkpoints],
    Return[Failure["CppBackend", <|"Detail" ->
      "transport endpoint observable identities and checkpoint identities must each be pairwise unique and retain request order"|>], Module]];
  If[normalized =!= {},
    requiredCompleteMax = Lookup[
      Lookup[normalized, "epsilon"], "required_complete_max"];
    If[IntegerQ[tokens["RequiredCompleteMax"]] &&
        Max[requiredCompleteMax] > tokens["RequiredCompleteMax"],
      Return[Failure["CppBackend", <|"Detail" ->
        "a transport endpoint required epsilon maximum exceeds the retained state's public target",
        "StateRequiredCompleteMax" -> tokens["RequiredCompleteMax"],
        "ObservableRequiredCompleteMax" -> Max[requiredCompleteMax]|>],
        Module]]];
  RunRequest[<|"schema" -> 2, "op" -> "transport.endpoint_batch",
    "session" -> tokens["Session"],
    "transport_state" -> tokens["TransportState"],
    "transport_state_checkpoint_identity" ->
      tokens["CheckpointIdentity"],
    "transport_state_provenance_identity" ->
      tokens["ProvenanceIdentity"],
    "checkpoint_policy" -> <|
      "schema" ->
        "diffexp2-deterministic-transport-endpoint-checkpoints-v1",
      "root" -> checkpointRoot|>,
    "observables" -> normalized|>]];

PersistentTransportArmStatistics[handle_Association] := Module[
  {tokens = persistentTransportArmHandles[handle]},
  If[FailureQ[tokens], Return[tokens, Module]];
  RunRequest[<|"schema" -> 2, "op" -> "transport.stats",
    "session" -> tokens["Session"],
    "transport_state" -> tokens["TransportState"]|>]];

ReleasePersistentTransportArm[handle_Association] := Module[
  {tokens = persistentTransportArmHandles[handle]},
  If[FailureQ[tokens], Return[tokens, Module]];
  RunRequest[<|"schema" -> 2, "op" -> "transport.release",
    "session" -> tokens["Session"],
    "transport_state" -> tokens["TransportState"]|>]];

persistentLineIntegralHandles[handle_Association] := Module[
  {session, line, checkpoint},
  session = Lookup[handle, "session", Lookup[handle, "Session", None]];
  line = Lookup[handle, "line", Lookup[handle, "Line", None]];
  checkpoint = Lookup[handle, "checkpoint_identity",
    Lookup[handle, "CheckpointIdentity", None]];
  If[!StringQ[session] || !StringQ[line] || !StringQ[checkpoint] ||
      StringLength[checkpoint] == 0,
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent line-integral handle requires exact session, line, and checkpoint tokens"|>],
      Module]];
  <|"Session" -> session, "Line" -> line,
    "CheckpointIdentity" -> checkpoint|>];

PersistentLineIntegralStatistics[handle_Association] := Module[
  {tokens = persistentLineIntegralHandles[handle]},
  If[FailureQ[tokens], Return[tokens, Module]];
  RunRequest[<|"schema" -> 2, "op" -> "integration.stats",
    "session" -> tokens["Session"], "line" -> tokens["Line"]|>]];

ExportPersistentLineIntegral[handle_Association,
    checkpointIdentity_String, outputDigits_:Automatic] := Module[
  {tokens = persistentLineIntegralHandles[handle], request},
  If[FailureQ[tokens], Return[tokens, Module]];
  If[StringLength[checkpointIdentity] == 0 ||
      (outputDigits =!= Automatic &&
       (!IntegerQ[outputDigits] || outputDigits < 1)),
    Return[Failure["CppBackend", <|"Detail" ->
      "line export requires a nonempty checkpoint identity and positive output digits"|>],
      Module]];
  request = <|"schema" -> 2, "op" -> "integration.export",
    "session" -> tokens["Session"], "line" -> tokens["Line"],
    "checkpoint_identity" -> checkpointIdentity|>;
  If[IntegerQ[outputDigits],
    request = Append[request, "output_digits" -> outputDigits]];
  RunRequest[request]];

ReleasePersistentLineIntegral[handle_Association] := Module[
  {tokens = persistentLineIntegralHandles[handle]},
  If[FailureQ[tokens], Return[tokens, Module]];
  RunRequest[<|"schema" -> 2, "op" -> "integration.release",
    "session" -> tokens["Session"], "line" -> tokens["Line"]|>]];

PersistentLocalStatistics[handle_Association] := Module[
  {tokens = persistentLocalHandles[handle]},
  If[FailureQ[tokens], Return[tokens, Module]];
  RunRequest[<|"schema" -> 2, "op" -> "local.stats",
    "session" -> tokens["Session"], "local" -> tokens["Local"]|>]];

ReleasePersistentLocal[handle_Association] := Module[
  {tokens = persistentLocalHandles[handle]},
  If[FailureQ[tokens], Return[tokens, Module]];
  RunRequest[<|"schema" -> 2, "op" -> "local.release",
    "session" -> tokens["Session"], "local" -> tokens["Local"]|>]];

persistentCheckpointSession[owner_String] := If[
  StringLength[owner] > 0, owner,
  Failure["CppBackend", <|"Detail" ->
    "persistent checkpoint session token must be nonempty"|>]];
persistentCheckpointSession[owner_Association] := Module[{session},
  session = Lookup[owner, "session", Lookup[owner, "Session", None]];
  If[StringQ[session] && StringLength[session] > 0, session,
    Failure["CppBackend", <|"Detail" ->
      "persistent checkpoint owner requires an exact session token"|>]]];
persistentCheckpointSession[_] := Failure["CppBackend", <|"Detail" ->
  "persistent checkpoint owner must be a session string or opaque persistent handle"|>];

SavePersistentCheckpoint[owner_, path_String, identity_String] := Module[
  {session = persistentCheckpointSession[owner], expanded},
  If[FailureQ[session], Return[session, Module]];
  If[StringLength[identity] == 0,
    Return[Failure["CppBackend", <|"Detail" ->
      "persistent checkpoint identity must be nonempty"|>], Module]];
  expanded = ExpandFileName[path];
  RunRequest[<|"schema" -> 2, "op" -> "checkpoint.save",
    "session" -> session, "path" -> expanded,
    "checkpoint_identity" -> identity|>]];

RestorePersistentCheckpoint[path_String, expectedIdentity_String] := Module[
  {expanded, response, session},
  If[StringLength[expectedIdentity] == 0,
    Return[Failure["CppBackend", <|"Detail" ->
      "expected persistent checkpoint identity must be nonempty"|>], Module]];
  expanded = ExpandFileName[path];
  response = RunRequest[<|"schema" -> 2, "op" -> "checkpoint.restore",
    "path" -> expanded, "expected_identity" -> expectedIdentity|>];
  If[persistentCommandOKQ[response],
    session = Lookup[response, "session", None];
    If[StringQ[session],
      AssociateTo[$persistentRestoredSessionHandles, session -> True]]];
  response];

ClosePersistentSession[owner_] := Module[
  {session = persistentCheckpointSession[owner], response},
  If[FailureQ[session], Return[session, Module]];
  response = RunRequest[<|"schema" -> 2, "op" -> "session.close",
    "session" -> session|>];
  If[persistentCommandOKQ[response], persistentForgetSessionHandle[session]];
  response];

ClearPersistentSessions[] := Module[{handles},
  handles = DeleteDuplicates@Join[
    (Lookup[#, "Handle"] & /@ Values[$persistentSessionCache]),
    Keys[$persistentRestoredSessionHandles]];
  Scan[Function[handle, Quiet[Check[RunRequest[<|"schema" -> 2,
      "op" -> "session.close", "session" -> handle|>], Null]]], handles];
  $persistentSessionCache = <||>;
  $persistentChartCache = <||>;
  $persistentSCCCache = <||>;
  $persistentPreparedTokenCache = <||>;
  $persistentRestoredSessionHandles = <||>;
  Null];

PersistentSessionInformation[] := Module[{handles},
  handles = DeleteDuplicates@Join[
    (Lookup[#, "Handle"] & /@ Values[$persistentSessionCache]),
    Keys[$persistentRestoredSessionHandles]];
  AssociationMap[RunRequest[<|"schema" -> 2, "op" -> "session.stats",
      "session" -> #|>] &, handles]];

End[];
EndPackage[];
