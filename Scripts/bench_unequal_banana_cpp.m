(* Direct DiffExp2/C++ benchmark for the legacy 15-master unequal-mass
   three-loop banana system.

   The Euclidean route starts from the independently frozen equal-mass
   checkpoint at p^2 = -1 and deforms only the squared masses

     {1,1,1,1} -> {2,3/2,4/3,1}.

   This avoids an artificial crossing of the p^2 = 0 differential-equation
   singularity.  The old slice exports are exact for this fixture: every
   order >= eps^2 is a structural-zero matrix, so M0 + eps M1 reconstructs
   the full rational system rather than an epsilon truncation.

   Environment knobs:
     UB_BACKEND=Cpp|Wolfram       recurrence backend (default Cpp)
     UB_WP=250                    working precision; WP500 requires a genuine
                                  high-precision seed (>=80 digits advised)
     UB_CHOP=34                   internal structural-zero/matching guard;
                                  this is not a claim of 34 output digits
     UB_REPORT_DIGITS=25          displayed central-value digits only;
                                  measured agreement is reported separately
     UB_EXPANSION_ORDER=50        Taylor order
     UB_EPSILON_ORDER=4           highest delivered epsilon order (1..4,
                                  fixed by the frozen boundary checkpoint)
     UB_DIVISION_ORDER=3          classic coupled matching fraction
     UB_ROC=10                    chart-coordinate rescaling
     UB_VALUE_TRANSPORT=0         classic basis-and-matching transport
     UB_PREFLIGHT_ONLY=0|1        print boundary/residual budget and exit
     DE2_CPP_THREADS=4            native homogeneous-column workers

   Run from the repository root:
     DE2_CPP_THREADS=4 wolframscript -file Scripts/bench_unequal_banana_cpp.m

   Output rows are prefixed UBENCH and use compact JSON. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];
scriptStart = AbsoluteTime[];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

envOr[name_String, default_String] := Module[{v = Environment[name]},
  If[StringQ[v] && StringLength[StringTrim[v]] > 0, v, default]];

wp = ToExpression[envOr["UB_WP", "250"]];
chop = ToExpression[envOr["UB_CHOP", "34"]];
reportDigits = ToExpression[envOr["UB_REPORT_DIGITS", "25"]];
expansionOrder = ToExpression[envOr["UB_EXPANSION_ORDER", "50"]];
epsilonOrder = ToExpression[envOr["UB_EPSILON_ORDER", "4"]];
divisionOrder = ToExpression[envOr["UB_DIVISION_ORDER", "3"]];
roc = ToExpression[envOr["UB_ROC", "10"]];
backend = envOr["UB_BACKEND", "Cpp"];
valueTransport = envOr["UB_VALUE_TRANSPORT", "0"] === "1";
preflightOnly = envOr["UB_PREFLIGHT_ONLY", "0"] === "1";

If[!IntegerQ[epsilonOrder] || !Between[epsilonOrder, {1, 4}],
  Print["UB_ERROR UB_EPSILON_ORDER must be an integer in 1..4"];
  Exit[2]];
If[!IntegerQ[reportDigits] || reportDigits < 1 || reportDigits > wp,
  Print["UB_ERROR UB_REPORT_DIGITS must be an integer in 1..UB_WP"];
  Exit[2]];

SetEnvironment["DE2_VALUE_TRANSPORT" -> If[valueTransport, "1", "0"]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Catch[expr, "DiffExp2Error"];

emit[phase_String, seconds_, extra_:<||>] := Print["UBENCH ", ExportString[
  Join[<|"Phase" -> phase, "Seconds" -> Round[N[seconds], 1.*^-6]|>, extra],
  "RawJSON", "Compact" -> True]];

(* Read a frozen checkpoint without JSON numericization.  Converting the
   Mathematica-shaped numeric substring directly preserves all 40 committed
   decimal digits; Import[...,"RawJSON"] would first round them to machine
   reals. *)
checkpointValues[file_String, selector_] := Module[
  {line, text, raw},
  line = SelectFirst[Import[file, "Lines"], selector, Missing["NotFound"]];
  If[MissingQ[line],
    Print["UB_ERROR checkpoint row not found in ", file]; Exit[2]];
  text = First@StringCases[line,
    RegularExpression["\\\"Values\\\":(\\[.*\\]),\\\"EpsOrders\\\""] :> "$1"];
  raw = ToExpression[StringReplace[text, {
    "[" -> "{", "]" -> "}",
    RegularExpression["([0-9.])e([+-]?[0-9]+)"] :> "$1*^$2"}]];
  Map[#[[1]] + I #[[2]] &, raw, {2}]];

equalMassOracle = FileNameJoin[
  {repoRoot, "Tests", "refs", "oracle_checkpoints_banana_equalmass.log"}];
equalMass4 = checkpointValues[equalMassOracle,
  StringContainsQ[#, "\"Line\":\"t_line_m1_to_10\",\"SegmentIndex\":1,\"Position\":\"Start\""] &];
equalMass4 = equalMass4[[All, 1 ;; epsilonOrder + 1]];

(* Legacy master ordering, frozen by Banana_example.m and the old checkpoint
   harness: six copies of equal-mass master 1, four of master 2, one of
   master 3, and four of master 4. *)
masterMap = {1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 4, 4, 4, 4};
boundary15 = equalMass4[[masterMap]];
preferredFile = FileNameJoin[{repoRoot, "Dependencies", "fire", "FIRE6",
  "examples", "bananaUnequal.preferred"}];
legacyMasterIndices = #[[2, 1 ;; 4]] & /@ Get[preferredFile];
If[Length[legacyMasterIndices] =!= 15 || legacyMasterIndices[[11]] =!= {1, 1, 1, 1},
  Print["UB_ERROR unexpected bananaUnequal.preferred master ordering"];
  Exit[2]];

(* Necessary input-accuracy gate, deliberately before the expensive matrix
   reconstruction and SegmentLine phases.  ResidTol is an independent
   recurrence certificate, while ChopPrecision is an internal
   structural-zero/matching guard rather than a promise of output digits.
   WP500 asks for 50
   residual digits and cannot honestly be certified from this frozen
   ~40-digit seed.  WP250 asks for 25 and retains measured input headroom.
   This gate is necessary, not a condition-number proof; a genuinely
   >=80-digit seed is recommended before selecting WP500 so the d=15
   propagation still has ample cancellation headroom. *)
ubBoundaryAccuracyPreflight[values_, wp_Integer, chop_Integer] := Module[
  {accs, minAcc, residTol, residDigits, required, exactInput, pass},
  accs = Cases[Flatten[values], z_?InexactNumberQ :> Accuracy[z]];
  accs = Select[accs, NumericQ[#] && # =!= Infinity &];
  exactInput = accs === {};
  minAcc = If[exactInput, Infinity, Min[accs]];
  residTol = DiffExp2`Tolerances`ResidTol[wp];
  residDigits = Floor[-Log[10, residTol]];
  required = Max[residDigits, chop] + DiffExp2`Tolerances`$SafetyDigits;
  pass = exactInput || TrueQ[minAcc >= required];
  <|"Pass" -> pass,
    "WorkingPrecision" -> wp,
    "ChopDigits" -> chop,
    "ResidTol" -> residTol,
    "ResidDigits" -> residDigits,
    "SafetyDigits" -> DiffExp2`Tolerances`$SafetyDigits,
    "RequiredInputAccuracy" -> required,
    "BoundaryMinAccuracy" -> If[exactInput, "Exact", N[minAcc, 8]],
    "InputHeadroomDigits" -> If[exactInput, "Exact", N[minAcc - required, 8]],
    "WP500RecommendedSeedDigits" -> 80|>];

boundaryPreflight = ubBoundaryAccuracyPreflight[boundary15, wp, chop];
Print["UBPREFLIGHT ", ExportString[boundaryPreflight,
  "RawJSON", "Compact" -> True]];
If[preflightOnly, Exit[If[TrueQ[boundaryPreflight["Pass"]], 0, 2]]];
If[!TrueQ[boundaryPreflight["Pass"]],
  Print["UB_ERROR frozen boundary cannot certify the requested residual/internal-guard " <>
    "budget; lower UB_WP consistently with the requested digits or provide a " <>
    "genuine >=80-digit seed for WP500"];
  Exit[2]];

matrixDir = FileNameJoin[{repoRoot, "Tests", "Banana_Matrices"}];
eps = Global`eps; x = Global`x;
psq = Global`psq; mm1 = Global`mm1; mm2 = Global`mm2;
mm3 = Global`mm3; mm4 = Global`mm4;

partialMatrix[var_Symbol, order_Integer] := Get[FileNameJoin[{matrixDir,
  "d" <> SymbolName[var] <> "_" <> ToString[order] <> ".m"}]];

(* Assert rather than assume that the stored slices really define an exact
   affine-in-epsilon system.  This is the analytic-regularization contract
   that makes reconstruction from legacy slices legitimate for this fixture. *)
sliceAuditTime = First@AbsoluteTiming[
  higherSlicesZero = And @@ Flatten@Table[
    partialMatrix[var, order] === ConstantArray[0, {15, 15}],
    {var, {psq, mm1, mm2, mm3, mm4}}, {order, 2, 4}]];
If[!TrueQ[higherSlicesZero],
  Print["UB_ERROR legacy matrices contain nonzero eps^2..eps^4 slices; " <>
    "an exact full export is required"];
  Exit[2]];
emit["SliceAudit", sliceAuditTime, <|"ExactAffineEpsilon" -> True|>];

matrixBuildTime = First@AbsoluteTiming[
  dmm1 = partialMatrix[mm1, 0] + eps partialMatrix[mm1, 1];
  dmm2 = partialMatrix[mm2, 0] + eps partialMatrix[mm2, 1];
  dmm3 = partialMatrix[mm3, 0] + eps partialMatrix[mm3, 1];
  (* Along the Euclidean mass deformation:
       psq=-1, mm1=1+x, mm2=1+x/2, mm3=1+x/3, mm4=1. *)
  massMatrix = Map[Cancel[Together[#]] &,
    (dmm1 + dmm2/2 + dmm3/3) /. {
      psq -> -1, mm1 -> 1 + x, mm2 -> 1 + x/2,
      mm3 -> 1 + x/3, mm4 -> 1}, {2}]];
Clear[dmm1, dmm2, dmm3];
emit["MatrixBuild", matrixBuildTime, <|"Dimension" -> 15|>];

configTime = First@AbsoluteTiming[
  configResult = catchDE2[DiffExp2`Config`LoadConfiguration[{
    "WorkingPrecision" -> wp,
    "ChopPrecision" -> chop,
    "LinearSolveChopPrecision" -> chop,
    "ExpansionOrder" -> expansionOrder,
    "EpsilonOrder" -> epsilonOrder,
    "DivisionOrder" -> divisionOrder,
    "StepDivisionOrder" -> divisionOrder,
    "RadiusOfConvergence" -> roc,
    "UsePade" -> True,
    "RecurrenceBackend" -> backend,
    "Variables" -> {}}]]];
If[FailureQ[configResult], Print["UB_ERROR configuration: ", configResult]; Exit[1]];
emit["Configure", configTime, <|"Backend" -> backend,
  "WorkingPrecision" -> wp, "ChopPrecision" -> chop,
  "ReportDigits" -> reportDigits,
  "ExpansionOrder" -> expansionOrder, "EpsilonOrder" -> epsilonOrder,
  "DivisionOrder" -> divisionOrder, "RadiusOfConvergence" -> roc,
  "ValueTransport" -> valueTransport,
  "CppThreads" -> If[backend === "Cpp",
    ToExpression[envOr["DE2_CPP_THREADS", "4"]], 0],
  "BoundarySource" -> "oracle_checkpoints_banana_equalmass.log:p2=-1",
  "BoundaryMinAccuracy" -> boundaryPreflight["BoundaryMinAccuracy"],
  "BoundaryInputHeadroomDigits" ->
    boundaryPreflight["InputHeadroomDigits"]|>];

loadTime = First@AbsoluteTiming[
  sys = catchDE2[DiffExp2`API`LoadSystem[
    <|"Matrix" -> massMatrix, "Variable" -> x|>]]];
If[FailureQ[sys], Print["UB_ERROR LoadSystem: ", sys]; Exit[1]];
emit["LoadSystem", loadTime, <|"SingularFactors" -> Length[sys["SingularFactors"]]|>];

planTime = First@AbsoluteTiming[
  plan = catchDE2[DiffExp2`Transport`SegmentLine[sys, {0, 1}]]];
If[FailureQ[plan], Print["UB_ERROR SegmentLine: ", plan]; Exit[1]];
emit["Plan", planTime, <|
  "Charts" -> plan["SegmentCount"],
  "SingularCharts" -> Count[plan["Charts"], c_ /; TrueQ[c["Singular"]]],
  "RealSingularities" -> Length[plan["Singularities", "Real"]],
  "ProjectedWaypoints" -> Length[plan["Singularities", "ProjectionWaypoints"]]|>];

transportTime = First@AbsoluteTiming[
  result = catchDE2[DiffExp2`Transport`TransportLine[sys, boundary15, plan]]];
If[FailureQ[result], Print["UB_ERROR TransportLine: ", result]; Exit[1]];
emit["Transport", transportTime, <|
  "Charts" -> result["SegmentCount"],
  "EndpointIsSingular" -> result["EndpointIsSingular"],
  "MaxErrorEstimate" -> N[Max[result["ErrorEstimate"]], 8]|>];

(* Emit compact endpoint coefficients for the cross-route FT comparison.
   Values are component vectors inside one honest EpsSeries. *)
endpointES = result["Value"];
endpoint = Transpose@Table[
  DiffExp2`EpsSeries`ESCoefficient[endpointES, k],
  {k, 0, epsilonOrder}];
(* The canonical top master is component 11.  Its normalization, already
   visible in the equal-mass boundary formula and independently pinned by
   the FT banana result, is

       J_11 = eps (1 + 3 eps) (1 + 4 eps) Exp[3 EulerGamma eps]
              I[1,1,1,1]_FT.

   In particular, the physical scalar integral's finite value is the eps^1
   coefficient of J_11.  Remove only a certified numerical-zero eps^0
   residue before displaying the longer de-normalized series. *)
topCanonical = endpoint[[11]];
topEps0Scale = Max[1, Max[Abs[N[topCanonical, Min[chop, 20]]]]];
topEps0Certified = TrueQ[Abs[N[First[topCanonical], Min[chop, 20]]] <
  10^-chop topEps0Scale];
If[topEps0Certified, topCanonical[[1]] = 0];
rawTopCoefficients = Table[SeriesCoefficient[
  Sum[topCanonical[[j + 1]] eps^j, {j, 0, epsilonOrder}]/
      (eps (1 + 3 eps) (1 + 4 eps) Exp[3 EulerGamma eps]), {eps, 0, k}],
    {k, 0, epsilonOrder - 1}];
(* Independent two-dimensional coordinate-space oracle at this Euclidean
   point:
     8 Integrate[r J0[r] K0[Sqrt[2]r] K0[Sqrt[3/2]r]
       K0[Sqrt[4/3]r] K0[r], {r,0,Infinity}].
   The same quadrature at unit masses reproduces the committed 8.268104...
   checkpoint, fixing conventions as well as the numerical value. *)
independentFiniteOracle =
  5.83402729266214946740741989567969814964058746213209;
oracleAbsDifference = Abs[rawTopCoefficients[[1]] - independentFiniteOracle];
oracleAgreementDigits = If[TrueQ[PossibleZeroQ[oracleAbsDifference]], chop,
  Max[0, -Log[10, N[oracleAbsDifference/independentFiniteOracle, 20]]]];
Print["UBRESULT ", ExportString[<|
  "Point" -> <|"psq" -> -1, "mm1" -> 2, "mm2" -> 3/2,
    "mm3" -> 4/3, "mm4" -> 1|>,
  "MasterOrdering" -> masterMap,
  "LegacyMasterIndices" -> legacyMasterIndices,
  "EpsOrders" -> Range[0, epsilonOrder],
  "Values" -> Map[{N[Re[#], reportDigits], N[Im[#], reportDigits]} &,
    endpoint, {2}],
  "ScalarIntegral" -> <|
    "LegacyMasterIndex" -> 11,
    "Normalization" ->
      "J11=eps(1+3eps)(1+4eps)Exp(3EulerGamma eps)I1111_FT",
    "CanonicalEps0CertifiedZero" -> topEps0Certified,
    "Finite" -> {N[Re[rawTopCoefficients[[1]]], reportDigits],
      N[Im[rawTopCoefficients[[1]]], reportDigits]},
    "IndependentBesselOracle" -> N[independentFiniteOracle, reportDigits],
    "OracleAbsDifference" -> N[oracleAbsDifference, Min[reportDigits, 12]],
    "OracleAgreementDigits" -> N[oracleAgreementDigits, 8],
    "ReportedDigitsAreCertified" -> False,
    "RawCoefficients" -> Map[
      {N[Re[#], reportDigits], N[Im[#], reportDigits]} &,
      rawTopCoefficients]|>,
  "ErrorEstimate" -> N[result["ErrorEstimate"], 8]|>,
  "RawJSON", "Compact" -> True]];

emit["TotalUnequalRoute", sliceAuditTime + matrixBuildTime + configTime +
    loadTime + planTime + transportTime,
  <|"IncludesMatrixReconstruction" -> True,
    "IncludesFrozenEqualMassSeedGeneration" -> False|>];
emit["WallTotal", AbsoluteTime[] - scriptStart,
  <|"IncludesPackageLoadAndOutput" -> True,
    "IncludesFrozenEqualMassSeedGeneration" -> False|>];

Exit[0];
