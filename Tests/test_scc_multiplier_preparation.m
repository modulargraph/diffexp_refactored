(* One direct parity check for the reusable rational/SCC multiplier seam. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

x = Global`x; t = Global`t; eps = Global`eps;
DiffExp2`Config`LoadConfiguration[{
  "WorkingPrecision" -> 80, "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 1, "Variables" -> {}}];

ls = <|"Center" -> 0, "ChartMap" -> <|"Center" -> 0, "Scale" -> 1|>,
  "Radius" -> 1,
  "Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
    "Coeffs" -> Table[{10 (k + 2) + n}, {k, -1, 1}, {n, 0, 3}]|>},
  "EpsWindow" -> <|"Min" -> -1, "CompleteMax" -> 1|>,
  "TWindow" -> <|"CompleteMax" -> 3|>,
  "ErrorEstimate" -> {0, 0, 0}, "Prescriptions" -> {}|>;

sys = <|"Variable" -> x, "Matrix" -> {
    {0, 0}, {(1 + eps x)/(eps x^2 (2 - x)), 0}}|>;
cs = DiffExp2`Solve`PrepareChart[sys, <|"ChartVar" -> t,
  "Center" -> 0, "Scale" -> 1, "Radius" -> 1, "LocalRadius" -> 1,
  "Name" -> "scc-multiplier-preparation", "Prescriptions" -> {},
  "UseSCCSkeleton" -> True|>];
coupling = cs["ThetaOriginal"][[2, 1]];
prepared = DiffExp2`SectorSeries`PrepareRationalMultiplier[ls, coupling, t];
product = DiffExp2`SectorSeries`MultiplyRational[ls, coupling, t];
zeroPrepared = DiffExp2`SectorSeries`PrepareRationalMultiplier[ls, 0, t];
inexactZero = 0``40;
inexactZeroPrepared = DiffExp2`SectorSeries`PrepareRationalMultiplier[
  ls, inexactZero, t];
inexactZeroProduct = DiffExp2`SectorSeries`MultiplyRational[
  ls, inexactZero, t];

kmin = ls["EpsWindow", "Min"];
kmax = ls["EpsWindow", "CompleteMax"];
jmin = prepared["EpsilonShift"];
kernels = prepared["TaylorKernels"];
arr = ls["Sectors"][[1, "Coeffs"]];
manual = Table[Sum[
    If[kmin <= k - jmin - j <= kmax,
      kernels[[j + 1, m + 1]]*
        arr[[k - jmin - j - kmin + 1, n - m + 1, c]], 0],
    {j, 0, Length[kernels] - 1}, {m, 0, n}],
  {k, kmin + jmin, kmax + jmin}, {n, 0, 3}, {c, 1}];

(* Signed epsilon shifts are serializer data.  Native scc.prepare retains
   them; its later solve work contract is responsible for halo sufficiency. *)
sparse = DiffExp2`Solve`PrepareSCCCouplingMatrix[cs, 1, 2, ls,
  <|"domain" -> "rational", "symbols" -> {}|>];
entry = First[sparse["entries"]];
parentRecords = DiffExp2`Solve`Private`sccParentExactRecords[cs];
originalCell = DiffExp2`Solve`Private`sccExactCellRecord[
  sys["Matrix"][[2, 1]], x];
thetaCell = DiffExp2`Solve`Private`sccExactCellRecord[coupling, t];
identityRecord = ImportString[prepared["ExactIdentity"], "RawJSON"];
groupIdentityRecord = ImportString[sparse["exact_identity"], "RawJSON"];
decodedKernels = Map[
  DiffExp2`CppBackend`DecodeScalar[#, 80] &,
  entry["multiplier", "kernels"], {2}];

(* The Acb preparation path deliberately grounds giant exact rationals once
   at 2x WP.  An exact Rational payload must instead retain the original Q
   coefficient; it may never inherit the Acb cache entry or reconstruct the
   value from that rounded Real.  This size is well above the production
   500-byte grounding threshold while remaining cheap to prepare. *)
largeRational = (10^2400 + 123456789)/(10^2400 + 987654321);
largeMultiplier = largeRational/(2 - t);
largeAcbPrepared = DiffExp2`SectorSeries`PrepareRationalMultiplier[
  ls, largeMultiplier, t, "SerializationDomain" -> "acb"];
largeRationalPrepared = DiffExp2`SectorSeries`PrepareRationalMultiplier[
  ls, largeMultiplier, t, "SerializationDomain" -> "rational"];
largeAcbReplay = DiffExp2`SectorSeries`PrepareRationalMultiplier[
  ls, largeMultiplier, t, "SerializationDomain" -> "acb"];
largeSolveAcb = Block[{DiffExp2`Solve`Private`$cppExactDomain = False},
  DiffExp2`Solve`Private`preparedEpsCoefficient[largeRational]];
largeSolveRational = Block[{
    DiffExp2`Solve`Private`$cppExactDomain = True},
  DiffExp2`Solve`Private`preparedEpsCoefficient[largeRational]];

ok = prepared["EpsilonShift"] === -1 &&
  prepared["CenterPoleOrder"] === 1 &&
  prepared["ProvenZero"] === False &&
  zeroPrepared["ProvenZero"] === True &&
  inexactZeroPrepared["ProvenZero"] === False &&
  inexactZeroProduct["EpsWindow"] === ls["EpsWindow"] &&
  !FreeQ[inexactZeroProduct["Sectors"][[All, "Coeffs"]],
    _?InexactNumberQ] &&
  StringQ[prepared["ExactIdentity"]] &&
  identityRecord["schema"] === "diffexp2-exact-expression-v1" &&
  identityRecord["variable"] === <|"node" -> "symbol",
    "context" -> "Global`", "name" -> "t"|> &&
  product["EpsWindow"] === <|"Min" -> -2, "CompleteMax" -> 0|> &&
  product["Sectors"][[1, "a"]] === -1 &&
  product["Sectors"][[1, "Coeffs"]] === manual &&
  sparse["domain"] === "rational" && sparse["symbols"] === {} &&
  {sparse["source_block"], sparse["target_block"],
    sparse["source_vertices"], sparse["target_vertices"]} ===
      {0, 1, {0}, {1}} &&
  {entry["row"], entry["column"], entry["source_vertex"],
    entry["target_vertex"]} === {0, 0, 0, 1} &&
  entry["exact_original_entry"] === originalCell["exact"] &&
  entry["exact_theta_entry"] === thetaCell["exact"] ===
    prepared["ExactIdentity"] &&
  originalCell["proven_zero"] === thetaCell["proven_zero"] === False &&
  Dimensions[parentRecords["exact_system_record"]] === {2, 2} &&
  Dimensions[parentRecords["exact_theta_record"]] === {2, 2} &&
  parentRecords["exact_system_record"][[2, 1]] === originalCell &&
  parentRecords["exact_theta_record"][[2, 1]] === thetaCell &&
  entry["multiplier", "epsilon_shift"] === -1 &&
  entry["multiplier", "exact_identity"] === prepared["ExactIdentity"] &&
  entry["multiplier", "proven_zero"] ===
    parentRecords["exact_system_record"][[2, 1]]["proven_zero"] === False &&
  groupIdentityRecord["schema"] === "diffexp2-scc-coupling-v1" &&
  groupIdentityRecord["serialization", "domain"] === "rational" &&
  groupIdentityRecord["serialization", "symbols"] === {} &&
  decodedKernels === kernels &&
  ByteCount[largeRational] > 500 &&
  !FreeQ[largeAcbPrepared["TaylorKernels"], _?InexactNumberQ] &&
  FreeQ[largeRationalPrepared["TaylorKernels"], _?InexactNumberQ] &&
  largeRationalPrepared["TaylorKernels"][[1, 1]] ===
    largeRational/2 &&
  largeRationalPrepared["ExactIdentity"] ===
    largeAcbPrepared["ExactIdentity"] &&
  largeAcbReplay === largeAcbPrepared &&
  InexactNumberQ[largeSolveAcb] &&
  largeSolveRational === largeRational;

Print[If[TrueQ[ok], "PASS", "FAIL"],
  ": reusable rational and sparse SCC multiplier preparation parity"];
If[!TrueQ[ok], Exit[1]];
