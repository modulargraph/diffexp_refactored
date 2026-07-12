(* The public rational-row handoff must retain a full analytic completion;
   finite Taylor kernels alone cannot certify the projected integration tail. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

x = Global`x; t = Global`t;
DiffExp2`Config`LoadConfiguration[{
  "WorkingPrecision" -> 80, "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 0, "Variables" -> {}}];

chart = <|"SystemSize" -> 1, "ChartVar" -> t, "Center" -> 0,
  "ChartMap" -> <|"Center" -> 0, "Scale" -> 1|>,
  "Radius" -> 2, "Prescriptions" -> {},
  "Name" -> "rational-row-tail-payload"|>;
shape = <|"EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 2|>,
  "TWindow" -> <|"CompleteMax" -> 4|>, "Dimension" -> 1|>;
row = DiffExp2`Solve`PrepareNativeRationalRow[
  chart, shape, {1/(1 - x/2)}, x,
  <|"domain" -> "rational", "symbols" -> {}|>];
multiplier = row["entries"][[1, "multiplier"]];
analytic = multiplier["analytic_coefficients"][[1]];
kernel = DiffExp2`CppBackend`DecodeScalar[#, 80] & /@
  multiplier["kernels"][[1]];
numerator = DiffExp2`CppBackend`DecodeScalar[#, 80] & /@
  analytic["numerator"];
denominator = DiffExp2`CppBackend`DecodeScalar[#, 80] & /@
  analytic["denominator"];
zeroAnalyticNumerators = Map[
  DiffExp2`CppBackend`DecodeScalar[#, 80] & /@ #["numerator"] &,
  Rest[multiplier["analytic_coefficients"]]];
replayed = ConstantArray[0, Length[kernel]];
Do[replayed[[n + 1]] = Together[(
    If[n + 1 <= Length[numerator], numerator[[n + 1]], 0] -
      Sum[If[lag + 1 <= Length[denominator], denominator[[lag + 1]], 0]*
        replayed[[n - lag + 1]], {lag, 1, n}]) / denominator[[1]]],
  {n, 0, Length[kernel] - 1}];

ok = row["schema"] === "diffexp2-prepared-rational-local-row-v1" &&
  TrueQ[DiffExp2`CppBackend`Private`persistentPreparedRationalRowQ[row]] &&
  Sort[Keys[multiplier]] === Sort[{"epsilon_shift", "center_pole_order",
    "kernels", "analytic_coefficients", "exact_identity", "proven_zero"}] &&
  Length[multiplier["analytic_coefficients"]] ===
    Length[multiplier["kernels"]] && denominator =!= {} &&
  denominator[[1]] =!= 0 && replayed === kernel &&
  zeroAnalyticNumerators === {{0}, {0}};

Print[If[TrueQ[ok], "PASS", "FAIL"],
  ": native rational-row analytic tail payload replays every finite kernel"];
If[!TrueQ[ok], Exit[1]];
