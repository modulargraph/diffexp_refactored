(* Closed-form integral pins for the full transport+integrate chain.
   These caught the ApplyCrossing per-column phase bug (the crossing
   phase is e^(sigma I Pi (a + n + b eps)) PER t-COLUMN, not per sector).

   System (su04 class): f1 = u^(-1+eps), f2 = -u^eps/eps with u = x
   (lower pin) or u = 1-x (reflected upper pin); J = eps*f is regular.
   Integral: Int_0^1 [u*f1 + f2] dx = (1 - 1/eps)/(1 + eps)
           = -1/eps + 2 - 2 eps + 2 eps^2 - ... *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

catchDE2[DiffExp2`Config`LoadConfiguration[{"WorkingPrecision" -> 100,
  "ExpansionOrder" -> 40, "EpsilonOrder" -> 6, "DivisionOrder" -> 4}]];
eps = Global`eps; x = Global`x;
anchor = 1/2;
want = <|-1 -> -1, 0 -> 2, 1 -> -2, 2 -> 2, 3 -> -2|>;

checkPin[label_, w_] := Module[{ok = True},
  If[FailureQ[w], assert[label, False],
    Do[Module[{got = DiffExp2`EpsSeries`ESCoefficient[w, k]},
      If[!(NumericQ[N[got]] && Abs[N[got - want[k], 20]] < 10^-10), ok = False]],
      {k, -1, 3}];
    assert[label, ok]]];

(* PIN B: singularity at the LOWER endpoint (positive-arm integration) *)
AB = (DiagonalMatrix[{-1 + eps, 2 eps}] + x*{{0, 0}, {1, 0}})/x;
sysB = catchDE2[DiffExp2`API`LoadSystem[<|"Matrix" -> AB, "Variable" -> x|>]];
JB1 = Normal[Series[eps*anchor^(-1 + eps), {eps, 0, 8}]];
JB2 = Normal[Series[-anchor^eps, {eps, 0, 8}]];
bvB = {Table[Coefficient[JB1, eps, k], {k, 0, 8}],
       Table[Coefficient[JB2, eps, k], {k, 0, 8}]};
wB = catchDE2[DiffExp2`API`LineIntegral[sysB, bvB, anchor, {0, 1},
  {x/eps, 1/eps}]];
checkPin["pinB_lower_endpoint_laurent_integral", wB];

(* PIN C: singularity at the UPPER endpoint (negative-arm + crossing) *)
AC = -(DiagonalMatrix[{-1 + eps, 2 eps}] + (1 - x)*{{0, 0}, {1, 0}})/(1 - x);
sysC = catchDE2[DiffExp2`API`LoadSystem[<|"Matrix" -> AC, "Variable" -> x|>]];
JC1 = Normal[Series[eps*(1 - anchor)^(-1 + eps), {eps, 0, 8}]];
JC2 = Normal[Series[-(1 - anchor)^eps, {eps, 0, 8}]];
bvC = {Table[Coefficient[JC1, eps, k], {k, 0, 8}],
       Table[Coefficient[JC2, eps, k], {k, 0, 8}]};
wC = catchDE2[DiffExp2`API`LineIntegral[sysC, bvC, anchor, {0, 1},
  {(1 - x)/eps, 1/eps}]];
checkPin["pinC_upper_endpoint_crossing_integral", wC];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
