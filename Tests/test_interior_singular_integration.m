(* Tests for definite integration over intervals that contain the local
   singular point (interior IBP poles, e.g. the box family at xx1 = 1/4):
   - regular powers of a meromorphic term must NOT be regulated (no
     artificial higher-epsilon tails);
   - the residue integrates to its principal value, real, with no tails;
   - an interval straddling zero must produce a REAL result (no complex
     Log leakage from the negative arm). *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];
Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];

passed = 0; failed = 0;
test[name_String, ok_] := If[TrueQ[ok],
  passed++; Print["  PASS: ", name],
  failed++; Print["  FAIL: ", name]
];

x = DiffExp`Symbols`x;
wp = 80;
DiffExp`State`DiffExpConfiguration["WorkingPrecision"] = wp;

A = SetPrecision[-1/5, wp];
B = SetPrecision[1/4, wp];

near[a_, b_, tol_:10^-20] := TrueQ[Abs[a - b] < tol];
imOf[z_] := Abs[Im[N[z, 30]]];

(* --- 1. pure regular content under x^-1: g = x over {A, B} --- *)
g0 = SeriesData[x, 0, {SetPrecision[0, wp], SetPrecision[1, wp]}, 0, 6, 1];
res = DiffExp`RegularizedIntegration`IntegrateSingularTermLaurent[
  -1, 0, -1, {g0, 0, 0}, {A, B}];
cs = res["Coefficients"]; mp = res["MinPower"];
val[p_] := Module[{i = p - mp + 1}, If[1 <= i <= Length[cs], cs[[i]], 0]];
test["regular part eps^-1 = B - A", near[val[-1], B - A]];
test["regular part eps^0 vanishes", near[val[0], 0]];
test["regular part eps^1 vanishes", near[val[1], 0]];
test["regular part is real", imOf[val[-1]] < 10^-25];

(* --- 2. residue + regular: g = 1 + x  ->  PV log(B/|A|) + (B - A) --- *)
g1 = SeriesData[x, 0, {SetPrecision[1, wp], SetPrecision[1, wp]}, 0, 6, 1];
res2 = DiffExp`RegularizedIntegration`IntegrateSingularTermLaurent[
  -1, 0, -1, {g1, 0, 0}, {A, B}];
cs2 = res2["Coefficients"]; mp2 = res2["MinPower"];
val2[p_] := Module[{i = p - mp2 + 1}, If[1 <= i <= Length[cs2], cs2[[i]], 0]];
test["PV residue eps^-1 = log(B/|A|) + (B - A)",
  near[val2[-1], Log[B/Abs[A]] + (B - A), 10^-18]];
test["PV residue eps^0 vanishes", near[val2[0], 0, 10^-18]];
test["PV residue eps^1 vanishes", near[val2[1], 0, 10^-18]];
test["PV residue is real", imOf[val2[-1]] < 10^-25];

(* --- 3. double pole + regular tail: x^-2 (c0 + c2 x^2) --- *)
(* finite part: c0 (-1/B - (-1/A)) + c2 (B - A); no tails *)
g2 = SeriesData[x, 0,
  {SetPrecision[1, wp], SetPrecision[0, wp], SetPrecision[1, wp]}, 0, 6, 1];
res3 = DiffExp`RegularizedIntegration`IntegrateSingularTermLaurent[
  -2, 0, -1, {g2, 0, 0}, {A, B}];
cs3 = res3["Coefficients"]; mp3 = res3["MinPower"];
val3[p_] := Module[{i = p - mp3 + 1}, If[1 <= i <= Length[cs3], cs3[[i]], 0]];
test["double pole finite part",
  near[val3[-1], (-1/B + 1/A) + (B - A), 10^-18]];
test["double pole eps^0 vanishes", near[val3[0], 0, 10^-18]];
test["double pole is real", imOf[val3[-1]] < 10^-25];

Print["\nResults: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
