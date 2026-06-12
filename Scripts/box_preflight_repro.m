(* Box preflight, phase C (DiffExp2 only): runtime confirmation of the two
   findings the static audit flags, on the vendored fixtures with generic
   boundary values (the geometry/crossing code paths do not depend on the
   physical boundary data):
   1. L2 [anchor -> 1] crosses the APPARENT singular chart at xx2 = 7/11
      whose local solution carries a t^(-1+eps) sector (b = 1): expected
      runtime E8 "crossing a multivalued singular chart without a derivable
      Im-sign" at the first far-side chart — the gate is syntactic on
      sector tags, so it fires regardless of the (physically ~zero) weight.
   2. L1 [anchor -> 0] crosses the REGULAR-indicial IBP chart at 1/4
      (sectors all (0,0,0)): expected to pass (sigma fallback 1, trivial
      phase) — plus per-chart timing for the runtime estimate.
   Also times L3 (d=1) and the L2 lo line for the estimate table. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

anchor = 11/23;
benchDir = FileNameJoin[{repoRoot, "Tests", "refs", "bench"}];
catch2[expr_] := Catch[expr, "DiffExp2Error"];
SetAttributes[catch2, HoldFirst];

runOne[level_Integer, target_] := Module[{fix, var, sys, sys2, bvals, kmax,
    t0, res, wp},
  fix = Get[FileNameJoin[{benchDir, "box_L" <> ToString[level] <> ".m"}]];
  var = fix["Variable"];
  wp = fix["FTSettings"]["WorkingPrecision"];
  catch2[DiffExp2`Config`LoadConfiguration[{
    "WorkingPrecision" -> wp,
    "ExpansionOrder" -> fix["FTSettings"]["ExpansionOrder"],
    "EpsilonOrder" -> fix["FTSettings"]["EpsilonOrder"],
    "DivisionOrder" -> fix["FTSettings"]["DivisionOrder"],
    "StepDivisionOrder" -> fix["FTSettings"]["StepDivisionOrder"],
    "Variables" -> {}, "Verbosity" -> 0}]];
  sys = catch2[DiffExp2`API`LoadSystem[<|"Matrix" -> fix["Matrix"],
    "Variable" -> var|>]];
  kmax = fix["FTSettings"]["EpsilonOrder"];
  (* generic real boundary values, d x (kmax+1), at working precision *)
  bvals = Table[N[1 + c/7 + (k - 1)/13, wp + 20],
    {c, Length[fix["Matrix"]]}, {k, kmax + 1}];
  t0 = SessionTime[];
  res = catch2[DiffExp2`API`TransportEndpoint[sys, bvals, anchor, target,
    "ExtraSingularFactors" -> fix["ExtraSingularFactors"]]];
  Print["TRANSPORT L", level, " -> ", target, ": ",
    If[FailureQ[res],
      (* the DE2Error one-line summary with the full payload is already in
         the log right above this line *)
      "FAIL (see DiffExp2 error line above)",
      "OK EndpointIsSingular=" <> ToString[res["EndpointIsSingular"]] <>
        " ncharts=" <> ToString[Length[res["Charts"]]]],
    "  dt=", N[SessionTime[] - t0, 4], "s"];
  res];

(* the blocker candidate first *)
r21 = runOne[2, 1];
(* the regular-indicial crossing + timing *)
r10 = runOne[1, 0];
r11 = runOne[1, 1];
r20 = runOne[2, 0];
r30 = runOne[3, 0];
r31 = runOne[3, 1];
Print["DONE t=", SessionTime[]];
Quit[0];
