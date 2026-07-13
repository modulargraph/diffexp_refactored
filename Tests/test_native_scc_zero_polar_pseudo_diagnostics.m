(* A CASE-P collision may have an identically zero exact polar gamma.  The
   Rational compensator then certifies the value but emits no compensation
   term.  The Wolfram summary validator must admit precisely that certified
   zero-polar shape without admitting unresolved or uncertified hits. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

record = <|
  "block" -> 2,
  "pseudo_hit_count" -> 1,
  "pseudo_compensation_count" -> 0,
  "max_pseudo_depth" -> 1,
  "pseudo_value_certified" -> True,
  "uncompensated_pseudo_hit_count" -> 0|>;

valid = DiffExp2`Solve`Private`sccNativePseudoDiagnosticsQ[
  <|"block_diagnostics" -> {record}|>];
unresolved = DiffExp2`Solve`Private`sccNativePseudoDiagnosticsQ[
  <|"block_diagnostics" -> {
    ReplacePart[record, "uncompensated_pseudo_hit_count" -> 1]}|>];
uncertified = DiffExp2`Solve`Private`sccNativePseudoDiagnosticsQ[
  <|"block_diagnostics" -> {
    ReplacePart[record, "pseudo_value_certified" -> False]}|>];
uninspected = DiffExp2`Solve`Private`sccNativePseudoDiagnosticsQ[
  <|"block_diagnostics" -> {
    ReplacePart[record, "max_pseudo_depth" -> 0]}|>];

ok = TrueQ[valid] && !TrueQ[unresolved] && !TrueQ[uncertified] &&
  !TrueQ[uninspected];

Print[If[ok, "PASS", "FAIL"],
  ": certified zero-polar CASE-P diagnostics contract"];
Exit[If[ok, 0, 1]];
