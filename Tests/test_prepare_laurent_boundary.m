(* Focused end-to-end contract for Laurent-valued regular-anchor data.
   The epsilon-dependent scalar system proves that transport budgets the
   private positive headroom consumed by the negative boundary power. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
Get[FileNameJoin[{repoRoot, "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

DiffExp2`LoadConfiguration[
  "RecurrenceBackend" -> "Wolfram",
  "WorkingPrecision" -> 60,
  "ExpansionOrder" -> 20,
  "EpsilonOrder" -> 2,
  "Verbosity" -> 0
];

system = DiffExp2`LoadSystem[<|
  "Matrix" -> {{Global`eps}},
  "Variable" -> Global`x
|>];
boundary = DiffExp2`PrepareLaurentBoundary[
  {1/Global`eps + 2},
  0,
  "EpsilonOrder" -> 2
];
result = DiffExp2`TransportEndpoint[system, boundary, 0, 1/2];
value = First[DiffExp2`EpsilonExpression[result, Global`eps]];
expected = Normal@Series[
  (1/Global`eps + 2) Exp[Global`eps/2],
  {Global`eps, 0, 2}
];

assert["transport retains the Laurent lower edge and requested complete top",
  DiffExp2`EpsilonWindow[result] ===
    <|"Min" -> -1, "CompleteMax" -> 2|>];
assert["transported Laurent point datum matches the exact scalar solution",
  Max[Abs@N[
    Coefficient[Expand[value - expected], Global`eps, #],
    50] & /@ Range[-1, 2]] < 10^-45];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Exit[1], Exit[0]];
