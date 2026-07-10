(* Export the physical all-unit level-0 integral directly from FTExamples.
   Unlike export_pysecdec_family_specs.m this never starts FIRE or builds the
   recursive FT ladder, so a numerical oracle can be minted before an
   expensive new topology is reduced. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["FeynmanTrick/FeynmanTrick.m"], {General::shdw, Symbol::shdw}];
Get[FileNameJoin[{repoRoot, "Scripts", "FTExamples.m"}]];

envOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]];

pythonString[expr_] := StringReplace[ToString[InputForm[expr]], {
  "^" -> "**", " " -> "", "FeynmanTrick`FTeps" -> "eps",
  "FTeps" -> "eps", "Global`" -> "", "Gamma[" -> "gamma(",
  "[" -> "(", "]" -> ")"}];

unitSimplexMap[n_Integer] := Module[{params, rescaled, jac},
  If[n <= 1, Return[{{1}, {}, 1}, Module]];
  params = Table[ToExpression["Global`z" <> ToString[i]], {i, n - 1}];
  rescaled = Table[Which[
    j === 1, Product[params[[i]], {i, 1, n - 1}],
    j <= n - 1, (1 - params[[j - 1]])*
      Product[params[[i]], {i, j, n - 1}],
    True, 1 - params[[n - 1]]], {j, n}];
  jac = FullSimplify[Abs[Det[D[rescaled[[1 ;; n - 1]], {params}]]],
    And @@ Thread[0 < params < 1]];
  {rescaled, params, jac}];

example = envOrDefault["FT_EXAMPLE", "kite"];
outFile = envOrDefault["PYSECDEC_SPEC_FILE",
  FileNameJoin[{$TemporaryDirectory, "pysecdec_" <> example <> "_l0.json"}]];
orders = Quiet[Check[
  ToExpression["{" <> envOrDefault["PYSECDEC_ORDERS", "0"] <> "}"],
  $Failed]];
If[!ListQ[orders] || !AllTrue[orders, IntegerQ],
  Print["PYSECDEC_ORDERS must be a comma-separated list of integers"];
  Exit[2]];

spec = FTExampleSpec[example];
If[!AssociationQ[spec], Exit[3]];
FeynmanTrick`SetFTOption["DimensionExpression", spec["Dimension"]];
props = spec["Propagators"];
loops = spec["LoopMomenta"];
repls = spec["Replacements"];
uf = FeynmanTrick`BoundaryConditions`ComputeSymanzikPolynomials[
  props, loops, repls];
If[uf === $Failed, Exit[4]];
{U, F, fvars} = uf;
n = Length[props];
loopCount = Length[loops];
v = n;
dExpr = spec["Dimension"];
gammaArg = Expand[v - loopCount*dExpr/2];
uPower = Expand[v - (loopCount + 1)*dExpr/2];
fPower = Expand[-gammaArg];
simplex = unitSimplexMap[n];
{rescaled, params, jac} = simplex;
rules = Thread[fvars -> rescaled];
Usub = FullSimplify[Together[U /. rules], And @@ Thread[0 < params < 1]];
Fsub = FullSimplify[Together[F /. rules], And @@ Thread[0 < params < 1]];

record = <|
  "Name" -> example <> "_L0_" <> StringRiffle[ConstantArray["1", n], "x"],
  "Example" -> example, "Level" -> 0, "Master" -> ConstantArray[1, n],
  "ActivePositions" -> Range[n], "ActivePowers" -> ConstantArray[1, n],
  "LoopCount" -> loopCount, "Variables" -> (SymbolName /@ params),
  "U" -> pythonString[Usub], "F" -> pythonString[Fsub],
  "UPower" -> pythonString[uPower], "FPower" -> pythonString[fPower],
  "Remainder" -> pythonString[jac],
  "Prefactor" -> pythonString[Gamma[gammaArg]],
  "RequestedOrders" -> orders|>;

Export[outFile, {record}, "RawJSON"];
Print["SpecFile=", outFile];
Print["Name=", record["Name"]];
Print["U=", record["U"]];
Print["F=", record["F"]];
