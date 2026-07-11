(* Planner-only profile for the exact 15-master unequal-mass banana line.
   This reconstructs the same Euclidean mass deformation as
   bench_unequal_banana_cpp.m, then stops immediately after SegmentLine. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

wp = ToExpression[Replace[Environment["UB_WP"],
  Except[_String?(StringLength[StringTrim[#]] > 0 &)] -> "250"]];
matrixDir = FileNameJoin[{repoRoot, "Tests", "Banana_Matrices"}];
eps = Global`eps; x = Global`x;
psq = Global`psq; mm1 = Global`mm1; mm2 = Global`mm2;
mm3 = Global`mm3; mm4 = Global`mm4;
partialMatrix[var_Symbol, order_Integer] := Get[FileNameJoin[{matrixDir,
  "d" <> SymbolName[var] <> "_" <> ToString[order] <> ".m"}]];

buildTime = First@AbsoluteTiming[
  dmm1 = partialMatrix[mm1, 0] + eps partialMatrix[mm1, 1];
  dmm2 = partialMatrix[mm2, 0] + eps partialMatrix[mm2, 1];
  dmm3 = partialMatrix[mm3, 0] + eps partialMatrix[mm3, 1];
  massMatrix = Map[Cancel[Together[#]] &,
    (dmm1 + dmm2/2 + dmm3/3) /. {
      psq -> -1, mm1 -> 1 + x, mm2 -> 1 + x/2,
      mm3 -> 1 + x/3, mm4 -> 1}, {2}]];
Clear[dmm1, dmm2, dmm3];
Print["DE2PLAN MatrixBuild ", N[buildTime, 8]];

DiffExp2`Config`LoadConfiguration[{
  "WorkingPrecision" -> wp, "ChopPrecision" -> 25,
  "ExpansionOrder" -> 50, "EpsilonOrder" -> 4,
  "DivisionOrder" -> 3, "RadiusOfConvergence" -> 10}];
sys = DiffExp2`API`LoadSystem[
  <|"Matrix" -> massMatrix, "Variable" -> x|>];
Print["DE2PLAN Factors ", Length[sys["SingularFactors"]], " ",
  ({Exponent[#, x], LeafCount[#]} & /@ sys["SingularFactors"])];
If[Environment["DE2_PLANNER_DUMP_FACTORS"] === "1",
  Print["DE2PLAN FactorExpressions ",
    ToString[sys["SingularFactors"], InputForm]]];
planTime = First@AbsoluteTiming[
  plan = DiffExp2`Transport`SegmentLine[sys, {0, 1}]];
Print["DE2PLAN Total ", N[planTime, 8], " Charts ", plan["SegmentCount"],
  " Real ", Length[plan["Singularities", "Real"]], " Projected ",
  Length[plan["Singularities", "ProjectionWaypoints"]]];
