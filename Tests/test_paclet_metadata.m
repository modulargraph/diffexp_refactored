(* Release metadata and source-tree paclet discovery smoke test. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

paclet = Get[FileNameJoin[{repoRoot, "PacletInfo.wl"}]];
citation = Import[FileNameJoin[{repoRoot, "CITATION.cff"}], "YAML"];
assert["paclet metadata object", Head[paclet] === PacletObject];
assert["paclet name and compatibility target",
  paclet["Name"] === "DiffExp2" && paclet["WolframVersion"] === "15.0+"];
assert["citation and paclet versions agree",
  AssociationQ[citation] && citation["version"] === paclet["Version"] &&
    citation["license"] === "GPL-3.0-or-later"];
assert["paclet and public API versions agree",
  Module[{loaded},
    Get[FileNameJoin[{repoRoot, "DiffExp2.m"}]];
    loaded = DiffExp2`$DiffExp2Version;
    loaded === paclet["Version"]]];
assert["paclet and Feynman-trick facade versions agree",
  Module[{loaded},
    Get[FileNameJoin[{repoRoot, "FeynmanTrick.m"}]];
    loaded = FeynmanTrick`$FeynmanTrickVersion;
    loaded === paclet["Version"]]];
assert["both root contexts are declared",
  ContainsAll[
    "Context" /. Rest[First[paclet["Extensions"]]],
    {"DiffExp2`", "FeynmanTrick`"}]];

PacletDirectoryLoad[repoRoot];
found = Select[PacletFind["DiffExp2"],
  ExpandFileName[# ["Location"]] === ExpandFileName[repoRoot] &];
assert["source tree is discoverable as a paclet", Length[found] === 1];
PacletDirectoryUnload[repoRoot];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Exit[1], Exit[0]];
