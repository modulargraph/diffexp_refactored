root = DirectoryName[DirectoryName[$InputFileName]];
SetDirectory[root];
Get["DiffExp2.m"];

eps = DiffExp2`Config`CanonicalEps[];
t = Unique["spectralFrameT$"];

makeChart[v_] := <|
  "ID" -> "spectral-frame-certificate-test",
  "SystemSize" -> 2, "ChartVar" -> t,
  "Center" -> 0, "ChartMap" -> <|"Center" -> 0, "Scale" -> 1|>,
  "Radius" -> 2, "Prescriptions" -> {},
  "V" -> v, "VInv" -> Map[Cancel[Together[#]] &, Inverse[v], {2}]|>;

goodV = {{1, 1}, {eps, 1 + eps}};
badV = {{1, 1}, {eps, -eps}};
tDependentV = {{1, t}, {0, 1}};

goodChart = makeChart[goodV];
good = DiffExp2`Solve`Private`sccSpectralFrameCertificate[goodChart];
bad = DiffExp2`Solve`Private`sccSpectralFrameCertificate[makeChart[badV]];
tDependent = DiffExp2`Solve`Private`sccSpectralFrameCertificate[
  makeChart[tDependentV]];

shape = <|
  "Center" -> 0, "ChartMap" -> <|"Center" -> 0, "Scale" -> 1|>,
  "Radius" -> 2, "Prescriptions" -> {},
  "EpsWindow" -> <|"Min" -> -2, "CompleteMax" -> 4|>,
  "TWindow" -> <|"CompleteMax" -> 3|>, "Dimension" -> 2|>;
transform = DiffExp2`Solve`Private`PrepareSCCSpectralSourceTransform[
  goodChart, shape, <|"domain" -> "rational", "symbols" -> {}|>];

entries = Lookup[transform, "entries", {}];
entryShifts = Lookup[Lookup[entries, "multiplier", <||>],
  "epsilon_shift", None];

ok = TrueQ[good["admissible"]] && !TrueQ[good["identity_v"]] &&
  good["det_epsilon_valuation"] === 0 &&
  !TrueQ[bad["admissible"]] && bad["det_epsilon_valuation"] === 1 &&
  StringContainsQ[bad["detail"], "val_eps(det(V)) == 0"] &&
  !TrueQ[tDependent["admissible"]] &&
  StringContainsQ[tDependent["detail"], "chart variable"] &&
  AssociationQ[transform] &&
  transform["schema"] ===
    "diffexp2-scc-spectral-source-transform-v1" &&
  transform["rows"] === 2 && transform["columns"] === 2 &&
  !TrueQ[transform["identity"]] &&
  TrueQ[transform["epsilon_unimodular"]] &&
  transform["det_epsilon_valuation"] === 0 &&
  transform["v_exact_identity"] === good["v_exact_identity"] &&
  transform["vinv_exact_identity"] === good["vinv_exact_identity"] &&
  Length[entries] === 4 && Sort[entryShifts] === {0, 0, 0, 1} &&
  AllTrue[entries,
    #["multiplier", "center_pole_order"] === 0 &&
      Length[#["multiplier", "kernels"]] === 7 &&
      AllTrue[#["multiplier", "kernels"], Length[#] === 4 &] &];

If[ok,
  Print["PASS: exact epsilon-unimodular SCC spectral-frame certificate"],
  Print["FAIL: exact epsilon-unimodular SCC spectral-frame certificate"];
  Print["good=", good];
  Print["bad=", bad];
  Print["tDependent=", tDependent];
  Print["transform=", transform];
  Exit[1]];
