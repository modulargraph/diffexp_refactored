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
  "IndicialData" -> <|"Reduction" -> <|
    "GaugeInverseCertified" -> True,
    "GaugeInverseCertificateSchema" ->
      "diffexp2-indicial-exact-gauge-inverse-v1"|>|>,
  "V" -> v, "VInv" -> Map[Cancel[Together[#]] &, Inverse[v], {2}]|>;

goodV = {{1, 1}, {eps, 1 + eps}};
laurentV = {{1, 1}, {eps, -eps}};
tDependentV = {{1, t}, {0, 1}};

goodChart = makeChart[goodV];
good = DiffExp2`Solve`Private`sccSpectralFrameCertificate[goodChart];
laurentChart = makeChart[laurentV];
laurent = DiffExp2`Solve`Private`sccSpectralFrameCertificate[laurentChart];
tDependent = DiffExp2`Solve`Private`sccSpectralFrameCertificate[
  makeChart[tDependentV]];

shape = <|
  "Center" -> 0, "ChartMap" -> <|"Center" -> 0, "Scale" -> 1|>,
  "Radius" -> 2, "Prescriptions" -> {},
  "EpsWindow" -> <|"Min" -> -2, "CompleteMax" -> 4|>,
  "TWindow" -> <|"CompleteMax" -> 3|>, "Dimension" -> 2|>;
transform = DiffExp2`Solve`Private`PrepareSCCSpectralSourceTransform[
  goodChart, shape, <|"domain" -> "rational", "symbols" -> {}|>];
laurentTransform = DiffExp2`Solve`Private`PrepareSCCSpectralSourceTransform[
  laurentChart, shape, <|"domain" -> "rational", "symbols" -> {}|>];

gauge = {{1, (1 + eps)/t}, {t, 2 + eps}};
gaugeChart = Join[goodChart, <|"Gauge" -> gauge,
  "GaugeInverse" -> Map[Cancel[Together[#]] &,
    Inverse[gauge], {2}]|>];
gaugeCertificate = DiffExp2`Solve`Private`sccGaugeFrameCertificate[
  gaugeChart];
toPhysical = DiffExp2`Solve`Private`PrepareSCCGaugeTransform[
  gaugeChart, "to_physical", shape,
  <|"domain" -> "rational", "symbols" -> {}|>];
toReduced = DiffExp2`Solve`Private`PrepareSCCGaugeTransform[
  gaugeChart, "to_reduced", shape,
  <|"domain" -> "rational", "symbols" -> {}|>];
gaugeEntries = Cases[Flatten[{gauge,
    gaugeChart["GaugeInverse"]}], value_ /; !TrueQ[PossibleZeroQ[value]]];
gaugePrepParity = AllTrue[gaugeEntries, Function[entry, Module[
    {identity = DiffExp2`SectorSeries`ExactExpressionIdentity[entry, t],
     dedicated, generic},
    dedicated = DiffExp2`Solve`Private`PrepareSCCGaugeMultiplier[
      shape, entry, t, identity];
    generic = DiffExp2`SectorSeries`PrepareRationalMultiplier[
      shape, entry, t];
    KeyTake[dedicated, {"EpsilonShift", "CenterPoleOrder",
        "TaylorKernels", "ExactIdentity", "ProvenZero"}] ===
      KeyTake[generic, {"EpsilonShift", "CenterPoleOrder",
        "TaylorKernels", "ExactIdentity", "ProvenZero"}]]]];

entries = Lookup[transform, "entries", {}];
entryShifts = Lookup[Lookup[entries, "multiplier", <||>],
  "epsilon_shift", None];

ok = TrueQ[good["admissible"]] && !TrueQ[good["identity_v"]] &&
  good["det_epsilon_valuation"] === 0 &&
  TrueQ[laurent["admissible"]] &&
  laurent["det_epsilon_valuation"] === 1 &&
  !TrueQ[tDependent["admissible"]] &&
  StringContainsQ[tDependent["detail"], "chart variable"] &&
  AssociationQ[transform] &&
  transform["schema"] ===
    "diffexp2-scc-spectral-source-transform-v1" &&
  transform["rows"] === 2 && transform["columns"] === 2 &&
  !TrueQ[transform["identity"]] &&
  TrueQ[transform["epsilon_unimodular"]] &&
  transform["det_epsilon_valuation"] === 0 &&
  laurentTransform["det_epsilon_valuation"] === 1 &&
  StringContainsQ[laurentTransform["exact_identity"],
    "\"det_epsilon_valuation\":1"] &&
  transform["v_exact_identity"] === good["v_exact_identity"] &&
  transform["vinv_exact_identity"] === good["vinv_exact_identity"] &&
  TrueQ[gaugeCertificate["admissible"]] &&
  !TrueQ[gaugeCertificate["identity_gauge"]] &&
  toPhysical["role"] === "to_physical" &&
  toReduced["role"] === "to_reduced" &&
  Max[Lookup[Lookup[toPhysical["entries"], "multiplier", <||>],
      "center_pole_order", 0]] === 1 &&
  StringContainsQ[toPhysical["exact_identity"],
    "diffexp2-scc-gauge-transform-identity-v1"] &&
  StringContainsQ[toReduced["exact_identity"],
    "\"role\":\"to_reduced\""] &&
  TrueQ[gaugePrepParity] &&
  AllTrue[Join[toPhysical["entries"], toReduced["entries"]],
    !KeyExistsQ[#["multiplier"], "kernels"] &&
      Length[#["multiplier", "analytic_coefficients"]] === 7 &] &&
  Length[entries] === 4 && Sort[entryShifts] === {0, 0, 0, 1} &&
  AllTrue[entries,
    #["multiplier", "center_pole_order"] === 0 &&
      !KeyExistsQ[#["multiplier"], "kernels"] &&
      Length[#["multiplier", "analytic_coefficients"]] === 7 &];

If[ok,
  Print["PASS: exact Laurent-unimodular SCC spectral-frame certificate"],
  Print["FAIL: exact Laurent-unimodular SCC spectral-frame certificate"];
  Print["good=", good];
  Print["laurent=", laurent];
  Print["tDependent=", tDependent];
  Print["transform=", transform];
  Exit[1]];
