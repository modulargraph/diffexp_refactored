(* Compute a literature family from Mathematica with the DiffExp 2.2 wrapper.
   The wrapper calls the same C++ executable through RunProcess and JSON.

   Build DiffExp first (build/diffexp), then evaluate this file. *)

here = DirectoryName[$InputFileName];
root = FileNameJoin[{here, "..", ".."}];                (* this checkout *)
Get[FileNameJoin[{root, "DiffExp.m"}]];
$DiffExpExecutable = FileNameJoin[{root, "build", "diffexp"}];

(* A family is an ordinary Association: the same keys as the JSON files. *)
family = Import[FileNameJoin[{here, "families", "hvp_chpt_e5e6_d2.json"}], "RawJSON"];

(* Change the kinematics or requested integrals directly, e.g. t = p^2/M^2 = -2: *)
(* family["external_gram"] = {{"-2"}}; *)

result = DiffExpFeynmanTrick[family, {
   "--ibp-provider", "ibp-solver",
   "--ibp-threads", ToString[$ProcessorCount],
   "--cache", FileNameJoin[{$HomeDirectory, ".cache", "diffexp-literature-2026"}]}];

result["coefficients"]  (* per integral: epsilon orders epsilon_low..epsilon_high *)
result["timings"]
