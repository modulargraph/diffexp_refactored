# Hygiene / refactor plan — diffexp_refactored support files

Survey scope: Scripts/, Tests/, Docs/, top-level loaders, AGENTS.md, README.md, Reference/, repo level. Static analysis only; nothing was executed or modified. Line numbers verified against current master (3b0d3c9).

---

## P0 — Test-battery integrity (cheap, high payoff, do first)

**P0.1 — bash/zsh array footgun** — `/Users/mhidding/Code/diffexp_refactored/Scripts/run_test_battery.sh:27`
`for t in $tests` iterates the whole array only under zsh; under `bash Scripts/run_test_battery.sh` it expands to the *first element only*, so the battery silently runs just `test_package_loading.m` and reports success.
**Action:** change line 27 to `for t in "${tests[@]}"; do` (valid in zsh and macOS bash 3.2; keeps the `#!/bin/zsh` shebang as belt-and-suspenders). Optionally add an early guard: `[ -n "$BASH_VERSION" ] && [ -z "$ZSH_VERSION" ]` is no longer needed once the quoting is fixed. **Risk: none** (pure expansion fix; zsh behavior identical).

**P0.2 — two battery tests can fail invisibly** — `Scripts/run_test_battery.sh:33` detects failure only via exit code or the literal `"Some tests FAILED"`. But:
- `Tests/test_symbols_namespacing.m:88` prints `"All tests: SOME FAILED"` and never calls `Exit`/`Quit` → wolframscript exits 0, grep misses it → battery reports PASS on failure.
- `Tests/test_resonant_2f1.m:122-131` prints `"TEST FAILED"` and never exits nonzero → same silent pass.
**Action:** standardize all 17 battery tests on the convention already used by the newest tests (`Tests/test_interior_singular_integration.m:66-68`): `Print["Some tests FAILED."]; Exit[1]` / `Exit[0]`. Minimum fix: add `Exit[1]`-on-failure tails to those two files. Defense in depth: broaden line 33's grep to `grep -Eq "Some tests FAILED|SOME FAILED|TEST FAILED"`. **Risk: none.**

**P0.3 — kernel path not overridable** — `Scripts/run_test_battery.sh:3` unconditionally overwrites `WolframKernel`.
**Action:** `export WolframKernel="${WolframKernel:-/Applications/Wolfram Engine.app/.../WolframKernel}"`. Keep the ONE-kernel "never parallelize" comment (line 2) — it is load-bearing. **Risk: none.**

---

## P1 — New tests for the campaign-revealed gaps (concrete sketches)

All four follow the existing harness pattern (`SetDirectory[repoRoot]`, `Quiet[Get["DiffExp.m"],{General::shdw,Symbol::shdw}]`, pass/fail counters, `Exit[0|1]`), and should be appended to the battery list in `Scripts/run_test_battery.sh`.

**P1.1 — `Tests/test_apparent_singularity_transport.m` (interior apparent singular point on a transport line).** Today nothing exercises this: the box L2 crossing at xx2=7/11 was the first ever and broke 5 ways (`Docs/FeynmanTrickBoxFamilyStatus.md:75-108`); the closest tests (`test_singular_recurrence.m`, `test_interior_singular_integration.m`) cover singular *endpoints* and interior *integration* poles, not an interior singular chart during transport.

Synthetic 2×2 with exact solution, apparent singularity at u=0, eps-sourced at every order:
- Fundamental solutions (by construction, verified analytically): `Y1 = {1 + eps*u^3, 3*eps*u^2}`, `Y2 = {u^2, 2*u}`. Wronskian `det F = u*(2 - eps*u^3)` has a simple zero at u=0 while every solution is analytic ⇒ apparent singularity for all eps.
- Matrix `A = F'.F^-1 = {{0, 1}, {6 eps u/(2-eps u^3), (2-4 eps u^3)/(u(2-eps u^3))}}`. Eps slices (write a new `Tests/ApparentSingularity_Matrices/` with `du_0.m … du_4.m`, same format as `Tests/Hypergeometric2F1_Matrices/dz_0.m`):
  - `du_0 = {{0,1},{0,1/u}}`
  - slice k≥1: `{{0,0},{3*u^(3k-2)/2^(k-1), -3*u^(3k-1)/2^k}}`, i.e. `du_1={{0,0},{3u, -3u^2/2}}`, `du_2={{0,0},{3u^4/2, -3u^5/4}}`, etc.
- Config as in `Tests/test_hypergeometric2f1.m:48-58` but `EpsilonOrder -> 4`, `UseMobius -> False`, WP 200.
- Boundary at u=-2/5, taking Y = Y1+Y2 (exact rationals): `BCs = {<|u -> -2/5|>, {{29/25, -8/125, 0, 0, 0}, {-4/5, 12/25, 0, 0, 0}}}`.
- `TransportTo[BCs, <|u -> 2/5|>, 1, True]` — the line crosses u=0 in the interior, forcing a chart centered on the apparent point (residue eigenvalues {0,1}).
- **Expected (exact closed form):** component 1 = `{29/25, +8/125, 0, 0, 0}`, component 2 = `{4/5, 12/25, 0, 0, 0}`; assert each order to `10^-30`. The eps≥2 orders are the poison detector: the matrix slices genuinely source them through the 1/u chart and they must cancel to zero by analyticity — precisely the mechanism that lost a constant at eps^2 in the box (`Docs/FeynmanTrickBoxFamilyStatus.md:88-97`).

**P1.2 — limitUpper/limitLower against a closed form.** `Tests/test_regularized_integration_edge_cases.m:609-657` unit-tests the inner `EvaluateEndpointLimitSectors`, but nothing tests the actual boundary path `EvaluateLimitFromTransport` (`FeynmanTrick/DiffExpIntegration.m:908`, dispatched at `:1596-1606`) end-to-end over a real transport result — the path whose collapsed-decomposition bug made the banana `{1,0,0,1}` boundary identically zero (`Docs/FeynmanTrickBananaStatus.md:339-348`).

Sketch (`Tests/test_limit_boundary_closed_form.m`): 2×2 *diagonal* system on (0,1) with masters `M1 = (1+x)(1-x)^(-2 eps)` and `M2 = (2-x) x^(3 eps)`, i.e. `A = diag[ 1/(1+x) + 2 eps/(1-x), -1/(2-x) + 3 eps/x ]`; write matrices in the `Level_*_Matrices` format `TransportLevel` consumes (mirror a dir produced by `test_feynmantrick_pipeline.m` for exact file naming). Boundary at x=1/2 is exact: `M1` tower `(3/2)(2 Log 2)^n/n!`, `M2` tower `(3/2)(-3 Log 2)^n/n!`, n=0..4. Run `TransportLevel`, then:
- `EvaluateLimitFromTransport[res, {1,1}, 1, 4, {0,0}, True]` (limitUpper): M1 is a pure b=-2 sector at x=1 ⇒ dropped wholesale; M2 is analytic at x=1 with value 1 ⇒ **expected `{1, 0, 0, 0, 0}`** exactly. The historical failure mode returns `{0,...}` (whole series dropped) — discriminating.
- `EvaluateLimitFromTransport[res, {1,1}, 0, 4, {0,0}, True]` (limitLower): M2 is a pure b=3 sector at x=0 ⇒ dropped; M1(0)=1 ⇒ **expected `{1, 0, 0, 0, 0}`**.
- Variant with x-dependent IBP coefficient `{1, 3-2x}` to catch coefficient-evaluation-at-endpoint bugs (limitUpper expected unchanged: `1*1 = 1`).

**P1.3 — transport slots are plain SeriesData (head-Plus type poisoning).** The class: inert one-past-the-end `SeriesCoefficient[...]` tails turn slot heads into `Plus`, silently disabling downstream `Head === SeriesData` probes (`Docs/FeynmanTrickBoxFamilyStatus.md:84-92`; normalization now at `DiffExp/IntegrationStrategies/Recurrence.m:558-585`). No test asserts the invariant.
Sketch: cheapest as extra asserts inside P1.1's test (it produces a saved transport crossing a singular chart — the worst case). For every `seg` in `result["SegmentData"]`: `series = DiffExp`RegularizedIntegration`Private`uncompressSeriesData[seg[[5]]]` (the pattern of `Scripts/eval_dump_generic.m` / `Scripts/check_transport_ode_residual.m:46`), then assert
- `AllTrue[Flatten[series, 1], MatchQ[#, _SeriesData] || TrueQ[PossibleZeroQ[#]] &]` (no `Plus` heads), and
- `FreeQ[series, HoldPattern[SeriesCoefficient[___]]]` (no inert tail requests survive into stored transport data).

**P1.4 — `ComputeSingularParticular` compound-source normalization unit test.** Zero test coverage: `ComputeSingularParticular` is referenced only in `DiffExp/IntegrationStrategies/Recurrence.m` (def `:550`, call `:852`) and the status doc.
Sketch (`Tests/test_singular_particular_sources.m`), calling `DiffExp`IntegrationStrategies`Private`ComputeSingularParticular` directly (tests already call privates, e.g. `test_integration_log_depth.m:34`), 1×1 system, `mode = "series"`, `eigenvalues = {1/3}` (non-resonant against s), `P = PInv = {{1}}`, `bHatCoeffs` from a regular part like `A_reg(x) = x`, minimal `ctx = <|"ExpansionOrder"->20, "WorkingPrecision"->80, "ChopPrecision"->40, "SystemSize"->1, "Label"->"test"|>`:
1. *Clean source*: `bClean = {SeriesData[x, 0, {3/7, 1, 2}, -1, 10, 1]}` (genuine x^-1 leading) → particular must lead at `x^0` (s = bLeadPow+1 = 0), head SeriesData.
2. *Poisoned source*: same series `+ SeriesCoefficient[SeriesData[x,0,{1},0,3,1], {x,0,5}]` (out-of-window inert tail; head becomes Plus — the upstream artifact verbatim). **Assert: result identical to (1) to working precision** — the normalization at `Recurrence.m:564-585` must zero the tail; pre-fix behavior fell back to leading power 0 and dropped the 3/7·x^-1 content.
3. *Noise-skip*: `{SeriesData[x, 0, {SetPrecision[10^-29,80], 3/7, 1, 2}, -2, 10, 1]}` (cancellation-residue junk below the true leading power, relative ~1e-29 ≪ the `10^-24` floor at `Recurrence.m:618`) → must produce the same particular as (1), not shift s onto -1.
4. *In-window inert tail* → must emit the "source entry is not a series" `PrintWarning` (`Recurrence.m:574-579`); capture and assert the message fires.
Correctness of the absolute value is carried by P1.1 end-to-end; this test pins the equality/normalization contract.

**Also (P1.5, dead test):** `Tests/test_transport_level.m:18` hardcodes a machine-/run-specific path `/private/var/folders/r9/.../FT_pipeline_test_32517/Level_3_Matrices` — it cannot run anywhere, has no pass/fail accounting, and is not in the battery. Delete it or rewrite it to generate its own matrices (its `TransportTo` step-through is the best starting scaffold for P1.2). **Risk: none** (already dead).

---

## P2 — Scripts/ dedup, headers, defaults

**P2.1 — shared prologue (one new file, ~40 lines saved, 4 files simplified).** Duplicated verbatim across the four .m scripts:
- `envOrDefault`: `Scripts/run_ft_stepwise.m:15-17`, `Scripts/export_pysecdec_family_specs.m:16-18`, `Scripts/check_transport_ode_residual.m:15-17`, `Scripts/eval_dump_generic.m:9-11`.
- repoRoot/SetDirectory/$Path/Quiet-Get boilerplate: `run_ft_stepwise.m:3-8`, `export_pysecdec_family_specs.m:3-8`, `check_transport_ode_residual.m:8-13`, `eval_dump_generic.m:3-7`.
- Identical FT option block (Threads/FThreads/Verbosity/ReductionCache/FIRETimeoutSeconds): `run_ft_stepwise.m:30-35` vs `export_pysecdec_family_specs.m:10-14`.
**Action:** add `Scripts/ScriptPrologue.m` exporting `envOrDefault`, `LoadRepoPackages[{"DiffExp","FeynmanTrick"}]`, `SetStandardFTOptions[]`; each script reduces to one `Get`. **Risk: low** — pure mechanical extraction; verify each script still runs once (the NormalizeLogPower incident shows cross-file symbol visibility is a real trap here, so keep the prologue context-free, plain Global` definitions).

**P2.2 — stale `/tmp` defaults (already env-overridable, but defaults point at dead artifacts):**
- `Scripts/eval_dump_generic.m:13` default `/tmp/diffexp_banana_l1_dumps/laurent_integral_0006.m` — the banana doc says the known-good dumps are in `/tmp/diffexp_banana_l1_dumps_v2final/` (`Docs/FeynmanTrickBananaStatus.md:206`). Make `DUMP_FILE` *required* (fail with a usage message) instead of defaulting to an ephemeral path.
- `Scripts/check_transport_ode_residual.m:19-21` defaults `TRANSPORT_FILE=/tmp/ft_transport_save/transport_level_1.m`, `SEG_INDEX=33` (banana-era segment count). Same treatment.
- `Scripts/pysecdec_family_driver.py:176` `--output-root` default `/tmp/pysecdec-ft-family` is a CLI flag — fine as is.
**Risk: none** (diagnostic tools only).

**P2.3 — missing usage headers (env-var contracts undocumented in-file):**
- `Scripts/run_ft_stepwise.m:1` — one-line header for a script driven by ~12 env vars (`FT_EXAMPLES`, `FT_WORKING_PRECISION`, `FT_EPS_ORDER`, `FT_EXPANSION_ORDER`, `FT_DIVISION_ORDER`, `FT_BOUNDARY_EXTRA_ORDER`, `FT_STOP_AFTER_BOUNDARY_LEVEL`, `FT_SAVE_TRANSPORT_DIR`, `FT_TRANSPORT_VERBOSITY`, `FT_POLE_ALLOWANCE`, `DEBUG_FUCHS_CHECK`, `DEBUG_BLOCK_RESID`). Add a usage comment block listing each with defaults.
- `Scripts/export_pysecdec_family_specs.m:1` — same for `FT_EXAMPLE`, `PYSECDEC_SPEC_FILE`, `INCLUDE_NEEDED`, `FT_FIXED_VALUE`.
- `Scripts/eval_dump_generic.m:1` — `DUMP_FILE`, `PRINT_SEGMENTS`.
- `Scripts/check_transport_ode_residual.m:1-6` — good prose header, but add the env list (`TRANSPORT_FILE`, `SEG_INDEX`, `CHECK_X_ORDER`, `MAX_EPS_CHECK`).
- Python scripts are fine (`compare_stepwise_log.py` docstring, `pysecdec_family_driver.py` argparse).
**Risk: none.**

**P2.4 — dead scripts: none.** All 8 Scripts/ files are alive and cross-referenced from the two status docs. Don't prune here.

---

## P3 — Docs staleness

**P3.1 — `Docs/FeynmanTrickBoxFamilyStatus.md` contradicts itself (highest doc priority — it's the "current" doc).** Commit 3b0d3c9 added the "RESOLVED: the eps^0 deficit" section (line 69, "Box L0 now reproduces the pin to 11 significant digits") but left behind:
- line 3: "Updated: 2026-06-10" (actual last update 2026-06-11);
- line 35: bisection table still says "L0 eps^0 | -47.9756 vs pin -41.2842 (deficit +6.69)";
- lines 150-152: "finite parts are NOT [trustworthy] until the endpoint-series defect is fixed. Each run prints loud warnings…" — both the endpoint-series defect (line 51) and the eps^0 deficit (line 69) are marked RESOLVED in the same file.
**Action:** refresh the table row to EXACT/11-digits, rewrite the practical-guidance bullet (box finite part now trusted to ~1e-11; remaining open item per latest commits is whatever survives of the seg12 follow-up list at lines 133-146). **Risk: none** (doc-only).

**P3.2 — `Docs/FeynmanTrickBananaStatus.md:462-469`:** "Remaining Notes / Next Steps" item 1 ("Fix `EvaluateLimitFromTransport`…") was completed and documented earlier in the same file (line 339, "Level-0 Fix: Sector-Aware Boundary Limits"). Delete or mark item 1 done. Rest of the doc is current (RESOLVED header is accurate). **Risk: none.**

**P3.3 — `FeynmanTrick/IMPLEMENTATION_STATUS.md` (last touched 2026-01-25) is badly stale:** "End-to-end pipeline (currently fails at boundary)", "Next steps: run pipeline test", a dead repo path `/Users/mhidding/Desktop/diffexp_refactored/Tests`, and `WOLFRAMSCRIPT_KERNELPATH` instead of the `WolframKernel` variable everything else uses. **Action:** either delete (the two Status docs + battery supersede it) or rewrite its "Status/Next Steps/How to Run" sections; its paper-equation summary (top half) is still useful. **Risk: none.**

**P3.4 — `Tests/SINGULARITY_DECOMPOSITION_PROBLEM.md` (2026-01-19):** a design spec whose closing TODO ("Implement … `DecomposeSingularity`") shipped months ago. Archive into Docs/ history or delete. **Risk: none.**

**P3.5 — `Docs/RegularizedIntegration.md` (2026-01-20):** not wrong (all six documented functions still exist in `DiffExp/RegularizedIntegration.m`) but it predates the entire campaign and omits the Laurent/sector machinery where everything now happens (`DefiniteIntegralWithPrefactorLaurent`, `IntegrateSingularTermLaurent`, `FitResidualEndpointSectors`, `EvaluateEndpointLimitSectors`). Add a "Laurent path" section or a banner pointing to the two Status docs. `Docs/recurrence_solver_prompt.md` is a one-off LLM prompt artifact — archive. The Jan-23 module references (CoreModules/Infrastructure/SeriesAndEvaluation/Transformations/SingularityDecomposition) and the June TransportAndStrategies.md spot-checked clean. **Risk: none.**

**P3.6 — README.md / AGENTS.md don't mention the battery or the example runner at all** (verified: zero hits for `Scripts|battery|run_ft` in both). `AGENTS.md:27-36` and `README.md:33-42` list only the six focused tests. **Action:** add to both: `Scripts/run_test_battery.sh` as the canonical gate (with the ONE-kernel / never-parallelize note, currently documented only in the script comment), and `Scripts/run_ft_stepwise.m` (`FT_EXAMPLES=bubble,sunrise,banana`) as the canonical end-to-end example run. **Risk: none.**

---

## P4 — Repo level

**P4.1 — `Reference/DiffExp_original.m`: keep.** Zero references anywhere in code (repo-wide grep: no `Get`/mention outside Reference/ itself). The `Reference/Examples/*.m` load `ParentDirectory[scriptDir]/DiffExp.m` (e.g. `Reference/Examples/Banana_example.m:10`) — i.e. `Reference/DiffExp.m`, which does not exist — so the whole Reference/ tree is archival-only and cannot run in place. It is the refactor baseline and cheap (120 KB); keep as-is, optionally add a one-line `Reference/README` note saying "archival; examples reference the original layout and do not run here."

**P4.2 — .gitignore** (`/Users/mhidding/Code/diffexp_refactored/.gitignore`):
- Vestigial entries: `kira/` (no such dir) and `FrontEnd/*` — harmless, optionally prune.
- Gap: `Reference/Examples/Tmp/` is a tracked output dir (`.gitkeep`) that `FivePointNonPlanar_example.m:88` writes into — add `Reference/Examples/Tmp/*` + `!Reference/Examples/Tmp/.gitkeep`.
- FIRE temp dirs and Laurent dumps need **no** ignore entries: FIRE runs go to `$TemporaryDirectory/FeynmanTrick` (`FeynmanTrick/FeynmanTrick.m:31`, run dirs deleted at `FIREInterface.m:438`), stepwise/pysecdec output to `$TemporaryDirectory` (`run_ft_stepwise.m:123`, `export_pysecdec_family_specs.m:149`), dumps to user-chosen `/tmp` paths. Nothing lands in the worktree.
- Leftover artifacts: none — `temp/` was fully removed in 6130b17; only ignored `.DS_Store` strays remain; worktree is clean.

**P4.3 — inventory (no judgment, per instruction):** `Papers/FeynmanTrick/` = paper source (main.tex, jheplike.sty, utphys.bst, 7 .eps figures, `anc/`); `Dependencies/fire/FIRE6/` = a full FIRE6 checkout (FIRE6.m, Makefile, configure, bin/, documentation, …) that `FeynmanTrick.m:27-29` points at via the `FIREPath` default. Both correctly gitignored as external clones.

**P4.4 — battery membership (optional):** 13 of 30 test files are outside the battery. Benchmarks (`test_recurrence_speedup.m`, `test_unequal_mass_banana_parity_speed.m`, `test_unequal_mass_full.m`, `test_banana_refactored.m`, `test_five_point_nonplanar.m`) are reasonably excluded for runtime; but `test_decomposition.m`, `test_topiecewise.m`, `test_rational_recurrence.m`, `test_local_series_solver.m`, `test_multiple_polylogarithms.m`, `test_hypergeometric2f1.m` are small and cover unique surface — consider adding them (after giving them the P0.2 exit-code convention; none of them currently calls `Exit`).

---

## What NOT to touch (validated numerics — frozen)

- `DiffExp/RegularizedIntegration.m` sector fitter / Laurent integration, `DiffExp/IntegrationStrategies/Recurrence.m`, `ResonantRecurrence.m`, `Dispatch.m` — all just validated against pySecDec (banana exact, box L0 to 11 digits, battery 17/17). In particular the `10^-24` relative floor at `Recurrence.m:618` is documented as load-bearing — do not "clean it up".
- Reference constants: `Scripts/compare_stepwise_log.py:18-29` (pinned pySecDec banana values), the kinematics in `Scripts/FTExamples.m`, pinned values inside `Docs/FeynmanTrickBoxFamilyStatus.md:14-25`, and the `Tests/*_Matrices/` directories.
- The "Combine into position 1, left-to-right" sequence convention (`Scripts/FTExamples.m:149-154`) and the negative-index refusal in `export_pysecdec_family_specs.m:73-82` (it prevents the silently-wrong-spec bug documented in the banana doc).
- The sequential, single-kernel structure of the battery (license limit) — fix its quoting, not its serialism.

Suggested order: P0 (minutes), P1.5 + P1.1 + P1.3 (one new matrices dir + one test file), P1.4, P1.2, then P2/P3 doc-and-script hygiene in one sweep, P4 last.