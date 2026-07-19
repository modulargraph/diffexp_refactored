# Four-loop banana debugging ledger

This file is the durable record for four-loop banana investigations. Update a
row before starting any long run, then fill in its duration, result, and
conclusion when the run ends.

Rules:

- Do not repeat a configuration already marked disproven unless a new change
  invalidates its conclusion.
- Record the source commit (and whether the tree is dirty), the exact runtime
  configuration, and the complete log/checkpoint path.
- A small local matching residual is not a correctness result. Compare the
  final Laurent coefficients with the oracle and record the scale of any
  discrepancy.
- Prefer checkpoint-producing runs. A long run which cannot be replayed must
  justify why no replayable route exists.
- Stop a run once it can no longer answer its stated hypothesis.

## Target and current status

- Equal-mass four-loop banana finite-part oracle:
  `39.655526834297652529992823046933581156446060218710...`
- Accurate incoming level-1 boundary:
  `/private/tmp/diffexp2-banana4-fresh-boundary-20260718/checkpoints/banana4_level1_boundary.mx`
- Current branch: `codex/all-examples-recovery`
- Current recovery commit: `3f0dfce`.
- Historical dirty-tree diff fingerprint for B4-11 through B4-13:
  `71955139205bf5cdef753235ad78a88e0963af7f1d25cfb5367a2196772192fb`.
- The incoming boundary is accurate. The known failure is in the final
  singular chart: locally accepted matching data can be amplified by roughly
  60 orders at the endpoint.
- No four-loop banana result has yet passed the oracle check.

## Experiment ledger

| ID | Hypothesis and change | Exact configuration | Duration | Result | Conclusion / next action |
|---|---|---|---:|---|---|
| B4-01 | Private high epsilon orders were changing public Laurent decisions because `ESTrim` chose its scale from the entire series. Added `ESTrimThrough` and used the public downstream top order. | Wide private profile `{20,15,4,0}`. Log: `/private/tmp/diffexp2-banana4-prefix-fix-20260719/exact-final.log` | ~42.5 min | Correct public coefficients survived early levels, but terminal result was wrong: pole `-6.8339e27`, finite `-1.53465e46`. Epsilon-series tests passed 31/31. | This was a real prefix-trimming bug and the fix stays. It is not the terminal banana cause. Do not rerun the old `ESTrim` behavior. |
| B4-02 | Raising terminal matching digits might stabilize the normalized singular frame. | `DE2_DIAGNOSTIC_TERMINAL_MATCH_DIGITS=75`, wide profile, stop after boundary 1. Log: `/private/tmp/diffexp2-banana4-prefix-fix-20260719/terminal75-boundary1.log` | Boundary-only | Boundary errors worsened to approximately `1e-3` through `1e-1`. | More digits alone do not repair the normalized-frame formulation. |
| B4-03 | Matching in a physical-frame fallback might avoid the normalized-frame reservoir deficit. | Strict terminal normalized frame disabled; minimal profile `{3,15,3,0}`. Log: `/private/tmp/diffexp2-banana4-prefix-fix-20260719/physical-fallback-halo3-final.log` | ~25.7 min | Wrong pole `-7.171e19`, finite `-1.473e36`. This version still transformed functionals through the exact-right frame, so it was not a genuinely direct physical contraction. | Physical-frame matching with the old consumer is insufficient. |
| B4-04 | Verify independently that the FIRE/earlier ladder is not corrupting the incoming terminal boundary. | WP 200, match digits 8, source expansion order 50, epsilon order 0, boundary extra 4, halo 0. Log: `/private/tmp/diffexp2-banana4-fresh-boundary-20260718/run.log`; checkpoint: `/private/tmp/diffexp2-banana4-fresh-boundary-20260718/checkpoints/banana4_level1_boundary.mx` | ~63 s | Produced accurate rational-looking level-1 boundary coefficients. | Treat this checkpoint as the trusted terminal input. Do not rerun FIRE for terminal experiments. |
| B4-05 | Use the Wolfram recurrence backend as an independent terminal control. | Accurate boundary replay, boundary extra 3, requested 4. EO50 log: `/private/tmp/diffexp2-banana4-level1-wolfram-replay.log`; EO25 log: `/private/tmp/diffexp2-banana4-level1-wolfram-eo25-replay.log` | >5 min each, stopped | Both remained on the lower arm and did not produce a result. | This control is too slow for quick diagnosis and supplied no correctness evidence. |
| B4-06 | The terminal functional convolution cap was computed from `basis shift + weight shift`, which can cancel and hide required orders. Changed it to depend on the minimum weight power plus primitive loss. | Focused native checkpoint/cancellation tests. | Test-scale | Tests passed. A basis at `+N` and a weight at `-N` now correctly demand `N` additional functional orders. | This was a real truncation bug. In the transformed path it also exposes the genuinely large order request; it does not alone fix banana. |
| B4-07 | A direct physical endpoint consumer might eliminate artificial exact-right-frame order demand. | Accurate boundary, EO25, WP200, terminal direct-physical diagnostic, no halo. Log: `/private/tmp/diffexp2-banana4-direct-physical-eo25-20260719.log` | ~11 min | Failed at terminal match 17 with a structured retry requesting 14 additional orders: strict normalized/exact-right residual had no common complete window. | Direct contraction cannot help while certification is performed in the incompatible transformed frame. |
| B4-08 | Combine physical-frame matching with a truly direct physical contraction `L(F) * w`. | Accurate boundary; EO25; WP200; seeded halo 3 using `{1->3,2->0,3->0,4->0}`; accepted plan `required={4,5,6,7}`. Log: `/private/tmp/diffexp2-banana4-direct-physical-halo3-eo25-20260719.log` | 1454.14 s (~24.2 min) | Locally accepted residual near `1e-8`, but wrong pole `3.1089124389272879e29` and finite `8.339242179006811e51`; final pole guard rejected it. | Direct physical weights, consumer, and residual frame are still insufficient. Endpoint evaluation amplifies locally tolerated errors enormously. Do not rerun this configuration unchanged. |
| B4-09 | Save the completed halo-3 native state so the terminal failure can be inspected without another full march. | Same as B4-08 with `FT_SAVE_NATIVE_TRANSPORT_CHECKPOINT=1`. The runner used `RunNativeTransportObservableBatchOwned`, which already forces `"ConsumeReceivingBases" -> True`. Log: `/private/tmp/diffexp2-banana4-terminal-debug-checkpoint-run.log`; requested directory: `/private/tmp/diffexp2-banana4-terminal-debug-checkpoints` | ~24 min | Arms completed, then atomic save failed: an unconsumed Rational-shadow SCC basis has a live-only exact witness and cannot be serialized. No checkpoint was published. | Initial interpretation (“the runner forgot consuming mode”) was wrong. The terminal factorized match deliberately retains its physical basis after registry consumption so later observable contraction can use it; that strong ownership keeps the witness-bearing locals in the checkpoint closure. Fix this specific post-match ownership/serialization case. |
| B4-10 | Cache source-operator fingerprints/references to avoid repeated hashing and copying of multi-megabyte identities. | Focused persistent transport round-trip tests. | Test-scale | Focused test passed. Full banana timing effect has not been isolated. Stack samples from B4-08/B4-09 remained dominated by exact Rational and parallel Acb recurrence work. | Keep as a performance improvement, but do not claim it solves the 24-minute recurrence. |
| B4-11 | Make the consumed terminal Rational-shadow factorization checkpointable without allowing a public/pre-match shadow basis to masquerade as replayable. The snapshot now permits witness-free local serialization only for private basis owners of a terminal factorized match whose public basis handles and match handle have both been consumed. | Dirty tree based on `90c2e3e`; pre-change diff fingerprint `71955139205bf5cdef753235ad78a88e0963af7f1d25cfb5367a2196772192fb`. | Test-scale | Full native build succeeded. `cpp_persistent_checkpoint`, `cpp_persistent_checkpoint_deferred_state`, and `cpp_persistent_transport_run_arms` passed. The Wolfram one-block CASE-P test still passed, including its assertion that checkpointing an unconsumed Rational-shadow basis fails closed. | The ownership distinction is now encoded and ordinary terminal round-trip plus pre-match rejection are covered. Add/obtain one consumed CASE-P terminal round-trip before spending 24 minutes on banana. |
| B4-12 | First consumed CASE-P terminal checkpoint fixture, forcing the same Rational-shadow specialization used by banana. | Public target ε⁰, source window `[-2,0]`, one lower singular terminal hop. | <1 s | Correctly failed before matching with `common_complete_max=-3`, `required_complete_max=0`, and structured retry `required_additional_epsilon_orders=3`; all four downstream assertions were therefore skipped/failed. | Widening the source through ε³ made the match complete, but this particular multi-block CASE-P system selected a direct physical materialization and therefore did not retain a terminal exact-right factorization. It cannot exercise the post-match checkpoint closure. Keep the separate CASE-P schedule coverage and use a deterministic scalar Rational-shadow terminal frame for the ownership test. |
| B4-13 | Prove the exact post-match ownership case with a small deterministic fixture: prepare a regular two-block SCC plan using Rational-shadow imports, release the eager public witness handles, stream fresh shadow bases through consuming hops, checkpoint both retained terminal factorizations, restore, and export. | Public target ε⁰, source `[-2,3]`, two one-hop arms, WP80/EO10. Test: `Tests/test_native_casep_terminal_checkpoint.m`. | ~2 s | 4/4 passed. Both states report `terminal_factorized_match=true`; atomic v2 checkpoint save succeeds; restore repeats zero arm marches; exported epsilon series is byte-identical; both sessions release cleanly. The separate one-block singular test still proves an unconsumed shadow basis fails closed. | The narrow serializer exception now has direct positive and negative lifecycle coverage. Safe to commit it before one replayable banana run. |
| B4-14 | Reuse an already captured terminal post-hop state to inspect endpoint amplification before paying for another march. This is a restore/contract experiment only; it does not validate the new serializer because the artifact predates `3f0dfce`. | Artifact: `.diagnostics/banana4-upper-balanced-md8-posthop.de2cp` (SHA-256 `8dd4699bece0de5f90c053fccf7b306fa7ee8be8dae0a8cb91f5a6e3fe4d180e`, 893 MiB); manifest SHA-256 `f79763e0ab67c03e5d4870b46ac1e34740ee9c06965ecb8d7d29319570069365`; EO25, match digits 8, upper arm, ε window `[-1,13]`, required through 5. Replay logs: `/private/tmp/diffexp2-banana4-b4-14-replay.log`, `/private/tmp/diffexp2-banana4-b4-14-factorized-replay.log`, `/private/tmp/diffexp2-banana4-b4-14-adjoint-replay.log`, and `/private/tmp/diffexp2-banana4-b4-14-state-replay.log`. | ~15–20 s per restore | Restore and contraction succeeded. Physical, factorized-route, and adjoint-route exports were byte-identical: ε⁰ `-1.474379443216296664084279224013345692996e36`, with spurious ε⁻¹ `-7.1759727083684920389e19`. State statistics prove this old v4 artifact has `terminal_factorized_match=false`; therefore route selection never entered the factorized terminal consumer. | The artifact is useful as a fast historical wrong-result replay, but it cannot reveal the current retained terminal match, modes, or weights. Do not mistake its route equality for evidence about the current factorization. A fresh current checkpoint-producing run is still required. |

### Implementation notes (not numerical experiments)

- While adding B4-14 state diagnostics, the first compile used
  `material_sector` outside its `local_detail` namespace and called a
  nonexistent `EpsilonFrame::width()` accessor. The build failed immediately;
  the implementation now uses `local_detail::material_sector` and
  `coefficients().size()`. This did not run or alter a banana calculation.

## Disproven configurations: do not repeat unchanged

- Increasing terminal match digits in the normalized frame.
- Physical-frame fallback while still consuming transformed exact-right
  functionals.
- Direct physical contraction while retaining transformed-frame certification.
- Combined physical certification and direct contraction with EO25/halo 3 as
  currently formulated: it is locally residual-valid but globally wrong.
- Re-running checkpoint mode without changing terminal-factorized-match
  ownership. Receiving-basis registry consumption is already enabled; the
  terminal match itself retains the witness-bearing physical basis.
- Restarting FIRE or the earlier ladder when investigating the final singular
  chart; use the trusted level-1 boundary checkpoint.

## Next controlled experiment

1. Add a checkpoint-safe post-match representation for the terminal
   factorization. It must retain the physical Acb basis/weights needed for
   endpoint and line contraction without making a pre-match Rational-shadow
   basis replayable as though its live exact witness still existed.
2. Prove with focused tests that a terminal factorized state round-trips, while
   the existing pre-match Rational-shadow checkpoint test continues to fail
   closed.
3. Run exactly one EO25/halo-3 terminal replay with checkpoint saving enabled.
4. Restore the checkpoint and inspect the materialized terminal local sector
   tags, physical weights, and per-sector contribution scales.
5. Use those data to implement an endpoint-regularity certificate/projection
   for forbidden Fuchsian modes. A local match residual by itself is not a
   sufficient terminal acceptance criterion.

## Additional historical evidence

Earlier terminal variants were also globally wrong:

- exact-right structural minimum, EO25: finite scale `5.915852757e28` after
  about 1870 seconds;
- strict residual variant: finite scale `6.225439972e38`;
- exact-right-frame variant: finite scale `9.22622045e28`;
- adjoint halo-8 variant: pole `-2.1855e23`.

These results reinforce that merely increasing halo/order or moving the same
residual check between closely related frames is not a systematic fix.
