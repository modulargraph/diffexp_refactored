# Four-loop banana investigation status

**State:** intentionally paused on 2026-07-26 while other integral families
are developed.

The four-loop banana is not yet an accepted example. No run has completed the
final publication and passed the independent Bessel-oracle comparison. This
does not invalidate the fixes already merged into `master`: they have focused
native regressions, and the complete native suite passes 49/49 at master
commit `6e0d761`.

The detailed, append-only evidence ledger is
[`banana4-debug-log.md`](banana4-debug-log.md). This page is the concise
handoff and should be updated whenever the investigation is resumed or paused.

## What is established

- FIRE preparation and the earlier ladder are not the active correctness
  blocker. A current-source level-1 boundary was produced by B4-241.
- Levels 3 and 2 complete with the demand-driven retained-work policy. B4-241
  wrote the level-1 boundary after about 914 seconds; the level-2 singular
  basis returned from the rejected fixed-halo timing of 308 seconds to about
  91 seconds.
- Both final-level algebraic crossings and both endpoint bases have completed
  in production-shaped runs. The final endpoint computations are still much
  too slow, but they are not known mathematical blockers.
- The B4-243 private/public Taylor-width mismatch is fixed. Exact-shadow
  publication now consumes the certified common public prefix without
  discarding the wider private tail reservoir.
- B4-245 reached the authoritative final upper residual. Its failure was
  localized: a small, correlated tail in the reduced/Fuchsian frame was
  converted to independent physical-coordinate balls and propagated through
  an ill-conditioned gauge. The resulting ordinary enclosure inflated to
  roughly `1e99999`, making all 105 required residual coefficients
  inconclusive. This is not evidence that more public epsilon orders are
  needed.
- The replacement implemented in B4-247 certifies the retained tail at a
  common rational clearance seed, keeps its correlations through the exact
  normal-frame map, and evaluates the incoming local at the same physical
  point. A production-shaped Rational-to-Acb lifecycle fixture passes
  matching, materialization, checkpoint save, and restore without building
  the singular-to-ordinary order-800 bridge.
- The native residual now distinguishes a definite mismatch from a typed,
  retryable singular-tail reservoir shortage. A Taylor-order retry must only
  follow the latter classification.

## What remains unknown

The B4-247 correlated-seed path has not reached the terminal match in a real
banana4 run. B4-248 was deliberately stopped during second-atlas preparation
at the user's request. It therefore supplied neither a pass nor a numerical
failure.

The next production question is narrow:

> Does the final upper match report a certified
> `exact_shadow_correlated_tail_seed_status` and avoid
> `singular-ordinary-bridge-start`?

If it does, the complete result must still pass the final Laurent audit and
the equal-mass finite-part oracle

```text
39.655526834297652529992823046933581156446060218710...
```

to at least eight decimal digits. Until then, banana4 must not be described as
working.

## Performance state

The best comparable final-level evidence before the pause is B4-245:

- atlas preparation: about 1072 seconds;
- algebraic singular owners: 392.5 and 412.8 seconds;
- lower endpoint basis: 2532.6 seconds;
- upper endpoint basis: 2080.7 seconds;
- old upper ordinary bridges: 1741.3 seconds;
- old exact proposal: 1836.0 seconds.

These timings are unacceptable for a normal example. The B4-247 path is
expected to remove the old order-800 bridge and proposal work, but that saving
has not yet been measured in production.

B4-248 is not a timing result. A competing exact Wolfram job and severe memory
pressure produced about 2.65 GiB of live swap; its second algebraic owner took
4874.5 seconds of kernel time instead of the previous 412.8 seconds. Ordinary
owners returned to roughly two seconds after pressure subsided.

Profiling implicates eager exact-Rational canonicalization, repeated block
movement, and FLINT `fmpq`/`fmpz_gcd` work in endpoint construction. A tested
alias-safe `fmpq_addmul`/zero-shift CASE-P micro-optimization was flat
(approximately 0.997--1.004x) and was discarded. Do not reintroduce it as a
claimed speedup.

## Durable local artifacts

These artifacts are diagnostic data and are intentionally not tracked by Git:

- trusted current-source level-1 boundary:
  `diffexp_refactored-master-recovery/.diagnostics/`
  `b4-241-clean-full-work-reservoir/checkpoints/`
  `banana4_level1_boundary.mx`;
- B4-245 production log:
  `diffexp_refactored-master-recovery/.diagnostics/`
  `b4-245-width-prefix-resume.log`;
- paused B4-248 log:
  `diffexp_refactored-master-recovery/.diagnostics/`
  `b4-248-correlated-upper-only.log`;
- prepared FIRE snapshot:
  `diffexp_refactored-master-recovery/.diagnostics/cache/prepared/`
  `banana4_2059e788bd7d9c1b25c6d059deb6ba4df79d96c1986722378a63e5ca588b548c.mx`.

Before relying on an artifact, confirm it still exists and retain the
checkpoint configuration/provenance validation. Do not weaken a failed
restore merely to reuse an obsolete chart plan.

## Resume protocol

1. Use `master` at or after `6e0d761`, make a Release build, and require the
   complete native suite to pass.
2. Run with quiet memory, no competing Wolfram kernel, and `caffeinate -i`.
   Do not use `timeout`, `alarm`, or another wall-clock kill around a
   production completion attempt.
3. Restore the B4-241 level-1 boundary. Do not rerun FIRE or levels 4--2 for
   the first terminal test.
4. Keep the established configuration: WP500, EO50, match digits 25, public
   epsilon order 0, boundary extra order 4, division order 4, zero automatic
   singular halo, ten C++ workers, value-aware planning, and basis-only native
   hop execution.
5. Run the upper arm first and inspect the correlated-seed diagnostics. Require
   the old singular-to-ordinary bridge to be absent.
6. If and only if the result is a typed
   `retryable_singular_tail_reservoir`, make one bounded upper-only Taylor
   retry at EO100 and require strict residual improvement. Do not increase
   public epsilon order or matching digits blindly.
7. After the upper terminal match passes, run both arms, the final Laurent
   publication audit, and the independent Bessel oracle.
8. Only after a correct checkpoint-assisted result should a clean-from-scratch
   timing/provenance acceptance run be attempted.

## Closed directions

Do not repeat these without a new theorem or implementation change:

- increasing terminal match digits alone;
- physical-frame fallback with the old transformed consumer;
- accepting a small one-point residual as a final result;
- blindly increasing public epsilon order or Taylor order;
- the fixed 32-row singular halo;
- B4-130's one-digit producer-profile scan;
- weakening checkpoint geometry/provenance validation;
- the discarded CASE-P add-multiply micro-optimization;
- any production run with a wall-clock alarm.

## Effect on other work

Banana4 is an opt-in release investigation and does not block work on other
families. The merged native fixes remain protected by fast deterministic
tests. Pentagon and planar double-box have their own opt-in correctness and
timing gates in [`PerformanceRegression.md`](PerformanceRegression.md), so
they can be developed and qualified independently while banana4 is paused.
