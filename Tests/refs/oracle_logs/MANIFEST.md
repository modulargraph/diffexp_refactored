# Old-core STEPWISE oracle logs (vendored from /tmp, 2026-06-11)

Campaign records of `Scripts/run_ft_stepwise.m` against the FROZEN old core
(code state f48cd94 — the N-root fitter commit; later commits before
vendoring were docs-only).  Single-license kernel, machine-local runs of
2026-06-11.  STEPWISE row format: JSON per boundary master with Laurent
`Coefficients` pairs [k, c_k] at ~50 digits, FT normalization with prefactor
included; FINAL row = the L0 leading coefficient.  Compare with
`Scripts/compare_stepwise_log.py` (see Tests/PINS.md section 5).

The exact commands/env for every log are in the three vendored run scripts:
`run_ladder2.sh` (l2_*), `run_final2.sh` (f_*), `m0_kernel_oracles.sh`
(../oracle_dumps/m0_oracle_*).  Approximate runtimes: bubble+sunrise ~5 min,
box ~25 min, banana ~20 min, box_bubble ~30 min, pentagon ~85 min,
box_triangle ~45 min per pass.

## PASS records (green campaign state — parity targets for M4/M5)

| log | example(s) | settings (env) | verdict |
|-----|------------|----------------|---------|
| `l2_box.log` | box, FT_EPS_ORDER=1 | WP=300 EXP=40 DIV=4 EXTRA=12 POLE=8 | FINAL eps^-2 = 12 exact, eps^-1 within 1.6e-14 of pin, eps^0 = -41.28416739754647 (11 digits vs pin, inside the reference's 4.7e-10 error).  Also contains the 50-digit L1 closed-form tower rows (PINS.md section 3.2).  NOTE: warning chatter in the log (singular-recurrence head-Plus, DeltaPrescriptions, munfl) is from probing segments and did not affect the validated values — the run IS the pass record. |
| `l2_banana.log` | banana, FT_EPS_ORDER=0 | WP=300 (canonical doc settings were 500/50/4/8) | comparator `failures: 0`; L0 = 8.2681045358689688 vs ref 8.26810451329511583 (2.7e-9 rel); all L1/L2/L3 rows match the hardcoded refs.  Contains honest-window warnings about orders >= 6 — those orders were not requested. |
| `l2_bubsun.log` | bubble + sunrise, FT_EPS_ORDER=2 | WP=200 | finals 0.86081788192800808 / 2.2367927002126465 match the campaign reference values (Docs/FeynmanTrickBananaStatus.md:478-481).  These rows are OLD-CORE values, i.e. parity oracles; independent pySecDec refs for bubble/sunrise are LOST (PINS.md section 6). |

Companion PASS record stored under `../oracle_dumps/`:
`m0_oracle_boxbubble.log` (box_bubble all 4 orders vs -1 x loop_package pin;
byte-identical to the campaign's ft_box_bubble_run3.log).  The /tmp
`l2_boxbubble.log` was truncated mid-run (no FINAL) and was deliberately NOT
vendored.

## FAIL baselines (defect records — M5 must beat these, R8)

| log | example | settings | failure signature |
|-----|---------|----------|-------------------|
| `f_pentagon.log` | pentagon | WP=300 EXP=40 DIV=4 EXTRA=12 POLE=8 | FAILS pin at ALL orders incl. leading pole: FINAL eps^-2 = 0.2521398 vs expected +13/42 = 0.3095238; eps^0 has a spurious imaginary part.  Log shows the mechanism: repeated "resonant endpoint coefficient with zero epsilon regulator; dropping contribution" (the D2 drop rule) + 2x "current point is not recognized as a branch point... add DeltaPrescriptions" (config gap, M0 task 16).  Byte-identical to /tmp/ft_pentagon_run3.log and /tmp/nroot_pentagon.log — the N-root fitter did not change the outcome. |
| `f_box_triangle.log` | box_triangle | WP=300 EXP=40 DIV=4 EXTRA=14 POLE=8 | FAILS pin: FINAL eps^-4 = 153.3649 vs pin 2.25; complex from eps^-3.  Mechanism in-log: 20+ zero-regulator drops, then "omitted 23 endpoint coefficient(s)... affects the reported Laurent window starting at epsilon order -3".  Confirmed D2-class; interior IBP pole at xx3 = 12167/12651 ~ 0.9617 near the L3 endpoint is the geometric trigger. |
| `ft_box_bubble_run2_budget_starved.log` | box_bubble | pre-budget-fix (commit c0b24f3 absent), EXTRA too small | the §3.4 BUDGET PIN record: lines 29/32 "WARNING level 2 transport capped at eps order 2 (needs 9)" / "level 1 ... (needs 11)"; resulting FINAL eps^-1/eps^0 are wrong/empty.  The plan's static budget formula must reproduce the 9 and 11. |

## Run scripts (provenance, exact env)

- `run_ladder2.sh` — produced the l2_* set (battery + box + banana +
  bubble/sunrise + box_bubble + pentagon + box_triangle ladder).
- `run_final2.sh` — produced f_pentagon.log / f_box_triangle.log.
- `m0_kernel_oracles.sh` — produced ../oracle_dumps/ artifacts
  (dispatch + Laurent dumps and the two m0_oracle_*.log records).
