# Old-core oracle artifacts (M0 pre-generation)

Generated at master=1130c72 (old DiffExp frozen as parity oracle) on
2026-06-11 with the single-license kernel. These let M3/M4 parity gates
run with ZERO old-code kernel time.

- dispatch_bubsun.tar.gz — 809 per-call solver-input dumps
  (DEBUG_DUMP_DISPATCH_DIR: each file holds <|"ctx","bVec","epsord"|> at
  DispatchStrategy entry) from FT_EXAMPLES=bubble,sunrise FT_EPS_ORDER=2
  FT_WORKING_PRECISION=200. Replay any dump through DiffExp2 Solve.m and
  compare against the old solver's output for M3 regular/singular-chart
  parity. Extract: tar -xzf dispatch_bubsun.tar.gz -C /tmp
- laurent_bubsun.tar.gz — 4 integral-level dumps
  (DIFFEXP_DUMP_LAURENT_DIR; replay with Scripts/eval_dump_generic.m).
- m0_oracle_bubsun.log — STEPWISE rows, comparator-verified (failures: 0).
- m0_oracle_boxbubble.log — complete box_bubble STEPWISE record
  (FINAL finite = 0.5, matches pin -1 x -0.5; replaces the power-cut-
  truncated l2_boxbubble.log).
