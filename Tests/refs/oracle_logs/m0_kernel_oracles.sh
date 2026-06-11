#!/bin/zsh
WK='/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/WolframKernel'
cd /Users/mhidding/Code/diffexp_refactored
echo "[K1] bubble+sunrise oracle dumps $(date +%T)"
rm -rf /tmp/m0_dispatch_bubsun /tmp/m0_laurent_bubsun
mkdir -p /tmp/m0_dispatch_bubsun /tmp/m0_laurent_bubsun
env WolframKernel=$WK DEBUG_DUMP_DISPATCH_DIR=/tmp/m0_dispatch_bubsun \
  DIFFEXP_DUMP_LAURENT_DIR=/tmp/m0_laurent_bubsun \
  FT_EXAMPLES=bubble,sunrise FT_EPS_ORDER=2 FT_WORKING_PRECISION=200 \
  wolframscript -file Scripts/run_ft_stepwise.m > /tmp/m0_oracle_bubsun.log 2>&1
echo "[K1] exit $? ; dispatch dumps: $(ls /tmp/m0_dispatch_bubsun | wc -l) ; laurent dumps: $(ls /tmp/m0_laurent_bubsun | wc -l)"
echo "[K2] box_bubble STEPWISE oracle $(date +%T)"
env WolframKernel=$WK FT_EXAMPLES=box_bubble FT_EPS_ORDER=0 FT_WORKING_PRECISION=300 \
  FT_EXPANSION_ORDER=40 FT_DIVISION_ORDER=4 FT_BOUNDARY_EXTRA_ORDER=16 FT_POLE_ALLOWANCE=8 \
  wolframscript -file Scripts/run_ft_stepwise.m > /tmp/m0_oracle_boxbubble.log 2>&1
echo "[K2] exit $? ; STEPWISE rows: $(grep -c STEPWISE /tmp/m0_oracle_boxbubble.log)"
echo "[K] DONE $(date +%T)"
