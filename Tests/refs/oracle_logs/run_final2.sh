#!/bin/zsh
WK='/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/WolframKernel'
cd /Users/mhidding/Code/diffexp_refactored
for ex in pentagon box_triangle; do
  EXTRA=12; [ $ex = box_triangle ] && EXTRA=14
  echo "[F] START $ex $(date +%T)"
  env WolframKernel=$WK FT_EXAMPLES=$ex FT_EPS_ORDER=0 FT_WORKING_PRECISION=300 \
    FT_EXPANSION_ORDER=40 FT_DIVISION_ORDER=4 FT_BOUNDARY_EXTRA_ORDER=$EXTRA \
    FT_POLE_ALLOWANCE=8 wolframscript -file Scripts/run_ft_stepwise.m \
    > /tmp/f_${ex}.log 2>&1
  echo "[F] DONE $ex exit $? $(date +%T)"
done
echo "[F] ALL DONE"
