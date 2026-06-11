#!/bin/zsh
WK='/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/WolframKernel'
cd /Users/mhidding/Code/diffexp_refactored
echo "[L2] battery $(date +%T)"
zsh Scripts/run_test_battery.sh > /tmp/l2_battery.log 2>&1
echo "[L2] battery exit $?"
echo "[L2] box $(date +%T)"
env WolframKernel=$WK FT_EXAMPLES=box FT_EPS_ORDER=1 FT_WORKING_PRECISION=300 FT_EXPANSION_ORDER=40 FT_DIVISION_ORDER=4 FT_BOUNDARY_EXTRA_ORDER=12 FT_POLE_ALLOWANCE=8 wolframscript -file Scripts/run_ft_stepwise.m > /tmp/l2_box.log 2>&1
echo "[L2] box exit $?"
echo "[L2] banana $(date +%T)"
env WolframKernel=$WK FT_EXAMPLES=banana FT_EPS_ORDER=0 FT_WORKING_PRECISION=300 wolframscript -file Scripts/run_ft_stepwise.m > /tmp/l2_banana.log 2>&1
echo "[L2] banana exit $?"
python3 Scripts/compare_stepwise_log.py /tmp/l2_banana.log 2>/dev/null | tail -1
echo "[L2] bubble+sunrise $(date +%T)"
env WolframKernel=$WK FT_EXAMPLES=bubble,sunrise FT_EPS_ORDER=2 FT_WORKING_PRECISION=200 wolframscript -file Scripts/run_ft_stepwise.m > /tmp/l2_bubsun.log 2>&1
python3 Scripts/compare_stepwise_log.py /tmp/l2_bubsun.log 2>/dev/null | tail -1
echo "[L2] box_bubble $(date +%T)"
env WolframKernel=$WK FT_EXAMPLES=box_bubble FT_EPS_ORDER=0 FT_WORKING_PRECISION=300 FT_EXPANSION_ORDER=40 FT_DIVISION_ORDER=4 FT_BOUNDARY_EXTRA_ORDER=16 FT_POLE_ALLOWANCE=8 wolframscript -file Scripts/run_ft_stepwise.m > /tmp/l2_boxbubble.log 2>&1
echo "[L2] box_bubble exit $?"
echo "[L2] pentagon $(date +%T)"
env WolframKernel=$WK FT_EXAMPLES=pentagon FT_EPS_ORDER=0 FT_WORKING_PRECISION=300 FT_EXPANSION_ORDER=40 FT_DIVISION_ORDER=4 FT_BOUNDARY_EXTRA_ORDER=12 FT_POLE_ALLOWANCE=8 wolframscript -file Scripts/run_ft_stepwise.m > /tmp/l2_pentagon.log 2>&1
echo "[L2] pentagon exit $?"
echo "[L2] box_triangle $(date +%T)"
env WolframKernel=$WK FT_EXAMPLES=box_triangle FT_EPS_ORDER=0 FT_WORKING_PRECISION=300 FT_EXPANSION_ORDER=40 FT_DIVISION_ORDER=4 FT_BOUNDARY_EXTRA_ORDER=14 FT_POLE_ALLOWANCE=8 wolframscript -file Scripts/run_ft_stepwise.m > /tmp/l2_boxtri.log 2>&1
echo "[L2] box_triangle exit $?"
echo "[L2] ALL DONE $(date +%T)"
