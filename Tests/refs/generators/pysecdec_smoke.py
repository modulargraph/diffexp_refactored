#!/usr/bin/env python3
from pathlib import Path

from pySecDec import LoopIntegralFromPropagators, loop_package

work = Path("/tmp/pysecdec-smoke")
work.mkdir(parents=True, exist_ok=True)
import os
os.chdir(work)

li = LoopIntegralFromPropagators(
    propagators=["k**2", "(k-p)**2"],
    loop_momenta=["k"],
    replacement_rules=[("p*p", "-1")],
    dimensionality="4-2*eps",
)

loop_package(
    name="bubble_smoke",
    loop_integral=li,
    requested_orders=[0],
    contour_deformation=False,
    form_threads=1,
    processes=1,
)

print(work / "bubble_smoke")
