#!/usr/bin/env python3
"""Pentagon loop_package reference: 1-loop massless 5-point,
{s12,s23,s34,s45,s15} = {-1,-2,-3,-5,-7}, dot products fixed numerically
to match Scripts/FTExamples.m (FT props are Euclidean -(...)^2: relative
sign (-1)^nprops = -1 for 5 props)."""
import os, shutil, subprocess
from pathlib import Path
from pySecDec import LoopIntegralFromPropagators, loop_package
from pySecDec.integral_interface import IntegralLibrary

name = "pentagon_ft"
root = Path("/tmp/pysecdec-ft-pentagon")
work = root / name
pkg = work / name
work.mkdir(parents=True, exist_ok=True)
os.chdir(work)

if not (pkg / f"{name}_pylink.so").exists():
    li = LoopIntegralFromPropagators(
        propagators=[
            "l1**2", "(l1+p1)**2", "(l1+p1+p2)**2",
            "(l1+p1+p2+p3)**2", "(l1+p1+p2+p3+p4)**2",
        ],
        loop_momenta=["l1"],
        external_momenta=["p1", "p2", "p3", "p4"],
        replacement_rules=[
            ("p1*p1", "0"), ("p2*p2", "0"), ("p3*p3", "0"), ("p4*p4", "0"),
            ("p1*p2", "-1/2"), ("p2*p3", "-1"), ("p3*p4", "-3/2"),
            ("p1*p3", "-1"), ("p2*p4", "-1"), ("p1*p4", "5"),
        ],
        dimensionality="4-2*eps",
    )
    loop_package(name=name, loop_integral=li, requested_orders=[1],
                 real_parameters=[], contour_deformation=False,
                 form_threads=1, processes=1)
    contrib_bin = Path("/tmp/pysecdec-venv/lib/python3.14/site-packages/pySecDecContrib/bin")
    env = os.environ.copy()
    env["PATH"] = f"{contrib_bin}:{env['PATH']}"
    subprocess.run(["make", "-j3"], cwd=pkg, env=env, check=True)

lib = IntegralLibrary(str(pkg / f"{name}_pylink.so"))
lib.use_Cuhre(epsrel=1e-8, epsabs=1e-12, maxeval=3000000)
res = lib(epsrel=1e-8, epsabs=1e-12, maxeval=3000000)
print("PENTAGON RESULT:")
print(res[2])
