#!/usr/bin/env python3
"""1-loop massless on-shell box via loop_package, to pin the convention
factor between pySecDec's loop_package route and the FT family normalization.
FT box values are independently verified (analytic) in the repo test suite."""
import json
import os
import shutil
import subprocess
from pathlib import Path

from pySecDec import LoopIntegralFromPropagators, loop_package
from pySecDec.integral_interface import IntegralLibrary

ROOT = Path("/tmp/pysecdec-ft-boxes")
NAME = "box1L_pin"


def main():
    work = ROOT / NAME
    pkg = work / NAME
    if not (pkg / f"{NAME}_pylink.so").exists():
        if work.exists():
            shutil.rmtree(work)
        work.mkdir(parents=True, exist_ok=True)
        os.chdir(work)
        li = LoopIntegralFromPropagators(
            propagators=["l1**2", "(l1+p1)**2", "(l1+p1+p2)**2", "(l1+p1+p2+p3)**2"],
            loop_momenta=["l1"],
            external_momenta=["p1", "p2", "p3"],
            replacement_rules=[
                ("p1*p1", "0"), ("p2*p2", "0"), ("p3*p3", "0"),
                ("p1*p2", "s/2"), ("p2*p3", "t/2"), ("p1*p3", "(-s-t)/2"),
            ],
            dimensionality="4-2*eps",
        )
        loop_package(
            name=NAME,
            loop_integral=li,
            requested_orders=[1],
            real_parameters=["s", "t"],
            contour_deformation=False,
            form_threads=1,
            processes=1,
        )
        contrib = subprocess.run(
            ["/tmp/pysecdec-venv/bin/python3", "-c",
             "import pySecDecContrib, os; print(os.path.join(os.path.dirname(pySecDecContrib.__file__), 'bin'))"],
            capture_output=True, text=True, check=True).stdout.strip()
        env = os.environ.copy()
        env["PATH"] = f"{contrib}:{env['PATH']}"
        subprocess.run(["make", "-j4"], cwd=pkg, env=env, check=True)

    lib = IntegralLibrary(str(pkg / f"{NAME}_pylink.so"))
    lib.use_Cuhre(epsrel=1e-7, epsabs=1e-12, maxeval=5_000_000)
    result = lib(real_parameters=[-1.0, -1.0 / 3.0], format="json")
    print("PYSECDEC_BOX1L_PIN", json.dumps(result, default=str))


if __name__ == "__main__":
    main()
