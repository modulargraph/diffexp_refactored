#!/usr/bin/env python3
import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path

from pySecDec import LoopIntegralFromPropagators, loop_package
from pySecDec.integral_interface import IntegralLibrary


PROPAGATORS = [
    "l1**2",
    "(l1+p1)**2",
    "(l1+p1+p2)**2",
    "(l1-l2)**2",
    "l2**2",
    "(l2+p1+p2)**2",
    "(l2+p1+p2+p3)**2",
]

EXAMPLES = {
    "box_bubble": [PROPAGATORS[i] for i in [2, 3, 4, 5, 6]],
    "box_triangle": [PROPAGATORS[i] for i in [1, 2, 3, 4, 5, 6]],
    "double_box_planar": PROPAGATORS,
}

REPLACEMENTS = [
    ("p1*p1", "0"),
    ("p2*p2", "0"),
    ("p3*p3", "0"),
    ("p1*p2", "s/2"),
    ("p2*p3", "t/2"),
    ("p1*p3", "(-s-t)/2"),
]


def package_dir(root: Path, name: str) -> Path:
    return root / name / name


def generate(name: str, root: Path, force: bool) -> Path:
    work = root / name
    pkg = work / name
    if force and work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True, exist_ok=True)
    if (pkg / f"{name}_pylink.so").exists():
        return pkg

    os.chdir(work)
    loop_integral = LoopIntegralFromPropagators(
        propagators=EXAMPLES[name],
        loop_momenta=["l1", "l2"],
        external_momenta=["p1", "p2", "p3"],
        replacement_rules=REPLACEMENTS,
        dimensionality="4-2*eps",
    )
    loop_package(
        name=name,
        loop_integral=loop_integral,
        requested_orders=[0],
        real_parameters=["s", "t"],
        contour_deformation=False,
        form_threads=1,
        processes=1,
    )
    return pkg


def compile_package(pkg: Path, jobs: int) -> None:
    if (pkg / f"{pkg.name}_pylink.so").exists():
        return
    contrib_bin = Path("/tmp/pysecdec-venv/lib/python3.14/site-packages/pySecDecContrib/bin")
    env = os.environ.copy()
    env["PATH"] = f"{contrib_bin}:{env['PATH']}"
    subprocess.run(["make", f"-j{jobs}"], cwd=pkg, env=env, check=True)


def integrate(pkg: Path, maxeval: int, epsrel: float, epsabs: float) -> object:
    lib = IntegralLibrary(str(pkg / f"{pkg.name}_pylink.so"))
    lib.use_Cuhre(epsrel=epsrel, epsabs=epsabs, maxeval=maxeval)
    return lib(
        real_parameters=[-1.0, -1.0 / 3.0],
        epsrel=epsrel,
        epsabs=epsabs,
        maxeval=maxeval,
        format="json",
    )


def jsonable(value):
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [jsonable(item) for item in value]
    if isinstance(value, list):
        return [jsonable(item) for item in value]
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("name", choices=sorted(EXAMPLES))
    parser.add_argument("--root", default="/tmp/pysecdec-ft-boxes")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--no-integrate", action="store_true")
    parser.add_argument("--maxeval", type=int, default=200000)
    parser.add_argument("--epsrel", type=float, default=1e-3)
    parser.add_argument("--epsabs", type=float, default=1e-7)
    args = parser.parse_args()

    root = Path(args.root)
    pkg = generate(args.name, root, args.force)
    compile_package(pkg, args.jobs)
    output = {"name": args.name, "package": str(pkg)}
    if not args.no_integrate:
        output["integral"] = jsonable(
            integrate(pkg, args.maxeval, args.epsrel, args.epsabs)
        )
    print("PYSECDEC_JSON", json.dumps(output, default=str, separators=(",", ":")))


if __name__ == "__main__":
    main()
