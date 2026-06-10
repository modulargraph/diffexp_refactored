#!/usr/bin/env python3
"""Generate, build, and run pySecDec packages from exported FT family specs."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import sympy as sp
from pySecDec import make_package


def is_one(expr: str) -> bool:
    return expr.strip() in {"1", "1.0", "(1)"}


def is_zero(expr: str) -> bool:
    return expr.strip() in {"0", "0.0", "(0)"}


def sympy_to_string(expr) -> str:
    return str(sp.simplify(expr)).replace(" ", "")


def factor_variable_monomial(expr, variables):
    if not variables:
        return [sp.Integer(0)] * len(variables), sp.sympify(expr)
    symbols = sp.symbols(" ".join(variables))
    if len(variables) == 1:
        symbols = (symbols,)
    expr_sp = sp.sympify(expr)
    poly = sp.Poly(sp.expand(expr_sp), *symbols, domain="EX")
    mins = [min(monom[i] for monom, _ in poly.terms()) for i in range(len(symbols))]
    monomial = sp.prod(symbol ** power for symbol, power in zip(symbols, mins))
    return mins, sp.simplify(expr_sp / monomial)


def factored_polynomials(spec: dict) -> tuple[list[str], str]:
    variables = list(spec["Variables"])
    eps = sp.symbols("eps")
    upower = sp.sympify(spec["UPower"])
    fpower = sp.sympify(spec["FPower"])
    u_exps, u_reg = factor_variable_monomial(spec["U"], variables)
    f_exps, f_reg = factor_variable_monomial(spec["F"], variables)
    r_exps, r_reg = factor_variable_monomial(spec["Remainder"], variables)

    symbols = sp.symbols(" ".join(variables)) if variables else ()
    if variables and len(variables) == 1:
        symbols = (symbols,)

    polynomials: list[str] = []
    for symbol, u_exp, f_exp, r_exp in zip(symbols, u_exps, f_exps, r_exps):
        exponent = sp.simplify(u_exp * upower + f_exp * fpower + r_exp)
        if exponent != 0:
            polynomials.append(f"({symbol})**({sympy_to_string(exponent)})")

    if sp.simplify(u_reg - 1) != 0 and not is_zero(spec["UPower"]):
        polynomials.append(f"({sympy_to_string(u_reg)})**({spec['UPower']})")
    if sp.simplify(f_reg - 1) != 0 and not is_zero(spec["FPower"]):
        polynomials.append(f"({sympy_to_string(f_reg)})**({spec['FPower']})")
    if not polynomials:
        polynomials = ["1"]
    return polynomials, sympy_to_string(r_reg)


def package_name(name: str) -> str:
    safe = re.sub(r"[^0-9A-Za-z_]", "_", name)
    if safe and safe[0].isdigit():
        safe = "p_" + safe
    return safe


def generate(
    spec: dict,
    output_root: Path,
    force: bool,
    decomposition_method: str,
    factor_monomials: bool,
    split: bool,
) -> Path:
    name = package_name(spec["Name"])
    package_dir = output_root / name
    if package_dir.exists():
        if force:
            shutil.rmtree(package_dir)
        else:
            return package_dir

    output_root.mkdir(parents=True, exist_ok=True)
    old_cwd = Path.cwd()
    os.chdir(output_root)
    try:
        variables = list(spec["Variables"])
        remainder_expression = spec["Remainder"]
        if not variables:
            variables = ["zdummy"]
            remainder_expression = f"({remainder_expression})"

        if factor_monomials:
            polynomials, remainder_expression = factored_polynomials(spec)
        else:
            polynomials = []
            if not is_one(spec["U"]) and not is_zero(spec["UPower"]):
                polynomials.append(f"({spec['U']})**({spec['UPower']})")
            if not is_one(spec["F"]) and not is_zero(spec["FPower"]):
                polynomials.append(f"({spec['F']})**({spec['FPower']})")
            if not polynomials:
                polynomials = ["1"]

        make_package(
            name=name,
            integration_variables=variables,
            regulators=["eps"],
            requested_orders=spec.get("RequestedOrders", [0]),
            polynomials_to_decompose=polynomials,
            prefactor=spec["Prefactor"],
            remainder_expression=remainder_expression,
            decomposition_method=decomposition_method,
            split=split,
            form_threads=1,
            processes=1,
            use_dreadnaut=False,
            use_Pak=False,
            use_light_Pak=False,
            use_iterative_sort=True,
        )
    finally:
        os.chdir(old_cwd)
    return package_dir


def run_cmd(cmd: list[str], cwd: Path, log_file: Path) -> str:
    proc = subprocess.run(
        cmd,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    log_file.write_text(proc.stdout)
    if proc.returncode != 0:
        raise SystemExit(
            f"command failed ({proc.returncode}): {' '.join(cmd)}\n"
            f"log: {log_file}"
        )
    return proc.stdout


def build_and_run(package_dir: Path, name: str) -> str:
    make_log = package_dir / "_build.log"
    run_log = package_dir / "_run.log"
    run_cmd(["make", "-j2"], package_dir, make_log)
    target = f"integrate_{name}"
    run_cmd(["make", target, "-j2"], package_dir, make_log)
    return run_cmd([str(package_dir / target)], package_dir, run_log)


def summarize_output(output: str) -> str:
    lines = [
        line.strip()
        for line in output.splitlines()
        if line.strip().startswith("amplitude")
    ]
    return lines[-1] if lines else output.strip().splitlines()[-1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("spec_file", type=Path)
    parser.add_argument("--output-root", type=Path, default=Path("/tmp/pysecdec-ft-family"))
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--only", nargs="*", default=None)
    parser.add_argument("--decomposition-method", default="geometric_no_primary")
    parser.add_argument("--factor-monomials", action="store_true")
    parser.add_argument("--split", action="store_true")
    args = parser.parse_args()

    specs = json.loads(args.spec_file.read_text())
    results = []
    for spec in specs:
        name = package_name(spec["Name"])
        if args.only and name not in args.only and spec["Name"] not in args.only:
            continue
        print(f"=== {name} ===", flush=True)
        package_dir = generate(
            spec,
            args.output_root,
            args.force,
            args.decomposition_method,
            args.factor_monomials,
            args.split,
        )
        output = build_and_run(package_dir, name)
        summary = summarize_output(output)
        result = {
            "name": name,
            "source_name": spec["Name"],
            "package_dir": str(package_dir),
            "summary": summary,
            "spec": spec,
        }
        results.append(result)
        print(summary, flush=True)

    results_file = args.output_root / "results.json"
    args.output_root.mkdir(parents=True, exist_ok=True)
    results_file.write_text(json.dumps(results, indent=2))
    print(f"Results={results_file}")


if __name__ == "__main__":
    main()
