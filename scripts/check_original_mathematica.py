#!/usr/bin/env python3
"""Run the original notebook computations through the native process wrapper.

Requires fetched original ancillary data and an explicitly selected licensed
WolframKernel. Calls run sequentially because installations may allow one kernel.
No Feynman-trick reconstruction is started by this runner.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
CASES = {
    "mpl": ("original_mpl.wl", {}),
    "1loop": ("original.wl", {"DIFFEXP_ORIGINAL_FAMILY": "1loop"}),
    "zmz": ("original.wl", {"DIFFEXP_ORIGINAL_FAMILY": "zmz"}),
    "mzz": ("original.wl", {"DIFFEXP_ORIGINAL_FAMILY": "mzz"}),
    "zzz": ("original.wl", {"DIFFEXP_ORIGINAL_FAMILY": "zzz"}),
    "henn": ("original.wl", {"DIFFEXP_ORIGINAL_FAMILY": "henn"}),
    "banana": ("original_banana.wl", {}),
    "banana-routes": ("original_banana_routes.wl", {}),
    "banana-plot": ("original_banana_plot.wl", {}),
    "zzz-high": ("original.wl", {"DIFFEXP_ORIGINAL_FAMILY": "zzz", "DIFFEXP_HIGH_PRECISION": "1"}),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kernel", required=True, type=Path)
    parser.add_argument("--case", action="append", choices=CASES)
    parser.add_argument("--include-high-precision", action="store_true")
    parser.add_argument("--timeout", type=int, help="Maximum seconds per workflow (default: 1800; 3600 for zzz-high)")
    args = parser.parse_args()
    if (args.timeout is not None and args.timeout <= 0) or not args.kernel.is_file():
        parser.error("Select an existing WolframKernel and a positive timeout")
    cases = args.case or [c for c in CASES if c != "zzz-high"]
    if args.include_high_precision and "zzz-high" not in cases:
        cases.append("zzz-high")
    if "banana-plot" in cases and not (ROOT / "build-reference/banana-wrapper-equal-saved.m").exists() and "banana-routes" not in cases:
        cases.insert(cases.index("banana-plot"), "banana-routes")
    if "banana-routes" in cases and not (ROOT / "build-reference/banana-wrapper-minus-one.m").exists() and "banana" not in cases:
        cases.insert(cases.index("banana-routes"), "banana")
    folder = ROOT / "build-reference/mathematica-workflows"
    folder.mkdir(parents=True, exist_ok=True)
    reports = []
    for case in cases:
        script, settings = CASES[case]
        environment = dict(os.environ)
        environment.pop("DIFFEXP_HIGH_PRECISION", None)
        environment.pop("DIFFEXP_ORIGINAL_FAMILY", None)
        environment.pop("DIFFEXP_REQUEST_ONLY", None)
        environment.pop("DIFFEXP_RESPONSE_FILE", None)
        environment.update(settings)
        timeout = args.timeout if args.timeout is not None else (3600 if case == "zzz-high" else 1800)
        start = time.monotonic()
        with (folder / f"{case}.log").open("w") as log:
            process = subprocess.Popen(
                [str(args.kernel.resolve()), "-noprompt", "-script", str(ROOT / "tests/mathematica" / script)],
                cwd=ROOT, env=environment, stdout=log, stderr=subprocess.STDOUT,
                start_new_session=True,
            )
            try:
                code = process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                # Include the child native process in timeout cleanup.
                import signal
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                code = "timeout"
        report = {"case": case, "exit": code, "wall_seconds": time.monotonic() - start, "time_limit_seconds": timeout}
        reports.append(report)
        (folder / "runner.json").write_text(json.dumps(reports, indent=2) + "\n")
        print(json.dumps(report), flush=True)
        if code != 0:
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
