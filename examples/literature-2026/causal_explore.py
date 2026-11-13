#!/usr/bin/env python3
"""Try the causal (+i0) side choices that DiffExp cannot derive by itself.

At a physical point DiffExp derives the side of every real singular point
from the Symanzik polynomial where it can.  On multi-loop levels it may stop
with an error naming one point:

    the causal side of the real singular point x = 0.2727... of level 2 could
    not be derived ...; supply an explicit `causal` object ...

For each such family this script reruns the calculation with that point set to
side +1 and to side -1.  If a run stops at a further point, it branches again.
Every combination that completes is kept as ``causal/results/<family>__<tag>.json``
(tag e.g. ``L2p+_L1p-``) together with the family file it ran, so that
``compare_references.py --causal`` can tell which choice reproduces the paper.
Choosing the side by matching published numbers is a fit, not a derivation;
the README says so for every family resolved this way.

    python3 causal_explore.py FAMILY --threads 8

FAMILY is a family in ``families/`` at a physical point (none of the shipped
families needs this; the README lists the physical variants still open).
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "causal"
POINT = re.compile(r"real singular point x = ([0-9.eE+-]+) of level (\d+) could not be derived")


def run(diffexp, family_path, log_path, args):
    command = [diffexp, "ft", str(family_path), "--ibp-provider", "ibp-solver",
               "--ibp-threads", str(args.threads), "--cache", str(args.cache),
               "--fire-seconds", "1800", "--level-seconds", "5400", "--total-seconds", "9000", "--json"]
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads))
    started = time.monotonic()
    with open(log_path, "w") as log:
        try:
            done = subprocess.run(command, stdout=subprocess.PIPE, stderr=log, text=True,
                                  env=env, timeout=args.timeout)
            stdout = done.stdout
        except subprocess.TimeoutExpired:
            return {"status": "timeout"}, time.monotonic() - started
    try:
        response = json.loads(stdout)
    except json.JSONDecodeError:
        response = {"status": "crashed"}
    response.pop("ibp_statistics", None)
    return response, time.monotonic() - started


def explore(diffexp, name, args):
    base = json.loads((HERE / "families" / f"{name}.json").read_text())
    (OUT / "families").mkdir(parents=True, exist_ok=True)
    (OUT / "results").mkdir(parents=True, exist_ok=True)
    queue, runs, summary = deque([[]]), 0, []
    # Reuse the error of a previous plain run instead of repeating it.
    previous = HERE / "results" / f"{name}.json"
    if previous.exists():
        match = POINT.search(json.loads(previous.read_text()).get("error") or "")
        if match:
            point = {"level": int(match.group(2)), "x": match.group(1)[:16]}
            queue = deque([[dict(point, side=1)], [dict(point, side=-1)]])
    while queue and runs < args.max_runs:
        points = queue.popleft()
        tag = "_".join(f"L{p['level']}{'p+' if p['side'] > 0 else 'p-'}{i}" for i, p in enumerate(points)) or "none"
        family = json.loads(json.dumps(base))
        family["name"] = f"{name}__{tag}"
        causal = dict(family.get("causal", {}))
        if points:
            causal["points"] = points
        if causal:
            family["causal"] = causal
        family_path = OUT / "families" / f"{name}__{tag}.json"
        family_path.write_text(json.dumps(family, indent=1) + "\n")
        response, seconds = run(diffexp, family_path, OUT / "results" / f"{name}__{tag}.log", args)
        runs += 1
        status = response.get("status")
        message = response.get("message", "")
        record = {"family": name, "tag": tag, "points": points, "status": status,
                  "wall_seconds": round(seconds, 1), "response": response}
        (OUT / "results" / f"{name}__{tag}.json").write_text(json.dumps(record, indent=1) + "\n")
        summary.append(f"{name} [{tag}]: {status} in {seconds:.0f} s {message[:120]}")
        print(summary[-1], flush=True)
        match = POINT.search(message) if status == "error" else None
        if match:
            x, level = match.group(1), int(match.group(2))
            for side in (1, -1):
                queue.append(points + [{"level": level, "x": x[:16], "side": side}])
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("families", nargs="+")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--jobs", type=int, default=1, help="families explored at the same time")
    parser.add_argument("--timeout", type=int, default=10800)
    parser.add_argument("--max-runs", type=int, default=7, help="runs per family (1 + 2 + 4 for two unknown points)")
    parser.add_argument("--cache", type=Path,
                        default=Path(os.environ.get("DIFFEXP_CACHE", Path.home() / ".cache" / "diffexp-literature-2026")))
    args = parser.parse_args()
    repo_build = HERE.parent.parent / "build" / "diffexp"
    diffexp = os.environ.get("DIFFEXP") or (str(repo_build) if repo_build.exists() else "diffexp")
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        list(pool.map(lambda n: explore(diffexp, n, args), args.families))


if __name__ == "__main__":
    main()
