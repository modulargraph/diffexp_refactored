#!/usr/bin/env python3
"""Run DiffExp 2.2 Feynman-trick calculations for the 2026 literature families.

Every family is a plain ``DiffExp.FeynmanFamily/v1`` JSON file in ``families/``
with a sibling ``<name>.meta.json`` describing its source paper.  This driver
only calls the public command line:

    diffexp ft FAMILY.json --ibp-provider ibp-solver --cache CACHE --json

and records the response, wall time and peak memory of each run in
``results/<name>.json``.  ``make_table.py`` turns those records into the
README table.

Usage (from anywhere; the executable defaults to ``build/diffexp`` of this
checkout, else ``$DIFFEXP``, else ``diffexp`` on the PATH):

    python3 run_families.py                          # all families, one after another
    python3 run_families.py hvp_chpt_e5e6_d2         # one family
    python3 run_families.py --threads 8 --jobs 4     # four families side by side

Only the Python standard library is used.
"""

import argparse
import json
import os
import resource
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
FAMILIES = HERE / "families"
RESULTS = HERE / "results"


def find_diffexp(explicit):
    repo_build = HERE.parent.parent / "build" / "diffexp"
    for candidate in (explicit, os.environ.get("DIFFEXP"), str(repo_build), shutil.which("diffexp")):
        if candidate and Path(candidate).is_file() and os.access(candidate, os.X_OK):
            return str(Path(candidate).resolve())
    sys.exit("diffexp executable not found: pass --diffexp or set DIFFEXP")


def family_names(selected):
    names = sorted(p.stem for p in FAMILIES.glob("*.json") if not p.name.endswith(".meta.json"))
    if selected:
        missing = [n for n in selected if n not in names]
        if missing:
            sys.exit("unknown families: " + ", ".join(missing))
        return selected
    return names


def children_maxrss_bytes():
    rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return rss if sys.platform == "darwin" else rss * 1024  # Linux reports KiB


def run_one(diffexp, name, args):
    family = FAMILIES / f"{name}.json"
    log_path = RESULTS / f"{name}.log"
    command = [
        diffexp, "ft", str(family),
        "--ibp-provider", "ibp-solver",
        "--ibp-threads", str(args.threads),
        "--cache", str(args.cache),
        "--fire-seconds", str(args.fire_seconds),
        "--level-seconds", str(args.level_seconds),
        "--total-seconds", str(args.prep_seconds),
        "--json",
    ] + args.extra
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads))
    record = {"family": name, "command": command, "threads": args.threads}
    started = time.monotonic()
    with open(log_path, "w") as log:
        try:
            done = subprocess.run(command, stdout=subprocess.PIPE, stderr=log, env=env,
                                  timeout=args.timeout, text=True)
            record["exit_code"] = done.returncode
            stdout = done.stdout
        except subprocess.TimeoutExpired as expired:
            record["exit_code"] = None
            record["status"] = "timeout"
            stdout = expired.stdout.decode() if isinstance(expired.stdout, bytes) else (expired.stdout or "")
    record["wall_seconds"] = round(time.monotonic() - started, 3)
    # RUSAGE_CHILDREN is a running maximum over all finished children; runs are
    # sequential, so this is an upper bound for this run.
    record["peak_rss_bytes_upper_bound"] = children_maxrss_bytes()
    try:
        response = json.loads(stdout) if stdout.strip() else None
    except json.JSONDecodeError:
        response = None
        record["stdout_tail"] = stdout[-2000:]
    if response is not None:
        response.pop("ibp_statistics", None)  # large; the full log keeps the details
        record["response"] = response
        record.setdefault("status", response.get("status", "unknown"))
        if response.get("status") == "error":
            record["error"] = response.get("message")
    else:
        record.setdefault("status", "crashed")
    if record["status"] != "completed":
        tail = log_path.read_text(errors="replace").splitlines()[-15:]
        record["log_tail"] = tail
    (RESULTS / f"{name}.json").write_text(json.dumps(record, indent=1) + "\n")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("families", nargs="*", help="family names (default: all)")
    parser.add_argument("--diffexp", help="path to the diffexp executable (or set DIFFEXP)")
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    parser.add_argument("--timeout", type=int, default=1800, help="wall-clock limit per family [s]")
    parser.add_argument("--prep-seconds", type=int, default=1800, help="IBP preparation budget [s]")
    parser.add_argument("--level-seconds", type=int, default=1200)
    parser.add_argument("--fire-seconds", type=int, default=600)
    parser.add_argument("--cache", type=Path,
                        default=Path(os.environ.get("DIFFEXP_CACHE", Path.home() / ".cache" / "diffexp-literature-2026")))
    parser.add_argument("--skip-done", action="store_true", help="skip families with a completed result")
    parser.add_argument("--jobs", type=int, default=1,
                        help="families run at the same time (each with --threads threads)")
    parser.add_argument("--extra", nargs=argparse.REMAINDER, default=[],
                        help="further diffexp arguments, e.g. --extra --epsilon-strategy sample")
    args = parser.parse_args()

    diffexp = find_diffexp(args.diffexp)
    RESULTS.mkdir(exist_ok=True)
    args.cache.mkdir(parents=True, exist_ok=True)
    version = subprocess.run([diffexp, "--version"], capture_output=True, text=True).stdout.strip()
    print(f"{version}; {args.threads} threads; timeout {args.timeout} s", flush=True)

    pending = []
    for name in family_names(args.families):
        previous = RESULTS / f"{name}.json"
        if args.skip_done and previous.exists() and json.loads(previous.read_text()).get("status") == "completed":
            print(f"{name}: skipped (completed earlier)", flush=True)
            continue
        pending.append(name)

    def report(name):
        print(f"{name}: running ...", flush=True)
        record = run_one(diffexp, name, args)
        detail = record.get("error") or ""
        print(f"{name}: {record['status']} in {record['wall_seconds']} s {detail}", flush=True)

    if args.jobs <= 1:
        for name in pending:
            report(name)
    else:
        # Families run side by side, each with --threads workers.  The memory
        # column is then a shared upper bound, not a per-family figure.
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            list(pool.map(report, pending))


if __name__ == "__main__":
    main()
