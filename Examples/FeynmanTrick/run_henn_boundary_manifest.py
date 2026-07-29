#!/usr/bin/env python3
"""Execute validated Henn FT runner commands without nesting Wolfram kernels."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path


SCHEMA = "DiffExp2.HennDoublePentagonBoundaryRun/v1"


def fail(message: str) -> int:
    print(message, file=sys.stderr)
    return 2


def main() -> int:
    if len(sys.argv) != 2:
        return fail(f"usage: {Path(sys.argv[0]).name} manifest.json")

    manifest_path = Path(sys.argv[1]).expanduser().resolve()
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return fail(f"cannot read Henn FT manifest {manifest_path}: {exc}")

    if manifest.get("Schema") != SCHEMA:
        return fail(f"unexpected Henn FT manifest schema in {manifest_path}")
    runs = manifest.get("Runs")
    if not isinstance(runs, list) or not runs:
        return fail(f"Henn FT manifest has no runs: {manifest_path}")

    for expected_ordinal, run in enumerate(runs, start=1):
        if not isinstance(run, dict) or run.get("Ordinal") != expected_ordinal:
            return fail("Henn FT manifest run order is malformed")
        command = run.get("Command")
        workdir = run.get("WorkingDirectory")
        assignments = run.get("Environment")
        log_value = run.get("LogFile")
        if (
            not isinstance(command, list)
            or not command
            or not all(isinstance(item, str) and item for item in command)
            or not isinstance(workdir, str)
            or not isinstance(assignments, dict)
            or not isinstance(log_value, str)
            or not log_value
            or not all(
                isinstance(key, str) and isinstance(value, str)
                for key, value in assignments.items()
            )
        ):
            return fail(f"Henn FT run {expected_ordinal} is malformed")
        log_path = Path(log_value).expanduser().resolve()

        log_path.parent.mkdir(parents=True, exist_ok=True)
        environment = os.environ.copy()
        environment.update(assignments)
        print(
            f"=== Henn FT runner {expected_ordinal}/{len(runs)} "
            f"log={log_path}",
            flush=True,
        )
        process: subprocess.Popen[str] | None = None
        try:
            with log_path.open("w", encoding="utf-8") as log:
                process = subprocess.Popen(
                    command,
                    cwd=workdir,
                    env=environment,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                )
                assert process.stdout is not None
                for line in process.stdout:
                    sys.stdout.write(line)
                    sys.stdout.flush()
                    log.write(line)
                    log.flush()
                status = process.wait()
        except KeyboardInterrupt:
            if process is not None and process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            print(
                f"Henn FT run {expected_ordinal} interrupted; "
                f"see {log_path}",
                file=sys.stderr,
            )
            return 130
        except OSError as exc:
            return fail(f"could not execute Henn FT run {expected_ordinal}: {exc}")
        if status != 0:
            print(
                f"Henn FT run {expected_ordinal} failed with status {status}; "
                f"see {log_path}",
                file=sys.stderr,
            )
            return status if status > 0 else 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
