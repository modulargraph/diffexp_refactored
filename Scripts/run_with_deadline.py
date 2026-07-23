#!/usr/bin/env python3
"""Run a command with a hard wall-clock deadline and tee its combined log."""

from __future__ import annotations

import os
import pathlib
import selectors
import signal
import subprocess
import sys
import time


def usage() -> None:
    raise SystemExit(
        "usage: run_with_deadline.py SECONDS LOG_PATH -- COMMAND [ARG ...]"
    )


def terminate_process_group(process: subprocess.Popen[bytes]) -> None:
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=3)
        return
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait()


def main(argv: list[str]) -> int:
    if len(argv) < 5 or argv[3] != "--":
        usage()
    try:
        seconds = int(argv[1])
    except ValueError:
        usage()
    if seconds < 1:
        usage()
    log_path = pathlib.Path(argv[2])
    command = argv[4:]
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("wb") as log:
        process = subprocess.Popen(
            command,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            start_new_session=True,
            bufsize=0,
        )
        assert process.stdout is not None
        selector = selectors.DefaultSelector()
        selector.register(process.stdout, selectors.EVENT_READ)
        deadline = time.monotonic() + seconds
        timed_out = False
        try:
            while selector.get_map():
                remaining = deadline - time.monotonic()
                if remaining <= 0 and process.poll() is None:
                    timed_out = True
                    message = (
                        f"HARD DEADLINE EXCEEDED: {seconds}s; "
                        "terminating command process group\n"
                    ).encode()
                    os.write(sys.stderr.fileno(), message)
                    log.write(message)
                    log.flush()
                    terminate_process_group(process)
                events = selector.select(timeout=max(0.0, min(0.25, remaining)))
                for key, _ in events:
                    chunk = os.read(key.fd, 65536)
                    if chunk:
                        os.write(sys.stdout.fileno(), chunk)
                        log.write(chunk)
                        log.flush()
                    else:
                        selector.unregister(key.fileobj)
            if process.poll() is None:
                process.wait()
        except KeyboardInterrupt:
            terminate_process_group(process)
            return 130
        finally:
            selector.close()
        return 124 if timed_out else process.returncode


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
