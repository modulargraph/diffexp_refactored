#!/usr/bin/env python3
"""Fail when a Git checkout or source archive differs from the release allowlist."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

ROOT_FILES = {
    ".gitignore",
    "CITATION.cff",
    "CHANGELOG.md",
    "CMakeLists.txt",
    "DiffExp2.m",
    "FeynmanTrick.m",
    "LICENSE",
    "PacletInfo.wl",
    "README.md",
}

DIFFEXP2_FILES = {
    "DiffExp2/API.m",
    "DiffExp2/Config.m",
    "DiffExp2/CppBackend.m",
    "DiffExp2/DiffExp2.m",
    "DiffExp2/EpsSeries.m",
    "DiffExp2/Indicial.m",
    "DiffExp2/Integrate.m",
    "DiffExp2/PublicAPI.m",
    "DiffExp2/SectorSeries.m",
    "DiffExp2/Solve.m",
    "DiffExp2/Tolerances.m",
    "DiffExp2/Transport.m",
}

FEYNMAN_TRICK_FILES = {
    "FeynmanTrick/BoundaryConditions.m",
    "FeynmanTrick/DiffExp2Pipeline.m",
    "FeynmanTrick/EpsPrefactors.m",
    "FeynmanTrick/FIREInterface.m",
    "FeynmanTrick/FeynmanTrick.m",
    "FeynmanTrick/FeynmanTrickIteration.m",
    "FeynmanTrick/LevelReduction.m",
    "FeynmanTrick/MatrixExport.m",
    "FeynmanTrick/PropagatorAlgebra.m",
}

CPP_FILES = {
    "cpp/CMakeLists.txt",
    "cpp/include/diffexp2/json_codec.hpp",
    "cpp/include/diffexp2/recurrence.hpp",
    "cpp/include/diffexp2/scalar.hpp",
    "cpp/src/json_codec.cpp",
    "cpp/src/librarylink.cpp",
    "cpp/tests/test_recurrence.cpp",
}

EXAMPLE_FILES = {
    "Examples/Direct/MinimalTransport.wl",
    "Examples/Direct/SingularEndpointAndSegments.wl",
    "Examples/FeynmanTrick/BoxBubble.sh",
    "Examples/FeynmanTrick/Bubble.sh",
    "Examples/FeynmanTrick/FourLoopUnequalBanana.sh",
    "Examples/FeynmanTrick/Kite.sh",
    "Examples/FeynmanTrick/README.md",
    "Examples/FeynmanTrick/Sunrise.sh",
    "Examples/FeynmanTrick/UnequalBanana.sh",
}

DOC_FILES = {
    "Docs/API.md",
    "Docs/AnalyticContinuation.md",
    "Docs/Citation.md",
    "Docs/CppBackend.md",
    "Docs/DirectSolver.md",
    "Docs/FeynmanTrick.md",
    "Docs/Installation.md",
    "Docs/Migration.md",
    "Docs/QuickStart.md",
    "Docs/ReleaseManifest.md",
    "Docs/Results.md",
}

SCRIPT_FILES = {
    "Scripts/FTExamples.m",
    "Scripts/banana4_bessel_oracle.m",
    "Scripts/check_release_tree.py",
    "Scripts/run_ft_stepwise2.m",
    "Scripts/run_release_tests.sh",
}

TEST_FILES = {
    "Tests/test_api.m",
    "Tests/test_banana4_bessel_oracle.m",
    "Tests/test_config.m",
    "Tests/test_cpp_arm_batch.m",
    "Tests/test_cpp_backend.m",
    "Tests/test_eps_series.m",
    "Tests/test_feynmantrick_algebra.m",
    "Tests/test_fast_tadpole_boundary.m",
    "Tests/test_feynmantrick_failure_semantics.m",
    "Tests/test_fire_inmemory_reduction_cache.m",
    "Tests/test_fire_level_reduction_batch.m",
    "Tests/test_feynmantrick_fire.m",
    "Tests/test_feynmantrick_iteration.m",
    "Tests/test_ft_example_specs.m",
    "Tests/test_ft_pipeline_facade.m",
    "Tests/test_ft_pipeline_process.m",
    "Tests/test_ft_ladder_checkpointing.m",
    "Tests/test_indicial.m",
    "Tests/test_integrate.m",
    "Tests/test_path_planner_algebraic.m",
    "Tests/test_paclet_metadata.m",
    "Tests/test_public_api.m",
    "Tests/test_sectorseries.m",
    "Tests/test_solve.m",
    "Tests/test_tolerances.m",
    "Tests/test_transport.m",
}

FIXTURE_FILES = {
    "Tests/refs/bench/banana_L1.m",
    "Tests/refs/oracle_logs/MANIFEST.md",
    "Tests/refs/oracle_logs/l2_banana.log",
    "Tests/refs/oracle_logs/l2_bubsun.log",
}

EXACT_FILES = (
    ROOT_FILES
    | DIFFEXP2_FILES
    | FEYNMAN_TRICK_FILES
    | CPP_FILES
    | EXAMPLE_FILES
    | DOC_FILES
    | SCRIPT_FILES
    | TEST_FILES
    | FIXTURE_FILES
)
REQUIRED_FILES = EXACT_FILES

# A source archive has no Git index. These paths are never release payload;
# they may nevertheless exist after an out-of-source-style local build or
# after the user installs the external FIRE dependency beside the sources.
ARCHIVE_IGNORED_TOP_LEVEL = {".git", "build", "Dependencies"}
ARCHIVE_IGNORED_PARTS = {"__pycache__"}
ARCHIVE_IGNORED_SUFFIXES = {".pyc", ".pyo"}


def archive_files() -> set[str]:
    files: set[str] = set()
    for path in ROOT.rglob("*"):
        if not (path.is_file() or path.is_symlink()):
            continue
        relative = path.relative_to(ROOT)
        if (relative.parts[0] in ARCHIVE_IGNORED_TOP_LEVEL or
                any(part in ARCHIVE_IGNORED_PARTS for part in relative.parts) or
                path.suffix in ARCHIVE_IGNORED_SUFFIXES):
            continue
        files.add(relative.as_posix())
    return files


def tracked_files() -> set[str]:
    probe = subprocess.run(
        ["git", "rev-parse", "--is-inside-work-tree"],
        cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        check=False,
    )
    if probe.returncode == 0:
        raw = subprocess.check_output(
            ["git", "ls-files", "-z"], cwd=ROOT
        ).decode("utf-8")
        return {item for item in raw.split("\0") if item}
    return archive_files()


def allowed(path: str) -> bool:
    return path in EXACT_FILES


def main() -> int:
    tracked = tracked_files()
    unexpected = sorted(path for path in tracked if not allowed(path))
    missing = sorted(REQUIRED_FILES - tracked)

    if unexpected:
        print("Unexpected tracked release files:", file=sys.stderr)
        for path in unexpected:
            print(f"  {path}", file=sys.stderr)
    if missing:
        print("Missing required release files:", file=sys.stderr)
        for path in missing:
            print(f"  {path}", file=sys.stderr)

    if unexpected or missing:
        return 1
    print(f"Release tree allowlist passed ({len(tracked)} source files).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
