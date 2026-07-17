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
    "DiffExp2/NativeTransport.m",
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
    "FeynmanTrick/FamilySpec.m",
    "FeynmanTrick/FIREInterface.m",
    "FeynmanTrick/FIRE7Runner.m",
    "FeynmanTrick/FeynmanTrick.m",
    "FeynmanTrick/FeynmanTrickIteration.m",
    "FeynmanTrick/LevelReduction.m",
    "FeynmanTrick/MatrixExport.m",
    "FeynmanTrick/PipelineRequest.m",
    "FeynmanTrick/PropagatorAlgebra.m",
}

CPP_FILES = {
    "cpp/CMakeLists.txt",
    "cpp/include/diffexp2/checkpoint.hpp",
    "cpp/include/diffexp2/integration.hpp",
    "cpp/include/diffexp2/immutable_recursive_cache.hpp",
    "cpp/include/diffexp2/json_codec.hpp",
    "cpp/include/diffexp2/line_integration.hpp",
    "cpp/include/diffexp2/local_algebra.hpp",
    "cpp/include/diffexp2/local_solution.hpp",
    "cpp/include/diffexp2/matching.hpp",
    "cpp/include/diffexp2/path_planner.hpp",
    "cpp/include/diffexp2/physical_ode.hpp",
    "cpp/include/diffexp2/recurrence.hpp",
    "cpp/include/diffexp2/scc_completeness.hpp",
    "cpp/include/diffexp2/scalar.hpp",
    "cpp/include/diffexp2/series_types.hpp",
    "cpp/include/diffexp2/singular_indicial.hpp",
    "cpp/include/diffexp2/tail_majorant.hpp",
    "cpp/src/checkpoint.cpp",
    "cpp/src/json_codec.cpp",
    "cpp/src/json_codec_backend_info.cpp",
    "cpp/src/json_codec_capabilities.hpp",
    "cpp/src/json_codec_chart_scc_state.ipp",
    "cpp/src/json_codec_checkpoint_objects.ipp",
    "cpp/src/json_codec_checkpoint_restore.ipp",
    "cpp/src/json_codec_checkpoint_snapshot.ipp",
    "cpp/src/json_codec_checkpoint_state.ipp",
    "cpp/src/json_codec_commands_integration.ipp",
    "cpp/src/json_codec_commands_session_tile.ipp",
    "cpp/src/json_codec_commands_solve.ipp",
    "cpp/src/json_codec_commands_transport.ipp",
    "cpp/src/json_codec_endpoint_state.ipp",
    "cpp/src/json_codec_local_state.ipp",
    "cpp/src/json_codec_matching_restore.ipp",
    "cpp/src/json_codec_matching_runtime.ipp",
    "cpp/src/json_codec_matching_state.ipp",
    "cpp/src/json_codec_multiplier_state.ipp",
    "cpp/src/json_codec_recurrence.cpp",
    "cpp/src/json_codec_recurrence.hpp",
    "cpp/src/json_codec_scc.cpp",
    "cpp/src/json_codec_scc.hpp",
    "cpp/src/json_codec_session_state.ipp",
    "cpp/src/json_codec_support.cpp",
    "cpp/src/json_codec_support.hpp",
    "cpp/src/json_codec_transport_algorithms.ipp",
    "cpp/src/json_codec_transport_state.ipp",
    "cpp/src/librarylink.cpp",
    "cpp/tests/test_integration.cpp",
    "cpp/tests/test_json_codec.cpp",
    "cpp/tests/test_line_integration.cpp",
    "cpp/tests/test_local_algebra.cpp",
    "cpp/tests/test_local_solution.cpp",
    "cpp/tests/test_matching.cpp",
    "cpp/tests/test_path_planner.cpp",
    "cpp/tests/test_physical_ode_evolution.cpp",
    "cpp/tests/test_persistent_acb_local_match.cpp",
    "cpp/tests/test_persistent_algebraic_tile_plan.cpp",
    "cpp/tests/test_persistent_analytic_prescription.cpp",
    "cpp/tests/test_persistent_checkpoint.cpp",
    "cpp/tests/test_persistent_checkpoint_deferred_state.cpp",
    "cpp/tests/test_persistent_endpoint_limit.cpp",
    "cpp/tests/test_persistent_local.cpp",
    "cpp/tests/test_persistent_local_match.cpp",
    "cpp/tests/test_persistent_local_solve_batch.cpp",
    "cpp/tests/test_persistent_matched_tail_majorant.cpp",
    "cpp/tests/test_persistent_owner_bound_residual.cpp",
    "cpp/tests/test_persistent_regular_equation_owner.cpp",
    "cpp/tests/test_persistent_parallel_arms.cpp",
    "cpp/tests/test_persistent_phase_planned_match.cpp",
    "cpp/tests/test_persistent_plan_match_hop.cpp",
    "cpp/tests/test_persistent_planned_endpoint.cpp",
    "cpp/tests/test_persistent_rational_row.cpp",
    "cpp/tests/test_persistent_rational_row_tail.cpp",
    "cpp/tests/test_persistent_scc_owner_bound_residual.cpp",
    "cpp/tests/test_persistent_scc_plan_match.cpp",
    "cpp/tests/test_persistent_scc_prepare.cpp",
    "cpp/tests/test_persistent_scc_pseudo_compensation.cpp",
    "cpp/tests/test_persistent_scc_solve_block_column.cpp",
    "cpp/tests/test_persistent_scc_solve_column.cpp",
    "cpp/tests/test_persistent_scc_solve_singular_column.cpp",
    "cpp/tests/test_persistent_scc_spectral_frames.cpp",
    "cpp/tests/test_persistent_single_arm_tile_plan.cpp",
    "cpp/tests/test_persistent_singular_scc_arms.cpp",
    "cpp/tests/test_persistent_tail_majorant.cpp",
    "cpp/tests/test_persistent_tile_line.cpp",
    "cpp/tests/test_persistent_transport_contract_pair.cpp",
    "cpp/tests/test_persistent_transport_endpoint_batch.cpp",
    "cpp/tests/test_persistent_transport_run_arms.cpp",
    "cpp/tests/test_recurrence.cpp",
    "cpp/tests/test_scc_formal_completeness.cpp",
    "cpp/tests/test_session_solve_many.cpp",
    "cpp/tests/test_singular_indicial.cpp",
    "cpp/tests/test_tail_majorant.cpp",
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
    "Docs/AlgebraicSquareRootPlan.md",
    "Docs/AnalyticContinuation.md",
    "Docs/AsymptoticBoundaryPlan.md",
    "Docs/Citation.md",
    "Docs/CppBackend.md",
    "Docs/DirectSolver.md",
    "Docs/FeynmanTrick.md",
    "Docs/Installation.md",
    "Docs/Migration.md",
    "Docs/PerformanceRegression.md",
    "Docs/QuickStart.md",
    "Docs/ReleaseManifest.md",
    "Docs/Results.md",
}

SCRIPT_FILES = {
    "Scripts/FTExamples.m",
    "Scripts/banana4_bessel_oracle.m",
    "Scripts/check_release_tree.py",
    "Scripts/diagnose_epsilon_lattice.m",
    "Scripts/pentagon_massive_oracle.m",
    "Scripts/run_banana_timing_regression.sh",
    "Scripts/run_double_box_timing_regression.sh",
    "Scripts/run_pentagon_timing_regression.sh",
    "Scripts/run_ft_stepwise2.m",
    "Scripts/run_release_tests.sh",
    "Scripts/verify_double_box_planar.m",
}

TEST_FILES = {
    "Tests/benchmark_banana_singular.m",
    "Tests/test_api.m",
    "Tests/test_banana4_bessel_oracle.m",
    "Tests/test_config.m",
    "Tests/test_cpp_arm_batch.m",
    "Tests/test_cpp_backend.m",
    "Tests/test_double_box_planar_oracle.m",
    "Tests/test_eps_series.m",
    "Tests/test_epsilon_regular_banana.m",
    "Tests/test_feynmantrick_algebra.m",
    "Tests/test_fast_tadpole_boundary.m",
    "Tests/test_feynmantrick_failure_semantics.m",
    "Tests/test_fire7_classical_exact_output.m",
    "Tests/test_fire7_lock_cache_contract.m",
    "Tests/test_fire7_modular_contract.m",
    "Tests/test_fire7_pipeline_plumbing.m",
    "Tests/test_fire7_rational_reconstruction.m",
    "Tests/test_fire7_real_bubble_parity.m",
    "Tests/test_fire7_real_sunrise_parity.m",
    "Tests/test_fire_inmemory_reduction_cache.m",
    "Tests/test_fire_level_reduction_batch.m",
    "Tests/test_feynmantrick_fire.m",
    "Tests/test_feynmantrick_iteration.m",
    "Tests/test_feynmantrick_family_spec.m",
    "Tests/test_ft_example_specs.m",
    "Tests/test_ft_pipeline_facade.m",
    "Tests/test_ft_pipeline_family_request.m",
    "Tests/test_ft_pipeline_process.m",
    "Tests/test_ft_runner_family_request_contract.m",
    "Tests/test_ft_ladder_checkpointing.m",
    "Tests/test_ft_native_checkpoint_resume.m",
    "Tests/test_ft_native_epsilon_preplanner.m",
    "Tests/test_ft_native_observable_dispatch.m",
    "Tests/test_ft_prep_cache_contract.m",
    "Tests/test_indicial.m",
    "Tests/test_integrate.m",
    "Tests/test_native_algebraic_chart_bridge.m",
    "Tests/test_native_algebraic_scale_bridge.m",
    "Tests/test_native_regular_basis_batch.m",
    "Tests/test_native_regular_prepare_owner.m",
    "Tests/test_native_checkpoint_single_visibility.m",
    "Tests/test_native_mixed_regular_singular_atlas.m",
    "Tests/test_native_rational_row_tail_payload.m",
    "Tests/test_native_regular_independent_arms.m",
    "Tests/test_native_scc_acb_regular_column.m",
    "Tests/test_native_scc_acb_polar_shadow.m",
    "Tests/test_native_scc_cache_reservation.m",
    "Tests/test_native_scc_pseudo_compensation.m",
    "Tests/test_native_scc_single_block_singular.m",
    "Tests/test_native_scc_singular_classifier.m",
    "Tests/test_native_scc_singular_prepare.m",
    "Tests/test_native_scc_solve_block_column.m",
    "Tests/test_native_scc_solve_column.m",
    "Tests/test_native_scc_spectral_frame_certificate.m",
    "Tests/test_native_scc_zero_polar_pseudo_diagnostics.m",
    "Tests/test_native_transport_observable_batch.m",
    "Tests/test_path_planner_algebraic.m",
    "Tests/test_pentagon_massive_oracle.m",
    "Tests/test_paclet_metadata.m",
    "Tests/test_public_api.m",
    "Tests/test_public_native_line_dispatch.m",
    "Tests/test_persistent_transport_contract_bridge.m",
    "Tests/test_persistent_transport_contract_pair_bridge.m",
    "Tests/test_persistent_transport_endpoint_batch_bridge.m",
    "Tests/test_persistent_transport_run_arms_bridge.m",
    "Tests/test_persistent_transport_pair_stream_wolfram.m",
    "Tests/test_persistent_weighted_planned_endpoint.m",
    "Tests/test_physical_ode_payload.m",
    "Tests/test_scc_integration_sequence.m",
    "Tests/test_scc_multiplier_preparation.m",
    "Tests/test_scc_physical_ode_manifest.m",
    "Tests/test_scc_structure_cache.m",
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
