# Release File-Disposition Manifest

This is a proposed clean-master manifest.  It does not delete, rename, or move
any existing file in this prototype branch.  The full development snapshot
must be preserved first on a dedicated development repository or branch; only
then should a release branch be assembled from the allowlist below.

## Principles

1. The release repository is a user product: DiffExp 2, the C++ recurrence
   backend, Feynman-trick support, concise user documentation, runnable
   examples, and focused tests.
2. Benchmark campaigns, design deliberations, agent reviews, large oracle
   dumps, legacy implementations, and debugging scripts remain available in
   the development repository without appearing as the product surface.
3. Current subprocess runtime dependencies of the public Feynman-trick facade
   remain in the release until that implementation moves fully into the
   package.
4. No history is destroyed.  “Archive” below means preserved on the
   development branch/repository and optionally tagged, not discarded.

## Clean release allowlist

| Path | Disposition | Notes |
| --- | --- | --- |
| ``README.md`` | release | replace experimental landing page with the DiffExp 2 README in this prototype |
| ``CITATION.cff`` | release | machine-readable software citation metadata |
| ``CHANGELOG.md`` | release | user-facing major-release changes and migration pointer |
| ``LICENSE`` | release | GNU GPL version 3; source notices permit version 3 or any later version |
| ``PacletInfo.wl`` | release | source-tree paclet metadata for the two root contexts and Wolfram 15.0 compatibility |
| ``DiffExp2.m`` | release | public root loader; installs the release facade and injects the C++ recurrence default |
| ``DiffExp2/`` | release | recurrence solver, transport, integration, API, configuration, C++ bridge |
| ``cpp/`` | release | C++20/FLINT recurrence engine and native unit tests |
| ``CMakeLists.txt`` | release | top-level compiled-backend build |
| ``FeynmanTrick.m`` and ``FeynmanTrick/*.m`` | release | topology, FIRE, iteration, exact matrix, boundary, and bridge code |
| ``Examples/`` | release | curated direct and Feynman-trick examples from this prototype |
| selected ``Docs/*.md`` | release | the user-doc allowlist in the next table |
| ``Scripts/run_ft_stepwise2.m`` | release, transitional | runtime dependency of the public facade's current subprocess implementation |
| ``Scripts/FTExamples.m`` | release, transitional | built-in topology registry consumed by the facade runtime |
| ``Scripts/banana4_bessel_oracle.m`` | release example/oracle | small independent normalization check; move under ``Examples/Oracles`` eventually |
| ``Scripts/run_release_tests.sh`` | release | exact sequential release verification, including direct-example smoke tests |
| ``Scripts/check_release_tree.py`` | release | machine-checks the exact tracked-file allowlist on the clean branch |
| focused ``Tests/*.m`` | release | exact first-release test allowlist is below |
| ``Tests/refs/bench/banana_L1.m`` | release test fixture | sole stored matrix fixture required by the focused C++ test allowlist |
| ``Tests/refs/oracle_logs/{MANIFEST.md,l2_bubsun.log,l2_banana.log}`` | release result provenance | exact files linked from ``Docs/Results.md``; they are documentation dependencies, not test fixtures |
| ``.gitignore`` | release | remove only entries that no longer apply after cleanup |

### User-document allowlist

Keep these as the release documentation surface:

```text
Docs/Installation.md
Docs/QuickStart.md
Docs/DirectSolver.md
Docs/FeynmanTrick.md
Docs/AnalyticContinuation.md
Docs/Results.md
Docs/API.md
Docs/Migration.md
Docs/CppBackend.md
Docs/ReleaseManifest.md
Docs/Citation.md
```

``Docs/CppBackend.md`` should receive one final editorial pass to remove
development-only microbenchmark detail or move that detail to the development
site while retaining build, correctness, scope, and troubleshooting sections.

## Developer-only maintained material

These files remain active engineering material in the development repository,
but should not appear on clean master:

| Path or pattern | Why it is developer-only |
| --- | --- |
| ``AGENTS.md`` | repository-agent instructions, not user documentation |
| ``Docs/specs/`` | binding implementation specifications and milestone decisions |
| ``Docs/reviews/`` | agent audits, triage notes, and review artifacts |
| ``Docs/RewritePlan.md`` | rewrite campaign plan |
| ``Docs/LessonsLedger.md`` | internal numerical-wisdom ledger |
| ``Docs/PerfGapAnalysis.md`` and ``Docs/SpeedIdeas.md`` | optimization research |
| ``Docs/BoxPreflight.md`` | campaign preflight analysis |
| ``Docs/FeynmanTrick*Status.md`` | historical status/debug narratives; distill verified outcomes into ``Results.md`` |
| ``Docs/recurrence_solver_prompt.md`` | development prompt, not product documentation |
| ``Docs/CoreModules.md``, ``Infrastructure.md``, ``SeriesAndEvaluation.md``, ``Transformations.md``, ``TransportAndStrategies.md``, ``SingularityDecomposition.md``, ``RegularizedIntegration.md`` | documentation of the legacy/refactored old core; archive or relabel on dev |
| ``Scripts/bench_*.m`` and ``Scripts/profile_*.m`` | performance campaign harnesses |
| ``Scripts/box_preflight_*.m`` and ``Scripts/check_transport_ode_residual.m`` | targeted diagnostics |
| ``Scripts/compare_*.py``, ``dump_*.m``, ``eval_dump_generic.m`` | replay/debug tooling |
| ``Scripts/export_pysecdec_*.m``, ``pysecdec_*.py``, ``box_verdict.py`` | external-oracle generation campaign |
| ``Scripts/gen_bench_fixtures.m`` and grouped-transform/adaptive-frame benches | fixture generation and solver research |
| large or campaign-specific ``Tests/refs/`` content | development provenance; keep only minimal release pins needed by retained tests |
| ``Tests/PINS.md`` and ``Tests/SINGULARITY_DECOMPOSITION_PROBLEM.md`` | internal campaign records |
| experimental topology tests not represented as supported examples | retain on dev until promoted |

Developer-only does not mean unimportant.  In particular, the specs, review
findings, and oracle generators are needed to maintain the numerical contract;
they simply belong in the development repository rather than the user-facing
release tree.

## Legacy/archive material

Preserve these paths immutably on the development repository or a named
archive tag:

| Path | Archive reason |
| --- | --- |
| ``DiffExp/`` | modularized DiffExp 1/reference implementation; not part of the DiffExp 2 product |
| root ``DiffExp.m`` in its current form | currently loads the old package; replace with a DiffExp 2 loader on release |
| ``Reference/`` | original monolithic source and historical examples |
| ``Docs/ExportDisposition.md`` | exhaustive old-export migration ledger |
| legacy ``Scripts/run_ft_stepwise.m`` | old-core Feynman-trick runner |
| old matrix-slice fixtures used only for parity | exact full-matrix format is the release input contract |
| obsolete campaign logs and failed-run dumps | retain for provenance, not as shipped examples |

## Tests proposed for clean master

Keep this exact first-release suite rather than every campaign artifact.

Core DiffExp 2 and public-API tests:

```text
Tests/test_tolerances.m
Tests/test_config.m
Tests/test_eps_series.m
Tests/test_indicial.m
Tests/test_sectorseries.m
Tests/test_solve.m
Tests/test_transport.m
Tests/test_path_planner_algebraic.m
Tests/test_integrate.m
Tests/test_api.m
Tests/test_public_api.m
Tests/test_paclet_metadata.m
```

Compiled-backend tests (set ``DE2_REQUIRE_CPP=1`` in the release job):

```text
Tests/test_cpp_backend.m
Tests/test_cpp_arm_batch.m
```

Process-free or fake-process Feynman-trick and shipped-oracle tests:

```text
Tests/test_feynmantrick_algebra.m
Tests/test_ft_example_specs.m
Tests/test_feynmantrick_failure_semantics.m
Tests/test_fire_inmemory_reduction_cache.m
Tests/test_fire_level_reduction_batch.m
Tests/test_ft_pipeline_facade.m
Tests/test_ft_pipeline_process.m
Tests/test_ft_ladder_checkpointing.m
Tests/test_banana4_bessel_oracle.m
```

Optional FIRE-enabled integration tests, run only in a job with a configured
FIRE installation:

```text
Tests/test_feynmantrick_fire.m
Tests/test_feynmantrick_iteration.m
```

The corresponding stored-file allowlist is deliberately smaller than
``Tests/refs/``:

```text
Tests/refs/bench/banana_L1.m
Tests/refs/oracle_logs/MANIFEST.md
Tests/refs/oracle_logs/l2_bubsun.log
Tests/refs/oracle_logs/l2_banana.log
```

``banana_L1.m`` is the only fixture loaded by the focused tests (through
``test_cpp_backend.m``).  The three ``oracle_logs`` entries are retained
because ``Docs/Results.md`` links to them as result provenance; no selected
test loads them.  In particular, this does not allowlist the rest of
``Tests/refs/``, its top-level ``MANIFEST.md``, the large oracle dumps, or the
pySecDec generation environment.

The current ``Scripts/run_test_battery.sh`` is a development runner and still
includes legacy-core tests. ``Scripts/run_release_tests.sh`` is the clean,
sequential runner for the exact required group above, the two direct examples,
and shell syntax checks. Set ``DE2_RUN_FIRE_TESTS=1`` to add the optional FIRE
group.

## Transitional runtime files

The public facade is now implemented as ``FeynmanTrick`PipelinePlan``,
``FeynmanTrick`RunIntegrationPipeline``, and
``FeynmanTrick`ResumeIntegrationPipeline``.  Its current implementation starts
``Scripts/run_ft_stepwise2.m`` as a subprocess, and that runner loads
``Scripts/FTExamples.m``.  Both scripts are therefore facade runtime
dependencies and release files for the first cut, despite living under
``Scripts/``.  They may move into the package later without changing the
public facade.

The same rule applies to FIRE: its source is not vendored, but configuration
and error messages must make the external dependency explicit.

## Release readiness exposed by this prototype

Resolved in the prototype:

- The root ``DiffExp2.m`` public loader injects ``RecurrenceBackend -> "Cpp"``
  through the public configuration facade.  The low-level configuration
  schema intentionally keeps ``"Wolfram"`` as its developer/reference default;
  ``Tests/test_public_api.m`` covers the release-facing default.
- The first-release Feynman-trick interface is the public facade named above.
  Its two transitional script dependencies are explicitly allowlisted until
  their implementation moves into the package.

Remaining release blockers:

1. Replace or retire the root ``DiffExp.m`` loader so it cannot silently load
   DiffExp 1 in a DiffExp 2 release.
2. Set the final release version/date in ``CITATION.cff``, ``PacletInfo.wl``,
   ``$DiffExp2Version``, and ``$FeynmanTrickVersion``. The first release
   targets Wolfram Language 15.0 and retains DiffExp's GPL-3.0-or-later terms.
3. Give the package a clean-install test outside the source checkout, including
   compiled-library discovery.
4. Run ``Scripts/run_release_tests.sh`` and at least the small FIRE-enabled
   bubble example on the exact release commit.
5. Re-run static documentation link checks after the disposition moves; links
   to development-only validation records must be replaced by curated release
   pins or stable external URLs.
6. Keep the root-in-basis limitation phrased precisely: line-variable
   non-integer powers in the input matrix are rejected; algebraic singular
   points and local exponents are not categorically unsupported.

## Suggested branch/repository sequence

1. Tag and push the full known-good development snapshot.
2. Create the long-lived development repository/branch and verify all research
   artifacts are present there.
3. Create clean master from the same known-good commit.
4. Apply the release allowlist without rewriting development history.
5. Fix the blockers above and run the clean-install/test matrix.
6. Tag the release only after documentation claims and implementation defaults
   agree.

This prototype intentionally stops before step 4: it supplies the proposed
release files and the disposition decision, but performs no destructive tree
cleanup.
