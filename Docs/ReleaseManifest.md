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
3. Current runtime dependencies of the Feynman-trick CLI remain in the release
   until a stable API replaces them.
4. No history is destroyed.  “Archive” below means preserved on the
   development branch/repository and optionally tagged, not discarded.

## Clean release allowlist

| Path | Disposition | Notes |
| --- | --- | --- |
| ``README.md`` | release | replace experimental landing page with the DiffExp 2 README in this prototype |
| ``DiffExp2/`` | release | recurrence solver, transport, integration, API, configuration, C++ bridge |
| ``cpp/`` | release | C++20/FLINT recurrence engine and native unit tests |
| ``CMakeLists.txt`` | release | top-level compiled-backend build |
| ``FeynmanTrick.m`` and ``FeynmanTrick/*.m`` | release | topology, FIRE, iteration, exact matrix, boundary, and bridge code |
| ``Examples/`` | release | curated direct and Feynman-trick examples from this prototype |
| selected ``Docs/*.md`` | release | the user-doc allowlist in the next table |
| ``Scripts/run_ft_stepwise2.m`` | release, transitional | implemented sector-native Feynman-trick CLI |
| ``Scripts/FTExamples.m`` | release, transitional | exact built-in topology fixtures consumed by the CLI |
| ``Scripts/banana4_bessel_oracle.m`` | release example/oracle | small independent normalization check; move under ``Examples/Oracles`` eventually |
| focused ``Tests/*.m`` | release | package loading, config/API, recurrence backend, transport, integration, singular endpoints, and shipped-example smoke tests |
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

Keep a focused, runnable suite rather than every campaign artifact:

- package loading and configuration;
- epsilon-series and sector-series invariants;
- API mini integration;
- regular and singular recurrence, resonance/log sectors, and no fallback;
- transport geometry, analytic continuation, endpoint limits, and
  regularized integration;
- C++ backend availability/parity when ``DE2_REQUIRE_CPP=1``;
- one small Feynman-trick algebra/iteration test not requiring FIRE;
- optional integration tests clearly labeled as requiring FIRE;
- smoke tests for every shipped direct example and shell syntax checks for
  every Feynman-trick example.

Large benchmark fixtures and pySecDec generation environments remain on dev.
Curated numeric pins may be copied into small release fixtures with their
normalization and provenance recorded in ``Docs/Results.md``.

## Transitional runtime files

The cleanest long-term layout would move the current runner and example
registry behind a public package function.  Until that work is complete,
removing ``Scripts/run_ft_stepwise2.m`` or ``Scripts/FTExamples.m`` would break
the documented Feynman-trick workflow.  They are therefore release files for
the first cut, despite living under ``Scripts/``.

The same rule applies to FIRE: its source is not vendored, but configuration
and error messages must make the external dependency explicit.

## Release blockers exposed by this prototype

1. Change the configuration schema default for ``RecurrenceBackend`` from
   ``"Wolfram"`` to ``"Cpp"`` and update tests that intentionally exercise the
   reference backend.
2. Replace or retire the root ``DiffExp.m`` loader so it cannot silently load
   DiffExp 1 in a DiffExp 2 release.
3. Decide whether the first release officially supports the transitional CLI
   or must first expose a stable Wolfram Language Feynman-trick facade.
4. Add a license file, citation metadata, release version, and supported
   Wolfram-version statement.  None is safely inferable from the snapshot.
5. Give the package a clean-install test outside the source checkout, including
   compiled-library discovery.
6. Run every direct example and at least the small FIRE-enabled bubble example
   on the exact release commit.
7. Re-run static documentation link checks after the disposition moves; links
   to development-only validation records must be replaced by curated release
   pins or stable external URLs.
8. Keep the root-in-basis limitation phrased precisely: line-variable
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
