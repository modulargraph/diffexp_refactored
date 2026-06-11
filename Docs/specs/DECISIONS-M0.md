# M0 orchestrator decisions (binding for the spec-amendment pass and M1+)

These resolve every cross-cutting conflict raised by REVIEW-math.md (D1-D8
blockers + majors) and REVIEW-minimalism.md (4 blockers + majors) and the
agents' open questions. Per-spec amendment agents apply these verbatim; a
spec may not contradict a decision.

DEC-1 ERROR PRIMITIVE. One library-wide helper, exported by Tolerances.m
  (bottom of the dependency order): DE2Error[id_String, payload_Association]
  prints a one-line summary and Throw[Failure["DiffExp2", payload], 
  "DiffExp2Error"]; the catch sits at every API.m entry point. Every module
  uses ONLY this idiom. payload always carries "ID", "Module", and whatever
  chart/sector/order context $ESErrorContext-style scoping provides.

DEC-2 LAURENT LEAD TOLERANCE. laurentLeadTol = the campaign-verified
  RELATIVE test with universal floor: zero iff |c| <= Max[10^(-ChopPrecision/2),
  10^-24] * scale. The 10^-24 floor is a hard constant (not configurable);
  the ambiguity band (DEC-3) does NOT apply at the floor. Band width
  recalibrated to 4 decades (the two documented populations are ~1e-29
  residues vs O(1) signal).

DEC-3 ZERO PREDICATE. Exactly ONE predicate library-wide: Tolerances.m's
  NumericallyZeroQ (ternary: True / False / loud ambiguity error within the
  4-decade band, exempting the 10^-24 floor). EpsSeries' ESCoeffZeroQ is a
  thin wrapper adding window context to the error payload - same semantics.

DEC-4 INTEGRATE CANCELLATION GATE. Relative, per-eps-order: the scale is
  the max |coefficient| of the merged combination at that eps order; the
  b=0 divergence error fires only if the offending coefficient exceeds
  laurentLeadTol relative to that scale (REVIEW-math D3 amendment verbatim).

DEC-5 CONFIG DEFAULTS. Self-consistent: AccuracyGoalValidate default False;
  AccuracyGoal "?" allowed when Validate is False; the E8 cross-field error
  fires only on the inconsistent combination. LoadConfiguration[{}] must
  succeed (REVIEW-math D4).

DEC-6 REGULAR-CHART SECTORS. Indicial returns the full d-dimensional basis
  description for a regular chart: one (a=0,b=0,p=0) sector FAMILY with
  d-dimensional coefficient space (REVIEW-math D5 amendment); Solve's basis
  completeness and Transport's T-1 assert match this.

DEC-7 CHART SYSTEM OWNERSHIP. Solve.m exports PrepareChart (assembles
  ChartSystem from Indicial output: V, VInv, J, collision data; +~40 lines
  in Solve's budget). Indicial.m exports EpsDegenerateFamilies (per-chart
  families with same a mod Z, distinct b, colliding eps=0 eigenvectors) -
  Transport's RecombineBasis is tag-driven off it (REVIEW-math D6+D8).

DEC-8 MATCHWEIGHTS WINDOWS. Storage kmin (=-p for log chains) and the value
  window MinPower are distinct; matching normalizes by the VALUE window
  (REVIEW-math D7 amendment verbatim).

DEC-9 BLOCK SOLVE + TWINDOW. Ratified: block-sequential SolveChart; NO
  TWindow coupling-depth degradation in the new core (recursion matrices
  are exact polynomials; the old degradation was a numeric-matrix artifact).
  TWindow stays in the object as the truncation-order record only.
  LessonsLedger entry for MaxCouplingOrder becomes "subsumed: exact
  polynomial recursion + ErrorEstimate covers it" (minimalism defect 3).

DEC-10 EXTRA SINGULAR FACTORS. "ExtraSingularFactors" is an option on BOTH
  TransportTo and IntegrateOverLine, threaded to Transport segmentation
  (closes gap G1; the FT pipeline cannot place charts without it).

DEC-11 ACCURACYGOALVALIDATE. Demoted to validation-only (Transport E11);
  the old adaptive expansion-order search is waived in Config's table with
  the reason "replaced by exact recursion + ErrorEstimate gate".

DEC-12 PADE ORDER. The reviewer-corrected order formula (fixing old
  Pade.m:35's off-by-one) is normative; SectorSeries test t14 is the gate.

DEC-13 EPS-LAURENT VALUE SHAPE. EpsSeries.md owns the canonical shape
  <|"EpsWindow" -> <|"Min","CompleteMax"|>, "Coeffs" -> ...|>; all modules
  adopt its key names verbatim. No alternative shapes anywhere.

DEC-14 SNAPTOL vs RANKTOL. Separate derivations (snapTol stricter) so the
  E7 gray zone is non-empty; exact values per the minimalism amendment.
  User-input target snapping uses an INPUT-scaled tolerance (a 1e-16-off
  user target must snap, not create a degenerate chart).

DEC-15 ABORT-ON-CONTINUATION-FAIL. Transport's final-chart-only semantics
  is normative across Transport/API/Config (REVIEW-math D14; pentagon
  triage does not override).

DEC-16 PRESCRIPTION ERROR TRIGGER. The missing-prescription loud error
  fires when the chart is multivalued AT ALL: b != 0 OR p > 0 OR
  Denominator[a] > 1 (pentagon-triage finding: the Kallen charts are
  fractional-a with b=0, p=0). Also: prescription-factor dedup must be
  SIGN-AWARE ({-1+x,+1} vs {1-x,+1} flip the implied i-delta side - the
  old deltaPrescriptionsForFactors dedup is sign-blind; new code errors on
  conflicting normalizations).

DEC-17 2F1 TEST VALUES. The reviewer-recomputed sector exponents from the
  committed dz_0.m residue ({0,-3/2} class, not {0,-1/2}) are normative for
  SU-03 and API test 12.

DEC-18 MOBIUS. Dropped from the new core entirely. The M4 banana classic
  parity oracle is re-baselined with UseMobius -> False (RoC chart
  rescaling KEPT - it is an affine rescaling, not a Mobius map). If the old
  code cannot complete the banana classic line without Mobius, the banana
  classic parity gate falls back to bubble/sunrise/2f1 lines plus the
  banana FT lines (which already run UseMobius -> False).

DEC-19 DEAD CONFIG KEYS. LogFile, Crosscheck* (grep-verified no consumers)
  dropped with dedicated migration errors; Config table records the waiver.

DEC-20 PARITY KEY NAMES. Old names "SeriesValues"/"KinematicPoint" kept
  verbatim in outputs that parity harnesses compare.

DEC-21 PENTAGON PATCH. The FTExamples Kallen-prescription patch from
  pentagon_triage.md is DEFERRED to M5 (R8 freeze: no old-core/config
  changes until the DiffExp2 ladder reaches pentagon; the patch text lives
  in the triage doc).

DEC-22 PREFACTOR POWERS. v1 contract narrowing per the minimalism
  amendment: IntegrateOverLine prefactors are rational in (x) times
  integer powers of eps only (eps-DEPENDENT closed-form prefactors like
  Beta-function pins are handled at the FT layer, which already carries
  EpsPrefactors as integer shifts). API test 17 replaced accordingly.

DEC-23 LOADSYSTEM INPUTS. Files only in v1 (no in-memory matrices); unit
  tests write temp files. The canonical dlog d_1.m format is declared
  legacy-parity-only (classic pentagon example), per the API open question.

DEC-24 INDEFINITE INTEGRALS. API exposes IntegrateOverLine only; the 6
  IndefiniteIntegral test sites are retargeted/retired at M6 (disposition
  table updated then).

DEC-25 IN-MEMORY ORACLE GAPS. The pins flagged NOT FOUND (banana z0=2/5
  refs, box L2 J1/J2 pointwise 20-digit values, etc.) are regenerated at
  their consuming milestone (M2/M3) by the orchestrator's kernel queue, not
  blocked on now.
