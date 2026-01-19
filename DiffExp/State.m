(* ::Package:: *)

(* DiffExp State Subpackage *)
(* This package handles configuration and state management *)

BeginPackage["DiffExp`State`", {"DiffExp`Symbols`"}];

(* Configuration options - exported *)
(* Note: AccuracyGoal, Variables, and WorkingPrecision are System symbols
   and should not have usage declarations here to avoid shadowing.
   They are used as configuration keys via their System` context. *)
ChopPrecision::usage = "Indicates the number off zeros after the decimal point after which terms should be set to 0 in intermediate computations.";
DeltaPrescriptions::usage = "A list of polynomials in the kinematic invariants and masses, each of which should contain an explicit factor \[PlusMinus]I\[Delta]. The zeros of the polynomials should describe singularities such as physical threshold singularities, or branch points of square roots.";
DivisionOrder::usage = "This option determines the inverse distance to the nearest singularity at which the line segments are evaluated, when the predivision strategy is used.";
EpsilonOrder::usage = "An integer specifying the highest order in the dimensional regulator \[Epsilon] in which the integrals should be computed.";
ExpansionOrder::usage = "Specifies the maximum power of the line parameter that should be kept in intermediate series expansions.";
IntegrationStrategy::usage = "Determines how the differential equations are solved. The possible values are \"Default\", and \"VOP\".";
LineParameter::usage = "The line parameter used for parsing lines to DiffExp.";
LogFile::usage = "Location of a log file on which to write all output of the current session.";
MatrixDirectory::usage = "The location of a directory on the file system which contains the partial derivative matrices.";
RadiusOfConvergence::usage = "This option has the effect of rescaling the line parameter of each line segment.";
RationalizationTolerance::usage = "Tolerance for rationalizing exponents (a, b) in singularity decomposition. Default 10^-10.";
SegmentationStrategy::usage = "This option determines which segmentation strategy is used.";
UseMobius::usage = "This option determines whether the line segments are obtained by linear transformations or by Mobius transformations.";
UsePade::usage = "Determines whether Pade approximants are used while transporting boundary conditions.";
Verbosity::usage = "Determines the level of printed output.";

(* State accessors - exported for internal use *)
FEC::usage = "Shorthand for DiffExpConfiguration.";
DiffExpConfiguration::usage = "Main configuration storage.";
DefaultConfiguration::usage = "Default configuration values.";

(* Value accessors *)
ChopPrecisionVal::usage = "Current ChopPrecision value.";
LinearSolveChopPrecisionVal::usage = "Current LinearSolveChopPrecision value.";
CrosscheckChopPrecision::usage = "Crosscheck chop precision.";
ExternalScalesVal::usage = "Current external scales.";
LineParameterVal::usage = "Current line parameter value.";
MatrixDirectoryVal::usage = "Current matrix directory.";
EpsilonOrderVal::usage = "Current epsilon order value.";
FEAccuracyGoal::usage = "Current accuracy goal.";
FEWorkingPrecision::usage = "Current working precision.";
DeltaPrescriptionsVal::usage = "Current delta prescriptions.";
UseMobiusVal::usage = "Current UseMobius value.";
RadiusOfConvergenceVal::usage = "Current radius of convergence.";
DivisionOrderVal::usage = "Current division order.";
ExpansionOrderVal::usage = "Current expansion order value.";
MaxCouplingOrder::usage = "Maximum coupling order.";

(* State variables - exported for internal use *)
AnalyticContinuationFailed::usage = "Flag indicating if analytic continuation failed.";
AnalyticContinuationReplacements::usage = "Replacement rules for analytic continuation.";
AnalyticContinuationReplacementsAssociation::usage = "Association of analytic continuation replacements.";
BenchmarkData::usage = "Benchmark timing data.";
CurrentSingularityWasAddedFromSquareRoot::usage = "Flag for square root singularity.";
CurrentSingularityHasIDeltaPrescription::usage = "Flag for iDelta prescription.";
DEqnSquareRoots::usage = "Square roots from differential equations.";
MultivaluedFail::usage = "Flag for multivalued function failure.";
UserDeltaPrescriptions::usage = "User-specified delta prescriptions.";
UsingClosedFormMatrix::usage = "Flag for closed form matrix.";
DEqnMatricesFactored::usage = "Factored differential equation matrices.";
DEqnMatricesFactoredClosedForm::usage = "Closed form factored matrices.";
DEqnMatricesExpanded::usage = "Expanded differential equation matrices.";
NumIntegrals::usage = "Number of integrals.";
IntegrationSequence::usage = "Order of integration.";
ExpansionMatrices::usage = "Expansion matrices.";
ExpansionMatricesCanonical1::usage = "Canonical expansion matrices.";
ExpansionMatricesClosedForm::usage = "Closed form expansion matrices.";
AlphabetLogs::usage = "Alphabet logarithms.";
AlphabetLogRules::usage = "Alphabet log rules.";
AlphabetLogRulesFactored::usage = "Factored alphabet log rules.";
AlphabetLogRulesExpanded::usage = "Expanded alphabet log rules.";
MatricesIrreducibleFactors::usage = "Irreducible factors from matrices.";
BufferedData::usage = "Buffered computation data.";
MyWronsk::usage = "Wronskian storage.";
MyWronskDetInv::usage = "Inverse Wronskian determinant storage.";
CurrCrosscheckFlags::usage = "Current crosscheck flags.";
CrosscheckFlags::usage = "Crosscheck flag definitions.";
DiffExpExtensions::usage = "DiffExp extensions.";
LogStream::usage = "Log stream for output.";
LastErrorContext::usage = "Stores diagnostic context from the last error for debugging. Access with DiffExp`State`LastErrorContext.";

(* Internal constants *)
ISeriesChangeCoefficient::usage = "Internal constant.";
IMaxLogOrder::usage = "Maximum log order.";
IMaxLogOrderDefault::usage = "Default maximum log order.";
ICheckMultivaluedChop::usage = "Multivalued check precision.";
ICrossCheckPrintResultOrder::usage = "Crosscheck print order.";
ICrossCheckVerifyResultOrder::usage = "Crosscheck verify order.";
ISafetyDigits::usage = "Safety digits.";
ISafetyExpansionSubtract::usage = "Safety expansion subtract.";
IExpansionOrdersAveraging::usage = "Expansion orders averaging.";
IExpansionOrderIncrease::usage = "Expansion order increase.";
IExpansionOrderDecrease::usage = "Expansion order decrease.";
IExpansionOrderIncrease2::usage = "Expansion order increase 2.";
IDigitsSurplusDecreaseExpansionOrder::usage = "Digits surplus decrease.";
ICurrEvalErrorSeriesDecrease::usage = "Current eval error series decrease.";
IDecreaseOrderByErrorPrecise::usage = "Decrease order by error precise.";
IMinExpansionOrder::usage = "Minimum expansion order.";

(* Helper function *)
SquareRootPrescriptionsAdded::usage = "Returns prescriptions added from square roots.";

Begin["`Private`"];

(* Default configuration *)
(* Use System` symbols explicitly for AccuracyGoal, Variables, WorkingPrecision to avoid shadowing *)
DefaultConfiguration = {
  System`AccuracyGoal -> "?",
  "AccuracyGoalValidate" -> "Before",
  ChopPrecision -> 250,
  "CrosscheckLevel" -> 0,
  DeltaPrescriptions -> {},
  DivisionOrder -> 3,
  EpsilonOrder -> 4,
  "EstimateError" -> "Fast",
  ExpansionOrder -> 50,
  "IgnoreIndicialCheck" -> False,
  "InvWronskSolver" -> "Auto",
  "KeepMatrixExpansions" -> False,
  "LinearSolveChopPrecision" -> 250,
  LineParameter -> Global`x,
  MatrixDirectory -> "",
  RadiusOfConvergence -> 1,
  RationalizationTolerance -> 10^-10,  (* Tolerance for rationalizing exponents in singularity decomposition *)
  "SaveExpansionsCompress" -> False,
  "HomogeneousSolve" -> "Expand",
  "Parallel" -> False,
  SegmentationStrategy -> "Predivision",
  IntegrationStrategy -> "Default",
  UseMobius -> False,
  UsePade -> False,
  System`Variables -> {},
  Verbosity -> 1,
  "VerbosityDebug" -> 0,
  System`WorkingPrecision -> 500
} // Association;

DiffExpExtensions = {};

(* Initialize configuration if not already set *)
If[!ValueQ[DiffExpConfiguration],
  DiffExpConfiguration = DefaultConfiguration;
];
FEC := DiffExpConfiguration;

(* Value accessors using delayed evaluation *)
ChopPrecisionVal := FEC[ChopPrecision];
LinearSolveChopPrecisionVal := FEC["LinearSolveChopPrecision"];
CrosscheckChopPrecision := 30;
ExternalScalesVal := FEC[System`Variables];
LineParameterVal := FEC[LineParameter];

(* Set the internal x to the line parameter *)
DiffExp`Symbols`x = FEC[LineParameter];

MatrixDirectoryVal := FEC[MatrixDirectory];
EpsilonOrderVal := FEC[EpsilonOrder];
FEAccuracyGoal := FEC[System`AccuracyGoal];
FEWorkingPrecision := FEC[System`WorkingPrecision];
DeltaPrescriptionsVal := FEC[DeltaPrescriptions];
UseMobiusVal := FEC[UseMobius];
RadiusOfConvergenceVal := FEC[RadiusOfConvergence];
DivisionOrderVal := FEC[DivisionOrder];

(* Initialize mutable state variables *)
If[!ValueQ[ExpansionOrderVal], ExpansionOrderVal = 50;];
If[!ValueQ[MaxCouplingOrder], MaxCouplingOrder = 1;];

If[!ValueQ[DEqnMatricesFactored],
  DEqnMatricesFactored = Association[{}];
];
If[!ValueQ[DEqnMatricesFactoredClosedForm],
  DEqnMatricesFactoredClosedForm = Association[{}];
];
If[!ValueQ[DEqnMatricesExpanded],
  DEqnMatricesExpanded = Association[{}];
];
If[!ValueQ[NumIntegrals],
  NumIntegrals = 0;
];

(* State variables initialization *)
AnalyticContinuationFailed = False;
AnalyticContinuationReplacements = {};
AnalyticContinuationReplacementsAssociation = Association[{}];
BenchmarkData = Association[];
CurrentSingularityWasAddedFromSquareRoot = False;
CurrentSingularityHasIDeltaPrescription = False;
DEqnSquareRoots = {};
MultivaluedFail = False;
UserDeltaPrescriptions = {};
UsingClosedFormMatrix = False;
LastErrorContext = {};

(* Crosscheck flags *)
CrosscheckFlags = {
  "FrobeniusSolutions" -> 1,
  "MatrixDelta" -> 1,
  "Wronskians" -> 1,
  "WronskInv" -> 0,
  "PeriodMatrix" -> 1,
  "GeneralSolutionMatrix" -> 2,
  "GeneralSolution" -> 1,
  "VariationOfParameters" -> 1,
  "SingularityCheck" -> 0
};
CurrCrosscheckFlags = {};

(* Internal constants *)
ISeriesChangeCoefficient = 2;
IMaxLogOrder = IMaxLogOrderDefault = 1;
ICheckMultivaluedChop = 5;
ICrossCheckPrintResultOrder = 5;
ICrossCheckVerifyResultOrder = 5;
ISafetyDigits = 2;
ISafetyExpansionSubtract = 5;
IExpansionOrdersAveraging = 3;
IExpansionOrderIncrease = 10;
IExpansionOrderDecrease = 10;
IExpansionOrderIncrease2 = 25;
IDigitsSurplusDecreaseExpansionOrder = 3;
ICurrEvalErrorSeriesDecrease := Ceiling[0.7(MaxCouplingOrder)] + 2;
IDecreaseOrderByErrorPrecise = MaxCouplingOrder;
IMinExpansionOrder = 10;

(* Helper function for square root prescriptions *)
SquareRootPrescriptionsAdded[] := {#, 1} & /@ Complement[
  Expand[DEqnSquareRoots],
  Flatten@Expand@{UserDeltaPrescriptions[[All, 1]], -UserDeltaPrescriptions[[All, 1]]}
];

End[];

EndPackage[];
