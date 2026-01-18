# DiffExp Function Assignments

This document maps every function and variable in DiffExp.m to its designated subpackage.

## Subpackage Overview

1. **Symbols.m** - Core variables used across all packages
2. **State.m** - Configuration, state variables, and related functions
3. **Utilities.m** - Small helper functions
4. **SeriesOps.m** - Series manipulation operations
5. **Integration.m** - DiffExpIntegrate and integration rules
6. **Pade.m** - Pade approximant functions
7. **Mobius.m** - Mobius transformation functions
8. **LineSegmentation.m** - Line segmentation and interval functions
9. **MatrixLoading.m** - Matrix loading and preparation
10. **Frobenius.m** - Frobenius solutions
11. **Wronskian.m** - Wronskian and F-matrix methods
12. **AnalyticContinuation.m** - Analytic continuation logic
13. **ErrorEstimates.m** - Error estimation functions
14. **Transport.m** - IntegrateSystem and TransportTo

---

## 1. Symbols.m

### Variables to define:
- `x` (line parameter, defaults to Global`x)
- `eps` / `\[Epsilon]` (dimensional regulator)
- `Logx` (symbol for Log[x])
- `\[Theta]p` (HeavisideTheta[x])
- `\[Theta]m` (HeavisideTheta[-x])

### Rule:
```mathematica
eps := \[Epsilon];
```

### Key implementation note:
The line parameter `x` is set to `FEC[LineParameter]` which defaults to `Global`x`.
This allows users to change it via UpdateConfiguration.

---

## 2. State.m

### Configuration Variables:
- `DefaultConfiguration` (Association)
- `DiffExpConfiguration` (main config storage)
- `FEC` (shorthand for DiffExpConfiguration)
- `DiffExpExtensions` (list)

### Configuration Value Accessors:
- `ChopPrecisionVal`
- `LinearSolveChopPrecisionVal`
- `CrosscheckChopPrecision`
- `ExternalScalesVal`
- `LineParameterVal`
- `MatrixDirectoryVal`
- `EpsilonOrderVal`
- `FEAccuracyGoal`
- `FEWorkingPrecision`
- `DeltaPrescriptionsVal`
- `UseMobiusVal`
- `RadiusOfConvergenceVal`
- `DivisionOrderVal`
- `ExpansionOrderVal`
- `MaxCouplingOrder`

### State Variables:
- `AnalyticContinuationFailed`
- `AnalyticContinuationReplacements`
- `AnalyticContinuationReplacementsAssociation`
- `BenchmarkData`
- `CurrentSingularityWasAddedFromSquareRoot`
- `CurrentSingularityHasIDeltaPrescription`
- `DEqnSquareRoots`
- `MultivaluedFail`
- `UserDeltaPrescriptions`
- `UsingClosedFormMatrix`
- `DEqnMatricesFactored`
- `DEqnMatricesFactoredClosedForm`
- `DEqnMatricesExpanded`
- `NumIntegrals`
- `IntegrationSequence`
- `ExpansionMatrices`
- `ExpansionMatricesCanonical1`
- `ExpansionMatricesClosedForm`
- `AlphabetLogs`
- `AlphabetLogRules`
- `AlphabetLogRulesFactored`
- `AlphabetLogRulesExpanded`
- `MatricesIrreducibleFactors`
- `BufferedData`
- `MyWronsk`
- `MyWronskDetInv`
- `CurrCrosscheckFlags`
- `CrosscheckFlags`
- `LogStream`

### Internal Constants:
- `ISeriesChangeCoefficient`
- `IMaxLogOrder`
- `IMaxLogOrderDefault`
- `ICheckMultivaluedChop`
- `ICrossCheckPrintResultOrder`
- `ICrossCheckVerifyResultOrder`
- `ISafetyDigits`
- `ISafetyExpansionSubtract`
- `IExpansionOrdersAveraging`
- `IExpansionOrderIncrease`
- `IExpansionOrderDecrease`
- `IExpansionOrderIncrease2`
- `IDigitsSurplusDecreaseExpansionOrder`
- `ICurrEvalErrorSeriesDecrease`
- `IDecreaseOrderByErrorPrecise`
- `IMinExpansionOrder`

### Functions:
- `CurrentConfiguration[]` (line 128)
- `LoadConfiguration[a__]` (line 491)
- `UpdateConfiguration[a__Rule]` (line 496)
- `UpdateConfiguration[l_List]` (line 497)
- `UpdateConfiguration[assoc_Association]` (line 498)
- `SquareRootPrescriptionsAdded[]` (line 177)

---

## 3. Utilities.m

### Attributes to set:
```mathematica
SetAttributes[{SApply,SN,SExpand,SN,SSN,SMultiply,SEval,SEval1,SEval2,DecreaseSeriesOrderBy,SeriesCoefficientMinus,SplitTimes,ApplyAnalyticContinuation},Listable];
```

### Functions:
- `PrintDebug[args__][lev_]` (line 212)
- `PrintInfo[args__][lev_]` (line 213)
- `PrintWarning[args__]` (line 214)
- `ReportError[mes__]` (line 215-216)
- `AllSameQ[l_,b_]` (line 220)
- `CA` (alias for ConstantArray, line 221)
- `GetCases[expr_,case_]` (line 222)
- `DependsQ[a_,b_]` (line 223)
- `ZeroQ[a_]` (line 224)
- `R` (alias for ReplaceAll, line 225)
- `FirstOrNull[l_]` (line 226)
- `FindPivots[Matrix_]` (line 227)
- `SplitTimes[Expr_]` (line 228)
- `SplitSum[Expr_]` (line 250)
- `PChop` (line 233)
- `LSPChop` (line 234)
- `CPChop` (line 235)
- `IsPoint[line_]` (line 236)
- `IsLine[line_]` (line 237)
- `IntervalOverlapQ[intv1_,intv2_]` (line 349)
- `IntervalIntersec[intv1_,intv2_]` (line 350)
- `IntervalContainsQ[intv_,point_]` (line 351)
- `ExactLineQ[line_Association]` (line 1041)
- `FactorOrTogether[line_Association]` (line 1046-1047)

---

## 4. SeriesOps.m

### Attributes:
```mathematica
SetAttributes[DiffExpSeries,Listable];
```

### Functions:
- `SApply[f_,0]` (line 238)
- `SApply[f_,a_SeriesData]` (line 239)
- `SApply[f_,a_]` (line 240)
- `SExpand[0]` (line 241)
- `SExpand[a_SeriesData]` (line 242)
- `SExpand[a_]` (line 243)
- `SN` (line 244)
- `SSN` (line 245)
- `SMultiply[a_,b_]` (line 246)
- `SeriesCoefficientMinus[a_SeriesData,k_]` (line 247)
- `SeriesCoefficientMinus[a_,k_]` (line 248)
- `ApplyAnalyticContinuation[s_SeriesData]` (line 249)
- `SafeReplaceSeries11[...]` (lines 251-256)
- `MaxLogxPower[ex_]` (line 229)
- `LogxCoeff[Ser_,Which_]` (line 230)
- `LogxCoeffNS[Ser_,Which_]` (line 231)
- `LogxCoeffList[Ser_]` (line 232)
- `MatrixMultiplySExpand[MatA_,MatB_]` (line 258)
- `MatrixPowerSExpand[a_,n_]` (line 296-297)
- `DiffExpSeries[Ser_,ord_]` (line 300)
- `DiffExpSeries[Ser_]` (line 317)
- `SeriesAlways[term_,{a_,b_,c_},ex_]` (line 301)
- `LeadingCoefficientSeries[Ser_,AddTo2_]` (line 319)
- `SeriesMinPower[Ser_]` (line 337)
- `SeriesMaxPower[Ser_]` (line 338)
- `DecreaseSeriesOrderBy[a_,k_]` (line 340)
- `SD[a_,b_]` (line 399) - Series derivative avoiding Log[x]
- `SD[a_,b__]` (line 405)
- `SD[a_]` (line 406)

---

## 5. Integration.m

### Attributes:
```mathematica
SetAttributes[DiffExpIntegrate1,Listable];
```

### Variables:
- `IntReps` (line 410)

### Functions:
- `DiffExpIntegrate[a__]` (line 411)
- `DiffExpIntegrate1[a_]` (line 420)
- `DiffExpIntegrate1[exp0_SeriesData,var_]` (line 421)
- `DiffExpIntegrate1[exp0_/;NumericQ[exp0],var_]` (line 428)
- `DiffExpIntegrate1[exp0_,var_]` (line 429)
- `UpdateIntReps[MaxOrd_]` (line 430)

---

## 6. Pade.m

### Functions:
- `GetPade[0]` (line 361)
- `GetPade[a_?NumericQ]` (line 362)
- `GetPade[a_SeriesData]` (line 363)
- `SEval1[a_SeriesData]` (line 386)
- `SEval1[0]` (line 387)
- `SEval1[a_/;NumericQ[a]]` (line 388)
- `SEval2[a_,at_]` (line 389)
- `SEval[a_SeriesData,at_]` (line 392)
- `SEval[0,at_]` (line 395)
- `SEval[a_/;NumericQ[a],at_]` (line 396)
- `ToPiecewise[SavedData2_,Pade,Ord_Integer]` (line 444)

---

## 7. Mobius.m

### Functions:
- `GetMobius[{-\[Infinity],zmid_,\[Infinity]}]` (line 354)
- `GetMobius[{zmin_,zmid_,\[Infinity]}]` (line 355)
- `GetMobius[{-\[Infinity],zmid_,zmax_}]` (line 356)
- `GetMobius[{zmin_,zmid_,zmax_}]` (line 357)
- `GetMobius[{zmin_,\[Infinity],zmax_}]` (line 358)
- `GetLineRescaled[line_Association,at_,{signsproj_,signsim_},nomobius_]` (line 2075)
- `GetMobiusCPL[...]` (lines 2092-2094)
- `GetMobiusCPR[...]` (lines 2095-2097)
- `GetCPLRep[MyEq_]` (line 2099)
- `GetCPL[...]` (lines 2103-2117)
- `GetCPR[...]` (lines 2118-2131)
- `FindNextCenterPointL[xbc_,singsproj_]` (line 2134)
- `FindNextCenterPointR[xbc_,singsproj_]` (line 2139)

---

## 8. LineSegmentation.m

### Functions:
- `RelateLines[a2_Association,b2_Association,noerror_]` (line 1240)
- `RelateLinesPoint[a_Association,b_Association,pointb_]` (line 1277)
- `FindMatrixSingularities[line_,getcomplex_,{fixat_,to_}]` (line 1280)
- `PrintMobiusNormalized[a_]` (line 1315)
- `GetLargestTerm[line_]` (line 1323)
- `GetMatricesPrecisionDistance[line_Association]` (line 2146)
- `CheckBoundaryConditionsAndReparametrize[bcs3_,line_Association]` (line 2179)
- `GetMatchingPoint[line_Association,bcsline_]` (line 2267)

---

## 9. MatrixLoading.m

### Functions:
- `LoadMatrices[Folder_]` (line 594)
- `PrepareMatrices[line_Association]` (line 990)
- `PrepareMatricesFrom1[lineorig_Association,linenew_Association]` (line 998)
- `PrepareMatricesFrom[lineorig_Association,linenew_Association]` (line 1035)
- `PrepareMatricesFactored[line_Association]` (line 1049)
- `PrepareMatricesExpanded[line_Association]` (line 1075)
- `ClearMatrices[line_]` (line 973)
- `ClearMatrices[]` (line 981)
- `NullSpaceTryAgainOnFail[ex_,r___]` (line 1090)
- `CombineDifferentialEquationsHomogeneous[Amat_,topind_]` (line 1116)
- `InitializeIntegrationSequence[line_]` (line 851)

---

## 10. Frobenius.m

### Functions:
- `Frobenius1[DEqn_]` (line 1143)
- `FrobeniusSolutions[DEqn_]` (line 1195)

---

## 11. Wronskian.m

### Functions:
- `MatrixLogxInverse[Mat_]` (line 2304)

Note: Wronskian computation is embedded in IntegrateSystem. The FMat/FMatInv logic is complex and tightly coupled.

---

## 12. AnalyticContinuation.m

### Attributes:
```mathematica
SetAttributes[Project\[Theta]s,Listable];
```

### Functions:
- `PrepareAnalyticContinuation[Line_]` (line 786)
- `Project\[Theta]s[Expr_,f_]` (line 2291)

---

## 13. ErrorEstimates.m

### Functions (embedded in TransportTo but can be extracted):
- `UpdateMatrixExpansionError[]` (line 1640) - defined inside TransportTo
- `ComputeErrorsPerIndeterminate[aaa_,bbb_,ExpansionsIndeterminates_]` (line 1828)
- `PrintError[]` (line 1846) - defined inside TransportTo

---

## 14. Transport.m

### Functions:
- `PrepareBoundaryConditions[bcs_List,line2_Association|line2_List]` (line 871)
- `TransportTo[bcs2_List,line2_,to2_,SaveExpansions_,SampleAtList_List]` (line 1327)
- `IntegrateSystem[bcs2_,line2_,opts2_]` (line 2315)
- `EvaluateCurrPoint[]` (line 1651) - defined inside TransportTo
- `GiveMultivaluedError[]` (line 3129) - defined inside IntegrateSystem
- `TurnOffPade[]` (line 2341) - defined inside IntegrateSystem

---

## Public API (exported from main package)

These symbols should be exported:
- `\[Epsilon]`
- `eps`
- `Logx`
- `\[Theta]p`
- `\[Theta]m`
- `CurrentConfiguration`
- `IntegrateSystem`
- `LoadConfiguration`
- `UpdateConfiguration`
- `PrepareBoundaryConditions`
- `ToPiecewise`
- `TransportTo`
- Configuration options: `ChopPrecision`, `DeltaPrescriptions`, `DivisionOrder`, `EpsilonOrder`, `ExpansionOrder`, `IntegrationStrategy`, `LineParameter`, `LogFile`, `MatrixDirectory`, `RadiusOfConvergence`, `SegmentationStrategy`, `UseMobius`, `UsePade`, `Verbosity`, `WorkingPrecision`, `Variables`, `AccuracyGoal`

---

## Implementation Order

1. Symbols.m - Must be first, imported by all
2. State.m - Core configuration, needs Symbols
3. Utilities.m - Helper functions, needs State for PChop etc.
4. SeriesOps.m - Series operations, needs Utilities
5. Integration.m - Integration rules, needs SeriesOps
6. Pade.m - Pade functions, needs SeriesOps and State
7. Mobius.m - Mobius transforms, needs State
8. LineSegmentation.m - Line handling, needs Mobius
9. Frobenius.m - Frobenius solutions, needs SeriesOps
10. Wronskian.m - Wronskian inverse, needs SeriesOps
11. AnalyticContinuation.m - Analytic continuation, needs State
12. MatrixLoading.m - Matrix handling, needs most above
13. ErrorEstimates.m - Error functions, needs SeriesOps
14. Transport.m - Main transport functions, needs all above
