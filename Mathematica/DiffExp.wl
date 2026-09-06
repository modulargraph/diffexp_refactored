(* DiffExp 2.1: input preparation and a process wrapper for native C++ transport. *)
BeginPackage["DiffExp`"];
$DiffExpExecutable::usage = "$DiffExpExecutable is Automatic or the path to the diffexp executable.";
DiffExpRun::usage = "DiffExpRun[{arguments}, input] runs diffexp and returns its exit code, stdout and stderr. input defaults to an empty string. A failed process returns Failure.";
DiffExpSeries::usage = "DiffExpSeries[request] sends an Association describing a regular Taylor-series problem to the C++ backend and returns an Association of exact coefficient strings.";
DiffExpFeynmanTrick::usage = "DiffExpFeynmanTrick[configurationOrFile, {arguments}] runs a generic family configuration through native Feynman-trick recursion and returns JSON data. DiffExpFamilyTemplate[name] supplies an editable starting configuration.";
DiffExpBackendInfo::usage = "DiffExpBackendInfo[] returns information about the C++ backend.";

CurrentConfiguration::usage = "CurrentConfiguration[] returns the current transport configuration.";
LoadConfiguration::usage = "LoadConfiguration[rules] resets and loads differential-equation configuration.";
UpdateConfiguration::usage = "UpdateConfiguration[rules] updates the current transport configuration.";
PrepareBoundaryConditions::usage = "PrepareBoundaryConditions[values, pointOrLine] prepares epsilon coefficients or asymptotic boundary input.";
TransportTo::usage = "TransportTo[boundary, pointOrLine, to:1, save:False] runs native transport and returns {point, values, errors}. Saved output includes replay metadata.";
IntegrateSystem::usage = "IntegrateSystem[boundary, line] returns numerical functions backed by native transport. Symbolic free-constant solutions are not supported.";
ToPiecewise::usage = "ToPiecewise[saved, pade:False] returns numerical functions. Evaluation replays and caches native transport at the requested physical point; pade is accepted for call compatibility and does not change the native method.";
DiffExpLastTimings::usage = "DiffExpLastTimings[] returns the most recent wrapper and native timings.";
DiffExpFamilyTemplate::usage = "DiffExpFamilyTemplate[name] returns a complete editable family configuration for DiffExpFeynmanTrick.";
BasisPrefactors::usage = "BasisPrefactors is a list of algebraic normalization factors, one per integral. Native transport converts their continued square-root sheets to principal endpoint values.";
ChopPrecision::usage = "DiffExp configuration or boundary-input symbol.";
DeltaPrescriptions::usage = "DiffExp configuration or boundary-input symbol.";
DivisionOrder::usage = "DiffExp configuration or boundary-input symbol.";
EpsilonOrder::usage = "DiffExp configuration or boundary-input symbol.";
ExpansionOrder::usage = "DiffExp configuration or boundary-input symbol.";
LineParameter::usage = "DiffExp configuration or boundary-input symbol.";
LogFile::usage = "DiffExp configuration or boundary-input symbol.";
MatrixDirectory::usage = "DiffExp configuration or boundary-input symbol.";
Recurrence::usage = "Recurrence selects native transport: auto (default), taylor, series, or spectral. Spectral requires positive AccuracyGoal and an ordinary epsilon-linear system.";
RadiusOfConvergence::usage = "DiffExp configuration or boundary-input symbol.";
SegmentationStrategy::usage = "DiffExp configuration or boundary-input symbol.";
IntegrationStrategy::usage = "DiffExp configuration or boundary-input symbol.";
UseMobius::usage = "DiffExp configuration or boundary-input symbol.";
UsePade::usage = "DiffExp configuration or boundary-input symbol.";
Verbosity::usage = "DiffExp configuration or boundary-input symbol.";
\[Epsilon]::usage = "DiffExp configuration or boundary-input symbol.";
eps::usage = "DiffExp configuration or boundary-input symbol.";
Logx::usage = "DiffExp configuration or boundary-input symbol.";
Begin["`Private`"];
$wrapperDirectory = DirectoryName[$InputFileName];
If[!ValueQ[$DiffExpExecutable], $DiffExpExecutable = Automatic];
Options[DiffExpRun] = {"Executable" -> Automatic};
Options[DiffExpSeries] = Options[DiffExpRun];
Options[DiffExpFeynmanTrick] = Options[DiffExpRun];
Options[DiffExpBackendInfo] = Options[DiffExpRun];

executable[option_] := Module[{selected = option, environment, candidates},
  If[selected === Automatic, selected = $DiffExpExecutable];
  If[selected =!= Automatic, Return[selected]];
  environment = Environment["DIFFEXP_EXECUTABLE"];
  If[StringQ[environment] && environment =!= "", Return[environment]];
  candidates = {
    FileNameJoin[{$wrapperDirectory, "..", "build", "diffexp"}],
    FileNameJoin[{$wrapperDirectory, "..", "..", "..", "bin", "diffexp"}]
  };
  SelectFirst[candidates, FileExistsQ, "diffexp"]
];

DiffExpRun[arguments : {___String}, input_String : "", OptionsPattern[]] :=
 Module[{program = executable[OptionValue["Executable"]], result},
  If[!StringQ[program] || program === "",
    Return[Failure["Executable", <|"MessageTemplate" -> "Set a nonempty executable path."|>]]];
  (* RunProcess requires normalized executable paths on some platforms. *)
  If[FileExistsQ[program], program = ExpandFileName[program]];
  (* Omit stdin for empty input: a fast process may exit before a pipe write. *)
  result = Quiet[If[input === "",
    RunProcess[Prepend[arguments, program], All],
    RunProcess[Prepend[arguments, program], All, input]]];
  If[!AssociationQ[result], Return[Failure["ProcessLaunch", <|
    "MessageTemplate" -> "Could not launch the DiffExp executable.", "Executable" -> program|>]]];
  If[Lookup[result, "ExitCode", -1] =!= 0,
    Return[Failure["ProcessExit", Join[<|"MessageTemplate" -> "DiffExp reported an error."|>, result]]]];
  result
];

jsonRun[arguments_, input_, program_] := Module[{process, data, elapsed},
  {elapsed,process} = AbsoluteTiming[DiffExpRun[arguments, input, "Executable" -> program]];
  $lastTimings=<|"process_seconds"->elapsed|>;
  If[FailureQ[process], Return[process]];
  data = Quiet[Check[ImportString[process["StandardOutput"], "RawJSON"], $Failed]];
  If[!AssociationQ[data], Return[Failure["JSONResponse", Join[<|
    "MessageTemplate" -> "DiffExp did not return a JSON object."|>, process]]]];
  $lastTimings=Join[$lastTimings,Lookup[data,"timings",<||>]];
  data
];

DiffExpSeries[request_Association, OptionsPattern[]] := Module[{input},
  input = Quiet[Check[ExportString[request, "RawJSON", "Compact" -> True], $Failed]];
  If[!StringQ[input], Return[Failure["JSONRequest", <|
    "MessageTemplate" -> "The request must contain JSON-compatible values."|>]]];
  jsonRun[{"series", "-"}, input, OptionValue["Executable"]]
];
DiffExpFeynmanTrick[family_String, arguments : {___String} : {}, OptionsPattern[]] :=
  jsonRun[Join[{"ft", family, "--json"}, arguments], "", OptionValue["Executable"]];
DiffExpBackendInfo[OptionsPattern[]] :=
  jsonRun[{"backend-info"}, "", OptionValue["Executable"]];
(* Compatibility input preparation and native process dispatch. *)
eps := \[Epsilon];
$defaultConfiguration = <|AccuracyGoal -> "?", BasisPrefactors -> {}, ChopPrecision -> 250,
 DeltaPrescriptions -> {}, DivisionOrder -> 3, EpsilonOrder -> 4, ExpansionOrder -> 50,
 LineParameter -> Global`x, MatrixDirectory -> "", RadiusOfConvergence -> 1, Recurrence -> "auto",
 SegmentationStrategy -> "Predivision", IntegrationStrategy -> "Default", UseMobius -> False,
 UsePade -> False, Variables -> {}, Verbosity -> 1, WorkingPrecision -> 500|>;
If[!AssociationQ[$configuration], $configuration = $defaultConfiguration];
$matrixEntries = {}; $matrixDimension = 0; $lastTimings = <||>;
compatFailure[tag_, message_, data_:<||>] := Throw[Failure[tag, Join[<|"MessageTemplate" -> message|>, data]], "DiffExpCompatibility"];
configValue[key_] := Lookup[$configuration, key, Lookup[$defaultConfiguration, key, Missing["NotConfigured"]]];
CurrentConfiguration[] := $configuration;
DiffExpLastTimings[] := $lastTimings;
optionKeys = Association[Map[(SymbolName[#] -> #)&, Keys[$defaultConfiguration]]];
normalizeConfiguration[a_Association] := Association[KeyValueMap[(Lookup[optionKeys, #1, #1] -> #2)&, a]];
LoadConfiguration[args___] := ($configuration = $defaultConfiguration; $matrixEntries = {}; $matrixDimension = 0; UpdateConfiguration[args]);
UpdateConfiguration[rules__Rule] := UpdateConfiguration[{rules}];
UpdateConfiguration[] := CurrentConfiguration[];
UpdateConfiguration[rules_List] := UpdateConfiguration[Association[rules]];
UpdateConfiguration[configuration_Association] := Module[{previous=$configuration,previousEntries=$matrixEntries,previousDimension=$matrixDimension,answer},
 answer=Catch[Module[{updates = normalizeConfiguration[configuration], result},
 $configuration = Join[$configuration, updates];
 If[!IntegerQ[configValue[EpsilonOrder]] || configValue[EpsilonOrder]<0 || !IntegerQ[configValue[ExpansionOrder]] || configValue[ExpansionOrder]<2,
  compatFailure["Configuration", "EpsilonOrder must be nonnegative and ExpansionOrder must be at least two."]];
 If[!NumericQ[configValue[WorkingPrecision]] || configValue[WorkingPrecision]<16,
  compatFailure["Configuration", "WorkingPrecision must be at least 16 decimal digits."]];
 If[!MemberQ[{"auto","taylor","series","spectral"},configValue[Recurrence]], compatFailure["Configuration", "Recurrence must be auto, taylor, series, or spectral."]];
 If[MemberQ[configValue[Variables], configValue[LineParameter]], compatFailure["Configuration", "LineParameter cannot be a matrix variable."]];
 If[KeyExistsQ[updates, MatrixDirectory] || (KeyExistsQ[updates, EpsilonOrder] && configValue[MatrixDirectory] =!= ""), loadMatrices[]];
 CurrentConfiguration[]], "DiffExpCompatibility"];
 If[FailureQ[answer],$configuration=previous;$matrixEntries=previousEntries;$matrixDimension=previousDimension];answer];

inputString[value_] := ToString[value, InputForm];
numericString[value_] := Module[{number = N[value, configValue[WorkingPrecision]], real, imaginary},
 If[!NumberQ[number], compatFailure["BoundaryValue", "Boundary coefficient is not numerical.", <|"Value" -> inputString[value]|>]];
 (* Preserve input precision annotations; rationalizing rounded boundary data
    would falsely present those decimals as exact values to the native solver. *)
 If[MatchQ[value,_Integer|_Rational], inputString[value], inputString[number]]
];
lineAssociation[line_Association] := line;
lineAssociation[line_List] := Association[line];
lineAssociation[_] := compatFailure["Line", "A point or line must be an Association or a list of rules."];
(* Plot and numerical endpoints supply approximate coordinates. Freeze their
   nominal values as rationals, bounded by the configured precision, before
   exact native pullback; boundary uncertainties are handled separately. *)
exactPathInput[value_] := value /. q_Real :> Rationalize[N[q,configValue[WorkingPrecision]],0];
wirePaths[line_Association] := Association[KeyValueMap[(SymbolName[#1] -> inputString[exactPathInput[#2] /. configValue[LineParameter] -> Global`x])&, line]];

loadMatrices[] := Module[{folder = configValue[MatrixDirectory], files, entries = {}, dimension = 0, file, name, parts, variable, order, matrix, logs, weight, remainder, encountered = {}, vars = configValue[Variables]},
 If[!StringQ[folder] || !DirectoryQ[folder], compatFailure["MatrixDirectory", "MatrixDirectory is not an existing directory."]];
 files = Sort[FileNames["d*_*.m", folder]];
 If[files === {}, compatFailure["Matrices", "No dVARIABLE_ORDER.m or d_1.m matrices were found."]];
 Do[
  name = FileNameTake[file]; parts = StringCases[name, RegularExpression["^d(.*)_([0-9]+)\\.m$"] -> {"$1", "$2"}];
  If[parts === {}, Continue[]]; {variable, order} = First[parts]; order = FromDigits[order];
  matrix = Quiet[Check[Block[{$Context = "Global`", $ContextPath = {"DiffExp`", "System`", "Global`"}}, Get[file]], $Failed]];
  If[Head[matrix] === SparseArray, matrix = Normal[matrix]];
  If[!MatrixQ[matrix] || Length[matrix]===0 || Dimensions[matrix][[1]] =!= Dimensions[matrix][[2]],
   compatFailure["MatrixShape", "A matrix file must contain a nonempty square matrix.", <|"File" -> file|>]];
  encountered=Join[encountered,Cases[matrix,s_Symbol /; Context[s] === "Global`" && s =!= configValue[LineParameter],Infinity]];
  If[variable =!= "",AppendTo[encountered,Symbol["Global`"<>variable]]];
  If[dimension===0, dimension=Length[matrix]];
  If[Length[matrix] =!= dimension, compatFailure["MatrixShape", "Matrix dimensions disagree."]];
  If[variable === "" && order =!= 1, compatFailure["CanonicalMatrix", "Canonical logarithmic matrices must be named d_1.m."]];
  If[order > configValue[EpsilonOrder],Continue[]];
  Do[If[matrix[[row,column]] === 0, Continue[]];
   If[variable === "",
    logs=DeleteDuplicates[Cases[matrix[[row,column]], _Log, {0,Infinity}]]; remainder=Expand[matrix[[row,column]]];
    Do[weight=Coefficient[remainder, log];
     If[!MatchQ[weight, _Integer|_Rational], compatFailure["CanonicalMatrix", "Canonical log weights must be exact rational constants."]];
     If[weight =!= 0, AppendTo[entries,<|"row"->row-1,"column"->column-1,"epsilon"->1,"variable"->"dlog","expression"->inputString[log[[1]]],"coefficient"->inputString[weight]|>]];
     remainder=Expand[remainder-weight log], {log,logs}];
    If[remainder =!= 0, compatFailure["CanonicalMatrix", "Canonical entries must be rational linear combinations of logarithms."]],
    AppendTo[entries,<|"row"->row-1,"column"->column-1,"epsilon"->order,"variable"->variable,"expression"->inputString[matrix[[row,column]]]|>]
   ],{row,dimension},{column,dimension}],{file,files}];
 If[dimension===0, compatFailure["Matrices", "No matrix was available within the requested epsilon window."]];
 If[vars === {},
  vars=DeleteDuplicates[Flatten[Table[
    If[entry["variable"] === "dlog", {}, {Symbol["Global`"<>entry["variable"]]}],{entry,entries}]]];
  (* Include symbols occurring in canonical letters and ordinary matrix coefficients. *)
  vars=DeleteDuplicates[Join[vars, encountered]];
  $configuration[Variables]=vars];
 $matrixDimension=dimension; $matrixEntries=entries;
];

(* Parse only a decimal grammar; native responses are never evaluated as Wolfram code. *)
decimalValue[s_String] := Module[{parts, mantissa, exponent=0, sign=1, digits, decimalPlaces, precision},
 If[!StringMatchQ[s, RegularExpression["[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?"]],
  compatFailure["NumericResponse", "The backend returned a non-decimal numeric field."]];
 parts=StringSplit[s,RegularExpression["[eE]"]]; mantissa=First[parts];
 If[Length[parts]===2, exponent=If[StringStartsQ[parts[[2]],"-"],-1,1] FromDigits[StringReplace[parts[[2]],StartOfString~~("+"|"-")->""]]];
 If[StringStartsQ[mantissa,"-"],sign=-1]; mantissa=StringReplace[mantissa,StartOfString~~("+"|"-")->""];
 decimalPlaces=If[StringContainsQ[mantissa,"."],StringLength[Last[StringSplit[mantissa,".",All]]],0];
 digits=StringDelete[mantissa,"."]; If[FromDigits[digits]===0,Return[0]]; precision=Max[20,StringLength[digits]];
 N[sign FromDigits[digits] 10^(exponent-decimalPlaces),precision]
];
decimalValue[_] := compatFailure["NumericResponse", "Numeric response fields must be decimal strings."];
wireNumber[a_Association] := decimalValue[a["real_midpoint"]]+I decimalValue[a["imaginary_midpoint"]];
wireRadius[a_Association] := decimalValue[a["radius"]];
wireNumber[_] := compatFailure["NumericResponse", "Expected a numeric response object."];
responseMatrices[data_Association] := Module[{values, errors, radii},
 values=Map[wireNumber,data["values"],{2}]; radii=Map[wireRadius,data["values"],{2}];
 errors=If[KeyExistsQ[data,"errors"], Map[(Abs[wireNumber[#]]+wireRadius[#])&,data["errors"],{2}],radii];
 {values,MapThread[Max,{errors,radii},2]}
];

boundaryCoefficients[bcs_List] := Module[{high=configValue[EpsilonOrder], regulator=\[Epsilon], result,series},
 result=Map[Function[row, Which[row === "?", ConstantArray["?",high+1],
  ListQ[row], If[Length[row]<high+1,compatFailure["BoundaryWindow","Too few epsilon coefficients in the boundary conditions."]];Take[row,high+1],
  True, series=Quiet[Check[Series[row,{regulator,0,high}],$Failed]];
   If[Head[series]===SeriesData && series[[4]]<0,compatFailure["BoundaryWindow","Boundary conditions must start at finite epsilon order."]];
   Table[Quiet[Check[SeriesCoefficient[series,{regulator,0,k}],$Failed]],{k,0,high}]]],bcs];
 If[!FreeQ[result,$Failed|_SeriesCoefficient],compatFailure["BoundaryExpansion","Could not expand the boundary in epsilon."]];result
];
PrepareBoundaryConditions[bcs_List,line_] := Catch[Module[{point=lineAssociation[line], coefficients, lp=configValue[LineParameter]},
 If[$matrixDimension>0 && Length[bcs] =!= $matrixDimension,compatFailure["BoundaryShape","The boundary has the wrong number of integrals."]];
 coefficients=boundaryCoefficients[bcs]; coefficients=coefficients /. Normal[point];
 If[FreeQ[Values[point],lp] && !FreeQ[coefficients,lp],compatFailure["BoundaryLine","Boundary terms depend on LineParameter but the boundary point does not."]];
 {point,coefficients}],"DiffExpCompatibility"];

asymptoticData[coefficients_,lp_] := Module[{constraints={},cutoffs={},expression,normal,cutoff,terms,power,degree,value,rowCutoffs,series},
 Do[rowCutoffs={};
  Do[expression=coefficients[[row,k+1]];If[expression==="?",Continue[]];
   If[Head[expression]===SeriesData,cutoff=expression[[5]]/expression[[6]];normal=Normal[expression],
    normal=Expand[expression];cutoff=Automatic];
   normal=Expand[normal /. Logx->Log[lp]];
   terms=If[Head[normal]===Plus,List@@normal,{normal}];
   Do[If[term===0,Continue[]];
    power=Exponent[term/.Log[lp]->Unique["log"],lp];
    If[!MatchQ[power,_Integer|_Rational],compatFailure["AsymptoticPower","The asymptotic boundary must have exact rational powers."]];
    degree=Exponent[term,Log[lp]]; value=Cancel[term/(lp^power Log[lp]^degree)];
    If[!FreeQ[value,lp],compatFailure["AsymptoticForm","Supply explicit powers/logarithms or a truncated asymptotic SeriesData boundary."]];
    AppendTo[constraints,<|"row"->row-1,"epsilon"->k,"power"->inputString[power],"log_degree"->degree,"value"->numericString[value]|>],{term,terms}];
   If[cutoff===Automatic,cutoff=If[normal===0,Infinity,1+Max[Exponent[#/.Log[lp]->Unique["log"],lp]&/@terms]]];AppendTo[rowCutoffs,cutoff],
   {k,0,configValue[EpsilonOrder]}];
  If[rowCutoffs=!={},AppendTo[cutoffs,<|"row"->row-1,"power"->inputString[Replace[Min[rowCutoffs],Infinity->1]]|>]],{row,Length[coefficients]}];
 <|"constraints"->constraints,"cutoffs"->cutoffs|>
];
baseRequest[paths_Association,save_] := Module[{request,prefactors=configValue[BasisPrefactors],aliases=<||>,renames,used=SymbolName/@Keys[paths],candidate,n=0,entries=$matrixEntries,mappedPaths},
 (* Kinematic x is legal when LineParameter is different. Only the wire
    representation must reserve x for the native path parameter. *)
 Do[If[name==="x" || name==="I" || name==="form" || StringMatchQ[name,RegularExpression["r[0-9]+"]],
  candidate="diffexpVariable"<>ToString[n++];While[MemberQ[used,candidate],candidate="diffexpVariable"<>ToString[n++]];
  AssociateTo[aliases,name->candidate];AppendTo[used,candidate]],{name,SymbolName/@Keys[paths]}];
 renames=KeyValueMap[(Symbol["Global`"<>#1]->Symbol["Global`"<>#2])&,aliases];
 mappedPaths=Association[KeyValueMap[(Symbol["Global`"<>Lookup[aliases,SymbolName[#1],SymbolName[#1]]]->#2)&,paths]];
 If[Length[renames]>0,
  entries=Map[Function[e,Join[e,<|"expression"->inputString[ToExpression[e["expression"],InputForm]/.renames],
   "variable"->Lookup[aliases,e["variable"],e["variable"]]|>]],entries];
  prefactors=prefactors/.renames];
 request=<|"schema"->"DiffExp.Transport/v1","dimension"->$matrixDimension,
 "epsilon_order"->configValue[EpsilonOrder],"taylor_order"->configValue[ExpansionOrder],
 "working_bits"->Ceiling[configValue[WorkingPrecision] Log[2,10]],
 "accuracy_goal"->If[NumberQ[configValue[AccuracyGoal]],configValue[AccuracyGoal],0],
 "recurrence"->configValue[Recurrence],"division_order"->configValue[DivisionOrder],"save_segments"->TrueQ[save],"paths"->wirePaths[mappedPaths],"entries"->entries|>;
 If[prefactors=!={},
  If[!ListQ[prefactors] || Length[prefactors]=!=$matrixDimension,compatFailure["BasisPrefactors","Supply one algebraic prefactor per integral."]];
  request["basis_prefactors"]=inputString/@prefactors];request];
nativeTransport[request_Association] := Module[{data,elapsed},
 $lastRequest=request;
 {elapsed,data}=AbsoluteTiming[jsonRun[{"transport","-"},ExportString[request,"RawJSON","Compact"->True],Automatic]];
 $lastTimings=Join[<|"process_seconds"->elapsed|>,If[AssociationQ[data],Lookup[data,"timings",<||>],<||>]];
 If[FailureQ[data],Throw[data,"DiffExpCompatibility"]];data
];

(* Coordinate deformation follows the supplied signed infinitesimal prescription. *)
deformedPaths[start_Association,finish_Association,base_:Automatic] := Module[{lp=configValue[LineParameter],paths,prescriptions=configValue[DeltaPrescriptions],polynomial,sign,derivative,candidates,crossings,atStart,atFinish,chosen,amplitude,delta=Global`\[Delta]},
 paths=If[base===Automatic,Association[KeyValueMap[(#1->(start[#1]+lp(#2-start[#1])))&,finish]],base];
 Do[candidates={};crossings={};
  Do[If[ListQ[prescription] && Length[prescription]===2,{polynomial,sign}=prescription,
    polynomial=prescription/.delta->0;sign=Coefficient[prescription,delta]/I];
   derivative=D[polynomial,variable];
   If[!FreeQ[derivative,Alternatives@@Keys[finish]] || !NumericQ[derivative] || derivative===0 || !NumericQ[sign],Continue[]];
   chosen=Sign[sign/derivative];If[!MemberQ[{-1,1},chosen],Continue[]];AppendTo[candidates,chosen];
   atStart=N[polynomial/.Normal[start]];atFinish=N[polynomial/.Normal[finish]];
   If[NumberQ[atStart] && NumberQ[atFinish] && TrueQ[Im[atStart]===0] && TrueQ[Im[atFinish]===0] && TrueQ[atStart atFinish<=0],AppendTo[crossings,chosen]],{prescription,prescriptions}];
  If[candidates=!={},chosen=First[If[crossings==={},candidates,crossings]];
   amplitude=Max[1,Abs[Rationalize[finish[variable]-start[variable],0]]]/4;
   paths[variable]=paths[variable]+I chosen amplitude lp(1-lp)],{variable,Keys[finish]}];paths
];

lineStartParameter[line_Association,start_Association,lp_] := Module[{dependent,candidates={},solutions,candidate,differences},
 dependent=Select[Keys[line],!FreeQ[line[#],lp]&];
 If[dependent==={},Return[0]];
 solutions=Quiet[Solve[line[First[dependent]]==start[First[dependent]],lp]];
 If[ListQ[solutions],candidates=lp/.solutions];
 candidate=SelectFirst[candidates,Function[q,
  differences=N[Values[line/.lp->q]-Lookup[start,Keys[line]],Min[40,configValue[WorkingPrecision]]];
  VectorQ[differences,NumericQ] && TrueQ[Max[Abs[differences]]<10^-25]],Missing["NotOnLine"]];
 If[MissingQ[candidate],compatFailure["BoundaryLine","The requested line does not pass through the numerical boundary point."]];candidate
];

TransportTo[bcs_List,line_,to_:1,save_:False] := Catch[Module[{began=AbsoluteTime[],boundary=bcs,target=Map[exactPathInput,lineAssociation[line]],toParameter=exactPathInput[to],lp=configValue[LineParameter],start,finish,coefficients,errors,paths,request,data,matrices,result,original=bcs,configuration=$configuration,explicitLine,parameter,from,basePaths,curve,nativeCalls={}},
 If[$matrixDimension===0,compatFailure["Matrices","LoadConfiguration must load differential-equation matrices first."]];
 If[Length[boundary]===2 && ListQ[boundary[[1]]] && AssociationQ[boundary[[1,1]]],boundary=boundary[[1]]];
 If[Length[boundary]<2 || !AssociationQ[boundary[[1]]],compatFailure["BoundaryShape","Use prepared boundary conditions or a previous transport result."]];
 start=Map[exactPathInput,boundary[[1]]];coefficients=boundary[[2]];
 If[Length[coefficients]=!=$matrixDimension || !And@@(ListQ[#]&&Length[#]>=configValue[EpsilonOrder]+1&/@coefficients),compatFailure["BoundaryShape","Boundary dimensions do not match the loaded system and epsilon window."]];
 coefficients=Take[#,configValue[EpsilonOrder]+1]&/@coefficients;
 If[!FreeQ[Values[start],lp],
  request=Join[baseRequest[start,False],<|"initial_only"->True,"asymptotic"->asymptoticData[coefficients,lp]|>];data=nativeTransport[request];AppendTo[nativeCalls,$lastTimings];
  parameter=decimalValue[data["parameter"]]; matrices=responseMatrices[data];
  start=Map[(#/.lp->Rationalize[parameter,0])&,start];coefficients=matrices[[1]];errors=matrices[[2]],
  errors=If[Length[boundary]>=3 && MatrixQ[boundary[[3]]],Take[#,configValue[EpsilonOrder]+1]&/@boundary[[3]],ConstantArray[0,Dimensions[coefficients]]]];
 explicitLine=!FreeQ[Values[target],lp];finish=If[explicitLine,Map[(#/.lp->toParameter)&,target],target];
 If[!And@@(KeyExistsQ[start,#]&/@Keys[finish]) || !And@@(KeyExistsQ[finish,#]&/@configValue[Variables]),compatFailure["LineVariables","The line must fix all matrix variables and include the boundary coordinates."]];
 If[explicitLine,
  from=lineStartParameter[target,start,lp];
  basePaths=Map[(#/.lp->(from+(toParameter-from)lp))&,target];curve=target,
  basePaths=Association[KeyValueMap[(#1->(start[#1]+lp(#2-start[#1])))&,finish]];curve=basePaths];
 paths=deformedPaths[start,finish,basePaths];
 request=Join[baseRequest[paths,False],<|"boundary"->Map[numericString,coefficients,{2}],"boundary_errors"->Map[numericString,errors,{2}]|>];
 data=nativeTransport[request];AppendTo[nativeCalls,$lastTimings];matrices=responseMatrices[data];result={finish,matrices[[1]],matrices[[2]]};
 $lastTimings=Join[$lastTimings,<|"wrapper_total_seconds"->AbsoluteTime[]-began,"native_total_seconds"->Total[Lookup[nativeCalls,"total_seconds",0]],"native_calls"->nativeCalls,"charts"->Lookup[data,"charts",Missing["NotAvailable"]]|>];
 If[TrueQ[configValue[Verbosity]>0],Print["DiffExp native transport: ",$lastTimings["wrapper_total_seconds"]," seconds."]];
 If[TrueQ[save],{result,{<|"schema"->"DiffExp.SavedTransport/v1","boundary"->original,"line"->curve,"to"->toParameter,"configuration"->configuration,"dimension"->$matrixDimension,"entries"->$matrixEntries,"segments"->Lookup[data,"segments",{}],"replay_at_physical_point"->True,"omitted_tails_certified"->False|>}},result]
],"DiffExpCompatibility"];

If[!AssociationQ[$savedMetadata],$savedMetadata=<||>];
If[!AssociationQ[$savedSamples],$savedSamples=<||>];
savedEvaluation[token_String,point_?NumericQ] := savedEvaluation[token,point] = Module[{metadata=$savedMetadata[token],samples,seed,result},
 Block[{$configuration=metadata["configuration"],$matrixEntries=metadata["entries"],$matrixDimension=metadata["dimension"]},
 $configuration[Verbosity]=0;seed=metadata["boundary"];
 samples=Lookup[$savedSamples,token,{}];
 If[TrueQ[Lookup[metadata,"linear_replay",False]] && samples=!={},
  seed=First[MinimalBy[samples,Abs[#[[1]]-point]&]][[2]]];
 result=TransportTo[seed,metadata["line"],point,False];
 If[FailureQ[result] && seed=!=metadata["boundary"],result=TransportTo[metadata["boundary"],metadata["line"],point,False]];
 If[!FailureQ[result],$savedSamples[token]=Append[samples,{point,result}]];result]];
savedComponent[evaluator_Symbol,row_,column_,point_?NumericQ] :=
 Module[{result=evaluator[point]},If[FailureQ[result],Indeterminate,result[[2,row,column]]]];
ToPiecewise[saved_List,pade_:False,order_:Null] := Catch[Module[{metadata,data=saved,token,endpoint=None,lp,line,dependent,samples={},boundary,from,evaluator},
 If[Length[data]===2 && ListQ[data[[1]]] && AssociationQ[data[[1,1]]],endpoint=data[[1]];data=data[[2]]];
 If[!MatchQ[data,{_Association}] || Lookup[data[[1]],"schema",""]=!="DiffExp.SavedTransport/v1",compatFailure["SavedSegments","Use TransportTo[..., True] to obtain saved transport data."]];metadata=First[data];
 lp=Lookup[metadata["configuration"],LineParameter,Global`x];line=metadata["line"];
 dependent=Select[Values[line],!FreeQ[#,lp]&];
 metadata["linear_replay"]=Length[dependent]===1 && PolynomialQ[First[dependent],lp] && Exponent[First[dependent],lp]===1;
 If[ListQ[endpoint],AppendTo[samples,{metadata["to"],endpoint}]];
 boundary=metadata["boundary"];
 If[TrueQ[metadata["linear_replay"]] && Length[boundary]>=2 && AssociationQ[boundary[[1]]] && FreeQ[Values[boundary[[1]]],lp],
  from=Block[{$configuration=metadata["configuration"]},lineStartParameter[line,boundary[[1]],lp]];AppendTo[samples,{from,boundary}]];
 (* Keep the line symbol in a private numeric callback. Plot substitutions
    cannot rewrite it, and dynamic plotting assignments are localized there. *)
 token=CreateUUID[];$savedMetadata[token]=metadata;$savedSamples[token]=samples;evaluator=Unique["savedEvaluator$"];
 With[{id=token,locals={lp}},evaluator[point_?NumericQ]:=Block[locals,savedEvaluation[id,point]]];
 Table[With[{evaluate=evaluator,r=row,c=column},Function[point,savedComponent[evaluate,r,c,point]]],{row,metadata["dimension"]},{column,Lookup[metadata["configuration"],EpsilonOrder]+1}]
],"DiffExpCompatibility"];
IntegrateSystem[bcs_List,line_] := Module[{saved=TransportTo[bcs,line,1,True]},If[FailureQ[saved],saved,ToPiecewise[saved]]];
IntegrateSystem[line_] := Failure["GeneralSolution",<|"MessageTemplate"->"A numerical boundary is required; symbolic free-constant general solutions are not available in the native interface."|>];
DiffExpFeynmanTrick[configuration_Association,arguments:{___String}:{},OptionsPattern[]] :=
 jsonRun[Join[{"ft","-","--json"},arguments],ExportString[configuration,"RawJSON","Compact"->True],OptionValue["Executable"]];
Options[DiffExpFamilyTemplate]=Options[DiffExpRun];
DiffExpFamilyTemplate[name_String,OptionsPattern[]] := jsonRun[{"family-template",name},"",OptionValue["Executable"]];

End[];
EndPackage[];
