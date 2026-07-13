(* Process-free contract for scoped native SCC-composite owner capacity. *)

repo = DirectoryName[DirectoryName[$InputFileName]];
Get[FileNameJoin[{repo, "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expression_] := Quiet[Catch[expression, "DiffExp2Error"]];

fakeCache[count_Integer?NonNegative] := AssociationThread[
  ("owner-" <> ToString[#] & /@ Range[count]),
  (<|"SyntheticOwner" -> #|> & /@ Range[count])];

cache = DiffExp2`Solve`Private`$nativeSCCCompositeCache;
cacheMax = DiffExp2`Solve`Private`$nativeSCCCompositeCacheMax;
reserved = DiffExp2`Solve`Private`$nativeSCCCompositeReservedCapacity;

DiffExp2`Solve`Private`$nativeSCCCompositeCacheMax = 32;
DiffExp2`Solve`Private`$nativeSCCCompositeReservedCapacity = 0;
DiffExp2`Solve`Private`$nativeSCCCompositeCache = fakeCache[32];

directFailure = catchDE2[
  DiffExp2`Solve`Private`nativeSCCCompositeRequireAdmission[
    <|"Center" -> "unreserved"|>]];
assert["native_scc_direct_unreserved_capacity_guard_remains_strict",
  FailureQ[directFailure] && directFailure["ID"] === "E6" &&
    directFailure["Capacity"] === 32 &&
    directFailure["DefaultCapacity"] === 32 &&
    directFailure["ReservedCapacity"] === 0 &&
    Length[DiffExp2`Solve`Private`$nativeSCCCompositeCache] === 32];

reservedProbe =
  DiffExp2`Solve`WithNativeSCCCompositeCacheReservation[2,
    beforeCapacity =
      DiffExp2`Solve`Private`nativeSCCCompositeEffectiveCapacity[];
    firstAdmission = catchDE2[
      DiffExp2`Solve`Private`nativeSCCCompositeRequireAdmission[
        <|"Center" -> "reserved-33"|>]];
    AssociateTo[DiffExp2`Solve`Private`$nativeSCCCompositeCache,
      "owner-33" -> <|"SyntheticOwner" -> 33|>];
    secondAdmission = catchDE2[
      DiffExp2`Solve`Private`nativeSCCCompositeRequireAdmission[
        <|"Center" -> "reserved-34"|>]];
    AssociateTo[DiffExp2`Solve`Private`$nativeSCCCompositeCache,
      "owner-34" -> <|"SyntheticOwner" -> 34|>];
    overReservation = catchDE2[
      DiffExp2`Solve`Private`nativeSCCCompositeRequireAdmission[
        <|"Center" -> "reserved-35"|>]];
    {beforeCapacity, firstAdmission, secondAdmission, overReservation,
      DiffExp2`Solve`Private`$nativeSCCCompositeReservedCapacity}];

assert["native_scc_scope_reserves_exactly_two_owners_beyond_32",
  reservedProbe[[1]] === 34 && reservedProbe[[2]] === Null &&
    reservedProbe[[3]] === Null && FailureQ[reservedProbe[[4]]] &&
    reservedProbe[[4]]["Capacity"] === 34 &&
    reservedProbe[[5]] === 34 &&
    Length[DiffExp2`Solve`Private`$nativeSCCCompositeCache] === 34 &&
    AllTrue["owner-" <> ToString[#] & /@ Range[34],
      KeyExistsQ[DiffExp2`Solve`Private`$nativeSCCCompositeCache, #] &]];

secondScope =
  DiffExp2`Solve`WithNativeSCCCompositeCacheReservation[2,
    secondCapacity =
      DiffExp2`Solve`Private`nativeSCCCompositeEffectiveCapacity[];
    secondFirst = catchDE2[
      DiffExp2`Solve`Private`nativeSCCCompositeRequireAdmission[
        <|"Center" -> "second-scope-35"|>]];
    AssociateTo[DiffExp2`Solve`Private`$nativeSCCCompositeCache,
      "owner-35" -> <|"SyntheticOwner" -> 35|>];
    secondSecond = catchDE2[
      DiffExp2`Solve`Private`nativeSCCCompositeRequireAdmission[
        <|"Center" -> "second-scope-36"|>]];
    AssociateTo[DiffExp2`Solve`Private`$nativeSCCCompositeCache,
      "owner-36" -> <|"SyntheticOwner" -> 36|>];
    secondOver = catchDE2[
      DiffExp2`Solve`Private`nativeSCCCompositeRequireAdmission[
        <|"Center" -> "second-scope-37"|>]];
    {secondCapacity, secondFirst, secondSecond, secondOver}];
assert["native_scc_back_to_back_scope_reserves_from_live_occupancy",
  secondScope[[1]] === 36 && secondScope[[2]] === Null &&
    secondScope[[3]] === Null && FailureQ[secondScope[[4]]] &&
    secondScope[[4]]["Capacity"] === 36 &&
    Length[DiffExp2`Solve`Private`$nativeSCCCompositeCache] === 36];

assert["native_scc_scope_restores_default_without_evicting_live_owners",
  DiffExp2`Solve`Private`$nativeSCCCompositeReservedCapacity === 0 &&
    DiffExp2`Solve`Private`nativeSCCCompositeEffectiveCapacity[] === 32 &&
    Length[DiffExp2`Solve`Private`$nativeSCCCompositeCache] === 36 &&
    FailureQ[catchDE2[
      DiffExp2`Solve`Private`nativeSCCCompositeRequireAdmission[
        <|"Center" -> "post-scope-direct"|>]]]];

regularMonolithic = <|"IndicialData" -> <|"Regular" -> True|>,
  "IntegrationSequence" -> <|"Components" -> {{1, 2}}|>|>;
regularDecomposed = <|"IndicialData" -> <|"Regular" -> True|>,
  "IntegrationSequence" -> <|"Components" -> {{1}, {2}}|>|>;
singular = <|"IndicialData" -> <|"Regular" -> False|>,
  "IntegrationSequence" -> <|"Components" -> {{1, 2}}|>|>;
assert["native_scc_owner_count_matches_actual_receiving_dispatch",
  !TrueQ[DiffExp2`NativeTransport`Private`nativeReceivingSystemUsesSCCCompositeQ[
      regularMonolithic]] &&
    TrueQ[DiffExp2`NativeTransport`Private`nativeReceivingSystemUsesSCCCompositeQ[
      regularDecomposed]] &&
    TrueQ[DiffExp2`NativeTransport`Private`nativeReceivingSystemUsesSCCCompositeQ[
      singular]]];

escaped = Catch[
  DiffExp2`Solve`WithNativeSCCCompositeCacheReservation[5,
    Throw["escaped", "reservation-scope"]], "reservation-scope"];
assert["native_scc_reservation_restores_after_nonlocal_exit",
  escaped === "escaped" &&
    DiffExp2`Solve`Private`$nativeSCCCompositeReservedCapacity === 0];

DiffExp2`Solve`Private`$nativeSCCCompositeReservedCapacity = 99;
DiffExp2`Solve`DropWolframPreparationCaches[];
assert["native_scc_drop_wolfram_caches_resets_reservation",
  DiffExp2`Solve`Private`$nativeSCCCompositeCache === <||> &&
    DiffExp2`Solve`Private`$nativeSCCCompositeReservedCapacity === 0 &&
    DiffExp2`Solve`Private`nativeSCCCompositeEffectiveCapacity[] === 32];

DiffExp2`Solve`Private`$nativeSCCCompositeCache = fakeCache[1];
DiffExp2`Solve`Private`$nativeSCCCompositeReservedCapacity = 99;
DiffExp2`Solve`ClearSolveCaches[];
assert["native_scc_clear_solve_caches_resets_reservation",
  DiffExp2`Solve`Private`$nativeSCCCompositeCache === <||> &&
    DiffExp2`Solve`Private`$nativeSCCCompositeReservedCapacity === 0 &&
    DiffExp2`Solve`Private`nativeSCCCompositeEffectiveCapacity[] === 32];

DiffExp2`Solve`Private`$nativeSCCCompositeCache = cache;
DiffExp2`Solve`Private`$nativeSCCCompositeCacheMax = cacheMax;
DiffExp2`Solve`Private`$nativeSCCCompositeReservedCapacity = reserved;

Print["Native SCC cache reservation tests: ", passed, " passed, ",
  failed, " failed."];
If[failed > 0, Exit[1]];
