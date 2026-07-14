#pragma once

#include <cstdint>

namespace diffexp2::json_codec_detail {

inline constexpr char kRegularTailMajorantCapability[] =
    "retained-regular-homogeneous-gronwall-cauchy-tail-v1";
inline constexpr char kOwnerBoundPhysicalResidualCapability[] =
    "owner-bound-physical-homogeneous-residual-v2";
inline constexpr char kRetainedEndpointLimitCapability[] =
    "retained-native-endpoint-sector-limit-v1";
inline constexpr char kRetainedPlannedEndpointLimitCapability[] =
    "retained-native-plan-bound-endpoint-sector-limit-v1";
inline constexpr char kRetainedTransportEndpointBatchCapability[] =
    "retained-native-transport-endpoint-batch-v1";
inline constexpr char kRetainedRationalRowCapability[] =
    "retained-native-rational-row-local-application-v1";
inline constexpr char kExactRegularLocalMatchCapability[] =
    "exact-rational-regular-local-match-v1";
inline constexpr char kRefinedAcbLocalMatchCapability[] =
    "exact-lattice-guided-acb-local-match-v1";
inline constexpr char kExactEvaluatedLatticeSchema[] =
    "diffexp2-exact-evaluated-epsilon-lattice-v1";
inline constexpr char kNativeUnitSaturationRequestSchema[] =
    "diffexp2-native-acb-unit-leading-saturation-request-v1";
inline constexpr char kNativeUnitSaturationCompactRequestSchema[] =
    "diffexp2-native-acb-unit-leading-saturation-request-v2";
inline constexpr char kNativeUnitSaturationProofSchema[] =
    "diffexp2-native-acb-unit-leading-saturation-proof-v1";
inline constexpr char kNativeSingularSCCSaturationRequestSchema[] =
    "diffexp2-native-acb-singular-scc-valuation-zero-saturation-request-v1";
inline constexpr char kNativeSingularSCCSaturationCompactRequestSchema[] =
    "diffexp2-native-acb-singular-scc-valuation-zero-saturation-request-v2";
inline constexpr char kNativeSingularSCCSaturationProofSchema[] =
    "diffexp2-native-acb-singular-scc-valuation-zero-saturation-proof-v1";
inline constexpr char kAcbSingularScalarSCCColumnCapability[] =
    "acb-regular-singular-scalar-block-dag-column-v1";
inline constexpr char kAcbSingularJordanSCCColumnCapability[] =
    "acb-regular-singular-jordan-block-dag-column-v1";

inline constexpr char kRetainedTilePlanCapability[] =
    "retained-exact-independent-arm-tile-plan-v1";
inline constexpr char kRetainedSingleArmTilePlanCapability[] =
    "retained-exact-single-arm-tile-plan-v1";
inline constexpr char kRetainedSingleArmTilePlanCheckpointSchema[] =
    "diffexp2-retained-single-arm-tile-plan-v1";
inline constexpr char kRetainedSingleArmTilePlanProvenanceSchema[] =
    "diffexp2-retained-exact-single-arm-tile-plan-v1";
inline constexpr char kRetainedPlannedMatchHopCapability[] =
    "retained-exact-plan-driven-local-match-hop-v1";
inline constexpr char kRetainedPlannedMatchMaterializationCapability[] =
    "retained-native-plan-match-local-materialization-v1";
inline constexpr char kRetainedStoredLineCapability[] =
    "retained-native-stored-truncation-physical-tile-integral-v1";
inline constexpr char kRetainedCertifiedLineCapability[] =
    "retained-native-certified-full-local-physical-tile-integral-v1";
inline constexpr char kRetainedParallelArmCapability[] =
    "retained-native-concurrent-two-arm-march-v1";
inline constexpr char kRetainedParallelTransportArmStateCapability[] =
    "retained-native-concurrent-two-arm-transport-states-v1";
inline constexpr char kRetainedTransportArmStateCapability[] =
    "retained-native-transport-arm-state-v1";
inline constexpr char kRetainedTransportArmContractionCapability[] =
    "retained-native-transport-arm-contraction-v1";
inline constexpr char kRetainedTransportPairContractionCapability[] =
    "retained-native-transport-pair-contraction-v1";
inline constexpr char kRetainedTransportPairStreamCapability[] =
    "retained-native-transport-pair-tile-stream-v1";
inline constexpr char kRetainedLineAggregateCapability[] =
    "retained-native-line-aggregate-v1";

inline constexpr char kCheckpointFormat[] =
    "diffexp2-persistent-native-session";
inline constexpr std::uint32_t kCheckpointPayloadSchema = 8;

}  // namespace diffexp2::json_codec_detail
