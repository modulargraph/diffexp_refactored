#include "diffexp/kernel/json_codec.hpp"

#include "json_codec_capabilities.hpp"

#include <boost/json.hpp>
#include <flint/flint.h>

namespace diffexp::kernel {
namespace json = boost::json;

using namespace json_codec_detail;

std::string backend_info_json() {
  return json::serialize(json::object{
      {"schema", 1},
      {"schemas", json::array{1, 2}},
      {"persistent_sessions", true},
      {"persistent_local_solutions", true},
      {"persistent_frame_independent_regular_equation_owners", true},
      {"persistent_frame_independent_regular_equation_owner_capability",
       kFrameIndependentRegularEquationOwnerCapability},
      {"persistent_frame_independent_regular_equation_owner_checkpoint",
       false},
      {"persistent_local_residual_certification", true},
      {"persistent_scc_prepare", true},
      {"persistent_scc_execute", false},
      {"persistent_scc_scalar_block_dag_column", true},
      {"persistent_scc_regular_block_dag_column", true},
      {"persistent_exact_regular_local_match", true},
      {"persistent_endpoint_limits", true},
      {"persistent_endpoint_limit_capability",
       kRetainedEndpointLimitCapability},
      {"persistent_plan_bound_endpoint_limits", true},
      {"persistent_plan_bound_endpoint_limit_capability",
       kRetainedPlannedEndpointLimitCapability},
      {"persistent_symbolic_endpoint_limits", false},
      {"persistent_exact_tile_plans", true},
      {"persistent_exact_tile_plan_capability", kRetainedTilePlanCapability},
      {"persistent_exact_single_arm_tile_plans", true},
      {"persistent_exact_single_arm_tile_plan_capability",
       kRetainedSingleArmTilePlanCapability},
      {"persistent_plan_driven_match_hop", true},
      {"persistent_plan_driven_match_hop_capability",
       kRetainedPlannedMatchHopCapability},
      {"persistent_plan_match_local_materialization", true},
      {"persistent_plan_match_local_materialization_capability",
       kRetainedPlannedMatchMaterializationCapability},
      {"persistent_stored_line_integration", true},
      {"persistent_stored_line_integration_capability",
       kRetainedStoredLineCapability},
      {"persistent_parallel_arm_march", true},
      {"persistent_parallel_arm_march_capability",
       kRetainedParallelArmCapability},
      {"persistent_parallel_transport_arm_states", true},
      {"persistent_parallel_transport_arm_states_capability",
       kRetainedParallelTransportArmStateCapability},
      {"persistent_transport_arm_state", true},
      {"persistent_transport_arm_state_capability",
       kRetainedTransportArmStateCapability},
      {"persistent_transport_arm_contraction", true},
      {"persistent_transport_arm_contraction_capability",
       kRetainedTransportArmContractionCapability},
      {"persistent_transport_pair_contraction", true},
      {"persistent_transport_pair_contraction_capability",
       kRetainedTransportPairContractionCapability},
      {"persistent_transport_pair_tile_stream", true},
      {"persistent_transport_pair_tile_stream_capability",
       kRetainedTransportPairStreamCapability},
      {"persistent_transport_endpoint_batch", true},
      {"persistent_transport_endpoint_batch_capability",
       kRetainedTransportEndpointBatchCapability},
      {"persistent_certified_tail_majorant", true},
      {"persistent_certified_tail_majorant_capability",
       kRegularTailMajorantCapability},
      {"persistent_certified_line_integration_capability",
       kRetainedCertifiedLineCapability},
      {"persistent_symbolic_line_integration", false},
      {"persistent_acb_local_match", true},
      {"persistent_acb_local_match_capability",
       kRefinedAcbLocalMatchCapability},
      {"persistent_checkpoint", true},
      {"persistent_checkpoint_schema", kCheckpointPayloadSchema},
      {"persistent_checkpoint_handle_scope",
       "complete-retained-native-ownership-closure"},
      {"persistent_scc_regular_singular_scalar_block_dag_column", true},
      {"persistent_scc_regular_singular_jordan_block_dag_column", true},
      {"persistent_scc_pseudo_compensation", false},
      {"backend", "DiffExp3 C++"},
      {"flint", flint_version},
      {"librarylink", false},
      {"standalone_cpp", true}});
}

}  // namespace diffexp::kernel
