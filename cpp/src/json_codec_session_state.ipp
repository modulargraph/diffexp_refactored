class TransportPairObservableStream;

struct SolverSession {
  std::string handle;
  std::string domain;
  slong precision_bits = 256;
  int output_digits = 50;
  std::vector<std::string> symbols;
  std::string analytic_identity;
  std::size_t chart_capacity = 256;
  std::size_t local_capacity = 1024;
  std::size_t scc_capacity = 128;
  std::size_t match_capacity = 1024;
  std::size_t endpoint_capacity = 1024;
  std::size_t tile_plan_capacity = 256;
  std::size_t transport_state_capacity = 256;
  std::size_t line_result_capacity = 2048;
  std::uint64_t next_chart = 1;
  std::uint64_t next_local = 1;
  std::uint64_t next_scc = 1;
  std::uint64_t next_match = 1;
  std::uint64_t next_endpoint = 1;
  std::uint64_t next_tile_plan = 1;
  std::uint64_t next_transport_state = 1;
  std::uint64_t next_line_result = 1;
  std::uint64_t next_transport_pair_stream = 1;
  std::size_t pending_local_solves = 0;
  std::size_t pending_matches = 0;
  std::size_t pending_endpoint_limits = 0;
  std::size_t pending_tile_plans = 0;
  std::size_t pending_transport_states = 0;
  std::size_t pending_line_integrations = 0;
  std::uint64_t total_local_solves = 0;
  std::uint64_t total_scc_column_solves = 0;
  std::uint64_t total_local_matches = 0;
  std::uint64_t total_endpoint_limits = 0;
  std::uint64_t total_endpoint_exports = 0;
  std::uint64_t total_tile_plans = 0;
  std::uint64_t total_transport_arm_marches = 0;
  std::uint64_t total_transport_contractions = 0;
  std::uint64_t total_transport_observables = 0;
  std::uint64_t total_transport_pair_contractions = 0;
  std::uint64_t total_transport_pair_observables = 0;
  std::uint64_t total_transport_endpoint_batches = 0;
  std::uint64_t total_transport_endpoint_rows = 0;
  std::uint64_t total_line_integrations = 0;
  std::uint64_t total_line_exports = 0;
  double total_local_run_parse_ms = 0.0;
  double total_local_kernel_ms = 0.0;
  double total_local_match_ms = 0.0;
  std::uint64_t checkpoint_generation = 0;
  std::uint64_t checkpoint_restore_count = 0;
  std::string restored_from_checkpoint_identity;
  double total_endpoint_limit_ms = 0.0;
  double total_endpoint_export_ms = 0.0;
  double total_tile_plan_ms = 0.0;
  double total_transport_arm_ms = 0.0;
  double total_transport_contraction_ms = 0.0;
  double total_transport_pair_contraction_ms = 0.0;
  double total_transport_endpoint_batch_ms = 0.0;
  double total_line_integration_ms = 0.0;
  double total_line_export_ms = 0.0;
  bool closed = false;
  mutable std::mutex mutex;
  std::unordered_map<std::string, std::shared_ptr<PreparedChartBase>> charts;
  std::unordered_map<std::string, std::string> handles_by_key;
  std::unordered_map<std::string, std::shared_ptr<StoredLocalBase>> locals;
  std::unordered_map<std::string, std::shared_ptr<StoredMatchBase>> matches;
  std::unordered_map<std::string, std::shared_ptr<StoredEndpointResult>>
      endpoints;
  std::unordered_map<std::string, std::shared_ptr<StoredTilePlan>> tile_plans;
  std::unordered_map<std::string, std::shared_ptr<StoredTransportArmState>>
      transport_states;
  std::unordered_map<std::string, std::shared_ptr<StoredLineResult>>
      line_results;
  // Ephemeral, never checkpointed.  A live stream owns one pending line
  // reservation and the two exact transport-state owners until finish,
  // abort, or session close.
  std::unordered_map<std::string,
                     std::shared_ptr<TransportPairObservableStream>>
      transport_pair_streams;
  std::unordered_map<std::string, std::shared_ptr<CompositeSCCChartBase>> sccs;
  std::unordered_map<std::string, std::string> scc_handles_by_key;
};

#include "json_codec_transport_algorithms.ipp"
std::vector<std::uint32_t> parse_index_vector(
    const json::value& raw, std::uint32_t bound, const char* label) {
  std::vector<std::uint32_t> values;
  for (const auto& value : as_array(raw, label)) {
    const auto index = as_u32(value, label);
    if (index >= bound)
      throw std::invalid_argument(std::string(label) + " is outside range");
    values.push_back(index);
  }
  return values;
}

struct ExactParentMatrix {
  struct Cell {
    std::string exact;
    bool proven_zero = false;
  };
  std::uint32_t dimension = 0;
  std::vector<Cell> cells;  // row (target), then column (source)
  std::string canonical_record;

  const Cell& at(std::uint32_t row, std::uint32_t column) const {
    return cells.at(static_cast<std::size_t>(row) * dimension + column);
  }
};

ExactParentMatrix parse_exact_parent_matrix(
    const json::value& raw, std::uint32_t dimension, const char* label) {
  const auto& rows = as_array(raw, label);
  if (rows.size() != dimension)
    throw std::invalid_argument(
        std::string(label) + " must have parent-dimension rows");
  ExactParentMatrix result;
  result.dimension = dimension;
  result.cells.reserve(static_cast<std::size_t>(dimension) * dimension);
  json::array canonical_rows;
  canonical_rows.reserve(dimension);
  for (std::uint32_t row = 0; row < dimension; ++row) {
    const auto& columns = as_array(rows[row], label);
    if (columns.size() != dimension)
      throw std::invalid_argument(
          std::string(label) + " must be a square parent-dimension matrix");
    json::array canonical_columns;
    canonical_columns.reserve(dimension);
    for (std::uint32_t column = 0; column < dimension; ++column) {
      const auto& cell = as_object(columns[column], "exact parent matrix cell");
      const auto exact = required_string(cell, "exact");
      const auto proven_zero = cell.at("proven_zero").as_bool();
      result.cells.push_back({exact, proven_zero});
      canonical_columns.push_back(json::object{
          {"exact", exact}, {"proven_zero", proven_zero}});
    }
    canonical_rows.push_back(std::move(canonical_columns));
  }
  result.canonical_record = json::serialize(canonical_rows);
  return result;
}

std::string canonical_exact_submatrix(
    const ExactParentMatrix& matrix,
    const std::vector<std::uint32_t>& vertices) {
  json::array rows;
  rows.reserve(vertices.size());
  for (const auto target : vertices) {
    json::array columns;
    columns.reserve(vertices.size());
    for (const auto source : vertices) {
      const auto& cell = matrix.at(target, source);
      columns.push_back(json::object{
          {"exact", cell.exact}, {"proven_zero", cell.proven_zero}});
    }
    rows.push_back(std::move(columns));
  }
  return json::serialize(rows);
}

std::string derived_coupling_identity(
    std::uint32_t source_block, std::uint32_t target_block,
    const std::vector<std::uint32_t>& source_vertices,
    const std::vector<std::uint32_t>& target_vertices,
    const ExactParentMatrix& original, const ExactParentMatrix& theta) {
  auto rectangular_record = [&](const ExactParentMatrix& matrix) {
    json::array rows;
    rows.reserve(target_vertices.size());
    for (const auto target : target_vertices) {
      json::array columns;
      columns.reserve(source_vertices.size());
      for (const auto source : source_vertices) {
        const auto& cell = matrix.at(target, source);
        columns.push_back(json::object{
            {"exact", cell.exact}, {"proven_zero", cell.proven_zero}});
      }
      rows.push_back(std::move(columns));
    }
    return rows;
  };
  return json::serialize(json::object{
      {"schema", "diffexp2-native-scc-coupling-v1"},
      {"source_block", source_block}, {"target_block", target_block},
      {"source_vertices", encode_indices(source_vertices)},
      {"target_vertices", encode_indices(target_vertices)},
      {"exact_original", rectangular_record(original)},
      {"exact_theta", rectangular_record(theta)}});
}

template <typename Scalar>
std::shared_ptr<CompositeSCCChartBase> parse_composite_scc_chart(
    const std::shared_ptr<SolverSession>& session, const json::object& root,
    const std::string& handle, const std::string& key,
    const std::string& exact_identity, std::string signature,
    const std::vector<std::shared_ptr<PreparedChartBase>>& erased_charts) {
  // Composite geometry is retained as an exact-rational Acb ball even when
  // recurrence coefficients themselves are exact.  Parse it under the same
  // guarded precision later used by block LocalSolutions.
  AcbPrecisionLease acb_lease(session->precision_bits);
  ComplexBall::set_precision(session->precision_bits);
  std::unique_lock<std::recursive_mutex> symbolic_lock;
  if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
    symbolic_lock =
        std::unique_lock<std::recursive_mutex>(symbolic_run_mutex());
    SymbolicRational::configure(session->symbols);
  }

  const auto& parent = as_object(root.at("parent"), "SCC parent manifest");
  const auto dimension = as_u32(parent.at("dimension"), "parent dimension");
  if (dimension == 0)
    throw std::invalid_argument("SCC parent dimension must be positive");
  std::shared_ptr<const PreparedPhysicalClearedODE<Scalar>>
      physical_equation;
  if (const auto* raw_physical = root.if_contains("physical_ode"))
    physical_equation = parse_prepared_physical_ode<Scalar>(
        *raw_physical, dimension, false);
  const auto exact_system = parse_exact_parent_matrix(
      parent.at("exact_system_record"), dimension,
      "exact parent system record");
  const auto exact_theta = parse_exact_parent_matrix(
      parent.at("exact_theta_record"), dimension,
      "exact parent theta record");
  validate_first_slice_rational_geometry(parent.at("chart"));
  const auto geometry_record = canonical_chart_geometry_record(
      parent.at("chart"));
  auto retained_geometry = parse_retained_composite_geometry(
      parent.at("chart"));
  auto graph = validate_scc_certificate(parent.at("scc"), dimension);
  const std::set<std::pair<std::uint32_t, std::uint32_t>> structural_edges(
      graph.structural_edges.begin(), graph.structural_edges.end());
  for (std::uint32_t target = 0; target < dimension; ++target) {
    for (std::uint32_t source = 0; source < dimension; ++source) {
      const auto original_zero = exact_system.at(target, source).proven_zero;
      const auto theta_zero = exact_theta.at(target, source).proven_zero;
      if (original_zero != theta_zero)
        throw std::invalid_argument(
            "parent original/theta structural-zero facts disagree");
      const auto graph_zero = !structural_edges.contains({source, target});
      if (original_zero != graph_zero)
        throw std::invalid_argument(
            "parent exact matrix zero facts do not reproduce the structural graph");
    }
  }

  const auto& execution = as_object(
      parent.at("execution"), "SCC execution contract");
  if (required_string(execution, "mode") != "BlockSequentialStrict")
    throw std::invalid_argument(
        "native SCC preparation only supports BlockSequentialStrict");
  CompositeWorkContract work;
  work.work_t_order = as_u32(
      execution.at("work_t_order"), "SCC work Taylor order");
  const auto& work_object = as_object(
      parent.at("work_contract"), "SCC work contract");
  work.work_min = as_i32(
      work_object.at("work_min"), "work epsilon minimum");
  if (const auto* audit = work_object.if_contains(
          "cancellation_audit_base"); audit != nullptr && !audit->is_null())
    work.cancellation_audit_base =
        as_i32(*audit, "work cancellation audit base");
  work.requested_min = as_i32(
      work_object.at("requested_min"), "requested epsilon minimum");
  work.requested_max = as_i32(
      work_object.at("requested_max"), "requested epsilon maximum");
  work.work_complete_max = as_i32(
      work_object.at("work_complete_max"), "work epsilon maximum");
  work.public_t_order = as_u32(
      work_object.at("public_t_order"), "public Taylor order");
  work.wolfram_coupling_depth = as_u32(
      work_object.at("wolfram_coupling_depth"),
      "Wolfram coupling depth");
  if (work.work_min > work.requested_min ||
      work.requested_min > work.requested_max ||
      work.requested_max > work.work_complete_max)
    throw std::invalid_argument(
        "SCC epsilon work contract has inconsistent ordered bounds");
  if (work.cancellation_audit_base.has_value() &&
      !(work.work_min < *work.cancellation_audit_base &&
        *work.cancellation_audit_base <= work.requested_min))
    throw std::invalid_argument(
        "SCC cancellation audit base is outside its lower guard contract");
  if (work.wolfram_coupling_depth != graph.coupling_depth + 1)
    throw std::invalid_argument(
        "Wolfram coupling depth must equal native edge depth plus one");
  const auto expected_work_t_order =
      static_cast<std::uint64_t>(work.public_t_order) + 2 +
      2 * static_cast<std::uint64_t>(work.wolfram_coupling_depth);
  if (expected_work_t_order > std::numeric_limits<std::uint32_t>::max() ||
      work.work_t_order != expected_work_t_order)
    throw std::invalid_argument(
        "SCC work Taylor order does not match the exact depth budget");
  const auto frame_width_i64 =
      static_cast<std::int64_t>(work.work_complete_max) -
      work.work_min + 1;
  if (frame_width_i64 <= 0 ||
      frame_width_i64 > std::numeric_limits<std::uint32_t>::max())
    throw std::invalid_argument("SCC work frame width is invalid");
  const auto frame_width = static_cast<std::uint32_t>(frame_width_i64);

  const auto& raw_blocks = as_array(root.at("blocks"), "SCC blocks");
  if (raw_blocks.size() != graph.component_count ||
      erased_charts.size() != raw_blocks.size())
    throw std::invalid_argument(
        "SCC preparation requires exactly one chart per component");
  std::vector<CompositeSCCBlock<Scalar>> blocks(graph.component_count);
  std::vector<std::uint8_t> block_seen(graph.component_count, 0);
  std::set<std::string> chart_handles;
  for (std::size_t index = 0; index < raw_blocks.size(); ++index) {
    const auto& raw_block = as_object(raw_blocks[index], "SCC block");
    const auto block = as_u32(raw_block.at("block"), "SCC block index");
    if (block != index || block >= graph.component_count || block_seen[block])
      throw std::invalid_argument(
          "SCC blocks must be in deterministic component order");
    block_seen[block] = 1;
    const auto declared_capabilities =
        canonical_native_scc_capabilities(raw_block);
    const bool regular = raw_block.at("regular").as_bool();
    const bool no_pseudo = raw_block.at("no_pseudo").as_bool();
    const bool identity_v = raw_block.at("identity_v").as_bool();
    const bool epsilon_unimodular_v =
        raw_block.if_contains("epsilon_unimodular_v") != nullptr
            ? raw_block.at("epsilon_unimodular_v").as_bool()
            : identity_v;
    const bool identity_gauge = raw_block.at("identity_gauge").as_bool();
    const bool exact_gauge = raw_block.if_contains("exact_gauge") != nullptr
        ? raw_block.at("exact_gauge").as_bool()
        : identity_gauge;
    if (!exact_gauge)
      throw std::invalid_argument(
          "native SCC preparation requires an exact invertible gauge");
    if (!epsilon_unimodular_v)
      throw std::invalid_argument(
          "native SCC preparation requires an exact t-independent Laurent-unimodular spectral frame");
    // `no_pseudo` is retained as producer provenance, not an admission
    // decision.  Exact Rational execution reconstructs every T/P/R branch
    // from the retained affine Jordan certificate; a producer may therefore
    // truthfully advertise false here without disabling native execution.
    // Rational execution may compensate a revalidated CASE-P event.  Acb
    // execution below admits only the producer-proved no-collision class and
    // still reconstructs every submitted schedule from an exact certificate.

    auto vertices = parse_index_vector(
        raw_block.at("vertices"), dimension, "SCC block vertex");
    if (vertices != graph.components[block])
      throw std::invalid_argument(
          "SCC block vertices do not equal their parent component in exact order");

    // This is the only type-erasure boundary for a composite.  Every later
    // execution method owns typed pointers directly and performs no casts.
    auto chart = std::dynamic_pointer_cast<PreparedChart<Scalar>>(
        erased_charts[index]);
    if (!chart)
      throw std::invalid_argument(
          "SCC block chart scalar domain differs from its session");
    const auto source_handle = required_string(raw_block, "chart");
    if (chart->handle() != source_handle ||
        !chart_handles.insert(source_handle).second)
      throw std::invalid_argument(
          "SCC block chart handles must be exact and one-to-one");
    const auto principal_identity = required_string(
        raw_block, "principal_identity");
    if (chart->exact_identity() != principal_identity)
      throw std::invalid_argument(
          "SCC principal identity differs from its prepared chart identity");
    if (!chart->native_scc_capabilities().has_value())
      throw std::invalid_argument(
          "SCC block chart lacks analytic.native_scc_capabilities metadata");
    if (*chart->native_scc_capabilities() != declared_capabilities)
      throw std::invalid_argument(
          "SCC block capability claims differ from the retained chart metadata");
    if (!chart->has_assembly())
      throw std::invalid_argument(
          "SCC spectral-frame capability requires a retained assembly operator");
    if (chart->has_identity_assembly() != identity_v)
      throw std::invalid_argument(
          "SCC identity_v capability contradicts the retained assembly operator");
    if (!chart->geometry_record().has_value())
      throw std::invalid_argument(
          "SCC block chart lacks analytic.geometry preparation metadata");
    if (*chart->geometry_record() != geometry_record)
      throw std::invalid_argument(
          "SCC block chart geometry differs from the parent manifest");
    if (!chart->principal_matrix_record().has_value())
      throw std::invalid_argument(
          "SCC block chart lacks analytic.principal_matrix metadata");
    if (*chart->principal_matrix_record() !=
        canonical_exact_submatrix(exact_system, vertices))
      throw std::invalid_argument(
          "SCC principal chart matrix differs from the indexed parent submatrix");
    if (chart->dimension() != vertices.size())
      throw std::invalid_argument(
          "SCC block chart dimension differs from its vertex count");
    if (chart->frame_base() != work.work_min ||
        chart->frame_width() != frame_width)
      throw std::invalid_argument(
          "SCC block chart frame differs from the exact work contract");

    CompositeSpectralSourceTransform<Scalar> source_transform;
    source_transform.identity = true;
    source_transform.epsilon_unimodular = true;
    source_transform.matrix.rows = chart->dimension();
    source_transform.matrix.columns = chart->dimension();
    if (const auto* raw_transform =
            raw_block.if_contains("source_transform")) {
      source_transform = parse_composite_spectral_source_transform<Scalar>(
          *raw_transform, chart->dimension(), frame_width, work,
          session->domain, session->symbols, identity_v,
          chart->assembly_exact_identity());
    } else if (!identity_v) {
      throw std::invalid_argument(
          "nonidentity SCC spectral assembly requires a certified target-source VInv transform");
    }

    CompositeGaugeTransform<Scalar> to_physical, to_reduced;
    to_physical.identity = identity_gauge;
    to_reduced.identity = identity_gauge;
    to_physical.role = "to_physical";
    to_reduced.role = "to_reduced";
    to_physical.matrix.rows = to_reduced.matrix.rows = chart->dimension();
    to_physical.matrix.columns = to_reduced.matrix.columns = chart->dimension();
    if (!identity_gauge) {
      const auto* raw_to_physical = raw_block.if_contains("to_physical");
      const auto* raw_to_reduced = raw_block.if_contains("to_reduced");
      if (!raw_to_physical || !raw_to_reduced)
        throw std::invalid_argument(
            "nonidentity SCC gauge requires both exact directional transforms");
      to_physical = parse_composite_gauge_transform<Scalar>(
          *raw_to_physical, "to_physical", chart->dimension(), frame_width,
          work, session->domain, session->symbols);
      to_reduced = parse_composite_gauge_transform<Scalar>(
          *raw_to_reduced, "to_reduced", chart->dimension(), frame_width,
          work, session->domain, session->symbols);
      if (to_physical.identity || to_reduced.identity ||
          to_physical.gauge_exact_identity != to_reduced.gauge_exact_identity ||
          to_physical.gauge_inverse_exact_identity !=
              to_reduced.gauge_inverse_exact_identity ||
          to_physical.gauge_det_exact_identity !=
              to_reduced.gauge_det_exact_identity)
        throw std::invalid_argument(
            "directional SCC gauge transforms do not share one exact frame certificate");
    } else if (raw_block.if_contains("to_physical") ||
               raw_block.if_contains("to_reduced")) {
      throw std::invalid_argument(
          "identity SCC gauge must not carry redundant directional transforms");
    }

    std::optional<ExactJordanIndicialCertificate> exact_indicial;
    if (const auto* raw_indicial =
            raw_block.if_contains("exact_affine_jordan_indicial"))
      exact_indicial = parse_exact_jordan_indicial_record(
          *raw_indicial, chart->dimension());
    if (const auto& retained = chart->exact_jordan_indicial();
        retained.has_value()) {
      if (exact_indicial.has_value() &&
          !same_exact_jordan_indicial(*exact_indicial, *retained))
        throw std::invalid_argument(
            "SCC block exact indicial manifest differs from its retained Rational operator certificate");
      exact_indicial = *retained;
    }
    if (exact_indicial.has_value() &&
        !chart->jordan_partition_matches(*exact_indicial))
      throw std::invalid_argument(
          "SCC block exact indicial certificate differs from its prepared Jordan partition");

    const auto& local_graph = chart->scc();
    std::vector<std::uint32_t> local_vertices(vertices.size());
    for (std::uint32_t local = 0; local < vertices.size(); ++local)
      local_vertices[local] = local;
    if (local_graph.component_count != 1 ||
        local_graph.components.size() != 1 ||
        local_graph.components.front() != local_vertices ||
        local_graph.coupling_depth != 0)
      throw std::invalid_argument(
          "SCC principal chart must be exactly one retained component");
    std::vector<std::uint32_t> local_of(dimension,
        std::numeric_limits<std::uint32_t>::max());
    for (std::uint32_t local = 0; local < vertices.size(); ++local)
      local_of[vertices[local]] = local;
    std::set<std::pair<std::uint32_t, std::uint32_t>> expected_intra;
    for (const auto [source, target] : graph.structural_edges)
      if (graph.component_of[source] == block &&
          graph.component_of[target] == block)
        expected_intra.insert({local_of[source], local_of[target]});
    const std::set<std::pair<std::uint32_t, std::uint32_t>> supplied_intra(
        local_graph.structural_edges.begin(),
        local_graph.structural_edges.end());
    if (supplied_intra != expected_intra)
      throw std::invalid_argument(
          "SCC principal chart graph does not cover its exact intra-component edges");

    blocks[block] = CompositeSCCBlock<Scalar>{
        block, std::move(vertices), source_handle, principal_identity,
        regular, no_pseudo, std::move(source_transform),
        std::move(to_physical), std::move(to_reduced),
        std::move(exact_indicial), std::move(chart)};
  }

  std::set<std::pair<std::uint32_t, std::uint32_t>> expected_cross;
  for (const auto edge : graph.structural_edges)
    if (graph.component_of[edge.first] != graph.component_of[edge.second])
      expected_cross.insert(edge);
  const std::set<std::pair<std::uint32_t, std::uint32_t>>
      expected_condensation(graph.condensation_edges.begin(),
                            graph.condensation_edges.end());
  std::set<std::pair<std::uint32_t, std::uint32_t>> supplied_condensation;
  std::set<std::pair<std::uint32_t, std::uint32_t>> active_edges;
  std::set<std::pair<std::uint32_t, std::uint32_t>> prepared_pairs;
  std::vector<CompositeSCCCoupling<Scalar>> couplings;
  const auto& raw_couplings = as_array(
      root.at("couplings"), "SCC couplings");
  if (raw_couplings.size() != graph.condensation_edges.size())
    throw std::invalid_argument(
        "SCC coupling groups must equal the condensation edge count");
  couplings.reserve(raw_couplings.size());
  for (std::size_t coupling_index = 0;
       coupling_index < raw_couplings.size(); ++coupling_index) {
    const auto& raw_value = raw_couplings[coupling_index];
    const auto& raw = as_object(raw_value, "SCC coupling");
    const auto source_block = as_u32(
        raw.at("source_block"), "source SCC block");
    const auto target_block = as_u32(
        raw.at("target_block"), "target SCC block");
    if (source_block >= graph.component_count ||
        target_block >= graph.component_count || source_block == target_block)
      throw std::invalid_argument("invalid cross-SCC coupling block pair");
    const auto block_pair = std::make_pair(source_block, target_block);
    if (block_pair != graph.condensation_edges[coupling_index])
      throw std::invalid_argument(
          "SCC coupling groups must use deterministic condensation-edge order");
    if (!supplied_condensation.insert(block_pair).second ||
        !expected_condensation.contains(block_pair))
      throw std::invalid_argument(
          "SCC coupling groups must match condensation edges one-to-one");
    const auto source_vertices = parse_index_vector(
        raw.at("source_vertices"), dimension, "coupling source vertex");
    const auto target_vertices = parse_index_vector(
        raw.at("target_vertices"), dimension, "coupling target vertex");
    if (source_vertices != blocks[source_block].vertices ||
        target_vertices != blocks[target_block].vertices)
      throw std::invalid_argument(
          "SCC coupling vertex bases differ from their block manifests");
    const auto columns = as_u32(raw.at("columns"), "coupling columns");
    const auto rows = as_u32(raw.at("rows"), "coupling rows");
    if (columns != source_vertices.size() || rows != target_vertices.size())
      throw std::invalid_argument(
          "SCC coupling matrix dimensions disagree with its blocks");
    if (raw.if_contains("symbols") == nullptr ||
        required_string(raw, "domain") != session->domain ||
        parse_symbols(raw) != session->symbols)
      throw std::invalid_argument(
          "SCC coupling coefficient field differs from its session");

    CompositeSCCCoupling<Scalar> coupling;
    coupling.source_block = source_block;
    coupling.target_block = target_block;
    coupling.producer_identity = required_string(raw, "exact_identity");
    coupling.matrix.rows = rows;
    coupling.matrix.columns = columns;
    coupling.matrix.exact_identity = derived_coupling_identity(
        source_block, target_block, source_vertices, target_vertices,
        exact_system, exact_theta);
    std::optional<std::pair<std::uint32_t, std::uint32_t>> previous_local;
    for (const auto& raw_entry_value : as_array(
             raw.at("entries"), "SCC coupling entries")) {
      const auto& raw_entry = as_object(
          raw_entry_value, "SCC coupling entry");
      const auto row = as_u32(raw_entry.at("row"), "coupling row");
      const auto column = as_u32(
          raw_entry.at("column"), "coupling column");
      if (row >= rows || column >= columns)
        throw std::invalid_argument("SCC coupling entry is out of range");
      const auto local_pair = std::make_pair(row, column);
      if (previous_local.has_value() && *previous_local >= local_pair)
        throw std::invalid_argument(
            "SCC coupling entries must use deterministic row-column order");
      previous_local = local_pair;
      const auto source_vertex = as_u32(
          raw_entry.at("source_vertex"), "coupling source vertex");
      const auto target_vertex = as_u32(
          raw_entry.at("target_vertex"), "coupling target vertex");
      if (source_vertex != source_vertices[column] ||
          target_vertex != target_vertices[row] ||
          graph.component_of[source_vertex] != source_block ||
          graph.component_of[target_vertex] != target_block)
        throw std::invalid_argument(
            "SCC coupling local/global vertex identities disagree");
      const auto edge = std::make_pair(source_vertex, target_vertex);
      if (!prepared_pairs.insert(edge).second)
        throw std::invalid_argument(
            "SCC coupling contains a duplicate global matrix entry");
      const auto exact_original_entry = required_string(
          raw_entry, "exact_original_entry");
      const auto exact_theta_entry = required_string(
          raw_entry, "exact_theta_entry");
      const auto& parent_original = exact_system.at(
          target_vertex, source_vertex);
      const auto& parent_theta = exact_theta.at(
          target_vertex, source_vertex);
      if (exact_original_entry != parent_original.exact ||
          exact_theta_entry != parent_theta.exact)
        throw std::invalid_argument(
            "SCC coupling exact entries differ from their indexed parent cells");
      const auto& raw_multiplier = as_object(
          raw_entry.at("multiplier"), "prepared SCC multiplier");
      auto multiplier = parse_prepared_rational_taylor_multiplier<Scalar>(
          raw_multiplier, frame_width,
          static_cast<std::size_t>(work.work_t_order) + 1, true,
          "prepared SCC multiplier");
      const auto shifted_work_min = static_cast<std::int64_t>(work.work_min) +
          multiplier.epsilon_shift;
      const auto shifted_work_max =
          static_cast<std::int64_t>(work.work_complete_max) +
          multiplier.epsilon_shift;
      if (shifted_work_min < std::numeric_limits<std::int32_t>::min() ||
          shifted_work_min > std::numeric_limits<std::int32_t>::max() ||
          shifted_work_max < std::numeric_limits<std::int32_t>::min() ||
          shifted_work_max > std::numeric_limits<std::int32_t>::max())
        throw std::invalid_argument(
            "SCC multiplier shift overflows the retained epsilon frame");
      if (multiplier.exact_identity != exact_theta_entry)
        throw std::invalid_argument(
            "prepared multiplier identity differs from the exact theta entry");
      if (multiplier.proven_zero != parent_original.proven_zero)
        throw std::invalid_argument(
            "prepared multiplier structural-zero fact differs from its exact parent cell");
      if (!multiplier.proven_zero) {
        if (multiplier.kernels.size() != frame_width ||
            std::any_of(multiplier.kernels.begin(), multiplier.kernels.end(),
                [&](const auto& kernel) {
                  return kernel.size() !=
                      static_cast<std::size_t>(work.work_t_order) + 1;
                }))
          throw std::invalid_argument(
              "active SCC multiplier kernels do not cover the exact work rectangle");
        if (!expected_cross.contains(edge))
          throw std::invalid_argument(
              "active SCC multiplier lacks an exact parent structural edge");
        active_edges.insert(edge);
      } else if (expected_cross.contains(edge)) {
        throw std::invalid_argument(
            "a proven-zero SCC multiplier contradicts a parent structural edge");
      }
      coupling.identities.push_back(CompositeCouplingIdentity{
          source_vertex, target_vertex, exact_original_entry,
          exact_theta_entry, multiplier.proven_zero});
      coupling.matrix.entries.push_back(
          typename PreparedSparseLocalMultiplierMatrix<Scalar>::Entry{
              row, column, std::move(multiplier)});
    }
    couplings.push_back(std::move(coupling));
  }
  if (supplied_condensation != expected_condensation)
    throw std::invalid_argument(
        "SCC coupling groups do not cover the condensation graph exactly");
  if (active_edges != expected_cross)
    throw std::invalid_argument(
        "active SCC coupling entries do not cover cross-component structural edges exactly");

  return make_retained_typed_shared<Scalar, CompositeSCCChart<Scalar>>(
      handle, key, exact_identity, std::move(signature),
      root.if_contains("rational_shadow_identity") != nullptr
          ? required_string(root, "rational_shadow_identity")
          : exact_identity,
      dimension,
      std::move(graph), exact_system.canonical_record,
      exact_theta.canonical_record,
      geometry_record, std::move(retained_geometry), work,
      std::move(blocks), std::move(couplings),
      std::move(physical_equation));
}

struct SessionRegistry {
  std::mutex mutex;
  std::uint64_t next_session = 1;
  std::unordered_map<std::string, std::shared_ptr<SolverSession>> sessions;
};

SessionRegistry& session_registry() {
  // LibraryLink's uninitialize hook clears this registry.  Deliberately leak
  // the trivial holder itself so static destruction can never run after the
  // FLINT context used by retained SymbolicRational values.
  static auto* registry = new SessionRegistry;
  return *registry;
}

std::shared_ptr<SolverSession> find_session(const std::string& handle) {
  auto& registry = session_registry();
  std::lock_guard<std::mutex> lock(registry.mutex);
  const auto found = registry.sessions.find(handle);
  if (found == registry.sessions.end())
    throw std::invalid_argument("unknown or closed persistent solver session");
  return found->second;
}

bool acb_contains_exact_rational(const ComplexBall& ball,
                                 const Rational& value) {
  fmpq_t exact;
  fmpq_init(exact);
  const auto text = value.str();
  const int parsed = fmpq_set_str(exact, text.c_str(), 10);
  if (parsed == 0) fmpq_canonicalise(exact);
  const bool contains = parsed == 0 && acb_contains_fmpq(ball.raw(), exact);
  fmpq_clear(exact);
  return contains;
}

bool acb_contains_exact_epsilon_rational(
    const ExactEpsilonRational<ComplexBall>& numeric,
    const ExactEpsilonRational<Rational>& exact) {
  if (numeric.zero != exact.zero || numeric.valuation != exact.valuation ||
      numeric.numerator.size() != exact.numerator.size() ||
      numeric.denominator.size() != exact.denominator.size())
    return false;
  for (std::size_t index = 0; index < exact.numerator.size(); ++index)
    if (!acb_contains_exact_rational(
            numeric.numerator[index], exact.numerator[index]))
      return false;
  for (std::size_t index = 0; index < exact.denominator.size(); ++index)
    if (!acb_contains_exact_rational(
            numeric.denominator[index], exact.denominator[index]))
      return false;
  return true;
}

bool acb_physical_ode_contains_rational_shadow(
    const PreparedPhysicalClearedODE<ComplexBall>& numeric,
    const PreparedPhysicalClearedODE<Rational>& exact) {
  if (numeric.dimension != exact.dimension ||
      numeric.q_lags.size() != exact.q_lags.size() ||
      numeric.c_lags.size() != exact.c_lags.size())
    return false;
  for (std::size_t lag = 0; lag < exact.q_lags.size(); ++lag)
    if (!acb_contains_exact_epsilon_rational(
            numeric.q_lags[lag], exact.q_lags[lag]))
      return false;
  for (std::size_t lag = 0; lag < exact.c_lags.size(); ++lag) {
    if (numeric.c_lags[lag].size() != exact.c_lags[lag].size()) return false;
    for (std::size_t entry = 0; entry < exact.c_lags[lag].size(); ++entry) {
      const auto& n = numeric.c_lags[lag][entry];
      const auto& e = exact.c_lags[lag][entry];
      if (n.row != e.row || n.column != e.column ||
          !acb_contains_exact_epsilon_rational(n.value, e.value))
        return false;
    }
  }
  return true;
}

LocalSolution<ComplexBall> specialize_rational_local_to_acb(
    const LocalSolution<Rational>& source,
    const RetainedCompositeGeometry& target_geometry,
    const std::string& checkpoint_identity) {
  validate_local_solution(source, false);
  if (!source.error.empty())
    throw std::invalid_argument(
        "Rational-shadow specialization rejects a local error envelope");
  if (source.chart.center_exact != target_geometry.chart.center_exact ||
      source.chart.scale_exact != target_geometry.chart.scale_exact ||
      source.chart.infinite_radius !=
          target_geometry.chart.infinite_radius ||
      !local_algebra_detail::same_prescriptions(
          source.prescriptions, target_geometry.prescriptions))
    throw std::invalid_argument(
        "Rational-shadow local geometry or analytic prescriptions differ from the target Acb SCC");
  if (!target_geometry.chart.infinite_radius &&
      !acb_contains_exact_rational(
          source.chart.radius, Rational(target_geometry.radius_exact)))
    throw std::invalid_argument(
        "Rational-shadow local radius does not contain the target exact radius");

  LocalSolution<ComplexBall> output;
  output.chart = target_geometry.chart;
  output.epsilon = source.epsilon;
  output.taylor_complete_max = source.taylor_complete_max;
  output.dimension = source.dimension;
  output.prescriptions = target_geometry.prescriptions;
  output.checkpoint_identity = checkpoint_identity;
  output.sectors.reserve(source.sectors.size());
  for (const auto& sector : source.sectors) {
    if (sector.a.domain != ExactDomain::Rational ||
        sector.b.domain != ExactDomain::Rational)
      throw std::invalid_argument(
          "Rational-shadow specialization requires exact Rational sector tags");
    LocalSector<ComplexBall> converted;
    converted.a = ExactScalarDescriptor::rational(sector.a.canonical);
    converted.b = ExactScalarDescriptor::rational(sector.b.canonical);
    converted.log_power = sector.log_power;
    converted.coefficients.reserve(sector.coefficients.size());
    for (const auto& coefficient : sector.coefficients)
      converted.coefficients.push_back(
          ComplexBall::from_strings(coefficient.str()));
    output.sectors.push_back(std::move(converted));
  }
  validate_local_solution(output, false);
  return output;
}

json::object solve_prepared_chart_safe(
    const std::shared_ptr<PreparedChartBase>& chart,
    const json::value& raw_run, int output_digits,
    const std::string& session_handle) {
  try {
    auto result = chart->solve(
        as_object(raw_run, "persistent recurrence run"), output_digits);
    result["session"] = session_handle;
    result["chart"] = chart->handle();
    return result;
  } catch (const RecurrenceError& error) {
    return json::object{{"status", "error"}, {"id", error.id},
                        {"detail", error.what()},
                        {"frame_base", error.frame_base},
                        {"shift", error.shift},
                        {"session", session_handle},
                        {"chart", chart->handle()}};
  } catch (const std::exception& error) {
    return json::object{{"status", "error"}, {"id", "CPP"},
                        {"detail", error.what()},
                        {"session", session_handle},
                        {"chart", chart->handle()}};
  }
}
