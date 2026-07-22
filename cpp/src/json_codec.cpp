#include "diffexp2/json_codec.hpp"

#include "json_codec_capabilities.hpp"
#include "json_codec_recurrence.hpp"
#include "json_codec_scc.hpp"
#include "json_codec_support.hpp"

#include "diffexp2/checkpoint.hpp"
#include "diffexp2/immutable_recursive_cache.hpp"
#include "diffexp2/line_integration.hpp"
#include "diffexp2/local_algebra.hpp"
#include "diffexp2/local_solution.hpp"
#include "diffexp2/matching.hpp"
#include "diffexp2/path_planner.hpp"
#include "diffexp2/physical_ode.hpp"
#include "diffexp2/recurrence.hpp"
#include "diffexp2/scc_completeness.hpp"
#include "diffexp2/singular_indicial.hpp"
#include "diffexp2/tail_majorant.hpp"

#include <boost/json.hpp>

#include <array>
#include <chrono>
#include <atomic>
#include <algorithm>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>

namespace diffexp2 {
namespace json = boost::json;
namespace {

using namespace json_codec_detail;

#include "json_codec_local_state.ipp"
#include "json_codec_multiplier_state.ipp"
#include "json_codec_endpoint_state.ipp"
#include "json_codec_matching_state.ipp"
#include "json_codec_chart_scc_state.ipp"
#include "json_codec_transport_state.ipp"
#include "json_codec_session_state.ipp"
#include "json_codec_checkpoint_state.ipp"
#include "json_codec_commands_session_tile.ipp"
#include "json_codec_commands_transport.ipp"
#include "json_codec_commands_integration.ipp"
#include "json_codec_commands_solve.ipp"
json::object run_one(const json::object& root) {
  if (as_i64(root.at("schema"), "schema") != 1)
    throw std::invalid_argument("unsupported recurrence schema");
  const auto domain = std::string(root.at("domain").as_string());
  const int digits = root.if_contains("output_digits")
                         ? static_cast<int>(as_i64(root.at("output_digits"), "output digits"))
                         : 50;
  if (domain == "rational") return run_typed<Rational>(root, digits);
  if (domain == "acb") {
    const auto precision_bits = static_cast<slong>(
        as_i64(root.at("precision_bits"), "precision bits"));
    AcbPrecisionLease lease(precision_bits);
    ComplexBall::set_precision(precision_bits);
    return run_typed<ComplexBall>(root, digits);
  }
  if (domain == "symbolic") {
    std::lock_guard<std::recursive_mutex> lock(symbolic_run_mutex());
    std::vector<std::string> variables;
    for (const auto& value : as_array(root.at("symbols"), "symbolic variables")) {
      if (!value.is_string())
        throw std::invalid_argument("symbolic variable names must be strings");
      variables.emplace_back(value.as_string());
    }
    SymbolicRational::configure(variables);
    return run_typed<SymbolicRational>(root, digits);
  }
  throw std::invalid_argument("unsupported scalar domain: " + domain);
}

json::object run_one_safe(const json::object& root) {
  try {
    return run_one(root);
  } catch (const RecurrenceError& error) {
    return json::object{{"status", "error"}, {"id", error.id},
      {"detail", error.what()}, {"frame_base", error.frame_base},
      {"shift", error.shift}};
  } catch (const std::exception& error) {
    return json::object{{"status", "error"}, {"id", "CPP"},
                        {"detail", error.what()}};
  }
}

}  // namespace

std::string run_recurrence_json(std::string_view input) {
  try {
    const auto parsed = json::parse(input);
    const auto& root = as_object(parsed, "root");
    if (const auto* schema = root.if_contains("schema");
        schema != nullptr && as_i64(*schema, "schema") == 2)
      return json::serialize(run_session_command(root));
    if (const auto* raw_batch = root.if_contains("batch"); raw_batch != nullptr) {
      const auto& requests = as_array(*raw_batch, "batch requests");
      const auto requested_threads = root.if_contains("threads")
          ? as_u32(root.at("threads"), "batch threads") : 1;
      auto thread_count = std::max<std::uint32_t>(1,
          std::min<std::uint32_t>(requested_threads,
                                 static_cast<std::uint32_t>(requests.size())));
      // Acb precision is process-global in the light RAII wrapper. All batch
      // requests must therefore agree before workers start.
      std::optional<std::int64_t> acb_precision;
      std::optional<std::vector<std::string>> symbolic_variables;
      for (const auto& raw_request : requests) {
        const auto& request = as_object(raw_request, "batch request");
        const auto domain = std::string(request.at("domain").as_string());
        if (domain == "acb") {
          const auto bits = as_i64(request.at("precision_bits"), "precision bits");
          if (acb_precision.has_value() && *acb_precision != bits)
            throw std::invalid_argument("all Acb batch requests must use one precision");
          acb_precision = bits;
        } else if (domain == "symbolic") {
          std::vector<std::string> variables;
          for (const auto& value : as_array(
                   request.at("symbols"), "symbolic variables"))
            variables.emplace_back(value.as_string());
          if (symbolic_variables.has_value() &&
              *symbolic_variables != variables)
            throw std::invalid_argument(
                "all symbolic batch requests must use one coefficient field");
          symbolic_variables = std::move(variables);
        }
      }
      if (acb_precision.has_value())
        ComplexBall::set_precision(static_cast<slong>(*acb_precision));
      if (symbolic_variables.has_value()) {
        SymbolicRational::configure(*symbolic_variables);
        // FLINT multivariate contexts are shared read-only, but keep this
        // first symbolic milestone single-threaded until its allocator-level
        // parallel behavior is separately certified.
        thread_count = 1;
      }
      std::vector<json::object> results(requests.size());
      std::atomic<std::size_t> next{0};
      auto worker = [&]() {
        while (true) {
          const auto index = next.fetch_add(1);
          if (index >= requests.size()) return;
          results[index] = run_one_safe(
              as_object(requests[index], "batch request"));
        }
      };
      std::vector<std::jthread> workers;
      workers.reserve(thread_count);
      for (std::uint32_t i = 0; i < thread_count; ++i)
        workers.emplace_back(worker);
      for (auto& thread : workers) thread.join();
      json::array encoded;
      encoded.reserve(results.size());
      for (auto& result : results) encoded.push_back(std::move(result));
      return json::serialize(json::object{{"status", "ok"},
                                          {"results", std::move(encoded)}});
    }
    return json::serialize(run_one_safe(root));
  } catch (const NativeIntegrationError& error) {
    json::object result{{"status", "error"}, {"id", error.id},
                        {"detail", error.what()}};
    if (!error.absolute_power.empty()) {
      result["absolute_power"] = error.absolute_power;
      result["log_power"] = error.log_power;
      result["epsilon_power"] = error.epsilon_power;
      result["component"] = error.component;
    }
    return json::serialize(result);
  } catch (const RecurrenceError& error) {
    return json::serialize(json::object{{"status", "error"}, {"id", error.id},
      {"detail", error.what()}, {"frame_base", error.frame_base}, {"shift", error.shift}});
  } catch (const std::exception& error) {
    return json::serialize(json::object{{"status", "error"}, {"id", "CPP"},
                                        {"detail", error.what()}});
  }
}

void reset_solver_sessions() {
  auto& registry = session_registry();
  std::unordered_map<std::string, std::shared_ptr<SolverSession>> sessions;
  {
    std::lock_guard<std::mutex> lock(registry.mutex);
    sessions.swap(registry.sessions);
  }
  // Destroy retained FLINT/Arb objects while the library and coefficient
  // contexts are still alive.  In-flight calls keep shared ownership until
  // they leave the backend.
  sessions.clear();
}

}  // namespace diffexp2
