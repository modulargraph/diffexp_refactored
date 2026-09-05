#include "diffexp/kernel/json_codec.hpp"

#include <boost/json.hpp>

#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>

namespace json = boost::json;

namespace {

json::object request(const json::object &payload) {
  return json::parse(diffexp::kernel::run_recurrence_json(json::serialize(payload)))
      .as_object();
}

json::object run() {
  return json::parse(R"json({
    "nmax":0,"p":0,"has_initial":true,"adaptive_probe":false,
    "a_target":"0","b_target":"0","a_shift_min":0,
    "a_shifts":["0"],
    "schedule":[[{"case":"R","da":"0","db":"0"}]],
    "initial":["1","0","0","0","0","0","0","0"],
    "initial_validity":[5],"source":null,"return_u":false
  })json")
      .as_object();
}

json::object metadata(const std::string &checkpoint) {
  auto result = json::parse(R"json({
    "chart":{"center_exact":"0","scale_exact":"1",
      "radius":"2","infinite_radius":false},
    "tag":{"a":{"domain":"rational","canonical":"0"},
      "b":{"domain":"rational","canonical":"0"}},
    "prescriptions":[],"checkpoint_identity":"placeholder"
  })json")
                    .as_object();
  result["checkpoint_identity"] = checkpoint;
  return result;
}

bool has_coefficient_slab(const json::object &value) {
  return value.if_contains("assembled") || value.if_contains("coefficients") ||
         value.if_contains("u") || value.if_contains("validity");
}

} // namespace

int main() {
  try {
    const auto created = request(json::object{{"schema", 2},
                                              {"op", "session.create"},
                                              {"domain", "rational"},
                                              {"output_digits", 30},
                                              {"local_capacity", 2}});
    const auto session = std::string(created.at("session").as_string());
    auto prepared_request = json::parse(R"json({
    "schema":2,"op":"chart.prepare","session":"placeholder",
    "key":"batch@0[-2,8]","identity":"de2-operator-persistent-local-batch-v1",
    "analytic":{"geometry":{"center_exact":"0","scale_exact":"1",
      "radius_exact":"2","infinite_radius":false,"prescriptions":[]}},
    "scc":{"components":[[0]],"structural_edges":[],
      "condensation_edges":[],"topological_order":[0],"coupling_depth":0},
    "problem":{"domain":"rational","d":1,"fb":-2,"w":8,
      "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
      "nhat_lags":[{"poly":[],"rat":[],"val":[null]}],
      "d0_inverse":"1","blocks":[[0]],
      "assembly":{"identity":true,"poly":[],"rat":[],"val":[0]},
      "physical_ode":{"schema":"diffexp3-physical-cleared-ode-v1",
        "basis":"physical-original-master","theta_coordinate":"local-t",
        "owner_signature_identity":"de2-operator-persistent-local-batch-v1",
        "payload_identity":"de2-physical-ode-persistent-local-batch-v1",
        "q":[{"zero":false,"valuation":0,"numerator":["1"],
          "denominator":["1"]}],"c":[[]]},
      "chop_digits":10}
  })json")
                                .as_object();
    prepared_request["session"] = session;
    const auto prepared = request(prepared_request);
    if (prepared.at("status") != "ok") {
      std::cerr << "chart prepare: " << json::serialize(prepared) << '\n';
      return EXIT_FAILURE;
    }
    const auto chart = std::string(prepared.at("chart").as_string());
    const auto equation_owner = request(json::parse(std::string(R"json({
    "schema":2,"op":"regular_equation.prepare","session":")json") +
                                                    session + R"json(",
    "capability":"frame-independent-regular-physical-equation-owner-v1",
    "key":"regular-equation:batch-fixture",
    "identity":"de2-equation-batch-fixture","dimension":1,
    "relative_accuracy_max_exact":"1/1000",
    "geometry":{"center_exact":"0","scale_exact":"1",
      "radius_exact":"2","infinite_radius":false,"prescriptions":[]},
    "physical_ode":{"schema":"diffexp3-physical-cleared-ode-v1",
      "basis":"physical-original-master","theta_coordinate":"local-t",
      "owner_signature_identity":"de2-equation-batch-fixture",
      "payload_identity":"de2-physical-ode-persistent-local-batch-v1",
      "q":[{"zero":false,"valuation":0,"numerator":["1"],
        "denominator":["1"]}],"c":[[]]}
  })json")
                                            .as_object());
    if (equation_owner.at("status") != "ok") {
      std::cerr << "equation owner: " << json::serialize(equation_owner)
                << '\n';
      return EXIT_FAILURE;
    }
    const auto equation_handle =
        std::string(equation_owner.at("equation_owner").as_string());

    json::array runs;
    runs.push_back(run());
    runs.push_back(run());
    json::array metadata_records;
    metadata_records.push_back(metadata("batch-column-1"));
    metadata_records.push_back(metadata("batch-column-2"));
    const auto batch =
        request(json::object{{"schema", 2},
                             {"op", "local.solve_batch"},
                             {"session", session},
                             {"chart", chart},
                             {"threads", 2},
                             {"runs", runs},
                             {"metadata", metadata_records},
                             {"equation_owner", equation_handle}});
    const auto &results = batch.at("results").as_array();
    const auto &first = results[0].as_object();
    const auto &second = results[1].as_object();
    const auto stats_after_success = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});

    for (const auto &result : results) {
      const auto &local = result.as_object().at("local");
      (void)request(json::object{{"schema", 2},
                                 {"op", "local.release"},
                                 {"session", session},
                                 {"local", local}});
    }

    json::array bad_metadata;
    bad_metadata.push_back(metadata("discarded-success"));
    bad_metadata.push_back(json::object{});
    const auto failed = request(json::object{{"schema", 2},
                                             {"op", "local.solve_batch"},
                                             {"session", session},
                                             {"chart", chart},
                                             {"threads", 2},
                                             {"runs", runs},
                                             {"metadata", bad_metadata}});
    const auto stats_after_failure = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});

    const bool ok =
        created.at("status") == "ok" && prepared.at("status") == "ok" &&
        batch.at("status") == "ok" && batch.at("attempted") == 2 &&
        batch.at("succeeded") == 2 && batch.at("failed") == 0 &&
        batch.at("requested_threads") == 2 && batch.at("worker_threads") == 2 &&
        batch.at("atomic_retention") == true &&
        batch.at("json_coefficients") == 0 && results.size() == 2 &&
        equation_owner.at("status") == "ok" &&
        std::string(first.at("chart").as_string()) == equation_handle &&
        std::string(second.at("chart").as_string()) == equation_handle &&
        first.at("source_operator_identity") == "de2-equation-batch-fixture" &&
        second.at("source_operator_identity") == "de2-equation-batch-fixture" &&
        first.at("checkpoint_identity") == "batch-column-1" &&
        second.at("checkpoint_identity") == "batch-column-2" &&
        first.at("native_retained") == true &&
        second.at("native_retained") == true && !has_coefficient_slab(batch) &&
        !has_coefficient_slab(first) && !has_coefficient_slab(second) &&
        stats_after_success.at("locals") == 2 &&
        stats_after_success.at("local_solves") == 2 &&
        stats_after_success.at("pending_local_solves") == 0 &&
        failed.at("status") == "error" &&
        stats_after_failure.at("locals") == 0 &&
        stats_after_failure.at("local_solves") == 2 &&
        stats_after_failure.at("pending_local_solves") == 0;
    if (!ok) {
      std::cerr << "batch: " << json::serialize(batch) << '\n'
                << "stats after success: "
                << json::serialize(stats_after_success) << '\n'
                << "failed batch: " << json::serialize(failed) << '\n'
                << "stats after failure: "
                << json::serialize(stats_after_failure) << '\n';
    }

    (void)request(json::object{
        {"schema", 2}, {"op", "session.close"}, {"session", session}});
    std::cout << (ok ? "PASS" : "FAIL")
              << ": atomic retained local.solve_batch\n";
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception &error) {
    std::cerr << "persistent local batch fixture threw: " << error.what()
              << '\n';
    return EXIT_FAILURE;
  }
}
