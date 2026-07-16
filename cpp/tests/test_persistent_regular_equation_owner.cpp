#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#include <unistd.h>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp2::run_recurrence_json(json::serialize(value)))
      .as_object();
}

void require_ok(const json::object& response, const char* label) {
  if (response.if_contains("status") == nullptr ||
      response.at("status") != "ok")
    throw std::runtime_error(std::string(label) + ": " +
                             json::serialize(response));
}

bool is_error(const json::object& response) {
  return response.if_contains("status") != nullptr &&
         response.at("status") == "error";
}

json::object physical_ode(const std::string& owner,
                          const std::string& payload) {
  const auto one = json::object{
      {"zero", false}, {"valuation", 0},
      {"numerator", json::array{"1"}},
      {"denominator", json::array{"1"}}};
  json::array c;
  c.emplace_back(json::array{});
  return json::object{
      {"schema", "diffexp2-physical-cleared-ode-v1"},
      {"basis", "physical-original-master"},
      {"theta_coordinate", "local-t"},
      {"owner_signature_identity", owner},
      {"payload_identity", payload},
      {"q", json::array{one}}, {"c", std::move(c)}};
}

json::object geometry(const std::string& center = "0") {
  return json::object{
      {"center_exact", center}, {"scale_exact", "1"},
      {"radius_exact", "2"}, {"infinite_radius", false},
      {"prescriptions", json::array{}}};
}

json::object prepare_request(const std::string& session,
                             const std::string& key,
                             const std::string& identity,
                             const std::string& center = "0") {
  return json::object{
      {"schema", 2}, {"op", "regular_equation.prepare"},
      {"session", session},
      {"capability",
       "frame-independent-regular-physical-equation-owner-v1"},
      {"key", key}, {"identity", identity}, {"dimension", 1},
      {"geometry", geometry(center)},
      {"physical_ode",
       physical_ode(identity, "de2-physical-ode-owner-fixture")}};
}

}  // namespace

int main() {
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"},
        {"domain", "rational"}, {"output_digits", 20},
        {"chart_capacity", 4}});
    require_ok(created, "session.create");
    if (created.at("regular_equation_owner_capability") !=
        "frame-independent-regular-physical-equation-owner-v1")
      throw std::runtime_error(
          "session did not advertise the regular equation-owner capability");
    const auto session = std::string(created.at("session").as_string());
    const std::string identity = "de2-equation-owner-fixture";

    const auto first = request(prepare_request(
        session, "regular-equation:fixture", identity));
    require_ok(first, "regular_equation.prepare");
    const auto handle =
        std::string(first.at("equation_owner").as_string());
    if (!handle.starts_with("eq:") || first.at("reused") != false ||
        std::string(first.at("identity").as_string()) != identity ||
        std::string(first.at("owner_signature_identity").as_string()) !=
            identity ||
        first.at("physical_payload_identity") !=
            "de2-physical-ode-owner-fixture")
      throw std::runtime_error(
          "regular equation-owner preparation lost exact provenance");

    const auto reused = request(prepare_request(
        session, "regular-equation:fixture", identity));
    require_ok(reused, "regular_equation.prepare reuse");
    if (reused.at("reused") != true ||
        std::string(reused.at("equation_owner").as_string()) != handle)
      throw std::runtime_error(
          "identical regular equation-owner preparation was not reused");

    const auto collision = request(prepare_request(
        session, "regular-equation:fixture", identity, "1"));
    if (!is_error(collision))
      throw std::runtime_error(
          "unequal geometry reused one regular equation-owner cache key");

    auto nonregular = prepare_request(
        session, "regular-equation:nonregular",
        "de2-equation-owner-nonregular");
    const auto one = json::object{
        {"zero", false}, {"valuation", 0},
        {"numerator", json::array{"1"}},
        {"denominator", json::array{"1"}}};
    json::array c0;
    c0.emplace_back(json::object{
        {"r", 0}, {"c", 0}, {"v", one}});
    json::array c_lags;
    c_lags.emplace_back(std::move(c0));
    nonregular.at("physical_ode").as_object()["c"] =
        std::move(c_lags);
    const auto nonregular_rejected = request(std::move(nonregular));
    if (!is_error(nonregular_rejected) ||
        std::string(nonregular_rejected.at("detail").as_string())
                .find("requires an ordinary center") ==
            std::string::npos)
      throw std::runtime_error(
          "nonregular physical q/C payload acquired a regular equation owner");

    const auto counters = request(json::object{
        {"schema", 2}, {"op", "session.counters"},
        {"session", session}});
    require_ok(counters, "session.counters");
    if (counters.at("regular_equation_owners") != 1)
      throw std::runtime_error(
          "session counters lost the retained regular equation owner");
    const auto stats = request(json::object{
        {"schema", 2}, {"op", "session.stats"},
        {"session", session}});
    require_ok(stats, "session.stats");
    if (stats.at("regular_equation_owners") != 1)
      throw std::runtime_error(
          "session stats lost the retained regular equation owner");

    const auto path = (std::filesystem::temp_directory_path() /
        ("diffexp2-regular-equation-owner-" +
         std::to_string(static_cast<long long>(::getpid())) + ".json"))
        .string();
    std::filesystem::remove(path);
    const auto unsupported_checkpoint = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"},
        {"session", session}, {"path", path},
        {"checkpoint_identity", "regular-equation-owner-live-v1"}});
    if (!is_error(unsupported_checkpoint) ||
        std::string(unsupported_checkpoint.at("detail").as_string())
                .find("does not yet support frame-independent") ==
            std::string::npos ||
        std::filesystem::exists(path))
      throw std::runtime_error(
          "checkpoint did not reject the unsupported live owner explicitly");

    const auto released = request(json::object{
        {"schema", 2}, {"op", "regular_equation.release"},
        {"session", session}, {"equation_owner", handle}});
    require_ok(released, "regular_equation.release");
    const auto released_twice = request(json::object{
        {"schema", 2}, {"op", "regular_equation.release"},
        {"session", session}, {"equation_owner", handle}});
    if (!is_error(released_twice))
      throw std::runtime_error(
          "a released regular equation owner remained publicly visible");

    const auto after_release = request(json::object{
        {"schema", 2}, {"op", "session.counters"},
        {"session", session}});
    require_ok(after_release, "session.counters after release");
    if (after_release.at("regular_equation_owners") != 0)
      throw std::runtime_error(
          "released regular equation owner remained in session counters");

    const auto checkpoint = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"},
        {"session", session}, {"path", path},
        {"checkpoint_identity", "regular-equation-owner-released-v1"}});
    require_ok(checkpoint, "checkpoint.save after owner release");
    std::filesystem::remove(path);

    const auto second = request(prepare_request(
        session, "regular-equation:second",
        "de2-equation-owner-second"));
    require_ok(second, "second regular_equation.prepare");
    const auto closed = request(json::object{
        {"schema", 2}, {"op", "session.close"},
        {"session", session}});
    require_ok(closed, "session.close");
    if (closed.at("released_regular_equation_owners") != 1)
      throw std::runtime_error(
          "session close did not account for its regular equation owner");
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
