#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <cassert>
#include <string>

int main() {
  namespace json = boost::json;

  const auto parsed = json::parse(diffexp2::backend_info_json());
  assert(parsed.is_object());
  const auto& info = parsed.as_object();

  assert(info.at("schema").as_int64() == 1);
  const auto& schemas = info.at("schemas").as_array();
  assert(schemas.size() == 2);
  assert(schemas[0].as_int64() == 1);
  assert(schemas[1].as_int64() == 2);
  assert(info.at("persistent_sessions").as_bool());
  assert(info.at("persistent_checkpoint_schema").as_int64() == 8);
  assert(info.at("persistent_exact_tile_plan_capability").as_string() ==
         "retained-exact-independent-arm-tile-plan-v1");
  assert(info.at("persistent_transport_pair_tile_stream_capability")
             .as_string() ==
         "retained-native-transport-pair-tile-stream-v1");
  assert(info.at("backend").as_string() == "DiffExp2 C++");
  assert(!info.at("flint").as_string().empty());
  assert(info.at("librarylink").as_bool());
}
