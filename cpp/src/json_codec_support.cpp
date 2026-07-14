#include "json_codec_support.hpp"

#include <algorithm>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace diffexp2::json_codec_detail {
namespace json = boost::json;

const json::object& as_object(const json::value& value, const char* label) {
  if (!value.is_object())
    throw std::invalid_argument(std::string(label) + " must be an object");
  return value.as_object();
}

const json::array& as_array(const json::value& value, const char* label) {
  if (!value.is_array())
    throw std::invalid_argument(std::string(label) + " must be an array");
  return value.as_array();
}

std::int64_t as_i64(const json::value& value, const char* label) {
  if (value.is_int64()) return value.as_int64();
  if (value.is_uint64() && value.as_uint64() <=
      static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))
    return static_cast<std::int64_t>(value.as_uint64());
  throw std::invalid_argument(std::string(label) + " must be an integer");
}

std::uint32_t as_u32(const json::value& value, const char* label) {
  const auto number = as_i64(value, label);
  if (number < 0 || number > std::numeric_limits<std::uint32_t>::max())
    throw std::invalid_argument(std::string(label) +
                                " is outside uint32 range");
  return static_cast<std::uint32_t>(number);
}

std::uint64_t as_u64(const json::value& value, const char* label) {
  if (value.is_uint64()) return value.as_uint64();
  if (value.is_int64() && value.as_int64() >= 0)
    return static_cast<std::uint64_t>(value.as_int64());
  throw std::invalid_argument(std::string(label) +
                              " must be a nonnegative integer");
}

double as_double(const json::value& value, const char* label) {
  if (value.is_double()) return value.as_double();
  if (value.is_int64()) return static_cast<double>(value.as_int64());
  if (value.is_uint64()) return static_cast<double>(value.as_uint64());
  throw std::invalid_argument(std::string(label) + " must be numeric");
}

std::int32_t as_i32(const json::value& value, const char* label) {
  const auto number = as_i64(value, label);
  if (number < std::numeric_limits<std::int32_t>::min() ||
      number > std::numeric_limits<std::int32_t>::max())
    throw std::invalid_argument(std::string(label) +
                                " is outside int32 range");
  return static_cast<std::int32_t>(number);
}

std::string required_string(const json::object& object, const char* key) {
  const auto& value = object.at(key);
  if (!value.is_string() || value.as_string().empty())
    throw std::invalid_argument(std::string(key) +
                                " must be a nonempty string");
  return std::string(value.as_string());
}

void require_exact_keys(
    const json::object& object,
    std::initializer_list<std::string_view> expected,
    const char* label) {
  std::set<std::string> actual;
  for (const auto& entry : object) actual.emplace(entry.key());
  std::set<std::string> wanted;
  for (const auto key : expected) wanted.emplace(key);
  if (actual != wanted)
    throw std::invalid_argument(std::string(label) +
                                " has unknown or missing fields");
}

json::value canonical_json_value(const json::value& value) {
  if (value.is_array()) {
    json::array output;
    output.reserve(value.as_array().size());
    for (const auto& item : value.as_array())
      output.push_back(canonical_json_value(item));
    return output;
  }
  if (value.is_object()) {
    std::vector<std::string> keys;
    keys.reserve(value.as_object().size());
    for (const auto& item : value.as_object())
      keys.emplace_back(item.key());
    std::sort(keys.begin(), keys.end());
    json::object output;
    for (const auto& key : keys)
      output[key] = canonical_json_value(value.as_object().at(key));
    return output;
  }
  return value;
}

}  // namespace diffexp2::json_codec_detail
