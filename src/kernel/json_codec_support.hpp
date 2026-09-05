#pragma once

#include <boost/json.hpp>

#include <cstdint>
#include <initializer_list>
#include <string>
#include <string_view>

namespace diffexp::kernel::json_codec_detail {

const boost::json::object& as_object(const boost::json::value& value,
                                     const char* label);
const boost::json::array& as_array(const boost::json::value& value,
                                   const char* label);
std::int64_t as_i64(const boost::json::value& value, const char* label);
std::uint32_t as_u32(const boost::json::value& value, const char* label);
std::uint64_t as_u64(const boost::json::value& value, const char* label);
double as_double(const boost::json::value& value, const char* label);
std::int32_t as_i32(const boost::json::value& value, const char* label);
std::string required_string(const boost::json::object& object,
                            const char* key);
void require_exact_keys(
    const boost::json::object& object,
    std::initializer_list<std::string_view> expected,
    const char* label);
boost::json::value canonical_json_value(const boost::json::value& value);

}  // namespace diffexp::kernel::json_codec_detail
