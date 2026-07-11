#pragma once

#include <string>
#include <string_view>

namespace diffexp2 {

std::string run_recurrence_json(std::string_view input);
std::string backend_info_json();

}  // namespace diffexp2
