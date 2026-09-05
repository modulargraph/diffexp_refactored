#pragma once

#include <string>
#include <string_view>

namespace diffexp::kernel {

std::string run_recurrence_json(std::string_view input);
std::string backend_info_json();
void reset_solver_sessions();

}  // namespace diffexp::kernel
