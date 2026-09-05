#include "diffexp/kernel/generated_protocol_contract.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

namespace {

void require(const bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

}  // namespace

int main() {
  using diffexp::kernel::protocol::ResponseKind;
  using diffexp::kernel::protocol::enforce_response_budget;
  using diffexp::kernel::protocol::operation_contract;

  const auto* hot = operation_contract("local.specialize_rational_shadow");
  require(hot != nullptr, "hot operation missing from generated contract");
  require(hot->kind == ResponseKind::Bounded,
          "hot operation must have a bounded response");
  enforce_response_budget(hot->operation, hot->max_response_bytes);
  bool rejected = false;
  try {
    enforce_response_budget(hot->operation, hot->max_response_bytes + 1U);
  } catch (const std::length_error&) {
    rejected = true;
  }
  require(rejected, "oversized hot response was not rejected");

  const auto* diagnostic = operation_contract("session.stats");
  require(diagnostic != nullptr && diagnostic->paged,
          "session diagnostics must be explicitly paged");
  require(diffexp::kernel::protocol::kDiagnosticPageDefault <=
              diffexp::kernel::protocol::kDiagnosticPageMax,
          "diagnostic page defaults are invalid");

  const auto* bulk = operation_contract("endpoint.export");
  require(bulk != nullptr && bulk->kind == ResponseKind::ExplicitBulk,
          "large endpoint transfer must be explicit bulk output");
  enforce_response_budget(bulk->operation, 1024U * 1024U * 1024U);

  bool unknown_rejected = false;
  try {
    enforce_response_budget("undeclared.operation", 1U);
  } catch (const std::logic_error&) {
    unknown_rejected = true;
  }
  require(unknown_rejected, "undeclared operation was not rejected");

  std::cout << "Generated protocol contract tests passed.\n";
  return EXIT_SUCCESS;
}
