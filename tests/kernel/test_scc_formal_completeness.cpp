#include "diffexp/kernel/scc_completeness.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

using diffexp::kernel::ComplexBall;
using diffexp::kernel::ExactEpsilonRational;
using diffexp::kernel::ExactScalarDescriptor;
using diffexp::kernel::JordanBlock;
using diffexp::kernel::LocalSector;
using diffexp::kernel::LocalSolution;
using diffexp::kernel::PhysicalODEMatrixEntry;
using diffexp::kernel::PreparedMatrix;
using diffexp::kernel::PreparedPhysicalClearedODE;
using diffexp::kernel::PreparedRecurrenceOperator;
using diffexp::kernel::Rational;
using diffexp::kernel::RecurrenceError;
using diffexp::kernel::RecurrenceProblem;
using diffexp::kernel::RecurrenceResult;

ExactEpsilonRational<Rational> epsilon_rational(
    std::int32_t valuation, const char* numerator) {
  ExactEpsilonRational<Rational> result;
  result.zero = false;
  result.valuation = valuation;
  result.numerator = {Rational(numerator)};
  result.denominator = {Rational(1)};
  return result;
}

PreparedPhysicalClearedODE<Rational> exponential_equation(
    std::int32_t c_valuation = 1) {
  PreparedPhysicalClearedODE<Rational> equation;
  equation.dimension = 1;
  equation.q_lags = {epsilon_rational(0, "1")};
  equation.c_lags.resize(1);
  equation.c_lags.front().push_back(
      PhysicalODEMatrixEntry<Rational>{
          0, 0, epsilon_rational(c_valuation, "1")});
  equation.owner_signature_identity = "scc-formal-owner-v1";
  equation.payload_identity = "scc-formal-payload-v1";
  equation.exact_payload_record = "scc-formal-payload-record-v1";
  return equation;
}

LocalSolution<Rational> exponential_local() {
  LocalSolution<Rational> local;
  local.chart.center_exact = "0";
  local.chart.scale_exact = "1";
  local.chart.radius = ComplexBall::from_strings("2");
  local.epsilon = {0, 2};
  local.taylor_complete_max = 1;
  local.dimension = 1;
  local.checkpoint_identity = "scc-formal-local-v1";
  LocalSector<Rational> sector;
  sector.a = ExactScalarDescriptor::rational("0");
  sector.b = ExactScalarDescriptor::rational("1");
  sector.log_power = 0;
  sector.coefficients.assign(local.sector_size(), Rational(0));
  // F=t^eps.  The exp(eps log(t)) factor is represented by b=1 rather
  // than expanded into a finite log polynomial.
  sector.coefficients[0] = Rational(1);
  local.sectors.push_back(std::move(sector));
  return local;
}

void test_integer_shifted_tower_alignment() {
  PreparedPhysicalClearedODE<Rational> equation;
  equation.dimension = 2;
  equation.q_lags = {epsilon_rational(0, "1")};
  equation.c_lags.resize(2);
  equation.c_lags[1].push_back(
      PhysicalODEMatrixEntry<Rational>{
          0, 1, epsilon_rational(0, "1")});
  equation.owner_signature_identity = "shifted-tower-owner-v1";
  equation.payload_identity = "shifted-tower-payload-v1";
  equation.exact_payload_record = "shifted-tower-record-v1";

  LocalSolution<Rational> local;
  local.chart.center_exact = "0";
  local.chart.scale_exact = "1";
  local.chart.radius = ComplexBall::from_strings("2");
  local.epsilon = {0, 0};
  local.taylor_complete_max = 2;
  local.dimension = 2;
  local.checkpoint_identity = "shifted-tower-local-v1";

  // theta(y0)=t y1, theta(y1)=0 with y0=t and y1=1.  The gauge path
  // represents y0 in the a=-1 tower and y1 in the integer-equivalent a=0
  // tower.  They are coefficients of the same exact formal series after
  // aligning a+n, even though their retained sector tags differ.
  LocalSector<Rational> lower;
  lower.a = ExactScalarDescriptor::rational("-1");
  lower.b = ExactScalarDescriptor::rational("0");
  lower.log_power = 0;
  lower.coefficients.assign(local.sector_size(), Rational(0));
  lower.coefficients[4] = Rational(1);  // n=2, component=0: t^(-1+2)

  LocalSector<Rational> upper;
  upper.a = ExactScalarDescriptor::rational("0");
  upper.b = ExactScalarDescriptor::rational("0");
  upper.log_power = 0;
  upper.coefficients.assign(local.sector_size(), Rational(0));
  upper.coefficients[1] = Rational(1);  // n=0, component=1
  local.sectors = {std::move(lower), std::move(upper)};

  const auto certificate =
      diffexp::kernel::certify_scc_parent_exact_formal_residual(
          equation, local, {0, 0}, 2);
  if (certificate.exact_tag_count != 1 ||
      certificate.coefficient_rows != 6)
    throw std::runtime_error(
        "integer-shifted formal towers were not aligned to one base");

  local.sectors.front().coefficients[4] = Rational(2);
  bool corruption_rejected = false;
  try {
    (void)diffexp::kernel::certify_scc_parent_exact_formal_residual(
        equation, local, {0, 0}, 2);
  } catch (const std::domain_error& error) {
    corruption_rejected =
        std::string(error.what()).find("not exact zero") !=
        std::string::npos;
  }
  if (!corruption_rejected)
    throw std::runtime_error(
        "aligned integer-shifted towers hid a genuine residual");
}

void test_formal_parent_certificate() {
  const auto equation = exponential_equation();
  auto local = exponential_local();
  const auto certificate =
      diffexp::kernel::certify_scc_parent_exact_formal_residual(
          equation, local, {0, 2}, 1);
  if (certificate.epsilon.min_power != 0 ||
      certificate.epsilon.complete_max != 2 ||
      certificate.taylor_complete_max != 1 ||
      certificate.exact_tag_count != 1 ||
      certificate.coefficient_rows != 6)
    throw std::runtime_error(
        "formal parent certificate changed its claimed slab");

  // Add eps^2 t inside the same formal sector.  q theta - C then has the
  // exact nonzero coefficient eps^2 t, so no numerical tolerance can accept
  // it as completeness evidence.
  local.sectors.front().coefficients[5] = Rational(1);
  const auto prefix =
      diffexp::kernel::certify_scc_parent_exact_formal_residual(
          equation, local, {0, 2}, 1, true);
  if (prefix.epsilon.min_power != 0 ||
      prefix.epsilon.complete_max != 1 ||
      prefix.coefficient_rows != 4)
    throw std::runtime_error(
        "formal parent prefix certificate did not stop before the first exact residual");
  bool corruption_rejected = false;
  try {
    (void)diffexp::kernel::certify_scc_parent_exact_formal_residual(
        equation, local, {0, 2}, 1);
  } catch (const std::domain_error& error) {
    corruption_rejected =
        std::string(error.what()).find("not exact zero") !=
            std::string::npos &&
        std::string(error.what()).find("coefficient=1") !=
            std::string::npos;
  }
  if (!corruption_rejected)
    throw std::runtime_error(
        "formal parent certificate accepted a nonzero Taylor residual");

  // eps^-1 F at output eps^2 requires the unknown input eps^3.  Even a
  // stored all-zero slab cannot turn that missing reservoir into structural
  // zero.
  bool reservoir_rejected = false;
  try {
    (void)diffexp::kernel::certify_scc_parent_exact_formal_residual(
        exponential_equation(-1), exponential_local(),
        {0, 2}, 1);
  } catch (const std::domain_error& error) {
    reservoir_rejected =
        std::string(error.what()).find("above the retained epsilon reservoir") !=
        std::string::npos;
  }
  if (!reservoir_rejected)
    throw std::runtime_error(
        "formal parent certificate ignored an unknown upper epsilon input");
}

void test_transaction_local_candidate() {
  PreparedRecurrenceOperator<Rational> prepared;
  prepared.dimension = 1;
  prepared.frame_base = -1;
  prepared.frame_width = 3;
  prepared.blocks = {JordanBlock{{0}}};
  PreparedMatrix<Rational> assembly;
  assembly.identity = true;
  assembly.valuations = {0};
  prepared.assembly_matrix = std::move(assembly);

  RecurrenceProblem<Rational> problem;
  problem.dimension = 1;
  problem.nmax = 0;
  problem.log_max = 0;
  problem.frame_base = -1;
  problem.frame_width = 3;

  RecurrenceResult<Rational> recurrence;
  // [n=0][log=0,1][component=0][eps=-1,0,1].
  recurrence.u.assign(6, Rational(0));
  recurrence.u[1] = Rational(1);
  recurrence.validity = {-1, diffexp::kernel::kCompleteInfinity};
  recurrence.top_valid = -1;
  const auto original_validity = recurrence.validity;

  bool ordinary_rejected = false;
  try {
    (void)diffexp::kernel::assemble_recurrence(prepared, problem, recurrence);
  } catch (const RecurrenceError& error) {
    ordinary_rejected = error.id == "E6";
  }
  if (!ordinary_rejected)
    throw std::runtime_error(
        "ordinary recurrence assembly unexpectedly promoted completeness");

  auto candidate = diffexp::kernel::assemble_scc_recurrence_candidate(
      prepared, problem, recurrence);
  if (!candidate.requires_parent_certificate ||
      candidate.recurrence_top_valid != -1 ||
      candidate.coefficients.min_power != 0 ||
      candidate.coefficients.complete_max != 1 ||
      candidate.coefficients.coefficients.size() != 2 ||
      !(candidate.coefficients.coefficients.front() == Rational(1)) ||
      recurrence.validity != original_validity)
    throw std::runtime_error(
        "transaction-local SCC candidate lost its deferred evidence");

  PreparedRecurrenceOperator<ComplexBall> acb_prepared;
  acb_prepared.dimension = 1;
  acb_prepared.frame_base = -1;
  acb_prepared.frame_width = 3;
  acb_prepared.blocks = {JordanBlock{{0}}};
  PreparedMatrix<ComplexBall> acb_assembly;
  acb_assembly.identity = true;
  acb_assembly.valuations = {0};
  acb_prepared.assembly_matrix = std::move(acb_assembly);
  RecurrenceProblem<ComplexBall> acb_problem;
  acb_problem.dimension = 1;
  acb_problem.nmax = 0;
  acb_problem.log_max = 0;
  acb_problem.frame_base = -1;
  acb_problem.frame_width = 3;
  RecurrenceResult<ComplexBall> acb_recurrence;
  acb_recurrence.u.assign(6, ComplexBall(0));
  acb_recurrence.u[1] = ComplexBall(1);
  acb_recurrence.validity = original_validity;
  acb_recurrence.top_valid = -1;
  auto acb_candidate = diffexp::kernel::assemble_scc_recurrence_candidate(
      acb_prepared, acb_problem, acb_recurrence);
  if (!acb_candidate.requires_parent_certificate)
    throw std::runtime_error(
        "Acb assembly exhaustion did not produce a deferred candidate");
  bool shadow_requested = false;
  try {
    diffexp::kernel::require_exact_domain_for_deferred_scc_candidate(acb_candidate);
  } catch (const RecurrenceError& error) {
    shadow_requested = error.id == "E5" &&
        std::string(error.what()).find("requires the exact Rational shadow") !=
            std::string::npos;
  }
  if (!shadow_requested)
    throw std::runtime_error(
        "Acb deferred completeness did not request its exact Rational shadow");
}

}  // namespace

int main() {
  try {
    ComplexBall::set_precision(256);
    test_transaction_local_candidate();
    test_formal_parent_certificate();
    test_integer_shifted_tower_alignment();
    std::cout << "PASS: SCC formal parent completeness certificate\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
