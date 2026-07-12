#include "diffexp2/singular_indicial.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using diffexp2::BlockStep;
using diffexp2::JordanBlock;
using diffexp2::MatrixEntry;
using diffexp2::MatrixShift;
using diffexp2::PreparedRecurrenceOperator;
using diffexp2::Rational;
using diffexp2::RecurrenceError;
using diffexp2::StepCase;

namespace {

enum class Fault { None, Superdiagonal, RootMismatch, OffBlock };

int passed = 0;
int failed = 0;

void check(const std::string& label, bool condition) {
  if (condition) {
    ++passed;
    std::cout << "  PASS: " << label << '\n';
  } else {
    ++failed;
    std::cout << "  FAIL: " << label << '\n';
  }
}

void add(MatrixShift<Rational>& matrix, std::uint32_t row,
         std::uint32_t column, const char* value) {
  matrix.entries.push_back({row, column, Rational(value)});
}

PreparedRecurrenceOperator<Rational> prepared_operator(
    Fault fault = Fault::None) {
  PreparedRecurrenceOperator<Rational> prepared;
  prepared.dimension = 3;
  prepared.frame_base = -4;
  prepared.frame_width = 10;
  // d0 = 2 + eps.  Leaving d0_inverse_scalar empty exercises exact
  // cross-multiplication rather than a trusted scalar inverse.
  prepared.d_lags = {{{0, Rational(2)}, {1, Rational(1)}}};
  prepared.nhat_lags.resize(1);
  prepared.nhat_lags.front().valuations.assign(
      9, diffexp2::kCompleteInfinity);
  prepared.blocks = {JordanBlock{{0, 1}}, JordanBlock{{2}}};

  // In the declared spectral order,
  //   Nhat_0/d0 = Jordan(5/2 + 2 eps/3, 2) (+)
  //                 Jordan(1/2 + eps/3, 1).
  MatrixShift<Rational> eps0;
  eps0.shift = 0;
  add(eps0, 0, 0, "5");
  add(eps0, 1, 1, fault == Fault::RootMismatch ? "7" : "5");
  add(eps0, 0, 1, fault == Fault::Superdiagonal ? "4" : "2");
  add(eps0, 2, 2, "1");
  if (fault == Fault::OffBlock) add(eps0, 0, 2, "1");

  MatrixShift<Rational> eps1;
  eps1.shift = 1;
  add(eps1, 0, 0, "23/6");
  add(eps1, 1, 1, fault == Fault::RootMismatch ? "29/6" : "23/6");
  add(eps1, 0, 1, fault == Fault::Superdiagonal ? "2" : "1");
  add(eps1, 2, 2, "7/6");

  MatrixShift<Rational> eps2;
  eps2.shift = 2;
  add(eps2, 0, 0, "2/3");
  add(eps2, 1, 1, "2/3");
  add(eps2, 2, 2, "1/3");
  prepared.nhat_lags.front().polynomial =
      {std::move(eps0), std::move(eps1), std::move(eps2)};
  return prepared;
}

std::vector<std::vector<BlockStep<Rational>>> high_root_schedule() {
  return {{{StepCase::Resonant, Rational(0), Rational(0)},
           {StepCase::Taylor, Rational(2), Rational("1/3")}}};
}

std::vector<std::vector<BlockStep<Rational>>> low_root_schedule() {
  return {
      {{StepCase::Taylor, Rational(-2), Rational("-1/3")},
       {StepCase::Resonant, Rational(0), Rational(0)}},
      {{StepCase::Taylor, Rational(-1), Rational("-1/3")},
       {StepCase::Taylor, Rational(1), Rational(0)}},
      {{StepCase::Pseudo, Rational(0), Rational("-1/3")},
       {StepCase::Taylor, Rational(2), Rational(0)}}};
}

template <typename Function>
bool rejects(Function&& function, const std::string& fragment) {
  try {
    function();
  } catch (const RecurrenceError& error) {
    return error.id == "E5" &&
           std::string(error.what()).find(fragment) != std::string::npos;
  }
  return false;
}

}  // namespace

int main() {
  const auto prepared = prepared_operator();
  const auto indicial =
      diffexp2::certify_exact_affine_jordan_operator(prepared);
  check("size-two Jordan block and distinct affine roots are certified",
        indicial.dimension == 3 && indicial.blocks.size() == 2 &&
            indicial.blocks[0].columns == std::vector<std::uint32_t>({0, 1}) &&
            indicial.blocks[0].root.a == Rational("5/2") &&
            indicial.blocks[0].root.b == Rational("2/3") &&
            indicial.blocks[1].root.a == Rational("1/2") &&
            indicial.blocks[1].root.b == Rational("1/3"));

  const auto resonant = diffexp2::certify_exact_affine_jordan_schedule(
      indicial, Rational("5/2"), Rational("2/3"), high_root_schedule());
  check("valid exact T/R schedule retains the true Jordan log width",
        resonant.contains_true_resonance && !resonant.contains_pseudo &&
            resonant.max_resonant_jordan_chain_length == 2 &&
            resonant.steps[0][0].kind == StepCase::Resonant &&
            resonant.steps[0][0].resonant_jordan_chain_length == 2 &&
            resonant.steps[0][0].intrinsic_homogeneous_log_degree == 1 &&
            resonant.steps[0][1].kind == StepCase::Taylor);

  const auto pseudo = diffexp2::certify_exact_affine_jordan_schedule(
      indicial, Rational("1/2"), Rational("1/3"), low_root_schedule());
  check("P is classified exactly with the target Jordan pole depth",
        pseudo.contains_true_resonance && pseudo.contains_pseudo &&
            pseudo.max_pseudo_epsilon_pole_depth == 2 &&
            pseudo.steps[2][0].kind == StepCase::Pseudo &&
            pseudo.steps[2][0].pseudo_epsilon_pole_depth == 2 &&
            pseudo.steps[2][0].d_a.is_zero() &&
            pseudo.steps[2][0].d_b == Rational("-1/3"));

  check("malformed prepared Jordan superdiagonal is rejected exactly",
        rejects([] {
          (void)diffexp2::certify_exact_affine_jordan_operator(
              prepared_operator(Fault::Superdiagonal));
        }, "superdiagonal"));
  check("unequal roots inside one declared Jordan block are rejected",
        rejects([] {
          (void)diffexp2::certify_exact_affine_jordan_operator(
              prepared_operator(Fault::RootMismatch));
        }, "unequal affine roots"));
  check("nonzero entries outside the block-Jordan structure are rejected",
        rejects([] {
          (void)diffexp2::certify_exact_affine_jordan_operator(
              prepared_operator(Fault::OffBlock));
        }, "outside the declared"));

  auto wrong_offsets = high_root_schedule();
  wrong_offsets[0][0].d_a = Rational(1);
  check("submitted schedule offsets cannot override the retained root",
        rejects([&] {
          (void)diffexp2::certify_exact_affine_jordan_schedule(
              indicial, Rational("5/2"), Rational("2/3"), wrong_offsets);
        }, "offsets contradict"));
  auto wrong_case = low_root_schedule();
  wrong_case[2][0].kind = StepCase::Resonant;
  check("submitted R cannot disguise an exact pseudo collision",
        rejects([&] {
          (void)diffexp2::certify_exact_affine_jordan_schedule(
              indicial, Rational("1/2"), Rational("1/3"), wrong_case);
        }, "contradicts exact case P"));

  std::cout << passed << " passed, " << failed << " failed\n";
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
