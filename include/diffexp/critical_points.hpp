#pragma once
// The number of master integrals of a merged family, from its graph alone.
//
// With the Lee-Pomeransky polynomial G = U + F of a family, the number of
// critical points of G^nu0 * prod_i b_i^nu_i on the torus outside G = 0, for
// generic exponents, is the dimension of the family's space of integrals: its
// number of master integrals without symmetry relations. It is the vector
// space dimension of the likelihood ideal, which Singular computes over a
// finite field in seconds where a finite-field IBP elimination of the same
// family takes minutes. Measured against IBP counts: kite 16 = 16, planar
// double box 12 = 12, Higgs+jet family F 75 = 75 and 49 = 49, double
// pentagon 69 = 69.
//
// The count needs the `Singular` program on the PATH; without it every
// function here reports "unknown" and callers fall back to other means.
#include "diffexp/artifact_store.hpp"
#include "diffexp/feynman.hpp"
#include <flint/ulong_extras.h>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#if defined(__unix__) || defined(__APPLE__)
#include <cerrno>
#include <fcntl.h>
#include <signal.h>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>
extern char **environ;
#endif

namespace diffexp::critical_points {
inline constexpr unsigned long prime = 32003;
// The likelihood ideal of G in Singular's language. `lines[k][i]` is the
// coefficient of the original propagator i in the present line k; the lines
// `left` and `right` are merged at a generic point x.
inline std::string script(const feynman::Family &family, const std::vector<std::vector<Rational>> &lines, unsigned left, unsigned right) {
  const std::size_t originals = family.lines.size();
  if (lines.size() < 2 || left >= lines.size() || right >= lines.size() || left == right)
    throw std::invalid_argument("critical point count: merge positions");
  std::vector<std::string> names;
  for (std::size_t k = 0; k + 1 < lines.size(); ++k) names.push_back("b" + std::to_string(k));
  ExactField field(names);
  Exact sample(field);
  // A value with no reason to be special: at a singular point of the level's
  // differential equation the count drops (1/3 is one for family F).
  const Rational x("173/417"), y("244/417");
  std::vector<Exact> weights(originals, sample.constant(0));
  std::size_t variable = 1;
  for (std::size_t k = 0; k < lines.size(); ++k) {
    if (lines[k].size() != originals) throw std::invalid_argument("critical point count: line composition shape");
    if (k == right) continue;
    const auto b = sample.variable(k == left ? 0 : variable++);
    for (std::size_t i = 0; i < originals; ++i) {
      const Rational coefficient = k == left ? x * lines[left][i] + y * lines[right][i] : lines[k][i];
      if (!coefficient.is_zero()) weights[i] = weights[i] + b * sample.constant(coefficient);
    }
  }
  const auto geometry = feynman::symanzik(family, weights);
  const auto g = geometry.U + geometry.F;
  if (g.denominator_terms().size() != 1) throw std::domain_error("critical point count: U + F is not a polynomial");
  const auto residue = [](const Rational &q) {
    const unsigned long denominator = fmpz_fdiv_ui(fmpq_denref(q.raw()), prime);
    if (!denominator) throw std::domain_error("critical point count: coefficient denominator divisible by the prime");
    return n_mulmod2(fmpz_fdiv_ui(fmpq_numref(q.raw()), prime), n_invmod(denominator, prime), prime);
  };
  std::ostringstream out;
  out << "ring r=" << prime << ",(t";
  for (const auto &name : names) out << ',' << name;
  out << "),dp;\npoly G=";
  bool any = false;
  for (const auto &term : g.numerator_terms()) {
    const auto c = residue(term.coefficient);
    if (!c) continue;
    out << (any ? "+" : "") << c;
    for (std::size_t k = 0; k < names.size(); ++k)
      if (term.powers[k]) out << '*' << names[k] << '^' << term.powers[k];
    any = true;
  }
  if (!any) out << '0';
  // Fixed generic exponents: the count must not depend on them.
  unsigned long state = 88172645463325252ull;
  const auto exponent = [&] { state ^= state << 13; state ^= state >> 7; state ^= state << 17; return 2 + state % (prime - 4); };
  out << ";\nideal I=";
  const auto nu0 = exponent();
  for (const auto &name : names) out << name << '*' << nu0 << "*diff(G," << name << ")+" << exponent() << "*G,";
  out << "t";
  for (const auto &name : names) out << '*' << name;
  out << "*G-1;\nideal J=std(I);\nprint(dim(J));print(vdim(J));\nquit;\n";
  return out.str();
}
#if defined(__unix__) || defined(__APPLE__)
// Runs Singular on a script; nothing if it is missing, fails or is too slow.
inline std::optional<std::string> run(const std::string &text, const std::filesystem::path &directory, unsigned seconds) {
  std::error_code ignored;
  std::filesystem::create_directories(directory, ignored);
  const auto stem = directory / ("run-" + std::to_string(getpid()) + "-" + std::to_string(std::hash<std::thread::id>()(std::this_thread::get_id())));
  const auto input = stem.string() + ".sing", output = stem.string() + ".out";
  { std::ofstream file(input); file << text; if (!file) return std::nullopt; }
  posix_spawn_file_actions_t actions;
  posix_spawn_file_actions_init(&actions);
  posix_spawn_file_actions_addopen(&actions, STDIN_FILENO, input.c_str(), O_RDONLY, 0);
  posix_spawn_file_actions_addopen(&actions, STDOUT_FILENO, output.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0600);
  posix_spawn_file_actions_addopen(&actions, STDERR_FILENO, "/dev/null", O_WRONLY, 0);
  std::string program = "Singular", quiet = "-q";
  char *arguments[] = {program.data(), quiet.data(), nullptr};
  pid_t pid = 0;
  const int failed = posix_spawnp(&pid, program.c_str(), &actions, nullptr, arguments, environ);
  posix_spawn_file_actions_destroy(&actions);
  std::optional<std::string> result;
  if (!failed) {
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(seconds);
    int status = 0;
    for (;;) {
      const pid_t waited = waitpid(pid, &status, WNOHANG);
      if (waited == pid) break;
      if (waited < 0 && errno != EINTR) { status = -1; break; }
      if (std::chrono::steady_clock::now() >= deadline) {
        kill(pid, SIGKILL);
        while (waitpid(pid, &status, 0) < 0 && errno == EINTR) {}
        status = -1;
        break;
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    if (status != -1 && WIFEXITED(status) && WEXITSTATUS(status) == 0) {
      std::ifstream file(output);
      std::ostringstream all;
      all << file.rdbuf();
      result = all.str();
    }
  }
  std::filesystem::remove(input, ignored);
  std::filesystem::remove(output, ignored);
  return result;
}
inline bool available() {
  static const bool found = [] {
    const auto reply = run("print(1);quit;\n", std::filesystem::temp_directory_path() / "diffexp-critical-points", 20);
    return reply && reply->find('1') != std::string::npos;
  }();
  return found;
}
// The count for one merge, kept in `cache` under the hash of its script. So
// is the time within which it did not finish: a count that was given up is
// not waited for again by a later run unless that run allows more time.
inline std::optional<std::size_t> count(const feynman::Family &family, const std::vector<std::vector<Rational>> &lines, unsigned left, unsigned right,
                                        const std::filesystem::path &cache, unsigned seconds = 300, bool *gave_up = nullptr) {
  if (gave_up) *gave_up = false;
  if (!available()) return std::nullopt;
  std::string text;
  try { text = script(family, lines, left, right); } catch (const std::exception &) { return std::nullopt; }
  const auto path = cache / (artifacts::detail::sha256(text) + ".count");
  {
    std::ifstream stored(path);
    std::size_t value = 0;
    if (stored >> value) return value;
  }
  const auto given_up = cache / (artifacts::detail::sha256(text) + ".unfinished");
  {
    std::ifstream stored(given_up);
    unsigned tried = 0;
    if (stored >> tried && tried >= seconds) {
      if (gave_up) *gave_up = true;
      return std::nullopt;
    }
  }
  const auto began = std::chrono::steady_clock::now();
  const auto reply = run(text, cache, seconds);
  if (!reply) {
    if (std::chrono::steady_clock::now() - began >= std::chrono::seconds(seconds)) {
      std::ofstream stored(given_up);
      stored << seconds << '\n';
      if (gave_up) *gave_up = true;
    }
    return std::nullopt;
  }
  std::istringstream in(*reply);
  long dimension = -2, value = -1;
  if (!(in >> dimension >> value) || dimension != 0 || value < 0) return std::nullopt;
  std::ofstream stored(path);
  stored << value << '\n';
  return static_cast<std::size_t>(value);
}
#else
inline bool available() { return false; }
inline std::optional<std::size_t> count(const feynman::Family &, const std::vector<std::vector<Rational>> &, unsigned, unsigned, const std::filesystem::path &, unsigned = 300, bool * = nullptr) { return std::nullopt; }
#endif
} // namespace diffexp::critical_points
