#pragma once
// Laurent coefficients in epsilon from values at shifted dimensions.
//
// `diffexp ft --epsilon-shift Q --epsilon-order 0` evaluates an integral in
// d0-2Q dimensions. A function with a pole of order P at epsilon=0 times
// epsilon^P is analytic there, so the polynomial through
// (Q_i, Q_i^P f(Q_i)) at symmetric nodes approximates its Taylor series and
// its coefficient k is the Laurent coefficient of epsilon^(k-P). Samples are
// independent of each other: they run as separate processes of this
// executable on the prepared exact cache.
#include "diffexp/kernel/scalar.hpp"
#include <boost/json.hpp>
#include <chrono>
#include <cerrno>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <vector>
#if defined(__unix__) || defined(__APPLE__)
#include <fcntl.h>
#include <signal.h>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>
extern char **environ;
#endif

namespace diffexp::epsilon_sampling {
using B = kernel::ComplexBall;
using Rational = kernel::Rational;

// The nodes +step, -step, +2 step, -2 step, ...: innermost first, so that
// dropping the last two leaves a symmetric set again.
inline std::vector<Rational> nodes(const Rational &step, unsigned samples) {
  if (samples < 2 || samples % 2 || step.sign() <= 0)
    throw std::invalid_argument("epsilon sampling needs an even number of samples and a positive step");
  std::vector<Rational> out;
  for (unsigned k = 1; k <= samples / 2; ++k) {
    out.push_back(step * Rational(static_cast<long>(k)));
    out.push_back(-(step * Rational(static_cast<long>(k))));
  }
  return out;
}
// Coefficients of the polynomial through (x_i, x_i^poles y_i), lowest first:
// Newton's divided differences, then Horner on the Newton basis. Ball
// arithmetic carries the radii of the samples through the cancellations.
inline std::vector<B> interpolate(const std::vector<Rational> &x, const std::vector<B> &y, unsigned poles) {
  if (x.size() != y.size() || x.empty())
    throw std::invalid_argument("epsilon sampling interpolation shape");
  const auto ball = [](const Rational &q) { return B::from_strings(q.str()); };
  std::vector<B> table;
  for (std::size_t i = 0; i < x.size(); ++i) {
    B scaled = y[i];
    for (unsigned p = 0; p < poles; ++p) scaled = scaled * ball(x[i]);
    table.push_back(scaled);
  }
  for (std::size_t level = 1; level < x.size(); ++level)
    for (std::size_t i = x.size() - 1; i >= level; --i)
      table[i] = (table[i] - table[i - 1]) / ball(x[i] - x[i - level]);
  std::vector<B> coefficients(x.size(), B(0));
  for (std::size_t i = x.size(); i-- > 0;) {
    std::vector<B> next(x.size(), B(0));
    for (std::size_t k = 0; k < x.size(); ++k) {
      if (k) next[k] = coefficients[k - 1];
      next[k] = next[k] - ball(x[i]) * coefficients[k];
    }
    next[0] = next[0] + table[i];
    coefficients = std::move(next);
  }
  return coefficients;
}
inline double log2_magnitude(const B &value) {
  mag_t m;
  mag_init(m);
  acb_get_mag(m, value.raw());
  const double out = mag_is_zero(m) ? -std::numeric_limits<double>::infinity() : mag_get_d_log2_approx(m);
  mag_clear(m);
  return out;
}
// Order of the pole at zero from the growth between the outermost and the
// innermost positive node: f ~ c/Q^P gives |f(near)|/|f(far)| = (far/near)^P.
inline unsigned pole_order(const std::vector<Rational> &x, const std::vector<std::vector<B>> &values, unsigned bound) {
  if (x.size() < 4) return 0;
  const std::size_t near = 0, far = x.size() - 2;
  const auto number = [](const Rational &q) {
    const auto text = q.str();
    const auto slash = text.find('/');
    return slash == std::string::npos ? std::stod(text) : std::stod(text.substr(0, slash)) / std::stod(text.substr(slash + 1));
  };
  const double spread = std::log2(number(x[far]) / number(x[near]));
  double order = 0;
  for (const auto &component : values) {
    const double a = log2_magnitude(component[near]), b = log2_magnitude(component[far]);
    if (!std::isfinite(a) || !std::isfinite(b)) continue;
    order = std::max(order, (a - b) / spread);
  }
  return static_cast<unsigned>(std::clamp<long>(std::lround(order), 0, bound));
}
struct Fit {
  unsigned poles = 0;
  // By component, then epsilon order from -poles: the coefficient, and the
  // change it makes to drop the two outermost samples.
  std::vector<std::vector<B>> coefficients;
  std::vector<std::vector<double>> log2_change;
};
// The balls returned are widened by that change: an estimate of the
// interpolation error, not a bound.
inline Fit fit(const std::vector<Rational> &x, const std::vector<std::vector<B>> &values, unsigned poles, unsigned high) {
  const std::size_t wanted = static_cast<std::size_t>(poles) + high + 1;
  if (x.size() < wanted + 2)
    throw std::invalid_argument("epsilon sampling: " + std::to_string(x.size()) + " samples cannot resolve a pole of order " +
                                std::to_string(poles) + " through epsilon^" + std::to_string(high) + "; at least " +
                                std::to_string(wanted + 2 + (wanted % 2)) + " are needed");
  Fit out;
  out.poles = poles;
  const std::vector<Rational> inner(x.begin(), x.end() - 2);
  for (const auto &component : values) {
    const auto all = interpolate(x, component, poles);
    const auto fewer = interpolate(inner, std::vector<B>(component.begin(), component.end() - 2), poles);
    std::vector<B> row;
    std::vector<double> changes;
    for (std::size_t k = 0; k < wanted; ++k) {
      const B change = all[k] - fewer[k];
      mag_t m;
      mag_init(m);
      acb_get_mag(m, change.raw());
      changes.push_back(mag_is_zero(m) ? -std::numeric_limits<double>::infinity() : mag_get_d_log2_approx(m));
      B widened = all[k];
      acb_add_error_mag(widened.raw(), m);
      mag_clear(m);
      row.push_back(widened);
    }
    out.coefficients.push_back(std::move(row));
    out.log2_change.push_back(std::move(changes));
  }
  return out;
}

#if defined(__unix__) || defined(__APPLE__)
struct Sample {
  Rational epsilon;
  std::vector<B> values; // by component
  double seconds = 0;
};
// Runs `command` + {"--epsilon-shift", Q} for every node, at most `jobs` at a
// time, each with its standard output in `directory`. A sample that fails
// stops the run with its own message.
inline std::vector<Sample> run(const std::vector<std::string> &command, const std::vector<Rational> &points,
                               const std::filesystem::path &directory, unsigned jobs, unsigned threads_per_sample,
                               const std::function<void(const Sample &, std::size_t done, std::size_t total)> &progress = {}) {
  std::filesystem::create_directories(directory);
  std::vector<Sample> samples(points.size());
  std::map<pid_t, std::size_t> running;
  std::size_t next = 0, done = 0;
  const auto output = [&](std::size_t i, const char *suffix) {
    return directory / ("sample-" + std::to_string(getpid()) + "-" + std::to_string(i) + suffix);
  };
  const auto stop = [&] {
    for (const auto &[pid, index] : running) kill(pid, SIGTERM);
    for (const auto &[pid, index] : running) {
      int status = 0;
      while (waitpid(pid, &status, 0) < 0 && errno == EINTR) {}
    }
    running.clear();
  };
  const auto launch = [&](std::size_t i) {
    std::vector<std::string> arguments = command;
    arguments.push_back("--epsilon-shift");
    arguments.push_back(points[i].str());
    std::vector<char *> argv;
    for (auto &argument : arguments) argv.push_back(argument.data());
    argv.push_back(nullptr);
    std::vector<std::string> variables;
    bool threads_given = false;
    for (char **entry = environ; *entry; ++entry) {
      variables.emplace_back(*entry);
      if (variables.back().rfind("OMP_NUM_THREADS=", 0) == 0) threads_given = true;
    }
    if (!threads_given) variables.push_back("OMP_NUM_THREADS=" + std::to_string(std::max(1u, threads_per_sample)));
    std::vector<char *> envp;
    for (auto &variable : variables) envp.push_back(variable.data());
    envp.push_back(nullptr);
    posix_spawn_file_actions_t actions;
    posix_spawn_file_actions_init(&actions);
    const auto out = output(i, ".json").string(), err = output(i, ".log").string();
    posix_spawn_file_actions_addopen(&actions, STDOUT_FILENO, out.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0600);
    posix_spawn_file_actions_addopen(&actions, STDERR_FILENO, err.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0600);
    pid_t pid = 0;
    const int failed = posix_spawn(&pid, argv[0], &actions, nullptr, argv.data(), envp.data());
    posix_spawn_file_actions_destroy(&actions);
    if (failed) throw std::runtime_error("epsilon sampling: cannot start " + arguments[0]);
    running.emplace(pid, i);
  };
  try {
    while (done < points.size()) {
      while (next < points.size() && running.size() < std::max(1u, jobs)) launch(next++);
      int status = 0;
      const pid_t pid = waitpid(-1, &status, 0);
      if (pid < 0) {
        if (errno == EINTR) continue;
        throw std::runtime_error("epsilon sampling: lost its sample processes");
      }
      const auto found = running.find(pid);
      if (found == running.end()) continue;
      const std::size_t i = found->second;
      running.erase(found);
      std::ifstream in(output(i, ".json"));
      std::string line, last;
      while (std::getline(in, line))
        if (!line.empty()) last = line;
      boost::json::value parsed;
      try {
        parsed = boost::json::parse(last);
      } catch (const std::exception &) {
        throw std::runtime_error("epsilon sampling: the sample at epsilon=" + points[i].str() + " returned no result (log: " +
                                 output(i, ".log").string() + ")");
      }
      const auto &report = parsed.as_object();
      if (!report.contains("coefficients")) {
        const auto *message = report.if_contains("message");
        throw std::runtime_error("epsilon sampling: the sample at epsilon=" + points[i].str() + " failed: " +
                                 (message && message->is_string() ? std::string(message->as_string()) : "no message"));
      }
      const long low = report.at("epsilon_low").to_number<long>();
      Sample sample{points[i], {}, report.at("timings").as_object().at("numerical_seconds").to_number<double>()};
      for (const auto &row : report.at("coefficients").as_array()) {
        const auto &entry = row.as_array().at(static_cast<std::size_t>(-low)).as_object();
        sample.values.push_back(B::from_strings(std::string(entry.at("real").as_string()), std::string(entry.at("imaginary").as_string())));
      }
      samples[i] = std::move(sample);
      std::filesystem::remove(output(i, ".json"));
      std::filesystem::remove(output(i, ".log"));
      ++done;
      if (progress) progress(samples[i], done, points.size());
    }
  } catch (...) {
    stop();
    throw;
  }
  return samples;
}
#endif
} // namespace diffexp::epsilon_sampling
