#include "diffexp/kernel/immutable_recursive_cache.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdlib>
#include <iostream>
#include <latch>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using Cache = diffexp::kernel::detail::ImmutableRecursiveCache<int, int>;
using namespace std::chrono_literals;

class TimedGate {
 public:
  explicit TimedGate(std::size_t target) : target_(target) {}

  bool arrive_and_wait() {
    std::unique_lock<std::mutex> lock(mutex_);
    ++arrived_;
    changed_.notify_all();
    return changed_.wait_for(lock, 2s, [&] { return arrived_ >= target_; });
  }

 private:
  const std::size_t target_;
  std::size_t arrived_ = 0;
  std::mutex mutex_;
  std::condition_variable changed_;
};

void update_peak(std::atomic<int>& peak, int value) {
  auto observed = peak.load();
  while (observed < value &&
         !peak.compare_exchange_weak(observed, value)) {}
}

struct Checks {
  void require(bool condition, const std::string& name) {
    if (condition) {
      ++passed;
    } else {
      ++failed;
      std::cerr << "FAIL: " << name << '\n';
    }
  }

  int passed = 0;
  int failed = 0;
};

void test_basic_contract_and_stats(Checks& checks) {
  Cache cache;
  auto first = cache.get_or_build(1, "one", [] { return 11; });
  bool duplicate_builder_ran = false;
  auto second = cache.get_or_build(1, "one", [&] {
    duplicate_builder_ran = true;
    return 12;
  });
  bool contract_error = false;
  try {
    (void)cache.get_or_build(1, "different", [] { return 13; });
  } catch (const diffexp::kernel::detail::ImmutableCacheContractError&) {
    contract_error = true;
  }
  const auto stats = cache.stats();
  checks.require(first.built && !second.built &&
                     first.value == second.value && *second.value == 11 &&
                     !duplicate_builder_ran,
                 "ready lookup reuses one immutable value");
  checks.require(contract_error,
                 "ready lookup rejects a changed contract");
  checks.require(stats.entries == 1 && stats.builds == 1 && stats.hits == 1,
                 "basic cache statistics retain their public meaning");
}

void test_distinct_keys_build_concurrently(Checks& checks) {
  Cache cache;
  TimedGate builders(2);
  std::latch start(1);
  std::atomic<int> active{0};
  std::atomic<int> peak{0};
  std::atomic<int> failures{0};
  auto run = [&](int key) {
    start.wait();
    try {
      (void)cache.get_or_build(key, "key:" + std::to_string(key), [&] {
        const auto now = active.fetch_add(1) + 1;
        update_peak(peak, now);
        if (!builders.arrive_and_wait()) {
          active.fetch_sub(1);
          throw std::runtime_error("distinct-key builders were serialized");
        }
        active.fetch_sub(1);
        return key * 10;
      });
    } catch (...) {
      ++failures;
    }
  };
  std::jthread left(run, 1);
  std::jthread right(run, 2);
  start.count_down();
  left.join();
  right.join();
  const auto stats = cache.stats();
  checks.require(failures == 0 && peak == 2,
                 "distinct keys construct concurrently");
  checks.require(stats.entries == 2 && stats.builds == 2 && stats.hits == 0,
                 "distinct concurrent builds publish exactly once");
}

void test_same_key_single_flight(Checks& checks) {
  constexpr int kCallers = 8;
  Cache cache;
  std::latch ready(kCallers);
  std::latch start(1);
  std::atomic<int> attempted{0};
  std::atomic<int> builders{0};
  std::atomic<int> built_results{0};
  std::atomic<int> failures{0};
  std::vector<std::shared_ptr<const int>> values(kCallers);
  std::vector<std::jthread> threads;
  threads.reserve(kCallers);
  for (int index = 0; index < kCallers; ++index) {
    threads.emplace_back([&, index] {
      ready.count_down();
      start.wait();
      ++attempted;
      try {
        auto result = cache.get_or_build(7, "seven", [&] {
          ++builders;
          while (attempted.load() != kCallers) std::this_thread::yield();
          return 77;
        });
        values[index] = result.value;
        if (result.built) ++built_results;
      } catch (...) {
        ++failures;
      }
    });
  }
  ready.wait();
  start.count_down();
  for (auto& thread : threads) thread.join();
  const auto reference = values.front();
  const bool identical = reference &&
      std::all_of(values.begin(), values.end(), [&](const auto& value) {
        return value == reference && *value == 77;
      });
  const auto stats = cache.stats();
  checks.require(failures == 0 && builders == 1 && built_results == 1 &&
                     identical,
                 "same-key callers share one in-flight build");
  checks.require(stats.entries == 1 && stats.builds == 1 &&
                     stats.hits == kCallers - 1,
                 "same-key waiters are counted as cache hits");
}

void test_inflight_contract_validation(Checks& checks) {
  Cache cache;
  std::latch builder_started(1);
  std::latch release_builder(1);
  std::atomic<int> owner_failures{0};
  std::jthread owner([&] {
    try {
      (void)cache.get_or_build(3, "stable", [&] {
        builder_started.count_down();
        release_builder.wait();
        return 30;
      });
    } catch (...) {
      ++owner_failures;
    }
  });
  builder_started.wait();
  bool contract_error = false;
  try {
    (void)cache.get_or_build(3, "changed", [] { return 31; });
  } catch (const diffexp::kernel::detail::ImmutableCacheContractError&) {
    contract_error = true;
  }
  release_builder.count_down();
  owner.join();
  const auto stats = cache.stats();
  checks.require(contract_error && owner_failures == 0,
                 "in-flight lookup validates the retained contract");
  checks.require(stats.entries == 1 && stats.builds == 1 && stats.hits == 0,
                 "contract rejection does not perturb build statistics");
}

void test_same_thread_cycles_and_retry(Checks& checks) {
  Cache direct;
  bool direct_cycle = false;
  try {
    (void)direct.get_or_build(1, "one", [&]() -> int {
      return *direct.get_or_build(1, "one", [] { return 2; }).value;
    });
  } catch (const diffexp::kernel::detail::ImmutableCacheCycleError&) {
    direct_cycle = true;
  }
  const auto direct_after_failure = direct.stats();
  auto retried = direct.get_or_build(1, "one", [] { return 5; });
  checks.require(direct_cycle && direct_after_failure.entries == 0 &&
                     direct_after_failure.builds == 0 && retried.built &&
                     *retried.value == 5,
                 "direct recursive cycle fails cleanly and can retry");

  Cache indirect;
  bool indirect_cycle = false;
  try {
    (void)indirect.get_or_build(1, "one", [&]() -> int {
      return *indirect.get_or_build(2, "two", [&]() -> int {
        return *indirect.get_or_build(1, "one", [] { return 3; }).value;
      }).value;
    });
  } catch (const diffexp::kernel::detail::ImmutableCacheCycleError&) {
    indirect_cycle = true;
  }
  const auto indirect_stats = indirect.stats();
  checks.require(indirect_cycle && indirect_stats.entries == 0 &&
                     indirect_stats.builds == 0,
                 "indirect same-thread cycle unwinds every in-flight key");
}

void test_two_thread_wait_cycle(Checks& checks) {
  Cache cache;
  TimedGate builders(2);
  std::atomic<int> cycles{0};
  std::atomic<int> other_failures{0};
  auto run = [&](int own, int dependency) {
    try {
      (void)cache.get_or_build(own, "key:" + std::to_string(own), [&]() -> int {
        if (!builders.arrive_and_wait())
          throw std::runtime_error("cycle builders did not overlap");
        return *cache.get_or_build(
            dependency, "key:" + std::to_string(dependency),
            [dependency] { return dependency; }).value;
      });
    } catch (const diffexp::kernel::detail::ImmutableCacheCycleError&) {
      ++cycles;
    } catch (...) {
      ++other_failures;
    }
  };
  std::jthread left(run, 1, 2);
  std::jthread right(run, 2, 1);
  left.join();
  right.join();
  const auto stats = cache.stats();
  checks.require(cycles == 2 && other_failures == 0,
                 "two-thread wait cycle is detected and propagated");
  checks.require(stats.entries == 0 && stats.builds == 0,
                 "cross-thread cycle leaves no poisoned entries");
}

void test_three_thread_wait_cycle(Checks& checks) {
  Cache cache;
  TimedGate builders(3);
  std::atomic<int> cycles{0};
  std::atomic<int> other_failures{0};
  auto run = [&](int own, int dependency) {
    try {
      (void)cache.get_or_build(own, "key:" + std::to_string(own), [&]() -> int {
        if (!builders.arrive_and_wait())
          throw std::runtime_error("cycle builders did not overlap");
        return *cache.get_or_build(
            dependency, "key:" + std::to_string(dependency),
            [dependency] { return dependency; }).value;
      });
    } catch (const diffexp::kernel::detail::ImmutableCacheCycleError&) {
      ++cycles;
    } catch (...) {
      ++other_failures;
    }
  };
  std::jthread one(run, 1, 2);
  std::jthread two(run, 2, 3);
  std::jthread three(run, 3, 1);
  one.join();
  two.join();
  three.join();
  const auto stats = cache.stats();
  checks.require(cycles == 3 && other_failures == 0,
                 "three-thread wait cycle is detected and propagated");
  checks.require(stats.entries == 0 && stats.builds == 0,
                 "multi-thread cycle cleanup removes every generation");
}

void test_builder_failure_retry(Checks& checks) {
  Cache cache;
  bool failed = false;
  try {
    (void)cache.get_or_build(9, "nine", []() -> int {
      throw std::runtime_error("expected build failure");
    });
  } catch (const std::runtime_error& error) {
    failed = std::string(error.what()) == "expected build failure";
  }
  const auto failed_stats = cache.stats();
  auto retry = cache.get_or_build(9, "nine", [] { return 99; });
  const auto retry_stats = cache.stats();
  checks.require(failed && failed_stats.entries == 0 &&
                     failed_stats.builds == 0,
                 "failed builder is not retained or counted");
  checks.require(retry.built && *retry.value == 99 &&
                     retry_stats.entries == 1 && retry_stats.builds == 1,
                 "failed key can be rebuilt successfully");
}

}  // namespace

int main() {
  Checks checks;
  test_basic_contract_and_stats(checks);
  test_distinct_keys_build_concurrently(checks);
  test_same_key_single_flight(checks);
  test_inflight_contract_validation(checks);
  test_same_thread_cycles_and_retry(checks);
  test_two_thread_wait_cycle(checks);
  test_three_thread_wait_cycle(checks);
  test_builder_failure_retry(checks);
  std::cout << "Immutable recursive cache tests: " << checks.passed
            << " passed, " << checks.failed << " failed\n";
  return checks.failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
