#pragma once

#include <condition_variable>
#include <cstdint>
#include <exception>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>

namespace diffexp::kernel::detail {

class ImmutableCacheContractError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

class ImmutableCacheCycleError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

// Construct recursively derived immutable values once per key. Builders run
// outside the cache mutex, so independent keys may be constructed in
// parallel. Callers for the same key wait for its one in-flight builder.
//
// Recursive same-thread dependencies and wait cycles between builder threads
// are rejected before waiting. The wait graph is deliberately cache-local:
// builders may recurse through this cache, while dependencies spanning two
// different cache instances are outside this class's contract.
template <typename Key, typename Value, typename Compare = std::less<Key>>
class ImmutableRecursiveCache {
 public:
  struct Lookup {
    std::shared_ptr<const Value> value;
    bool built = false;
  };

  struct Stats {
    std::uint64_t entries = 0;
    std::uint64_t builds = 0;
    std::uint64_t hits = 0;
  };

  template <typename Builder>
  Lookup get_or_build(const Key& key, const std::string& contract,
                      Builder&& builder) {
    if (contract.empty())
      throw std::invalid_argument(
          "immutable recursive cache received an empty contract");

    const auto thread = std::this_thread::get_id();
    std::shared_ptr<BuildState> state;
    {
      std::unique_lock<std::mutex> lock(mutex_);
      if (const auto found = entries_.find(key); found != entries_.end()) {
        state = found->second;
        if (state->contract != contract)
          throw ImmutableCacheContractError(
              "immutable recursive cache contract changed for a retained key");
        if (state->phase == BuildPhase::Ready) {
          ++hits_;
          return {state->value, false};
        }
        if (state->phase != BuildPhase::Building)
          throw std::logic_error(
              "immutable recursive cache retained a failed build");
        if (state->owner == thread)
          throw ImmutableCacheCycleError(
              "immutable recursive cache dependency is cyclic");
        register_wait_locked(thread, state);
        try {
          state->changed.wait(lock, [&] {
            return state->phase != BuildPhase::Building;
          });
        } catch (...) {
          unregister_wait_locked(thread, state.get());
          state->waiters.erase(thread);
          throw;
        }
        // Completion normally removes this edge before notifying. Keep the
        // waiter-side removal as a defensive cleanup for exceptional waits.
        unregister_wait_locked(thread, state.get());
        state->waiters.erase(thread);
        if (state->phase == BuildPhase::Ready) {
          ++hits_;
          return {state->value, false};
        }
        if (!state->failure)
          throw std::logic_error(
              "immutable recursive cache failed without an exception");
        auto failure = state->failure;
        lock.unlock();
        std::rethrow_exception(failure);
      }

      state = std::make_shared<BuildState>(contract, thread);
      auto [stored, inserted] = entries_.emplace(key, state);
      if (!inserted)
        throw std::logic_error(
            "immutable recursive cache in-flight insertion failed");
    }

    std::shared_ptr<const Value> value;
    try {
      value = std::make_shared<const Value>(
          std::forward<Builder>(builder)());
    } catch (...) {
      auto failure = std::current_exception();
      {
        std::lock_guard<std::mutex> lock(mutex_);
        require_owned_build_locked(key, state, thread);
        state->phase = BuildPhase::Failed;
        state->failure = failure;
        state->owner = std::thread::id{};
        clear_waiters_locked(state);
        const auto found = entries_.find(key);
        if (found == entries_.end() || found->second != state)
          throw std::logic_error(
              "immutable recursive cache lost a failed build");
        entries_.erase(found);
      }
      state->changed.notify_all();
      std::rethrow_exception(failure);
    }

    {
      std::lock_guard<std::mutex> lock(mutex_);
      require_owned_build_locked(key, state, thread);
      state->value = std::move(value);
      state->phase = BuildPhase::Ready;
      state->owner = std::thread::id{};
      clear_waiters_locked(state);
      ++ready_entries_;
      ++builds_;
    }
    state->changed.notify_all();
    return {state->value, true};
  }

  Stats stats() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return {ready_entries_, builds_, hits_};
  }

 private:
  enum class BuildPhase { Building, Ready, Failed };

  struct BuildState {
    BuildState(std::string retained_contract, std::thread::id builder)
        : contract(std::move(retained_contract)), owner(builder) {}

    std::string contract;
    BuildPhase phase = BuildPhase::Building;
    std::thread::id owner;
    std::shared_ptr<const Value> value;
    std::exception_ptr failure;
    std::condition_variable changed;
    std::set<std::thread::id> waiters;
  };

  struct WaitEdge {
    std::thread::id owner;
    const BuildState* state = nullptr;
  };

  void register_wait_locked(
      const std::thread::id& waiter,
      const std::shared_ptr<BuildState>& state) {
    if (wait_for_.find(waiter) != wait_for_.end())
      throw std::logic_error(
          "immutable recursive cache thread already waits for a build");

    auto cursor = state->owner;
    std::set<std::thread::id> visited;
    while (true) {
      if (cursor == waiter)
        throw ImmutableCacheCycleError(
            "immutable recursive cache dependency is cyclic");
      if (!visited.insert(cursor).second)
        throw std::logic_error(
            "immutable recursive cache wait graph is already cyclic");
      const auto dependency = wait_for_.find(cursor);
      if (dependency == wait_for_.end()) break;
      cursor = dependency->second.owner;
    }

    const auto [edge, edge_inserted] = wait_for_.emplace(
        waiter, WaitEdge{state->owner, state.get()});
    if (!edge_inserted)
      throw std::logic_error(
          "immutable recursive cache wait registration failed");
    try {
      if (!state->waiters.insert(waiter).second) {
        wait_for_.erase(edge);
        throw std::logic_error(
            "immutable recursive cache waiter registration failed");
      }
    } catch (...) {
      const auto retained = wait_for_.find(waiter);
      if (retained != wait_for_.end() && retained->second.state == state.get())
        wait_for_.erase(retained);
      throw;
    }
  }

  void unregister_wait_locked(const std::thread::id& waiter,
                              const BuildState* state) {
    const auto found = wait_for_.find(waiter);
    if (found != wait_for_.end() && found->second.state == state)
      wait_for_.erase(found);
  }

  void clear_waiters_locked(const std::shared_ptr<BuildState>& state) {
    for (const auto& waiter : state->waiters)
      unregister_wait_locked(waiter, state.get());
    state->waiters.clear();
  }

  void require_owned_build_locked(
      const Key& key, const std::shared_ptr<BuildState>& state,
      const std::thread::id& thread) const {
    const auto found = entries_.find(key);
    if (found == entries_.end() || found->second != state ||
        state->phase != BuildPhase::Building || state->owner != thread)
      throw std::logic_error(
          "immutable recursive cache lost ownership of an in-flight build");
  }

  mutable std::mutex mutex_;
  std::map<Key, std::shared_ptr<BuildState>, Compare> entries_;
  std::map<std::thread::id, WaitEdge> wait_for_;
  std::uint64_t ready_entries_ = 0;
  std::uint64_t builds_ = 0;
  std::uint64_t hits_ = 0;
};

}  // namespace diffexp::kernel::detail
