#pragma once

#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace diffexp2::detail {

class ImmutableCacheContractError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

class ImmutableCacheCycleError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

// A small single-writer cache for recursively derived immutable values.
// Construction is serialized deliberately: a builder may request another key,
// and the explicit in-flight set turns A -> A or A -> B -> A into a loud error
// instead of a promise deadlock. Completed values are shared read-only.
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
    std::lock_guard<std::recursive_mutex> lock(mutex_);
    if (const auto found = entries_.find(key); found != entries_.end()) {
      if (found->second.contract != contract)
        throw ImmutableCacheContractError(
            "immutable recursive cache contract changed for a retained key");
      ++hits_;
      return {found->second.value, false};
    }
    if (!active_builds_.insert(key).second)
      throw ImmutableCacheCycleError(
          "immutable recursive cache dependency is cyclic");

    std::shared_ptr<const Value> value;
    try {
      value = std::make_shared<const Value>(
          std::forward<Builder>(builder)());
    } catch (...) {
      active_builds_.erase(key);
      throw;
    }
    active_builds_.erase(key);
    auto [stored, inserted] = entries_.emplace(
        key, Entry{contract, std::move(value)});
    if (!inserted)
      throw std::logic_error(
          "immutable recursive cache insertion failed");
    ++builds_;
    return {stored->second.value, true};
  }

  Stats stats() const {
    std::lock_guard<std::recursive_mutex> lock(mutex_);
    return {static_cast<std::uint64_t>(entries_.size()), builds_, hits_};
  }

 private:
  struct Entry {
    std::string contract;
    std::shared_ptr<const Value> value;
  };

  mutable std::recursive_mutex mutex_;
  std::map<Key, Entry, Compare> entries_;
  std::set<Key, Compare> active_builds_;
  std::uint64_t builds_ = 0;
  std::uint64_t hits_ = 0;
};

}  // namespace diffexp2::detail
