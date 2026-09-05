#include "diffexp/kernel/checkpoint.hpp"

#include <boost/crc.hpp>
#include <boost/json.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <charconv>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <limits>
#include <new>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

namespace diffexp::kernel::checkpoint {
namespace {

namespace json = boost::json;

constexpr std::array<unsigned char, 8> kMagic{
    'D', 'E', '2', 'C', 'P', '0', '0', '1'};
constexpr std::size_t kFixedHeaderBytes = 32;
// Darwin rejects a single write(2) whose byte count exceeds its internal
// syscall limit with EINVAL, even though size_t and ssize_t can represent the
// complete request.  Checkpoints may legitimately contain multi-gigabyte
// retained local slabs, so keep every read/write syscall comfortably below
// platform limits and stream the container sections incrementally.
constexpr std::size_t kMaximumIoChunkBytes =
    std::size_t{1} * 1024 * 1024;
constexpr std::uint64_t kMaximumPayloadBytes =
    std::uint64_t{8} * 1024 * 1024 * 1024;
constexpr std::string_view kCanonicalJsonDagCodec =
    "diffexp3-canonical-json-string-dag-v1";
constexpr std::size_t kMinimumInternedStringBytes = 128;
constexpr std::size_t kMaximumCanonicalJsonEntries = 1'000'000;
constexpr std::uint64_t kMinimumDecodedJsonBudget =
    std::uint64_t{64} * 1024 * 1024;
constexpr std::uint64_t kMaximumDecodedJsonBudget =
    std::uint64_t{8} * 1024 * 1024 * 1024;
// Large retained singular charts can legitimately exceed the ordinary
// decoded-string ceiling even though their canonical-DAG container remains
// below the 8 GiB payload limit.  Keep the normal untrusted-input limit, but
// permit an explicit diagnostic restore of a locally produced checkpoint.
// This is deliberately bounded so the override cannot disable expansion
// accounting altogether.
constexpr std::uint64_t kMaximumDiagnosticDecodedJsonBudget =
    std::uint64_t{32} * 1024 * 1024 * 1024;
constexpr std::uint64_t kDecodedToEncodedBudgetRatio = 1024;
constexpr std::size_t kMaximumCanonicalJsonReferenceDepth = 4096;
constexpr char kEncodedStringPrefix = '\0';
constexpr char kEncodedLiteralTag = 'L';
constexpr char kEncodedReferenceTag = 'R';

std::runtime_error system_error(const std::string& action,
                                const std::string& path) {
  return std::runtime_error(action + " " + path + ": " +
                            std::strerror(errno));
}

std::uint32_t checksum(std::string_view bytes) {
  boost::crc_32_type crc;
  crc.process_bytes(bytes.data(), bytes.size());
  return crc.checksum();
}

std::string dump_arb_exact(const arb_t value) {
  char* raw = arb_dump_str(value);
  if (raw == nullptr) throw std::bad_alloc();
  std::string output(raw);
  flint_free(raw);
  return output;
}

void load_arb_exact(arb_t output, const std::string& dump,
                    const char* component) {
  if (dump.empty() || arb_load_str(output, dump.c_str()) != 0)
    throw std::invalid_argument(
        std::string("invalid exact Arb checkpoint ") + component +
        " component");
}

void append_u32(std::vector<unsigned char>& output, std::uint32_t value) {
  for (int shift = 24; shift >= 0; shift -= 8)
    output.push_back(static_cast<unsigned char>((value >> shift) & 0xffu));
}

void append_u64(std::vector<unsigned char>& output, std::uint64_t value) {
  for (int shift = 56; shift >= 0; shift -= 8)
    output.push_back(static_cast<unsigned char>((value >> shift) & 0xffu));
}

std::uint32_t read_u32(const unsigned char* bytes) {
  std::uint32_t value = 0;
  for (int index = 0; index < 4; ++index)
    value = (value << 8) | bytes[index];
  return value;
}

std::uint64_t read_u64(const unsigned char* bytes) {
  std::uint64_t value = 0;
  for (int index = 0; index < 8; ++index)
    value = (value << 8) | bytes[index];
  return value;
}

void write_all(int descriptor, const unsigned char* data, std::size_t size,
               const std::string& path) {
  std::size_t written = 0;
  while (written < size) {
    const auto request =
        std::min(kMaximumIoChunkBytes, size - written);
    const auto amount = ::write(descriptor, data + written, request);
    if (amount < 0) {
      if (errno == EINTR) continue;
      throw system_error("cannot write checkpoint", path);
    }
    if (amount == 0)
      throw std::runtime_error("short write while writing checkpoint " + path);
    written += static_cast<std::size_t>(amount);
  }
}

void read_all(int descriptor, unsigned char* data, std::size_t size,
              const std::string& path) {
  std::size_t read_bytes = 0;
  while (read_bytes < size) {
    const auto request =
        std::min(kMaximumIoChunkBytes, size - read_bytes);
    const auto amount = ::read(descriptor, data + read_bytes, request);
    if (amount < 0) {
      if (errno == EINTR) continue;
      throw system_error("cannot read checkpoint", path);
    }
    if (amount == 0)
      throw std::runtime_error("truncated checkpoint " + path);
    read_bytes += static_cast<std::size_t>(amount);
  }
}

std::string temporary_path(const std::string& path) {
  static std::atomic<std::uint64_t> sequence{0};
  return path + ".tmp." + std::to_string(static_cast<long long>(::getpid())) +
         "." + std::to_string(sequence.fetch_add(1) + 1);
}

void close_checked(int descriptor, const std::string& path) {
  if (::close(descriptor) != 0)
    throw system_error("cannot close checkpoint", path);
}

void require_exact_keys(const json::object& object,
                        std::initializer_list<std::string_view> expected,
                        const char* label) {
  if (object.size() != expected.size())
    throw std::invalid_argument(std::string(label) +
                                " has unknown or missing fields");
  for (const auto key : expected)
    if (object.if_contains(key) == nullptr)
      throw std::invalid_argument(std::string(label) +
                                  " has unknown or missing fields");
}

class CanonicalJsonStringDagEncoder {
 public:
  json::value encode(const json::value& value) { return encode_value(value); }

  json::array release_entries() { return std::move(entries_); }

 private:
  json::value encode_value(const json::value& value) {
    if (value.is_string()) return encode_string(value.as_string());
    if (value.is_array()) {
      json::array encoded;
      encoded.reserve(value.as_array().size());
      for (const auto& item : value.as_array())
        encoded.push_back(encode_value(item));
      return encoded;
    }
    if (value.is_object()) {
      json::object encoded;
      for (const auto& item : value.as_object())
        encoded[item.key()] = encode_value(item.value());
      return encoded;
    }
    return value;
  }

  json::value encode_string(const json::string& raw) {
    const std::string text(raw.data(), raw.size());
    if (!text.empty() && text.front() == kEncodedStringPrefix)
      return json::value(json::string(
          std::string{kEncodedStringPrefix, kEncodedLiteralTag} + text));
    if (text.size() < kMinimumInternedStringBytes ||
        (text.front() != '{' && text.front() != '['))
      return json::value(json::string(text));
    if (const auto found = identities_.find(text);
        found != identities_.end())
      return json::value(json::string(reference(found->second)));

    json::value parsed;
    try {
      parsed = json::parse(text);
    } catch (const std::exception&) {
      return json::value(json::string(text));
    }
    if ((!parsed.is_object() && !parsed.is_array()) ||
        json::serialize(parsed) != text)
      return json::value(json::string(text));
    if (entries_.size() >= kMaximumCanonicalJsonEntries)
      throw std::invalid_argument(
          "checkpoint canonical-JSON identity table is too large");
    const auto index = entries_.size();
    identities_.emplace(text, index);
    entries_.emplace_back(nullptr);
    auto encoded = encode_value(parsed);
    entries_[index] = std::move(encoded);
    return json::value(json::string(reference(index)));
  }

  static std::string reference(std::size_t index) {
    return std::string{kEncodedStringPrefix, kEncodedReferenceTag} +
           std::to_string(index);
  }

  json::array entries_;
  std::unordered_map<std::string, std::size_t> identities_;
};

class CanonicalJsonStringDagDecoder {
 public:
  CanonicalJsonStringDagDecoder(const json::array& entries,
                                std::uint64_t decoded_budget)
      : entries_(entries), decoded_(entries.size()), state_(entries.size()),
        decoded_budget_(decoded_budget) {
    if (entries.size() > kMaximumCanonicalJsonEntries)
      throw std::invalid_argument(
          "checkpoint canonical-JSON identity table is too large");
    if (decoded_budget == 0 ||
        decoded_budget > kMaximumDiagnosticDecodedJsonBudget)
      throw std::invalid_argument(
          "checkpoint canonical-JSON decoded-size budget is invalid");
  }

  json::value decode(const json::value& value) { return decode_value(value); }

  void require_all_entries_used() const {
    if (std::any_of(state_.begin(), state_.end(),
                    [](unsigned char state) { return state != 2; }))
      throw std::invalid_argument(
          "checkpoint canonical-JSON identity table contains unused entries");
  }

 private:
  json::value decode_value(const json::value& value) {
    if (value.is_string()) return decode_string(value.as_string());
    if (value.is_array()) {
      json::array decoded;
      decoded.reserve(value.as_array().size());
      for (const auto& item : value.as_array())
        decoded.push_back(decode_value(item));
      return decoded;
    }
    if (value.is_object()) {
      json::object decoded;
      for (const auto& item : value.as_object())
        decoded[item.key()] = decode_value(item.value());
      return decoded;
    }
    return value;
  }

  json::value decode_string(const json::string& raw) {
    const std::string text(raw.data(), raw.size());
    if (text.empty() || text.front() != kEncodedStringPrefix) {
      charge(text.size());
      return json::value(json::string(text));
    }
    if (text.size() >= 2 && text[1] == kEncodedLiteralTag) {
      if (text.size() < 3 || text[2] != kEncodedStringPrefix)
        throw std::invalid_argument(
            "checkpoint canonical-JSON literal token is malformed");
      charge(text.size() - 2);
      return json::value(json::string(text.substr(2)));
    }
    if (text.size() < 3 || text[1] != kEncodedReferenceTag)
      throw std::invalid_argument(
          "checkpoint canonical-JSON identity token is malformed");
    std::size_t index = 0;
    const auto* first = text.data() + 2;
    const auto* last = text.data() + text.size();
    const auto parsed = std::from_chars(first, last, index);
    if (parsed.ec != std::errc{} || parsed.ptr != last ||
        index >= entries_.size() ||
        std::string_view(first, static_cast<std::size_t>(last - first)) !=
            std::to_string(index))
      throw std::invalid_argument(
          "checkpoint canonical-JSON identity reference is invalid");
    if (reference_depth_ >= kMaximumCanonicalJsonReferenceDepth)
      throw std::invalid_argument(
          "checkpoint canonical-JSON identity nesting is too deep");
    ++reference_depth_;
    try {
      const auto& identity = decoded_identity(index);
      --reference_depth_;
      charge(identity.size());
      return json::value(json::string(identity));
    } catch (...) {
      --reference_depth_;
      throw;
    }
  }

  const std::string& decoded_identity(std::size_t index) {
    if (decoded_[index].has_value()) return *decoded_[index];
    if (state_[index] != 0)
      throw std::invalid_argument(
          "checkpoint canonical-JSON identity table contains a cycle");
    state_[index] = 1;
    auto value = decode_value(entries_[index]);
    if (!value.is_object() && !value.is_array())
      throw std::invalid_argument(
          "checkpoint canonical-JSON identity table entry is not composite");
    auto identity = json::serialize(value);
    if (identity.size() < kMinimumInternedStringBytes)
      throw std::invalid_argument(
          "checkpoint canonical-JSON identity table contains a short entry");
    charge(identity.size());
    decoded_[index] = identity;
    state_[index] = 2;
    return *decoded_[index];
  }

  void charge(std::size_t bytes) {
    if (bytes > decoded_budget_ - decoded_work_)
      throw std::invalid_argument(
          "checkpoint canonical-JSON expansion exceeds its decoded-size "
          "budget (consumed=" + std::to_string(decoded_work_) +
          ", next=" + std::to_string(bytes) +
          ", budget=" + std::to_string(decoded_budget_) + ")");
    decoded_work_ += bytes;
  }

  const json::array& entries_;
  std::vector<std::optional<std::string>> decoded_;
  std::vector<unsigned char> state_;
  std::uint64_t decoded_budget_ = 0;
  std::uint64_t decoded_work_ = 0;
  std::size_t reference_depth_ = 0;
};

std::uint64_t decoded_json_budget(std::size_t header_bytes,
                                  std::size_t payload_bytes) {
  auto maximum = kMaximumDecodedJsonBudget;
  if (const auto* diagnostic =
          std::getenv("DIFFEXP_CHECKPOINT_DECODED_JSON_BUDGET_GIB");
      diagnostic != nullptr && *diagnostic != '\0') {
    std::uint64_t gib = 0;
    const auto* first = diagnostic;
    const auto* last = diagnostic + std::strlen(diagnostic);
    const auto parsed = std::from_chars(first, last, gib);
    if (parsed.ec != std::errc{} || parsed.ptr != last || gib < 8 ||
        gib > 32)
      throw std::invalid_argument(
          "DIFFEXP_CHECKPOINT_DECODED_JSON_BUDGET_GIB must be an integer "
          "between 8 and 32");
    maximum = gib * std::uint64_t{1024} * 1024 * 1024;
  }
  const auto encoded_bytes = static_cast<std::uint64_t>(header_bytes) +
      static_cast<std::uint64_t>(payload_bytes);
  const auto scaled = encoded_bytes >
          maximum / kDecodedToEncodedBudgetRatio
      ? maximum
      : encoded_bytes * kDecodedToEncodedBudgetRatio;
  return std::min(maximum,
                  std::max(kMinimumDecodedJsonBudget, scaled));
}

struct RawContainer {
  std::uint32_t schema = 0;
  std::string header_json;
  std::string payload_json;
};

void write_raw_atomic(const std::string& path, std::uint32_t schema,
                      std::string_view header_json,
                      std::string_view payload_json) {
  if (schema != kContainerSchema && schema != kLegacyContainerSchema)
    throw std::invalid_argument("unsupported checkpoint container schema");
  if (path.empty()) throw std::invalid_argument("checkpoint path is empty");
  if (header_json.size() > std::numeric_limits<std::uint32_t>::max())
    throw std::invalid_argument("checkpoint header is too large");
  if (payload_json.size() > kMaximumPayloadBytes)
    throw std::invalid_argument("checkpoint payload exceeds the 8 GiB limit");

  const auto total =
      static_cast<std::uint64_t>(kFixedHeaderBytes) +
      static_cast<std::uint64_t>(header_json.size()) +
      static_cast<std::uint64_t>(payload_json.size());
  if (total < header_json.size() || total < payload_json.size() ||
      total > std::numeric_limits<std::size_t>::max())
    throw std::overflow_error("checkpoint container size overflow");
  std::vector<unsigned char> fixed;
  fixed.reserve(kFixedHeaderBytes);
  fixed.insert(fixed.end(), kMagic.begin(), kMagic.end());
  append_u32(fixed, schema);
  append_u32(fixed, static_cast<std::uint32_t>(header_json.size()));
  append_u64(fixed, static_cast<std::uint64_t>(payload_json.size()));
  append_u32(fixed, checksum(header_json));
  append_u32(fixed, checksum(payload_json));
  if (fixed.size() != kFixedHeaderBytes)
    throw std::logic_error(
        "checkpoint fixed-header serialization changed size");

  const std::filesystem::path target(path);
  const auto directory = target.has_parent_path()
      ? target.parent_path() : std::filesystem::path(".");
  if (!std::filesystem::is_directory(directory))
    throw std::invalid_argument("checkpoint directory does not exist: " +
                                directory.string());
  const auto temporary = temporary_path(path);
  int descriptor = ::open(temporary.c_str(), O_WRONLY | O_CREAT | O_EXCL,
                          S_IRUSR | S_IWUSR);
  if (descriptor < 0) throw system_error("cannot create checkpoint", temporary);
  bool open = true;
  try {
    write_all(descriptor, fixed.data(), fixed.size(), temporary);
    write_all(
        descriptor,
        reinterpret_cast<const unsigned char*>(header_json.data()),
        header_json.size(), temporary);
    write_all(
        descriptor,
        reinterpret_cast<const unsigned char*>(payload_json.data()),
        payload_json.size(), temporary);
    if (::fsync(descriptor) != 0)
      throw system_error("cannot fsync checkpoint", temporary);
    close_checked(descriptor, temporary);
    open = false;
    if (::rename(temporary.c_str(), path.c_str()) != 0)
      throw system_error("cannot atomically replace checkpoint", path);
    const int directory_descriptor = ::open(directory.c_str(), O_RDONLY);
    if (directory_descriptor < 0)
      throw system_error("cannot open checkpoint directory", directory.string());
    if (::fsync(directory_descriptor) != 0) {
      const auto saved_errno = errno;
      ::close(directory_descriptor);
      errno = saved_errno;
      throw system_error("cannot fsync checkpoint directory", directory.string());
    }
    close_checked(directory_descriptor, directory.string());
  } catch (...) {
    if (open) ::close(descriptor);
    ::unlink(temporary.c_str());
    throw;
  }
}

RawContainer read_raw(const std::string& path) {
  if (path.empty()) throw std::invalid_argument("checkpoint path is empty");
  const int descriptor = ::open(path.c_str(), O_RDONLY);
  if (descriptor < 0) throw system_error("cannot open checkpoint", path);
  bool open = true;
  try {
    struct stat metadata {};
    if (::fstat(descriptor, &metadata) != 0)
      throw system_error("cannot stat checkpoint", path);
    if (metadata.st_size < static_cast<off_t>(kFixedHeaderBytes))
      throw std::runtime_error("checkpoint is shorter than its fixed header");
    std::array<unsigned char, kFixedHeaderBytes> fixed{};
    read_all(descriptor, fixed.data(), fixed.size(), path);
    if (!std::equal(kMagic.begin(), kMagic.end(), fixed.begin()))
      throw std::runtime_error("checkpoint magic or container version is invalid");
    const auto schema = read_u32(fixed.data() + 8);
    if (schema != kContainerSchema && schema != kLegacyContainerSchema)
      throw std::runtime_error("unsupported checkpoint container schema " +
                               std::to_string(schema));
    const auto header_size = read_u32(fixed.data() + 12);
    const auto payload_size = read_u64(fixed.data() + 16);
    const auto expected_header_checksum = read_u32(fixed.data() + 24);
    const auto expected_payload_checksum = read_u32(fixed.data() + 28);
    if (payload_size > kMaximumPayloadBytes ||
        payload_size > std::numeric_limits<std::size_t>::max())
      throw std::runtime_error("checkpoint payload length is unsupported");
    const auto expected_size = static_cast<std::uint64_t>(kFixedHeaderBytes) +
        header_size + payload_size;
    if (expected_size != static_cast<std::uint64_t>(metadata.st_size))
      throw std::runtime_error("checkpoint length does not match its header");
    RawContainer output;
    output.schema = schema;
    output.header_json.resize(header_size);
    output.payload_json.resize(static_cast<std::size_t>(payload_size));
    read_all(descriptor,
             reinterpret_cast<unsigned char*>(output.header_json.data()),
             output.header_json.size(), path);
    read_all(descriptor,
             reinterpret_cast<unsigned char*>(output.payload_json.data()),
             output.payload_json.size(), path);
    close_checked(descriptor, path);
    open = false;
    if (checksum(output.header_json) != expected_header_checksum)
      throw std::runtime_error("checkpoint JSON header checksum mismatch");
    if (checksum(output.payload_json) != expected_payload_checksum)
      throw std::runtime_error("checkpoint payload checksum mismatch");
    return output;
  } catch (...) {
    if (open) ::close(descriptor);
    throw;
  }
}

}  // namespace

ExactComplexBallDump dump_complex_ball_exact(const ComplexBall& value) {
  if (!value.is_finite())
    throw std::invalid_argument(
        "cannot checkpoint a non-finite complex ball");
  return {dump_arb_exact(acb_realref(value.raw())),
          dump_arb_exact(acb_imagref(value.raw()))};
}

ComplexBall load_complex_ball_exact(const ExactComplexBallDump& dump) {
  ComplexBall value;
  load_arb_exact(acb_realref(value.raw()), dump.real, "real");
  load_arb_exact(acb_imagref(value.raw()), dump.imaginary, "imaginary");
  if (!value.is_finite())
    throw std::invalid_argument(
        "exact Arb checkpoint decoded a non-finite complex ball");
  return value;
}

void write_atomic(const std::string& path, std::string_view header_json,
                  std::string_view payload_json) {
  write_atomic(path, json::parse(header_json), json::parse(payload_json));
}

void write_atomic(const std::string& path, const json::value& header,
                  const json::value& payload) {
  CanonicalJsonStringDagEncoder encoder;
  auto encoded_header = encoder.encode(header);
  auto encoded_payload = encoder.encode(payload);
  json::object header_envelope{
      {"codec", kCanonicalJsonDagCodec},
      {"root", std::move(encoded_header)}};
  json::object payload_envelope{
      {"codec", kCanonicalJsonDagCodec},
      {"canonical_json_strings", encoder.release_entries()},
      {"root", std::move(encoded_payload)}};
  write_raw_atomic(path, kContainerSchema, json::serialize(header_envelope),
                   json::serialize(payload_envelope));
}

JsonContainer read_json(const std::string& path) {
  auto raw = read_raw(path);
  if (raw.schema == kLegacyContainerSchema)
    return {json::parse(raw.header_json), json::parse(raw.payload_json)};

  const auto header_value = json::parse(raw.header_json);
  const auto payload_value = json::parse(raw.payload_json);
  if (!header_value.is_object() || !payload_value.is_object())
    throw std::invalid_argument(
        "checkpoint canonical-JSON DAG envelope is not an object");
  const auto& header = header_value.as_object();
  const auto& payload = payload_value.as_object();
  require_exact_keys(header, {"codec", "root"},
                     "checkpoint canonical-JSON header envelope");
  require_exact_keys(payload,
                     {"codec", "canonical_json_strings", "root"},
                     "checkpoint canonical-JSON payload envelope");
  if (!header.at("codec").is_string() ||
      !payload.at("codec").is_string() ||
      header.at("codec").as_string() != kCanonicalJsonDagCodec ||
      payload.at("codec").as_string() != kCanonicalJsonDagCodec)
    throw std::invalid_argument(
        "checkpoint canonical-JSON DAG codec identity is invalid");
  if (!payload.at("canonical_json_strings").is_array())
    throw std::invalid_argument(
        "checkpoint canonical-JSON identity table is not an array");
  CanonicalJsonStringDagDecoder decoder(
      payload.at("canonical_json_strings").as_array(),
      decoded_json_budget(raw.header_json.size(), raw.payload_json.size()));
  JsonContainer decoded{decoder.decode(header.at("root")),
                        decoder.decode(payload.at("root"))};
  decoder.require_all_entries_used();
  return decoded;
}

Container read(const std::string& path) {
  auto decoded = read_json(path);
  return {json::serialize(decoded.header), json::serialize(decoded.payload)};
}

}  // namespace diffexp::kernel::checkpoint
