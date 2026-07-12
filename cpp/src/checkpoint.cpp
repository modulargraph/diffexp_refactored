#include "diffexp2/checkpoint.hpp"

#include <boost/crc.hpp>

#include <array>
#include <atomic>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

namespace diffexp2::checkpoint {
namespace {

constexpr std::array<unsigned char, 8> kMagic{
    'D', 'E', '2', 'C', 'P', '0', '0', '1'};
constexpr std::size_t kFixedHeaderBytes = 32;
constexpr std::uint64_t kMaximumPayloadBytes =
    std::uint64_t{8} * 1024 * 1024 * 1024;

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
    const auto amount = ::write(descriptor, data + written, size - written);
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
    const auto amount = ::read(descriptor, data + read_bytes,
                               size - read_bytes);
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

}  // namespace

void write_atomic(const std::string& path, std::string_view header_json,
                  std::string_view payload_json) {
  if (path.empty()) throw std::invalid_argument("checkpoint path is empty");
  if (header_json.size() > std::numeric_limits<std::uint32_t>::max())
    throw std::invalid_argument("checkpoint header is too large");
  if (payload_json.size() > kMaximumPayloadBytes)
    throw std::invalid_argument("checkpoint payload exceeds the 8 GiB limit");

  std::vector<unsigned char> bytes;
  const auto total = kFixedHeaderBytes + header_json.size() +
                     payload_json.size();
  if (total < header_json.size() || total < payload_json.size())
    throw std::overflow_error("checkpoint container size overflow");
  bytes.reserve(total);
  bytes.insert(bytes.end(), kMagic.begin(), kMagic.end());
  append_u32(bytes, kContainerSchema);
  append_u32(bytes, static_cast<std::uint32_t>(header_json.size()));
  append_u64(bytes, static_cast<std::uint64_t>(payload_json.size()));
  append_u32(bytes, checksum(header_json));
  append_u32(bytes, checksum(payload_json));
  bytes.insert(bytes.end(), header_json.begin(), header_json.end());
  bytes.insert(bytes.end(), payload_json.begin(), payload_json.end());

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
    write_all(descriptor, bytes.data(), bytes.size(), temporary);
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

Container read(const std::string& path) {
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
    if (schema != kContainerSchema)
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
    Container output;
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

}  // namespace diffexp2::checkpoint
