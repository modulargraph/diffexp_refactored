#pragma once

#include "diffexp2/scalar.hpp"

#include <cstdint>
#include <string>
#include <string_view>

namespace diffexp2::checkpoint {

// Bump this identity whenever a change makes a previously emitted payload
// unsafe to reconstruct as the same retained native state.  It is deliberately
// independent of compiler timestamps so checkpoints remain usable after an
// identical source rebuild.
inline constexpr std::string_view kBuildIdentity =
    "diffexp2-persistent-checkpoint-core-v1";
inline constexpr std::uint32_t kContainerSchema = 1;

struct Container {
  std::string header_json;
  std::string payload_json;
};

// Arb's dump/load representation preserves the exact midpoint limbs and
// radius, unlike the deliberately presentation-oriented decimal bridge.  Live
// numeric checkpoint sections must use this codec so restore neither loses a
// bit nor manufactures an unjustified enclosure.
struct ExactComplexBallDump {
  std::string real;
  std::string imaginary;
};

ExactComplexBallDump dump_complex_ball_exact(const ComplexBall& value);
ComplexBall load_complex_ball_exact(const ExactComplexBallDump& dump);

// The container is endian-stable and checksums its JSON header and payload
// independently.  write_atomic performs temporary-file + fsync + rename and
// fsyncs the containing directory before returning.
void write_atomic(const std::string& path, std::string_view header_json,
                  std::string_view payload_json);
Container read(const std::string& path);

}  // namespace diffexp2::checkpoint
