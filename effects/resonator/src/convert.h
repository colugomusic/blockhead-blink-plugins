#pragma once

#include <cstdint>
#include <blink/quiet-include.hpp>
BLINK_QUIET_INCLUDE_BEGIN
#include <DSP/MLDSPOps.h>
BLINK_QUIET_INCLUDE_END

namespace resonator {
namespace convert {

[[nodiscard]] inline
auto linear_to_ratio(float x) {
	return std::pow(4.0f, x);
}

[[nodiscard]] inline
auto linear_to_ratio(const ml::DSPVector& x) {
	return ml::pow({4.0f}, x);
}

[[nodiscard]] inline
auto ratio_to_linear(float x) {
	return std::log(x) / std::log(4.0f);
};

} // namespace convert
} // namespace resonator
