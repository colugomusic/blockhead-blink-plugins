#pragma once

#include <cstdint>
#include <blink/const_math.hpp>
#include <blink/quiet-include.hpp>
#include <tweak/const-math.hpp>
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

[[nodiscard]] constexpr inline
auto ratio_to_linear(float x) {
	return ::tweak::const_math::log(x) / ::tweak::const_math::log(4.0f);
};

} // namespace convert
} // namespace resonator
