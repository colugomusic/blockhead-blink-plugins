#pragma once

#include <blink/math.hpp>

namespace eq::convert {

constexpr auto GAIN_MAX = 30.0f;

[[nodiscard]] inline auto gain_db_to_linear(float v) { return v / GAIN_MAX; }
[[nodiscard]] inline auto linear_to_gain_db(float v) { return v * GAIN_MAX; }

} // eq::convert
