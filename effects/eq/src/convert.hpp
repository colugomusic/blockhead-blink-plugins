#pragma once

#include <blink/math.hpp>

namespace eq::convert {

constexpr auto GAIN_MAX = 30.0f;

[[nodiscard]] inline auto gain_db_to_raw(float v) { return blink::math::convert::bi_to_uni(v / GAIN_MAX); }
[[nodiscard]] inline auto raw_to_gain_db(float v) { return blink::math::convert::uni_to_bi(v) * GAIN_MAX; }

} // eq::convert
