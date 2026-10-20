#pragma once

#include "convert.hpp"
#include <snd/filters.hpp>

namespace eq::filters {

[[nodiscard]] inline
auto freq_param_to_filter_hz(float v) -> float {
	return blink::math::convert::linear_to_filter_hz(v);
}

[[nodiscard]] inline
auto gain_param_to_band_shelf_db(float v) -> snd::filters::gain_db {
	return {convert::raw_to_gain_db(v)};
}

[[nodiscard]] inline
auto gain_param_to_band_pass_db(float) -> snd::filters::gain_db {
	return {30.0f};
}

[[nodiscard]] inline
auto gain_param_to_hi_shelf_db(float v) -> snd::filters::inverted_gain_db {
	return {-convert::raw_to_gain_db(v)};
}

[[nodiscard]] inline
auto gain_param_to_lo_shelf_db(float v) -> snd::filters::gain_db {
	return {convert::raw_to_gain_db(v)};
}

[[nodiscard]] inline
auto slope_param_to_filter_slope(float v) -> float {
	return (1.0f - v) * 2.0f;
}

[[nodiscard]] inline
auto bandwidth_param_to_filter_bw(float v) -> float {
	return std::lerp(0.2f, 4.0f, v);
}

[[nodiscard]] inline
auto make_coeffs_for_band_shelf(float freq, float gain01, float slope, float bw, float SR) -> snd::filters::shelving::band_shelf::coeffs {
	const auto filter_freq  = freq_param_to_filter_hz(freq);
	const auto filter_gain  = gain_param_to_band_shelf_db(gain01);
	const auto filter_slope = slope_param_to_filter_slope(slope);
	const auto filter_bw    = bandwidth_param_to_filter_bw(bw);
	return snd::filters::shelving::band_shelf::make_coeffs(filter_freq, filter_gain, filter_slope, filter_bw, SR);
}

[[nodiscard]] inline
auto make_coeffs_for_band_pass(float freq, float gain01, float slope, float bw, float SR) -> snd::filters::shelving::band_pass::coeffs {
	const auto filter_freq  = freq_param_to_filter_hz(freq);
	const auto filter_gain  = gain_param_to_band_pass_db(gain01);
	const auto filter_slope = slope_param_to_filter_slope(slope);
	const auto filter_bw    = bandwidth_param_to_filter_bw(bw);
	return snd::filters::shelving::band_pass::make_coeffs(filter_freq, filter_gain, filter_slope, filter_bw, SR);
}

[[nodiscard]] inline
auto make_coeffs_for_hi_shelf(float freq, float gain01, float slope, float SR) -> snd::filters::shelving::hi_shelf::coeffs {
	const auto filter_freq  = freq_param_to_filter_hz(freq);
	const auto filter_gain  = gain_param_to_hi_shelf_db(gain01);
	const auto filter_slope = slope_param_to_filter_slope(slope);
	return snd::filters::shelving::hi_shelf::make_coeffs(filter_freq, filter_gain, filter_slope, SR);
}

[[nodiscard]] inline
auto make_coeffs_for_lo_shelf(float freq, float gain01, float slope, float SR) -> snd::filters::shelving::lo_shelf::coeffs {
	const auto filter_freq  = freq_param_to_filter_hz(freq);
	const auto filter_gain  = gain_param_to_lo_shelf_db(gain01);
	const auto filter_slope = slope_param_to_filter_slope(slope);
	return snd::filters::shelving::lo_shelf::make_coeffs(filter_freq, filter_gain, filter_slope, SR);
}

} // eq::filters
