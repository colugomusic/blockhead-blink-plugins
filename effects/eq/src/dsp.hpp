#pragma once

#include "model.h"
#include <blink/search.hpp>
#include <DSP/MLDSPFunctional.h>
#include <algorithm>
#include <cmath>
#include <tuple>

namespace eq {
namespace dsp {

constexpr auto BAND_COUNT = 8u;
constexpr auto MAGNITUDE_MAX = 30.0f;

enum struct curve_type {
	shelf_lo,
	shelf_hi,
	pass_lo,
	pass_hi,
	bell,
};

struct AudioData {
	struct {
		std::array<blink::uniform::Option, BAND_COUNT> band_on;
		std::array<blink::uniform::Option, BAND_COUNT> band_curve;
	} option;
	struct {
		std::array<blink::uniform::SliderReal, BAND_COUNT> band_freq;
		std::array<blink::uniform::SliderReal, BAND_COUNT> band_mag;
		std::array<blink::uniform::SliderReal, BAND_COUNT> band_q;
	} slider;
};

[[nodiscard]]
auto linear_to_magnitude_db(float v) -> float {
	return v * MAGNITUDE_MAX;
}

[[nodiscard]]
auto make_audio_data(const Model& model, const blink_UniformParamData* param_data) -> AudioData {
	auto out = AudioData{};
	for (size_t i = 0; i < BAND_COUNT; ++i) {
		out.option.band_on[i]    = blink::make_option_data(model.plugin, param_data, model.params.option.band_on[i]);
		out.option.band_curve[i] = blink::make_option_data(model.plugin, param_data, model.params.option.band_curve[i]);
		out.slider.band_freq[i]  = blink::make_slider_real_data(model.plugin, param_data, model.params.slider.band_freq[i]);
		out.slider.band_mag[i]   = blink::make_slider_real_data(model.plugin, param_data, model.params.slider.band_mag[i]);
		out.slider.band_q[i]     = blink::make_slider_real_data(model.plugin, param_data, model.params.slider.band_q[i]);
	}
	return out;
}

auto process(Model* model, UnitDSP* unit_dsp, const blink_VaryingData& varying, const blink_UniformData& uniform, const float* in, float* out) -> blink_Error {
	unit_dsp->block_positions.add(varying.positions, BLINK_VECTOR_SIZE);
	const auto data = make_audio_data(*model, uniform.param_data);
	static const auto F_MAX = blink::math::convert::linear_to_filter_hz(1.0f);
	auto display_to_omega = [F_MAX](float display_x) -> float {
		return std::clamp(blink::math::convert::linear_to_filter_hz(display_x) / F_MAX, 0.0001f, 0.9999f);
	};

	for (size_t s = 0; s < kFloatsPerDSPVector; ++s) {
		auto left  = in[s];
		auto right = in[s + kFloatsPerDSPVector];
		for (size_t i = 0; i < BAND_COUNT; ++i) {
			if (data.option.band_on[i].value <= 0) {
				continue;
			}

			const auto omega = display_to_omega(data.slider.band_freq[i].value);
			const auto q = data.slider.band_q[i].value;
			const auto k = std::lerp(0.1f, 1.0f, q);
			const auto magnitude = data.slider.band_mag[i].value;
			const auto A = blink::math::convert::db_to_linear(linear_to_magnitude_db(magnitude) / 2.0f);

			switch (static_cast<curve_type>(data.option.band_curve[i].value)) {
				case curve_type::shelf_lo: {
					const auto coeffs = filters::shelf_lo::make_coeffs<float>(omega, k, A);
					std::tie(unit_dsp->shelf_lo_states[i][0], left) = filters::shelf_lo::process(unit_dsp->shelf_lo_states[i][0], left, coeffs);
					std::tie(unit_dsp->shelf_lo_states[i][1], right) = filters::shelf_lo::process(unit_dsp->shelf_lo_states[i][1], right, coeffs);
					break;
				}
				case curve_type::shelf_hi: {
					const auto coeffs = filters::shelf_hi::make_coeffs<float>(omega, k, A);
					std::tie(unit_dsp->shelf_hi_states[i][0], left) = filters::shelf_hi::process(unit_dsp->shelf_hi_states[i][0], left, coeffs);
					std::tie(unit_dsp->shelf_hi_states[i][1], right) = filters::shelf_hi::process(unit_dsp->shelf_hi_states[i][1], right, coeffs);
					break;
				}
				case curve_type::pass_lo: {
					const auto coeffs = filters::pass_lo::make_coeffs<float>(omega, k);
					std::tie(unit_dsp->pass_lo_states[i][0], left) = filters::pass_lo::process(unit_dsp->pass_lo_states[i][0], left, coeffs);
					std::tie(unit_dsp->pass_lo_states[i][1], right) = filters::pass_lo::process(unit_dsp->pass_lo_states[i][1], right, coeffs);
					break;
				}
				case curve_type::pass_hi: {
					const auto coeffs = filters::pass_hi::make_coeffs<float>(omega, k);
					std::tie(unit_dsp->pass_hi_states[i][0], left) = filters::pass_hi::process(unit_dsp->pass_hi_states[i][0], left, coeffs);
					std::tie(unit_dsp->pass_hi_states[i][1], right) = filters::pass_hi::process(unit_dsp->pass_hi_states[i][1], right, coeffs);
					break;
				}
				case curve_type::bell:
				default: {
					const auto coeffs = filters::bell::make_coeffs<float>(omega, k, A, 8);
					std::tie(unit_dsp->bell_states[i][0], left) = filters::bell::process(unit_dsp->bell_states[i][0], left, coeffs);
					std::tie(unit_dsp->bell_states[i][1], right) = filters::bell::process(unit_dsp->bell_states[i][1], right, coeffs);
					break;
				}
			}
		}
		out[s] = left;
		out[s + kFloatsPerDSPVector] = right;
	}

	return BLINK_OK;
}

auto reset(Model*, UnitDSP* unit_dsp) -> void {
	unit_dsp->shelf_lo_states = {};
	unit_dsp->shelf_hi_states = {};
	unit_dsp->pass_lo_states  = {};
	unit_dsp->pass_hi_states  = {};
	unit_dsp->bell_states     = {};
}

} // namespace dsp
} // namespace eq
