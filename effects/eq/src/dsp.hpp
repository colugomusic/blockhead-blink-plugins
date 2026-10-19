#pragma once

#include "curve-type.hpp"
#include "model.h"
#include <blink/search.hpp>
#include <DSP/MLDSPFunctional.h>

namespace eq {
namespace dsp {

struct AudioData {
	struct {
		std::array<blink::uniform::Option, 8> band_on;
		std::array<blink::uniform::Option, 8> band_curve;
	} option;
	struct {
		std::array<blink::uniform::SliderReal, 8> band_freq;
		std::array<blink::uniform::SliderReal, 8> band_gain;
		std::array<blink::uniform::SliderReal, 8> band_slope;
		std::array<blink::uniform::SliderReal, 8> band_bandwidth;
	} slider;
};

[[nodiscard]]
auto make_audio_data(const Model& model, const blink_UniformParamData* param_data) -> AudioData {
	auto out = AudioData{};
	for (size_t i = 0; i < BAND_COUNT; i++) {
		out.option.band_on[i]        = blink::make_option_data(model.plugin, param_data, model.params.option.band_on[i]);
		out.option.band_curve[i]     = blink::make_option_data(model.plugin, param_data, model.params.option.band_curve[i]);
		out.slider.band_freq[i]      = blink::make_slider_real_data(model.plugin, param_data, model.params.slider.band_freq[i]);
		out.slider.band_gain[i]      = blink::make_slider_real_data(model.plugin, param_data, model.params.slider.band_gain[i]);
		out.slider.band_slope[i]     = blink::make_slider_real_data(model.plugin, param_data, model.params.slider.band_slope[i]);
		out.slider.band_bandwidth[i] = blink::make_slider_real_data(model.plugin, param_data, model.params.slider.band_bandwidth[i]);
	}
	return out;
}

[[nodiscard]]
auto process(filters::shelving::band_shelf::state* band, ml::DSPVector signal, ml::DSPVector freq, ml::DSPVector gain01, ml::DSPVector slope, ml::DSPVector bw, float SR) -> ml::DSPVector {
	for (int i = 0; i < kFloatsPerDSPVector; i++) {
		auto coeffs = filters::shelving::band_shelf::make_coeffs(freq[i], gain01[i], slope[i], bw[i], SR);
		std::tie(*band, signal[i]) = filters::shelving::band_shelf::process(*band, signal[i], coeffs);
	}
	return signal;
}

[[nodiscard]]
auto process(filters::shelving::hi_shelf::state* band, ml::DSPVector signal, ml::DSPVector freq, ml::DSPVector gain01, ml::DSPVector slope, float SR) -> ml::DSPVector {
	for (int i = 0; i < kFloatsPerDSPVector; i++) {
		auto coeffs = filters::shelving::hi_shelf::make_coeffs(freq[i], gain01[i], slope[i], SR);
		std::tie(*band, signal[i]) = filters::shelving::hi_shelf::process(*band, signal[i], coeffs);
	}
	return signal;
}

[[nodiscard]]
auto process(filters::shelving::lo_shelf::state* band, ml::DSPVector signal, ml::DSPVector freq, ml::DSPVector gain01, ml::DSPVector slope, float SR) -> ml::DSPVector {
	for (int i = 0; i < kFloatsPerDSPVector; i++) {
		auto coeffs = filters::shelving::lo_shelf::make_coeffs(freq[i], gain01[i], slope[i], SR);
		std::tie(*band, signal[i]) = filters::shelving::lo_shelf::process(*band, signal[i], coeffs);
	}
	// auto test_coeffs = filters::shelving::lo_shelf::make_coeffs(freq[0], gain01[0], slope[0], 44100.0f);
	// auto test0 = filters::shelving::lo_shelf::process({}, 0.5f, test_coeffs);
	// auto test1 = filters::shelving::lo_shelf::process({}, 0.6f, test_coeffs);
	// auto test2 = filters::shelving::lo_shelf::process({}, 0.7f, test_coeffs);
	// auto test3 = filters::shelving::lo_shelf::process({}, 0.8f, test_coeffs);
	// auto test4 = filters::shelving::lo_shelf::process({}, 0.9f, test_coeffs);
	return signal;
}

auto process(Model* model, UnitDSP* unit_dsp, const blink_VaryingData& varying, const blink_UniformData& uniform, const float* in, float* out) -> blink_Error {
	unit_dsp->block_positions.add(varying.positions, BLINK_VECTOR_SIZE);
	const auto data = make_audio_data(*model, uniform.param_data);
	const auto SRf  = static_cast<float>(unit_dsp->SR.value);
	auto in_vec = ml::DSPVectorArray<2>{in};
	auto signal = in_vec;
	for (size_t ch = 0; ch < 2; ch++) {
		auto& channel = unit_dsp->channels[ch];
		for (size_t i = 0; i < BAND_COUNT; i++) {
			const auto on = data.option.band_on[i].value > 0;
			if (on) {
				const auto freq  = blink::math::convert::linear_to_filter_hz(data.slider.band_freq[i].value);
				const auto gain  = data.slider.band_gain[i].value;
				const auto slope = data.slider.band_slope[i].value * 2.0f;
				const auto bw    = std::lerp(0.2, 4.0f, data.slider.band_bandwidth[i].value);
				const auto type  = static_cast<curve_type>(data.option.band_curve[i].value);
				switch (type) {
					case curve_type::band_shelf: { signal.row(ch) = process(&channel.bands[i].band_shelf, signal.constRow(ch), freq, gain, slope, bw, SRf); break; }
					case curve_type::hi_shelf:   { signal.row(ch) = process(&channel.bands[i].hi_shelf, signal.constRow(ch), freq, gain, slope, SRf); break; }
					case curve_type::lo_shelf:   { signal.row(ch) = process(&channel.bands[i].lo_shelf, signal.constRow(ch), freq, gain, slope, SRf); break; }
					default:                     { break; }
				}
			}
		}
	}
	ml::storeAligned(signal.constRow(0), out);
	ml::storeAligned(signal.constRow(1), out + kFloatsPerDSPVector);
	return BLINK_OK;
}

auto reset(Model* model, UnitDSP* unit_dsp) -> void {
	for (auto& channel : unit_dsp->channels) {
		channel = {};
	}
}

} // namespace dsp
} // namespace eq
