#pragma once

#include "band-type.hpp"
#include "effects/eq/src/filters.hpp"
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

auto settle(band_shelf_filter* filter, float freq, float gain01, float slope, float bw, float SR) -> void {
	const auto coeffs  = filters::make_coeffs_for_band_shelf(freq, gain01, slope, bw, SR);
	filter->state      = snd::filters::shelving::band_shelf::settle(filter->state, coeffs);
	filter->settled    = true;
}

auto settle(band_pass_filter* filter, float freq, float gain01, float slope, float bw, float SR) -> void {
	const auto coeffs  = filters::make_coeffs_for_band_pass(freq, gain01, slope, bw, SR);
	filter->state      = snd::filters::shelving::band_pass::settle(filter->state, coeffs);
	filter->settled    = true;
}

auto settle(hi_shelf_filter* filter, float freq, float gain01, float slope, float SR) -> void {
	const auto coeffs  = filters::make_coeffs_for_hi_shelf(freq, gain01, slope, SR);
	filter->state      = snd::filters::shelving::hi_shelf::settle(filter->state, coeffs);
	filter->settled    = true;
}

auto settle(lo_shelf_filter* filter, float freq, float gain01, float slope, float SR) -> void {
	const auto coeffs  = filters::make_coeffs_for_lo_shelf(freq, gain01, slope, SR);
	filter->state      = snd::filters::shelving::lo_shelf::settle(filter->state, coeffs);
	filter->settled    = true;
}

[[nodiscard]]
auto process(band_shelf_filter* filter, ml::DSPVector signal, ml::DSPVector freq, ml::DSPVector gain01, ml::DSPVector slope, ml::DSPVector bw, float SR) -> ml::DSPVector {
	if (!filter->settled) {
		settle(filter, freq[0], gain01[0], slope[0], bw[0], SR);
	}
	for (int i = 0; i < kFloatsPerDSPVector; i++) {
		const auto coeffs = filters::make_coeffs_for_band_shelf(freq[i], gain01[i], slope[i], bw[i], SR);
		std::tie(filter->state, signal[i]) = snd::filters::shelving::band_shelf::process(filter->state, signal[i], coeffs);
	}
	return signal;
}

[[nodiscard]]
auto process(band_pass_filter* filter, ml::DSPVector signal, ml::DSPVector freq, ml::DSPVector gain01, ml::DSPVector slope, ml::DSPVector bw, float SR) -> ml::DSPVector {
	if (!filter->settled) {
		settle(filter, freq[0], gain01[0], slope[0], bw[0], SR);
	}
	for (int i = 0; i < kFloatsPerDSPVector; i++) {
		const auto coeffs = filters::make_coeffs_for_band_pass(freq[i], gain01[i], slope[i], bw[i], SR);
		std::tie(filter->state, signal[i]) = snd::filters::shelving::band_pass::process(filter->state, signal[i], coeffs);
	}
	return signal;
}

[[nodiscard]]
auto process(hi_shelf_filter* filter, ml::DSPVector signal, ml::DSPVector freq, ml::DSPVector gain01, ml::DSPVector slope, float SR) -> ml::DSPVector {
	if (!filter->settled) {
		settle(filter, freq[0], gain01[0], slope[0], SR);
	}
	for (int i = 0; i < kFloatsPerDSPVector; i++) {
		const auto coeffs = filters::make_coeffs_for_hi_shelf(freq[i], gain01[i], slope[i], SR);
		std::tie(filter->state, signal[i]) = snd::filters::shelving::hi_shelf::process(filter->state, signal[i], coeffs);
	}
	return signal;
}

[[nodiscard]]
auto process(lo_shelf_filter* filter, ml::DSPVector signal, ml::DSPVector freq, ml::DSPVector gain01, ml::DSPVector slope, float SR) -> ml::DSPVector {
	if (!filter->settled) {
		settle(filter, freq[0], gain01[0], slope[0], SR);
	}
	for (int i = 0; i < kFloatsPerDSPVector; i++) {
		const auto coeffs = filters::make_coeffs_for_lo_shelf(freq[i], gain01[i], slope[i], SR);
		std::tie(filter->state, signal[i]) = snd::filters::shelving::lo_shelf::process(filter->state, signal[i], coeffs);
	}
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
				const auto type  = static_cast<band_type>(data.option.band_curve[i].value);
				const auto freq  = data.slider.band_freq[i].value;
				const auto gain  = data.slider.band_gain[i].value;
				const auto slope = data.slider.band_slope[i].value;
				const auto bw    = data.slider.band_bandwidth[i].value;
				switch (type) {
					case band_type::band_shelf: { signal.row(ch) = process(&channel.bands[i].band_shelf, signal.constRow(ch), freq, gain, slope, bw, SRf); break; }
					case band_type::band_pass:  { signal.row(ch) = process(&channel.bands[i].band_pass, signal.constRow(ch), freq, gain, slope, bw, SRf); break; }
					case band_type::hi_shelf:   { signal.row(ch) = process(&channel.bands[i].hi_shelf, signal.constRow(ch), freq, gain, slope, SRf); break; }
					case band_type::lo_shelf:   { signal.row(ch) = process(&channel.bands[i].lo_shelf, signal.constRow(ch), freq, gain, slope, SRf); break; }
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
