#pragma once

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
		std::array<blink::uniform::SliderReal, 8> band_mag;
		std::array<blink::uniform::SliderReal, 8> band_q;
	} slider;
};

[[nodiscard]]
auto make_audio_data(const Model& model, const blink_UniformParamData* param_data) -> AudioData {
	auto out = AudioData{};
	return out;
}

auto process(Model* model, UnitDSP* unit_dsp, const blink_VaryingData& varying, const blink_UniformData& uniform, const float* in, float* out) -> blink_Error {
	unit_dsp->block_positions.add(varying.positions, BLINK_VECTOR_SIZE);
	const auto data = make_audio_data(*model, uniform.param_data);
	auto in_vec     = ml::DSPVectorArray<2>{in};
	auto out_vec    = ml::DSPVectorArray<2>{};
	ml::storeAligned(out_vec.constRow(0), out);
	ml::storeAligned(out_vec.constRow(1), out + kFloatsPerDSPVector);
	return BLINK_OK;
}

auto reset(Model* model, UnitDSP* unit_dsp) -> void {
}

} // namespace dsp
} // namespace eq
