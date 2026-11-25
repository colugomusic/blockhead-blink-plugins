#pragma once

#include "model.h"
#include <MLDSPOps.h>
#include <blink/search.hpp>
#include <block_positions.hpp>

namespace spring_reverb {
namespace dsp {

static constexpr auto DECAY_LO = 0.01f;
static constexpr auto DECAY_HI = 20.f;
static const auto UNITY_TO_DECAY = ml::Projection(ml::projections::intervalMap({0, 1}, {DECAY_LO, DECAY_HI}, ml::projections::log({DECAY_LO, DECAY_HI})));

struct AudioData {
	struct {
		blink::uniform::Env size;
		blink::uniform::Env decay;
		blink::uniform::Env mix;
	} env;
};

[[nodiscard]]
auto make_audio_data(const Model& model, const blink_UniformParamData* param_data) -> AudioData {
	AudioData out;
	out.env.decay = blink::make_env_data(model.plugin, param_data, model.params.env.decay);
	out.env.size  = blink::make_env_data(model.plugin, param_data, model.params.env.size);
	out.env.mix   = blink::make_env_data(model.plugin, param_data, model.params.env.mix);
	return out;
}

[[nodiscard]]
auto make_input_ramp(int frame_count, float SR) -> ml::DSPVector {
	const auto ramp_samples = 0.01f * SR;
	const auto ramp_begin   = ml::DSPVector{0.0f};
	const auto ramp_end     = ml::DSPVector{ramp_samples};
	const auto frame_ramp   = ml::columnIndex() + frame_count;
	return ml::clamp(ml::inverseLerp(ramp_begin, ramp_end, frame_ramp), ml::DSPVector{0.0f}, ml::DSPVector{1.0f});
}

// This is adapted from code taken from madronalib. Thank you Randy!
auto process(Model* model, UnitDSP* unit_dsp, const blink_VaryingData& varying, const blink_UniformData& uniform, const float* in, float* out) -> blink_Error {
	unit_dsp->block_positions.add(varying.positions, BLINK_VECTOR_SIZE);
	const auto data = make_audio_data(*model, uniform.param_data);
	unit_dsp->glide_size.setGlideTimeInSamples(0.1f * unit_dsp->SR.value);
	unit_dsp->glide_decay.setGlideTimeInSamples(0.1f * unit_dsp->SR.value);
	const auto size  = blink::search::one(data.env.size, unit_dsp->block_positions);
	const auto decay = blink::search::one(data.env.decay, unit_dsp->block_positions);
	const auto mix   = blink::search::vec(data.env.mix, unit_dsp->block_positions);
	const float sr = float(unit_dsp->SR.value);
	const float RT60const = 0.001f;
	float sizeU = ml::lerp(0.01f, 1.0f, size);
	float decayU = decay;
	// generate delay and feedback scalars
	float decayTime = UNITY_TO_DECAY(decayU);
	float decayIterations = decayTime / (sizeU * 0.5f);
	float feedback = (decayU < 1.0f) ? powf(RT60const, 1.0f / decayIterations) : 1.0f;
	if (unit_dsp->buffer_count++ < 1) {
		unit_dsp->glide_size.setValue(sizeU * 2.0f);
		unit_dsp->glide_decay.setValue(feedback);
	}
	// generate smoothed delay time and feedback gain vectors
	ml::DSPVector vSmoothDelay = unit_dsp->glide_size(sizeU * 2.0f);
	ml::DSPVector vSmoothFeedback = unit_dsp->glide_decay(feedback);
	// get the minimum possible delay in samples, which is the length of a DSPVector.
	ml::DSPVector vMin(kFloatsPerDSPVector);
	// get smoothed allpass times in samples
	ml::DSPVector delayParamInSamples = sr * vSmoothDelay;
	ml::DSPVector vt1 = max(0.00476f * delayParamInSamples, vMin);
	ml::DSPVector vt2 = max(0.00358f * delayParamInSamples, vMin);
	ml::DSPVector vt3 = max(0.00973f * delayParamInSamples, vMin);
	ml::DSPVector vt4 = max(0.00830f * delayParamInSamples, vMin);
	ml::DSPVector vt5 = max(0.029f * delayParamInSamples, vMin);
	ml::DSPVector vt6 = max(0.021f * delayParamInSamples, vMin);
	ml::DSPVector vt7 = max(0.078f * delayParamInSamples, vMin);
	ml::DSPVector vt8 = max(0.090f * delayParamInSamples, vMin);
	ml::DSPVector vt9 = max(0.111f * delayParamInSamples, vMin);
	ml::DSPVector vt10 = max(0.096f * delayParamInSamples, vMin);
	const auto in_vec    = ml::DSPVectorArray<2>{in};
	const auto in_ramp   = make_input_ramp(unit_dsp->frame_count, sr);
	const auto in_ramped = in_vec * ml::repeatRows<2>(in_ramp);
	const auto final_mix = mix * in_ramp;
	// sum stereo inputs and diffuse with four allpass filters in series
	ml::DSPVector monoInput = (in_ramped.constRow(0) + in_ramped.constRow(1));
	ml::DSPVector diffusedInput = unit_dsp->ap4(unit_dsp->ap3(unit_dsp->ap2(unit_dsp->ap1(monoInput, vt1), vt2), vt3), vt4);
	// get delay times in samples, subtracting the constant delay of one DSPVector and clamping to zero
	ml::DSPVector vDelayTimeL = max(0.0313f * delayParamInSamples - vMin, ml::DSPVector(0.f));
	ml::DSPVector vDelayTimeR = max(0.0371f * delayParamInSamples - vMin, ml::DSPVector(0.f));
	// sum diffused input with feedback, and apply late diffusion of two more allpass filters to each channel
	ml::DSPVector vTapL = unit_dsp->ap7(unit_dsp->ap5(diffusedInput + unit_dsp->delayL(unit_dsp->feedbackL, vDelayTimeL), vt5), vt7);
	ml::DSPVector vTapR = unit_dsp->ap8(unit_dsp->ap6(diffusedInput + unit_dsp->delayR(unit_dsp->feedbackR, vDelayTimeR), vt6), vt8);
	// apply final allpass filter and gain, and store the feedback
	unit_dsp->feedbackR = unit_dsp->ap9(vTapL, vt9) * vSmoothFeedback;
	unit_dsp->feedbackL = unit_dsp->ap10(vTapR, vt10) * vSmoothFeedback;
	// append the left and right taps and return the stereo output
	const auto out_vec = ml::lerp(in_vec, ml::concatRows(vTapL, vTapR), ml::repeatRows<2>(final_mix));
	ml::storeAligned(out_vec.constRow(0), out);
	ml::storeAligned(out_vec.constRow(1), out + kFloatsPerDSPVector);
	unit_dsp->frame_count += kFloatsPerDSPVector;
	return BLINK_OK;
}

auto init(Model* model, UnitDSP* unit_dsp) -> void {
	unit_dsp->ap1.mGain = 0.75f;
	unit_dsp->ap2.mGain = 0.70f;
	unit_dsp->ap3.mGain = 0.625f;
	unit_dsp->ap4.mGain = 0.625f;
	unit_dsp->ap5.mGain = unit_dsp->ap6.mGain = 0.7f;
	unit_dsp->ap7.mGain = unit_dsp->ap8.mGain = 0.6f;
	unit_dsp->ap9.mGain = unit_dsp->ap10.mGain = 0.5f;
	unit_dsp->ap1.setMaxDelayInSamples(500.f);
	unit_dsp->ap2.setMaxDelayInSamples(500.f);
	unit_dsp->ap3.setMaxDelayInSamples(1000.f);
	unit_dsp->ap4.setMaxDelayInSamples(1000.f);
	unit_dsp->ap5.setMaxDelayInSamples(2600.f);
	unit_dsp->ap6.setMaxDelayInSamples(2600.f);
	unit_dsp->ap7.setMaxDelayInSamples(8000.f);
	unit_dsp->ap8.setMaxDelayInSamples(8000.f);
	unit_dsp->ap9.setMaxDelayInSamples(10000.f);
	unit_dsp->ap10.setMaxDelayInSamples(10000.f);
	unit_dsp->delayL.setMaxDelayInSamples(3500.f);
	unit_dsp->delayR.setMaxDelayInSamples(3500.f);
}

auto reset(Model* model, UnitDSP* unit_dsp) -> void {
	unit_dsp->ap1.clear();
	unit_dsp->ap2.clear();
	unit_dsp->ap3.clear();
	unit_dsp->ap4.clear();
	unit_dsp->ap5.clear();
	unit_dsp->ap6.clear();
	unit_dsp->ap7.clear();
	unit_dsp->ap8.clear();
	unit_dsp->ap9.clear();
	unit_dsp->ap10.clear();
	unit_dsp->delayL.clear();
	unit_dsp->delayR.clear();
	unit_dsp->feedbackL = ml::DSPVector();
	unit_dsp->feedbackR = ml::DSPVector();
	unit_dsp->frame_count = 0;
}

} // namespace dsp
} // namespace spring_reverb
