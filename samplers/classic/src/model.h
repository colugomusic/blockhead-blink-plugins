#pragma once

#include <blink/plugin-impl.hpp>
#include <blink/transform/tape.hpp>
#include "shared/noise_generator.h"

namespace classic_sampler {

struct Params {
	struct {
		blink_ParamIdx amp;
		blink_ParamIdx noise_amount;
		blink_ParamIdx noise_color;
		blink_ParamIdx pan;
		blink_ParamIdx pitch;
	} env;
	struct {
		blink_ParamIdx loop;
		blink_ParamIdx noise_mode;
		blink_ParamIdx reverse_mode;
		blink_ParamIdx reverse_toggle;
	} option;
	struct {
		blink_ParamIdx samp_offs;
	} slider_int;
	struct {
		blink_ParamIdx amp;
		blink_ParamIdx noise_width;
		blink_ParamIdx pan;
		blink_ParamIdx pitch;
	} slider_real;
};

struct ReverseCorrection {
	struct Grain {
		bool on     = false;
		double beg  = 0.0;
		double pos  = 0.0;
		double ff   = 0.0;
		double vpos = 0.0;
		double vend = 0.0;
		double vff  = 0.0;
	} grain;
	blink::BlockPositions dry_positions;
};

struct DrawState {
	blink::transform::Tape tape_transformer;
};

struct UnitDSP {
	blink_SR SR;
	blink::BlockPositions block_positions;
	blink::transform::Tape tape_transformer;
	NoiseGenerator noise_gen;
	ReverseCorrection reverse_correction;
};

using Instance = blink::Instance<>;
using Unit     = blink::Unit<UnitDSP>;

struct Model {
	blink::Plugin plugin;
	blink::Entities<Instance, Unit> entities;
	DrawState draw;
	Params params;
};

} // namespace classic_sampler

