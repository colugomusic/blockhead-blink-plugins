#pragma once

#include "filters.hpp"
#include <blink/plugin-impl.hpp>

namespace eq {

struct Params {
	struct {
		std::array<blink_ParamIdx, 8> band_on;
		std::array<blink_ParamIdx, 8> band_curve;
	} option;
	struct {
		std::array<blink_ParamIdx, 8> band_freq;
		std::array<blink_ParamIdx, 8> band_mag;
		std::array<blink_ParamIdx, 8> band_q;
	} slider;
};

struct UnitDSP {
	blink_SR SR;
	blink::BlockPositions block_positions;
	std::array<std::array<filters::multi_state<float>, 2>, 8> shelf_lo_states;
	std::array<std::array<filters::multi_state<float>, 2>, 8> shelf_hi_states;
	std::array<std::array<filters::multi_state<float>, 2>, 8> pass_lo_states;
	std::array<std::array<filters::multi_state<float>, 2>, 8> pass_hi_states;
	std::array<std::array<filters::multi_state<float>, 2>, 8> bell_states;
};

using Instance = blink::Instance<>;
using Unit     = blink::Unit<UnitDSP>;

struct Model {
	blink::Plugin plugin;
	blink::Entities<Instance, Unit> entities;
	Params params;
};

} // namespace eq
