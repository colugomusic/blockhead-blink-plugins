#pragma once

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
		std::array<blink_ParamIdx, 8> band_slope;
		std::array<blink_ParamIdx, 8> band_bandwidth;
	} slider;
};

struct UnitDSP {
	blink_SR SR;
	blink::BlockPositions block_positions;
};

using Instance = blink::Instance<>;
using Unit     = blink::Unit<UnitDSP>;

struct Model {
	blink::Plugin plugin;
	blink::Entities<Instance, Unit> entities;
	Params params;
};

} // namespace eq
