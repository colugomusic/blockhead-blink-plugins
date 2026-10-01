#pragma once

#include <blink/plugin-impl.hpp>
#include <DSP/MLDSPFilters.h>

namespace filter_highpass {

struct Params {
	struct {
		blink_ParamIdx frequency;
		blink_ParamIdx resonance;
		blink_ParamIdx mix;
	} env;
};

struct UnitDSP {
	blink_SR SR;
	blink::BlockPositions block_positions;
	std::array<ml::Hipass, 2> filter;
};

using Instance = blink::Instance<>;
using Unit     = blink::Unit<UnitDSP>;

struct Model {
	blink::Plugin plugin;
	blink::Entities<Instance, Unit> entities;
	Params params;
};

} // namespace filter_highpass
