#pragma once

#include <blink/plugin_impl.hpp>
BLINK_QUIET_INCLUDE_BEGIN
#include <DSP/MLDSPFilters.h>
BLINK_QUIET_INCLUDE_END

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
	std::array<ml::Bandpass, 2> filter;
};

using Instance = blink::Instance<>;
using Unit     = blink::Unit<UnitDSP>;

struct Model {
	blink::Plugin plugin;
	blink::Entities<Instance, Unit> entities;
	Params params;
};
