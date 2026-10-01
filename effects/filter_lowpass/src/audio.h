#pragma once

#include <blink/effect_unit.hpp>
#include <blink/quiet-include.hpp>
#include <DSP/MLDSPFilters.h>
#include <DSP/MLDSPGens.h>

namespace lowpass {

class Plugin;
class Instance;

class Audio : public blink::EffectUnit
{
public:

	Audio(Instance* instance);

	blink_Error process(const blink_EffectBuffer& buffer, const blink_EffectUnitState& unit_state, const float* in, float* out) override;
	void reset() override;

private:
	
	const Plugin* plugin_;
	std::array<ml::Lopass, 2> lopass_;
};

} // lowpass
