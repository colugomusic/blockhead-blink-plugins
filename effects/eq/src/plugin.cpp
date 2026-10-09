#define BLINK_EXPORT

#ifndef _USE_MATH_DEFINES
#	define _USE_MATH_DEFINES
#endif

#include "dsp.hpp"
#include "filters.hpp"
#include "model.h"
#include <blink/tweak.hpp>
#include <blink_std.h>
#include <cmrc/cmrc.hpp>
#include <plugin-impl.hpp>
#include <ranges>

using namespace eq;

namespace eq {
namespace { // -----------------------------------------------------------------------------------------------

Model model;

constexpr auto PLUGIN_UUID = blink_UUID{"6fccd4a1-8da9-45c3-8174-8b5d7f8cf846"};
constexpr auto BAND_COUNT  = 8;

constexpr auto BAND_ON_UUID = std::array<blink_UUID, BAND_COUNT>{
	blink_UUID{"6ecdb955-1a41-4c44-83a9-5fe2b2c27bb4"},
	blink_UUID{"e28cb6db-29e4-4cbf-9623-5f1bd8f1c6fb"},
	blink_UUID{"b85e6352-9f57-46bb-bf2f-aaf829288aaa"},
	blink_UUID{"f2b5907e-40cc-4045-9e2a-7a58bbea3ece"},
	blink_UUID{"70357d43-2519-4cf8-a822-eece8af2e9ce"},
	blink_UUID{"2ccccf85-5546-45dc-a44d-f7626ecd98ef"},
	blink_UUID{"9587cfda-2fef-4dba-8704-b13ae4b571c4"},
	blink_UUID{"5dd691f3-243b-4420-8712-a7dd232831e1"},
};

constexpr auto BAND_FREQ_UUID = std::array<blink_UUID, BAND_COUNT>{
	blink_UUID{"77de5482-d145-4034-bb17-b8a93393a444"},
	blink_UUID{"a6c4025e-2ea0-45e7-8731-041f4cee9036"},
	blink_UUID{"0b72ef70-74c7-4314-9373-470adfbe15f0"},
	blink_UUID{"c7750543-95b0-40f6-8f7f-802d771bb367"},
	blink_UUID{"f4e5e0eb-18a0-4161-90c0-3804f979235a"},
	blink_UUID{"bdae68c4-e3a4-4ec1-a296-c7e35eb34be1"},
	blink_UUID{"02d79ea4-f438-47e5-8307-9fcf18d90ebb"},
	blink_UUID{"a1972b5d-b817-4a9a-820a-aea65c3b7584"},
};

constexpr auto BAND_MAGNITUDE_UUID = std::array<blink_UUID, BAND_COUNT>{
	blink_UUID{"9e926878-e152-46df-8677-ceecb8a490d3"},
	blink_UUID{"b7b1b6e4-88d7-4541-b291-5597650e35eb"},
	blink_UUID{"c402588e-d9a6-4104-aaac-f95e4518dfa8"},
	blink_UUID{"1398abd2-0128-49ab-b8fa-583f0ddda531"},
	blink_UUID{"365e453c-1cfd-4db2-bfee-372d3d4b8fba"},
	blink_UUID{"adc935fa-9f25-4dc6-97d2-4095d2da75fc"},
	blink_UUID{"338c8d86-48f5-4595-b801-cd13b3f79f67"},
	blink_UUID{"065ce3c5-1e39-4edb-9f4b-7828a4813897"},
};

constexpr auto BAND_Q_UUID = std::array<blink_UUID, BAND_COUNT>{
	blink_UUID{"f27589ef-6b0b-4d05-ae32-691d3d19430a"},
	blink_UUID{"24b26e1a-cc37-4a02-8b83-e28e103b433a"},
	blink_UUID{"be20bbd1-3d19-48f5-a732-14499a2aff38"},
	blink_UUID{"4871e651-69eb-4a39-8871-79a461377943"},
	blink_UUID{"d17a1a19-ee8d-48da-8ced-02f91137c42e"},
	blink_UUID{"54e99e42-efd2-447c-abf7-f55b7bc32c6a"},
	blink_UUID{"a703f46f-3d3e-4f72-a499-ae5f99f3b48d"},
	blink_UUID{"0aad6c3e-ffcd-4cdc-a85a-b76c820ccfe8"},
};

constexpr auto BAND_CURVE_UUID = std::array<blink_UUID, BAND_COUNT>{
	blink_UUID{"3a603ae9-db0c-4975-96d3-efad14a75162"},
	blink_UUID{"84135c5b-f7a4-47de-a3e1-847dab52e2a8"},
	blink_UUID{"72cfb8f0-3cd1-4617-a3bb-2e3d2f2266bf"},
	blink_UUID{"95129991-801b-4e81-b510-e21548a4399b"},
	blink_UUID{"89276999-337c-4b33-944b-b55f1f5414a1"},
	blink_UUID{"ab06f298-db34-47e4-b415-93aa80acd3ab"},
	blink_UUID{"5befd6de-f833-468d-a2ef-25d01ac98f58"},
	blink_UUID{"afe67063-4f7b-4570-a42b-cc8b18ec6243"},
};

enum struct curve_type {
	shelf_lo,
	shelf_hi,
	pass_lo,
	pass_hi,
	bell,
};

struct band_spec {
	std::string_view on_param_name;
	std::string_view freq_param_name;
	std::string_view mag_param_name;
	std::string_view q_param_name;
	std::string_view curve_param_name;
	bool enabled = false;
	float frequency = 1000.0f;
	float q         = 0.0f;
	float magnitude = 0.0f;
	curve_type curve = curve_type::bell;
};

static const auto BAND_SPECS = std::array<band_spec, BAND_COUNT>{
	band_spec{
		.on_param_name    = "Band 1 Enabled",
		.freq_param_name  = "Band 1 Frequency",
		.mag_param_name   = "Band 1 Magnitude",
		.q_param_name     = "Band 1 Q",
		.curve_param_name = "Band 1 Curve Type",
		.enabled          = false,
		.frequency        = blink::math::convert::filter_hz_to_linear(100.0f),
		.magnitude        = 0.0f,
		.curve            = curve_type::bell
	},
	band_spec{
		.on_param_name    = "Band 2 Enabled",
		.freq_param_name  = "Band 2 Frequency",
		.mag_param_name   = "Band 2 Magnitude",
		.q_param_name     = "Band 2 Q",
		.curve_param_name = "Band 2 Curve Type",
		.enabled          = false,
		.frequency        = blink::math::convert::filter_hz_to_linear(600.0f),
		.q                = 1.0f,
		.magnitude        = 0.5f,
		.curve            = curve_type::bell
	},
	band_spec{
		.on_param_name    = "Band 3 Enabled",
		.freq_param_name  = "Band 3 Frequency",
		.mag_param_name   = "Band 3 Magnitude",
		.q_param_name     = "Band 3 Q",
		.curve_param_name = "Band 3 Curve Type",
		.enabled          = true,
		.frequency        = blink::math::convert::filter_hz_to_linear(1000.0f),
		.q                = 0.5f,
		.magnitude        = 0.5f,
		.curve            = curve_type::bell
	},
	band_spec{
		.on_param_name    = "Band 4 Enabled",
		.freq_param_name  = "Band 4 Frequency",
		.mag_param_name   = "Band 4 Magnitude",
		.q_param_name     = "Band 4 Q",
		.curve_param_name = "Band 4 Curve Type",
		.enabled          = false,
		.frequency        = blink::math::convert::filter_hz_to_linear(2000.0f),
		.magnitude        = -1.0f,
		.curve            = curve_type::bell
	},
	band_spec{
		.on_param_name    = "Band 5 Enabled",
		.freq_param_name  = "Band 5 Frequency",
		.mag_param_name   = "Band 5 Magnitude",
		.q_param_name     = "Band 5 Q",
		.curve_param_name = "Band 5 Curve Type",
		.enabled          = false,
		.frequency        = blink::math::convert::filter_hz_to_linear(4000.0f),
		.magnitude        = -0.5f,
		.curve            = curve_type::bell
	},
	band_spec{
		.on_param_name    = "Band 6 Enabled",
		.freq_param_name  = "Band 6 Frequency",
		.mag_param_name   = "Band 6 Magnitude",
		.q_param_name     = "Band 6 Q",
		.curve_param_name = "Band 6 Curve Type",
	},
	band_spec{
		.on_param_name    = "Band 7 Enabled",
		.freq_param_name  = "Band 7 Frequency",
		.mag_param_name   = "Band 7 Magnitude",
		.q_param_name     = "Band 7 Q",
		.curve_param_name = "Band 7 Curve Type",
	},
	band_spec{
		.on_param_name    = "Band 8 Enabled",
		.freq_param_name  = "Band 8 Frequency",
		.mag_param_name   = "Band 8 Magnitude",
		.q_param_name     = "Band 8 Q",
		.curve_param_name = "Band 8 Curve Type",
	},
};

struct audio_data {
	struct {
		std::array<blink::uniform::Option, BAND_COUNT> band_on;
		std::array<blink::uniform::Option, BAND_COUNT> band_curve;
	} option;
	struct {
		std::array<blink::uniform::SliderReal, BAND_COUNT> band_freq;
		std::array<blink::uniform::SliderReal, BAND_COUNT> band_mag;
		std::array<blink::uniform::SliderReal, BAND_COUNT> band_q;
	} slider;
};

[[nodiscard]]
auto make_audio_data(const Model& model, const blink_UniformParamData* param_data) -> audio_data {
	auto out = audio_data{};
	for (size_t i = 0; i < BAND_COUNT; i++) {
		out.option.band_on[i]    = blink::make_option_data(model.plugin, param_data, model.params.option.band_on[i]);
		out.option.band_curve[i] = blink::make_option_data(model.plugin, param_data, model.params.option.band_curve[i]);
		out.slider.band_freq[i]  = blink::make_slider_real_data(model.plugin, param_data, model.params.slider.band_freq[i]);
		out.slider.band_mag[i]   = blink::make_slider_real_data(model.plugin, param_data, model.params.slider.band_mag[i]);
		out.slider.band_q[i]     = blink::make_slider_real_data(model.plugin, param_data, model.params.slider.band_q[i]);
	}
	return out;
}

constexpr auto MAGNITUDE_DEFAULT = 0.0f;
constexpr auto MAGNITUDE_MAX     = 30.0f;

auto magnitude_db_to_linear(float v) -> float                   { return v / MAGNITUDE_MAX; }
auto linear_to_magnitude_db(float v) -> float                   { return v * MAGNITUDE_MAX; }
auto magnitude_stepify(float v) -> float                        { return magnitude_db_to_linear(tweak::math::stepify<100>(linear_to_magnitude_db(v))); }
auto magnitude_constrain(float v) -> float                      { return std::clamp(v, -1.0f, 1.0f); }
auto magnitude_increment(float v, bool precise) -> float        { return magnitude_db_to_linear(tweak::increment<1, 10>(linear_to_magnitude_db(v), precise)); };
auto magnitude_decrement(float v, bool precise) -> float        { return magnitude_db_to_linear(tweak::decrement<1, 10>(linear_to_magnitude_db(v), precise)); };
auto magnitude_drag(float v, int amount, bool precise) -> float { return magnitude_db_to_linear(tweak::drag<float, 1, 10>(linear_to_magnitude_db(v), amount / 5, precise)); };

auto magnitude_from_string(const char* str, float* out) -> blink_Bool {
	if (const auto mag = tweak::find_number<float>(str)) {
		*out = magnitude_db_to_linear(*mag);
		return {true};
	}
	return {false};
}

auto magnitude_to_string(float v, char buffer[BLINK_STRING_MAX]) -> void {
	auto ss = std::stringstream{};
	ss << linear_to_magnitude_db(v) << " dB";
	blink::tweak::write_string(ss.str(), buffer);
}

auto magnitude_tweaker() -> blink_TweakerReal {
	auto out        = blink_TweakerReal{0};
	out.constrain   = magnitude_constrain;
	out.decrement   = magnitude_decrement;
	out.drag        = magnitude_drag;
	out.from_string = magnitude_from_string;
	out.increment   = magnitude_increment;
	out.stepify     = magnitude_stepify;
	out.to_string   = magnitude_to_string;
	return out;
}

auto add_magnitude_slider(const blink::Plugin& plugin, float default_value) -> blink_SliderRealIdx {
	const auto idx = blink::add::slider::empty_real(plugin.host);
	plugin.host.write_slider_real_default_value(plugin.host.usr, idx, default_value);
	plugin.host.write_slider_real_tweaker(plugin.host.usr, idx, magnitude_tweaker());
	return idx;
}

auto add_band_on_option_param(const blink::Plugin& plugin, blink_UUID uuid, const band_spec& spec) -> blink_ParamIdx {
	const auto param_idx  = blink::add::param::option(plugin, uuid);
	const auto flags      = blink_ParamFlags_IsToggle | blink_ParamFlags_MovesDisplay;
	blink::write::param::name(plugin, param_idx, {spec.on_param_name.data()});
	blink::write::param::option_default_value(plugin, param_idx, spec.enabled ? 1 : 0);
	blink::write::param::add_flags(plugin, param_idx, flags);
	return param_idx;
}

auto add_band_freq_slider_param(const blink::Plugin& plugin, blink_UUID uuid, const band_spec& spec) -> blink_ParamIdx {
	const auto param_idx  = blink::add::param::slider_real(plugin, uuid);
	const auto slider_idx = blink::add::slider::filter_frequency(plugin.host, spec.frequency);
	const auto flags      = blink_ParamFlags_MovesDisplay;
	blink::write::param::slider(plugin, param_idx, slider_idx);
	blink::write::param::add_flags(plugin, param_idx, flags);
	blink::write::param::name(plugin, param_idx, {spec.freq_param_name.data()});
	return param_idx;
}

auto add_band_mag_slider_param(const blink::Plugin& plugin, blink_UUID uuid, const band_spec& spec) -> blink_ParamIdx {
	const auto param_idx  = blink::add::param::slider_real(plugin, uuid);
    const auto slider_idx = add_magnitude_slider(plugin, spec.magnitude);
	const auto flags      = blink_ParamFlags_MovesDisplay;
    blink::write::param::slider(plugin, param_idx, slider_idx);
	blink::write::param::name(plugin, param_idx, {spec.mag_param_name.data()});
	blink::write::param::add_flags(plugin, param_idx, flags);
	return param_idx;
}

auto add_band_q_slider_param(const blink::Plugin& plugin, blink_UUID uuid, const band_spec& spec) -> blink_ParamIdx {
	// Create this by starting with the default filter resonaance slider parameter.
	// The UUID is overwritten with ours.
	const auto param_idx  = blink::add::param::slider_real(plugin, {BLINK_STD_UUID_FILTER_RESONANCE});
	const auto slider_idx = blink::read::slider_real(plugin, param_idx);
	const auto flags      = blink_ParamFlags_MovesDisplay;
	blink::write::slider::default_value(plugin, slider_idx, spec.q);
	blink::write::param::uuid(plugin, param_idx, uuid);
	blink::write::param::name(plugin, param_idx, {spec.q_param_name.data()});
	blink::write::param::add_flags(plugin, param_idx, flags);
	return param_idx;
}

auto add_band_curve_option_param(const blink::Plugin& plugin, blink_UUID uuid, const band_spec& spec) -> blink_ParamIdx {
	const auto param_idx = blink::add::param::option(plugin, uuid);
	const auto flags     = blink_ParamFlags_MovesDisplay;
	blink::write::param::name(plugin, param_idx, {spec.curve_param_name.data()});
	blink::write::param::add_flags(plugin, param_idx, flags);
	blink::write::param::strings(plugin, param_idx, {{"Shelf Low", "Shelf High", "Cut Low", "Cut High", "Bell"}});
	blink::write::param::option_default_value(plugin, param_idx, static_cast<int64_t>(spec.curve));
	return param_idx;
}

auto add_band_params(const blink::Plugin& plugin, auto add_fn) -> std::array<blink_ParamIdx, BAND_COUNT> {
	auto arr = std::array<blink_ParamIdx, BAND_COUNT>{};
	for (const auto i : std::views::iota(0, 8)) {
		arr[i] = add_fn(plugin, i);
	}
	return arr;
}

auto add_band_on_params(const blink::Plugin& plugin) -> std::array<blink_ParamIdx, BAND_COUNT> {
	auto fn_add = [](const blink::Plugin& plugin, size_t band_index) {
		return add_band_on_option_param(plugin, BAND_ON_UUID[band_index], BAND_SPECS[band_index]);
	};
	return add_band_params(plugin, fn_add);
}

auto add_band_freq_params(const blink::Plugin& plugin) -> std::array<blink_ParamIdx, BAND_COUNT> {
	auto fn_add = [](const blink::Plugin& plugin, size_t band_index) {
		return add_band_freq_slider_param(plugin, BAND_FREQ_UUID[band_index], BAND_SPECS[band_index]);
	};
	return add_band_params(plugin, fn_add);
}

auto add_band_mag_params(const blink::Plugin& plugin) -> std::array<blink_ParamIdx, BAND_COUNT> {
	auto fn_add = [](const blink::Plugin& plugin, size_t band_index) {
		return add_band_mag_slider_param(plugin, BAND_MAGNITUDE_UUID[band_index], BAND_SPECS[band_index]);
	};
	return add_band_params(plugin, fn_add);
}

auto add_band_q_params(const blink::Plugin& plugin) -> std::array<blink_ParamIdx, BAND_COUNT> {
	auto fn_add = [](const blink::Plugin& plugin, size_t band_index) {
		return add_band_q_slider_param(plugin, BAND_Q_UUID[band_index], BAND_SPECS[band_index]);
	};
	return add_band_params(plugin, fn_add);
}

auto add_band_curve_params(const blink::Plugin& plugin) -> std::array<blink_ParamIdx, BAND_COUNT> {
	auto fn_add = [](const blink::Plugin& plugin, size_t band_index) {
		return add_band_curve_option_param(plugin, BAND_CURVE_UUID[band_index], BAND_SPECS[band_index]);
	};
	return add_band_params(plugin, fn_add);
}

[[nodiscard]]
auto frequency_response(const eq::audio_data& audio_data, float x) -> float {
	// Display coords (filter_hz_to_linear) are log-pitch, not proportional to Hz.
	// Convert to actual Hz then normalise by the display's top frequency (~20 kHz)
	// so that filter bandwidths/slopes scale correctly in octaves. No sample rate needed.
	static const auto F_MAX   = blink::math::convert::linear_to_filter_hz(1.0f);
	auto display_to_omega = [F_MAX](float display_x) -> float {
		return std::clamp(blink::math::convert::linear_to_filter_hz(display_x) / F_MAX, 0.0001f, 0.9999f);
	};
	const auto omega_x = display_to_omega(x);
	auto y = float{1.0f};
	for (size_t i = 0; i < BAND_COUNT; i++) {
		if (audio_data.option.band_on[i].value > 0) {
			const auto omega     = display_to_omega(audio_data.slider.band_freq[i].value);
			const auto magnitude = audio_data.slider.band_mag[i].value;
			const auto q         = audio_data.slider.band_q[i].value;
			const auto k         = std::lerp(0.1f, 1.0f, q);
			const auto A         = blink::math::convert::db_to_linear(linear_to_magnitude_db(magnitude) / 2.0f);
			switch (static_cast<curve_type>(audio_data.option.band_curve[i].value)) {
				case curve_type::shelf_lo: {
					const auto coeffs = filters::shelf_lo::make_coeffs<float>(omega, k, A, 8);
					y *= filters::shelf_lo::transfer(coeffs, omega_x);
					break;
				}
				case curve_type::shelf_hi: {
					const auto coeffs = filters::shelf_hi::make_coeffs<float>(omega, k, A, 8);
					y *= filters::shelf_hi::transfer(coeffs, omega_x);
					break;
				}
				case curve_type::pass_lo: {
					const auto coeffs = filters::pass_lo::make_coeffs<float>(omega, k, 8);
					y *= filters::pass_lo::transfer(coeffs, omega_x);
					break;
				}
				case curve_type::pass_hi: {
					const auto coeffs = filters::pass_hi::make_coeffs<float>(omega, k, 8);
					y *= filters::pass_hi::transfer(coeffs, omega_x);
					break;
				}
				case curve_type::bell:
				default: {
					const auto coeffs = filters::bell::make_coeffs<float>(omega, k, A, 8);
					y *= filters::bell::transfer(coeffs, omega_x);
					break;
				}
			}
		}
	}
	return blink::math::convert::linear_to_db(y) / MAGNITUDE_MAX;
}

} // ---------------------------------------------------------------------------------------------------------
} // eq

auto blink_get_error_string(blink_Error error) -> blink_TempString {
	return {blink::get_std_error_string(static_cast<blink_StdError>(error))};
}

auto blink_effect_get_info(blink_InstanceIdx) -> blink_EffectInstanceInfo {
	return {-1, -1, -1, -1};
}

auto blink_get_plugin_info() -> blink_PluginInfo {
	blink_PluginInfo out = {0};
	out.uuid     = PLUGIN_UUID;
	out.name     = {"EQ"};
	out.category = {BLINK_STD_CATEGORY_FILTERS};
	out.version  = {PLUGIN_VERSION};
	out.has_icon = {true};
	return out;
}

auto blink_init(blink_PluginIdx plugin_idx, blink_HostFns host) -> blink_Error {
	blink::init(&model.plugin, plugin_idx, host);
	model.params.option.band_on    = add_band_on_params(model.plugin);
	model.params.slider.band_freq  = add_band_freq_params(model.plugin);
	model.params.slider.band_mag   = add_band_mag_params(model.plugin);
	model.params.slider.band_q     = add_band_q_params(model.plugin);
	model.params.option.band_curve = add_band_curve_params(model.plugin);
	auto fr_info                   = blink_FrequencyResponseInfo{};
	fr_info.band_count             = BAND_COUNT;
	fr_info.extra_count            = 0;
	fr_info.enabled                = model.params.option.band_on.data();
	fr_info.frequency              = model.params.slider.band_freq.data();
	fr_info.magnitude              = model.params.slider.band_mag.data();
	fr_info.mb_right_horizontal    = model.params.option.band_curve.data();
	fr_info.mb_right_vertical      = model.params.slider.band_q.data();
	blink::add::frequency_response(model.plugin, fr_info);
	return BLINK_OK;
}

auto blink_instance_destroy(blink_InstanceIdx instance_idx) -> blink_Error {
	return blink::destroy_instance(&model.entities, instance_idx);
}

auto blink_instance_make() -> blink_InstanceIdx {
	return blink::make_instance(&model.entities);
}

auto blink_instance_reset(blink_InstanceIdx) -> blink_Error {
	return BLINK_OK;
}

auto blink_instance_stream_init(blink_InstanceIdx, blink_SR) -> blink_Error {
	return BLINK_OK;
}

auto blink_effect_process(blink_UnitIdx unit_idx, const blink_VaryingData* varying, const blink_UniformData* uniform, const float* in, float* out) -> blink_Error {
	auto& unit_dsp = model.entities.unit.get<UnitDSP>(unit_idx.value);
	return dsp::process(&model, &unit_dsp, *varying, *uniform, in, out);
}

auto blink_terminate() -> blink_Error {
	return blink::terminate(&model.entities);
}

auto blink_unit_add(blink_InstanceIdx instance_idx) -> blink_UnitIdx {
	return blink::add_unit(&model.entities, instance_idx);
}

auto blink_unit_reset(blink_UnitIdx unit_idx) -> blink_Error {
	auto& unit_dsp = model.entities.unit.get<UnitDSP>(unit_idx.value);
	dsp::reset(&model, &unit_dsp);
	return BLINK_OK;
}

auto blink_unit_stream_init(blink_UnitIdx unit_idx, blink_SR SR) -> blink_Error {
	auto& unit_dsp = model.entities.unit.get<UnitDSP>(unit_idx.value);
	unit_dsp.SR = SR;
	dsp::reset(&model, &unit_dsp);
	return BLINK_OK;
}

auto blink_frequency_response(const blink_UniformParamData* param_data, blink_FrequencyResponseIdx, blink_FrameCount n, const float* in_x_01, float* out_y_01) -> blink_Error {
	const auto audio_data = make_audio_data(model, param_data);
	for (uint64_t i = 0; i < n.value; i++) {
		out_y_01[i] = frequency_response(audio_data, in_x_01[i]);
	}
	return BLINK_OK;
}

CMRC_DECLARE(plugin);

auto blink_get_resource_data(const char* path) -> blink_ResourceData {
	return blink::get_resource_data(&model.plugin, cmrc::plugin::get_filesystem(), path);
}
