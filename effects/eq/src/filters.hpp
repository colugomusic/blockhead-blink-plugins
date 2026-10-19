#pragma once

#include "blink/math.hpp"
#include <cassert>
#include <complex>
#include <span>
#include <numbers>
#include <snd/misc.hpp>

namespace eq::filters {

constexpr auto GAIN_MAX = 12.0f;

[[nodiscard]]
auto reciprocal(float v) -> float {
	return 1.0f / v;
}

[[nodiscard]]
auto bi_clip(float v, float max) -> std::pair<float,float> {
	if (v > max) { return {max, v}; }
	else         { return {v, max}; }
}

[[nodiscard]] auto gain_db_to_linear(float v) { return v / GAIN_MAX; }
[[nodiscard]] auto linear_to_gain_db(float v) { return v * GAIN_MAX; }

} // eq::filters

namespace eq::filters::shelving {

struct svf_result {
	std::complex<float> LP;
	std::complex<float> BPn;
	std::complex<float> HP;
};

struct zdf_signal {
	std::complex<float> state;
	float signal;
};

struct integrator_state {
	float fbk = 0.0f;
};

[[nodiscard]]
auto g2R(float g, float M) -> float {
	const auto a = ((2.0f - g) / g) * M;
	return sqrt(a + 2.0f) / 2.0f;
}

[[nodiscard]]
auto R2p1(float _2R, float M) -> float {
	const auto a = std::max(1.0f, _2R);
	const auto b = std::sqrt(std::pow(a, 2.0f) - 1.0f);
	const auto c = a + b;
	if (M >= 1.0f) { return c; }
	else           { return reciprocal(c); }
}

[[nodiscard]]
auto a22R(float v) -> float {
	return cos(v) * 2.0f;
}

[[nodiscard]]
auto approx_tan_rational(float v) -> float {
	auto a = (v * -0.0896638f) + 0.0388452f;
	auto b = (v * -0.430871f) + 0.0404318f;
	auto c = (a * v) + 1.00005f;
	auto d = (b * v) + 1.0f;
	return (c * v) / d;
}

[[nodiscard]]
auto compute_wdiv2(float f, float SR) -> float {
	f = std::max(f, 0.01f);
	f *= std::numbers::pi_v<float> / SR;
	f = std::min(f, 1.50845f);
	return approx_tan_rational(f);
}

[[nodiscard]]
auto integrator(integrator_state state, zdf_signal zdf, float wdiv2) -> std::tuple<integrator_state, zdf_signal> {
	const auto a = zdf.signal * wdiv2;
	const auto r = wdiv2 * zdf.state.real();
	const auto i = state.fbk + (zdf.state.imag() * wdiv2);
	const auto out = a + state.fbk;
	state.fbk = a + out;
	snd::flush_denormal_to_zero(state.fbk);
	zdf = {
		.state = {r, i},
		.signal = out,
	};
	return std::make_tuple(state, zdf);
}

[[nodiscard]]
auto svf_input_fbk_solve(float x, std::complex<float> zdf) -> std::tuple<float, float> {
	const auto a = x * zdf.real();
	const auto b = a + zdf.imag();
	const auto c = 1.0f - zdf.real();
	const auto d = 1.0f / c;
	return std::make_tuple(b * d, d);
}

[[nodiscard]]
auto svf_input_fbk(zdf_signal in) -> zdf_signal {
	const auto [y, _] = svf_input_fbk_solve(in.signal, in.state * -1.0f);
	return {
		.state  = {1.0f, 0.0f},
		.signal = in.signal + y
	};
}

[[nodiscard]]
auto gain(zdf_signal zdf, float v) -> zdf_signal {
	zdf.signal *= v;
	zdf.state  *= v;
	return zdf;
}

} // eq::shelving

namespace eq::filters::shelving::tf {

[[nodiscard]]
auto svf(std::complex<float> s, float _2R) -> svf_result {
	const auto spow2 = s * s;
	const auto a     = spow2 + 1.0f;
	const auto b     = s * _2R;
	const auto c     = a + b;
	const auto LP    = 1.0f / c;
	const auto BPn   = (LP * s) * _2R;
	const auto HP    = (LP * s) * s;
	return {
		.LP  = LP,
		.BPn = BPn,
		.HP  = HP
	};
}

} // eq::filters::shelving::tf

namespace eq::filters::shelving::band_shelf {

struct lp_to_bp_state {
	integrator_state integrator1;
	integrator_state integrator2;
	std::complex<float> fbk = {0.0f, 0.0f};
};

struct svf_state {
	lp_to_bp_state lp_to_bp_1;
	lp_to_bp_state lp_to_bp_2;
	std::complex<float> fbk = {0.0f, 0.0f};
};

struct tilt_state {
	svf_state svf;
};

struct state {
	state() = default;
	tilt_state tilt1;
	tilt_state tilt2;
};

struct tilt_coeffs {
	float damp;
	float M;
	float recip_M;
	float BW;
	float wMid;
};

struct coeffs {
	float F;
	float recip_Mpow2;
	float wdiv2;
	tilt_coeffs tilt1;
	tilt_coeffs tilt2;
};

[[nodiscard]]
auto make_coeffs(float f, float gain01, float slope, float bw, float SR) -> coeffs {
	const auto gain_db              = linear_to_gain_db(gain01);
	const auto slope_mod            = std::pow(2.0f, slope);
	const auto M                    = blink::math::convert::db_to_linear(gain_db / -2.0f);
	const auto Mpow2                = std::pow(M, 2.0f);
	const auto recip_Mpow2          = reciprocal(Mpow2);
	const auto [slope_lo, slope_hi] = bi_clip(slope_mod, 2.0f);
	const auto slope_hi_half        = slope_hi / 2.0f;
	const auto g2R_M                = Mpow2 + recip_Mpow2;
	const auto g2R1                 = g2R(slope_lo, g2R_M);
	const auto g2R2                 = g2R(slope_hi_half, g2R_M);
	const auto g2R1_clamp           = std::min(1.0f, g2R1);
	const auto g2R2_clamp           = std::min(1.0f, g2R2);
	const auto R1x                  = g2R1 / g2R1_clamp;
	const auto a1                   = std::acos(g2R1_clamp);
	const auto a2                   = std::acos(g2R2_clamp);
	const auto g2R2_R2p1            = R2p1(g2R2, M);
	const auto tilt1_M              = M * g2R2_R2p1;
	const auto tilt2_M              = M / g2R2_R2p1;
	const auto a2div2               = a2 / 2.0f;
	const auto tilt1_damp           = a22R(a1 + a2div2) * R1x;
	const auto tilt2_damp           = a22R(a1 - a2div2);
	const auto tilt1_wMid           = std::sqrt(tilt1_M);
	const auto tilt2_wMid           = std::sqrt(tilt2_M);
	const auto BWpow2               = std::pow(2.0f, bw);
	const auto tilt1_recip_M        = reciprocal(tilt1_M);
	const auto tilt2_recip_M        = reciprocal(tilt2_M);
	const auto tilt1_BW             = (BWpow2 - reciprocal(BWpow2)) * tilt1_wMid;
	const auto tilt2_BW             = (BWpow2 - reciprocal(BWpow2)) * tilt2_wMid;
	const auto wdiv2                = compute_wdiv2(f, SR);
	return {
		.F             = f,
		.recip_Mpow2   = recip_Mpow2,
		.wdiv2         = wdiv2,
		.tilt1 = {
			tilt1_damp,
			tilt1_M,
			tilt1_recip_M,
			tilt1_BW,
			tilt1_wMid,
		},
		.tilt2 = {
			tilt2_damp,
			tilt2_M,
			tilt2_recip_M,
			tilt2_BW,
			tilt2_wMid,
		},
	};
}

[[nodiscard]]
auto zdf_output(zdf_signal zdf, std::complex<float> ext) -> zdf_signal {
	zdf.state = {
		zdf.state.real() * ext.real(),
		zdf.state.imag() + (zdf.state.real() * ext.imag())
	};
	return zdf;
}

[[nodiscard]]
auto lp_to_bp_input_fbk(zdf_signal in, std::complex<float> fbk) -> std::tuple<zdf_signal, std::complex<float>> {
	fbk *= -1.0f;
	const auto [y, G] = svf_input_fbk_solve(in.signal, fbk);
	const auto out    = zdf_signal{{1.0f, 0.0f}, in.signal + y};
	const auto ext    = std::complex<float>{G + in.state.real(), (fbk.imag() + in.state.imag()) * G};
	return std::make_tuple(out, ext);
}

[[nodiscard]]
auto lp_to_bp(lp_to_bp_state state, zdf_signal in, float wdiv2, float BW) -> std::tuple<lp_to_bp_state, zdf_signal> {
	zdf_signal a, b, c;
	auto ext                       = std::complex<float>{};
	in                             = gain(in, BW);
	std::tie(a, ext)               = lp_to_bp_input_fbk(in, state.fbk);
	std::tie(state.integrator1, b) = integrator(state.integrator1, a, wdiv2);
	std::tie(state.integrator2, c) = integrator(state.integrator2, b, wdiv2);
	state.fbk                      = c.state;
	return std::make_tuple(state, zdf_output(b, ext));
}

[[nodiscard]]
auto svf(svf_state state, float in, float wdiv2, float damp, float BW) -> std::tuple<svf_state, float, float, float> {
	zdf_signal HP, BP, LP;
	HP                             = svf_input_fbk({state.fbk, in});
	std::tie(state.lp_to_bp_1, BP) = lp_to_bp(state.lp_to_bp_1, HP, wdiv2, BW);
	std::tie(state.lp_to_bp_2, LP) = lp_to_bp(state.lp_to_bp_2, BP, wdiv2, BW);
	const auto BPn = gain(BP, damp);
	state.fbk = LP.state + BPn.state;
	return std::make_tuple(state, LP.signal, HP.signal, BPn.signal);
}

[[nodiscard]]
auto tilt(tilt_state state, float in, float wdiv2, tilt_coeffs coeffs) -> std::tuple<tilt_state, float> {
	float LP, HP, BPn;
	std::tie(state.svf, LP, HP, BPn) = svf(state.svf, in, wdiv2, coeffs.damp, coeffs.BW);
	return std::make_tuple(state, (LP / coeffs.M) + (HP * coeffs.M) + BPn);
}

[[nodiscard]]
auto process(band_shelf::state state, float in, const band_shelf::coeffs& coeffs) -> std::tuple<band_shelf::state, float> {
	auto signal = in;
	std::tie(state.tilt1, signal) = tilt(state.tilt1, signal, coeffs.wdiv2, coeffs.tilt1);
	std::tie(state.tilt2, signal) = tilt(state.tilt2, signal, coeffs.wdiv2, coeffs.tilt2);
	signal *= coeffs.recip_Mpow2;
	return std::make_tuple(state, signal);
}

} // eq::filters::shelving::band_shelf

namespace eq::filters::shelving::band_shelf::tf {

[[nodiscard]]
auto convert(std::complex<float> v, float _2R) -> std::complex<float> {
	static constexpr auto EPSILON = 0.00000001f;
	const auto a = 1.0f / (v + EPSILON);
	const auto b = v + a;
	return b * (1.0f / _2R);
}

[[nodiscard]]
auto svf(float f, float F, float _2R, float BW) -> svf_result {
	return shelving::tf::svf(convert({0.0f, f / F}, BW), _2R);
}

[[nodiscard]]
auto tilt(std::complex<float> t, float f, float F, tilt_coeffs coeffs) -> std::complex<float> {
	const auto [LP, BPn, HP] = svf(f, F, coeffs.damp, coeffs.BW);
	const auto a = LP * coeffs.recip_M;
	const auto b = HP * coeffs.M;
	const auto c = a + b;
	const auto d = c + BPn;
	return t * d;
}

[[nodiscard]]
auto compute(const band_shelf::coeffs& coeffs, float f) -> float {
	auto s = std::complex<float>{1.0f, 0.0f};
	s = tilt(s, f, coeffs.F, coeffs.tilt1);
	s = tilt(s, f, coeffs.F, coeffs.tilt2);
	s *= coeffs.recip_Mpow2;
	return abs(s);
}

auto compute(const band_shelf::coeffs& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert (in_omega.size() == out_magnitude.size());
	auto fn_compute = [&coeffs](float omega) -> float { return compute(coeffs, omega); };
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(), fn_compute);
}

} // eq::filters::shelving::band_shelf::tf

namespace eq::filters::shelving::hilo {

struct tilt_coeffs {
	float damp;
	float M;
	float recip_M;
	float wMid;
};

struct coeffs {
	float F;
	float Mpow2;
	float recip_Mpow2;
    float wdiv2;
	tilt_coeffs tilt1;
	tilt_coeffs tilt2;
	tilt_coeffs tilt3;
	tilt_coeffs tilt4;
};

struct svf_state {
	integrator_state integrator1;
	integrator_state integrator2;
	std::complex<float> fbk = {0.0f, 0.0f};
};

struct tilt_state {
	svf_state svf;
};

struct state {
	state() = default;
	tilt_state tilt1;
	tilt_state tilt2;
	tilt_state tilt3;
	tilt_state tilt4;
};

[[nodiscard]]
auto make_coeffs(float f, float gain_db, float slope, float SR) -> hilo::coeffs {
	const auto slope_mod              = std::pow(2.0f, slope);
	const auto M                      = blink::math::convert::db_to_linear(gain_db / -2.0f);
	const auto Mpow2                  = std::pow(M, 2.0f);
	const auto recip_Mpow2            = reciprocal(Mpow2);
	const auto [slope1_lo, slope1_hi] = bi_clip(slope_mod, 2.0f);
	const auto slope1_hi_half         = slope1_hi / 2.0f;
	const auto [slope2_lo, slope2_hi] = bi_clip(slope1_hi_half, 2.0f);
	const auto slope2_hi_half         = slope2_hi / 2.0f;
	const auto g2R_M                  = Mpow2 + recip_Mpow2;
	const auto g2R1                   = g2R(slope1_lo, g2R_M);
	const auto g2R2                   = g2R(slope2_lo, g2R_M);
	const auto g2R3                   = g2R(slope2_hi_half, g2R_M);
	const auto g2R1_clamp             = std::min(1.0f, g2R1);
	const auto g2R2_clamp             = std::min(1.0f, g2R2);
	const auto g2R3_clamp             = std::min(1.0f, g2R3);
	const auto R1x                    = g2R1 / g2R1_clamp;
	const auto a1                     = std::acos(g2R1_clamp);
	const auto a2                     = std::acos(g2R2_clamp);
	const auto a3                     = std::acos(g2R3_clamp);
	const auto g2R2_R2p1              = R2p1(g2R2, M);
	const auto g2R3_R2p1              = R2p1(g2R3, M);
	const auto M1                     = std::sqrt(M * g2R3_R2p1);
	const auto tilt34_M               = std::sqrt(M / g2R3_R2p1);
	const auto tilt1_M                = M1 * g2R2_R2p1;
	const auto tilt2_M                = M1 / g2R2_R2p1;
	const auto a2div2                 = a2 / 2.0f;
	const auto a3div4                 = a3 / 4.0f;
	const auto tilt1_damp             = a22R(a1 + a2div2) * R1x;
	const auto tilt2_damp             = a22R(a1 - a2div2);
	const auto tilt3_damp             = a22R((a1 + a2div2) - a3div4);
	const auto tilt4_damp             = a22R((a1 - a2div2) - a3div4);
	const auto tilt1_recip_M          = reciprocal(tilt1_M);
	const auto tilt2_recip_M          = reciprocal(tilt2_M);
	const auto tilt34_recip_M         = reciprocal(tilt34_M);
	const auto tilt1_wMid             = std::sqrt(tilt1_M);
	const auto tilt2_wMid             = std::sqrt(tilt2_M);
	const auto tilt34_wMid            = std::sqrt(tilt34_M);
	const auto wdiv2                  = compute_wdiv2(f, SR);
	return {
		.F           = f,
		.Mpow2       = Mpow2,
		.recip_Mpow2 = recip_Mpow2,
		.wdiv2       = wdiv2,
		.tilt1 = {
			tilt1_damp,
			tilt1_M,
			tilt1_recip_M,
			tilt1_wMid,
		},
		.tilt2 = {
			tilt2_damp,
			tilt2_M,
			tilt2_recip_M,
			tilt2_wMid,
		},
		.tilt3 = {
			tilt3_damp,
			tilt34_M,
			tilt34_recip_M,
			tilt34_wMid,
		},
		.tilt4 = {
			tilt4_damp,
			tilt34_M,
			tilt34_recip_M,
			tilt34_wMid,
		},
	};
}

[[nodiscard]]
auto svf(svf_state state, float in, float wdiv2, float damp) -> std::tuple<svf_state, float, float, float> {
	zdf_signal HP, BP, LP;
	HP                              = svf_input_fbk({state.fbk, in});
	std::tie(state.integrator1, BP) = integrator(state.integrator1, HP, wdiv2);
	std::tie(state.integrator2, LP) = integrator(state.integrator2, BP, wdiv2);
	const auto BPn = gain(BP, damp);
	state.fbk = LP.state + BPn.state;
	return std::make_tuple(state, LP.signal, HP.signal, BPn.signal);
}

[[nodiscard]]
auto tilt(tilt_state state, float in, float wdiv2, tilt_coeffs coeffs) -> std::tuple<tilt_state, float> {
	wdiv2 *= std::sqrt(coeffs.M);
	float LP, HP, BPn;
	std::tie(state.svf, LP, HP, BPn) = svf(state.svf, in, wdiv2, coeffs.damp);
	return std::make_tuple(state, (LP / coeffs.M) + (HP * coeffs.M) + BPn);
}

[[nodiscard]]
auto process(hilo::state state, float in, const hilo::coeffs& coeffs) -> std::tuple<hilo::state, float> {
	auto signal = in;
	std::tie(state.tilt1, signal) = tilt(state.tilt1, signal, coeffs.wdiv2, coeffs.tilt1);
	std::tie(state.tilt2, signal) = tilt(state.tilt2, signal, coeffs.wdiv2, coeffs.tilt2);
	std::tie(state.tilt3, signal) = tilt(state.tilt3, signal, coeffs.wdiv2, coeffs.tilt3);
	std::tie(state.tilt4, signal) = tilt(state.tilt4, signal, coeffs.wdiv2, coeffs.tilt4);
	return std::make_tuple(state, signal);
}

} // eq::filters::shelving::hilo

namespace eq::filters::shelving::hilo::tf {

[[nodiscard]]
auto svf(float f, float F, float _2R) -> svf_result {
	return shelving::tf::svf({0.0f, f / F}, _2R);
}

[[nodiscard]]
auto tilt(std::complex<float> t, float f, float F, tilt_coeffs coeffs) -> std::complex<float> {
	const auto [LP, BPn, HP] = svf(f, F * coeffs.wMid, coeffs.damp);
	const auto a = LP * coeffs.recip_M;
	const auto b = HP * coeffs.M;
	const auto c = a + b;
	const auto d = c + BPn;
	return t * d;
}

[[nodiscard]]
auto compute(const hilo::coeffs& coeffs, float f) -> std::complex<float> {
	auto s = std::complex<float>{1.0f, 0.0f};
	s = hilo::tf::tilt(s, f, coeffs.F, coeffs.tilt1);
	s = hilo::tf::tilt(s, f, coeffs.F, coeffs.tilt2);
	s = hilo::tf::tilt(s, f, coeffs.F, coeffs.tilt3);
	s = hilo::tf::tilt(s, f, coeffs.F, coeffs.tilt4);
	return s;
}

} // eq::filters::shelving::hilo::tf

namespace eq::filters::shelving::hi_shelf {

struct coeffs { hilo::coeffs hilo; };
struct state { hilo::state hilo; };

[[nodiscard]]
auto make_coeffs(float f, float gain01, float slope, float SR) -> coeffs {
	return {.hilo = hilo::make_coeffs(f, linear_to_gain_db(0.0f-gain01), slope, SR)};
}

[[nodiscard]]
auto process(hi_shelf::state state, float in, const hi_shelf::coeffs& coeffs) -> std::tuple<hi_shelf::state, float> {
	auto signal = in;
	std::tie(state.hilo, signal) = hilo::process(state.hilo, signal, coeffs.hilo);
	signal *= coeffs.hilo.Mpow2;
	return std::make_tuple(state, signal);
}

} // eq::filters::shelving::hi_shelf

namespace eq::filters::shelving::hi_shelf::tf {

[[nodiscard]]
auto compute(const hi_shelf::coeffs& coeffs, float f) -> float {
	return abs(hilo::tf::compute(coeffs.hilo, f)) * coeffs.hilo.Mpow2;
}

auto compute(const hi_shelf::coeffs& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert (in_omega.size() == out_magnitude.size());
	auto fn_compute = [&coeffs](float omega) -> float { return compute(coeffs, omega); };
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(), fn_compute);
}

} // eq::filters::shelving::hi_shelf::tf

namespace eq::filters::shelving::lo_shelf {

struct coeffs { hilo::coeffs hilo; };
struct state { hilo::state hilo; };

[[nodiscard]]
auto make_coeffs(float f, float gain01, float slope, float SR) -> coeffs {
	return {.hilo = hilo::make_coeffs(f, linear_to_gain_db(gain01), slope, SR)};
}

[[nodiscard]]
auto process(lo_shelf::state state, float in, const lo_shelf::coeffs& coeffs) -> std::tuple<lo_shelf::state, float> {
	auto signal = in;
	std::tie(state.hilo, signal) = hilo::process(state.hilo, signal, coeffs.hilo);
	signal *= coeffs.hilo.recip_Mpow2;
	return std::make_tuple(state, signal);
}

} // eq::filters::shelving::lo_shelf

namespace eq::filters::shelving::lo_shelf::tf {

[[nodiscard]]
auto compute(const lo_shelf::coeffs& coeffs, float f) -> float {
	return abs(hilo::tf::compute(coeffs.hilo, f)) * coeffs.hilo.recip_Mpow2;
}

auto compute(const lo_shelf::coeffs& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert (in_omega.size() == out_magnitude.size());
	auto fn_compute = [&coeffs](float omega) -> float { return compute(coeffs, omega); };
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(), fn_compute);
}

} // eq::filters::shelving::lo_shelf::tf

