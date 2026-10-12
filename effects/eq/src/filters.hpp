#pragma once

#include "blink/math.hpp"
#include <cassert>
#include <complex>
#include <DSP/MLDSPOps.h>
#include <span>
#include <algorithm>
#include <numbers>

namespace eq::filters {

auto acos(float v) -> float         { return std::acos(v); }
auto cos(float v) -> float          { return std::cos(v); }
auto pow(float v, float n) -> float { return std::pow(v, n); }
auto sqrt(float v) -> float         { return std::sqrt(v); }
auto min(float a, float b) -> float { return std::min(a, b); }
auto max(float a, float b) -> float { return std::max(a, b); }

template <typename T> [[nodiscard]]
auto reciprocal(T v) -> T {
	return T(1) / v;
}

template <typename T> [[nodiscard]]
auto bi_clip(T v, T max) -> std::pair<T,T> {
	if (v > max)  { return {max, v}; }
	else          { return {v, max}; }
}

} // eq::filters

namespace eq::filters::shelving {

struct svf_result {
	std::complex<float> LP;
	std::complex<float> BPn;
	std::complex<float> HP;
};

template <typename T> [[nodiscard]]
auto g2R(T g, T M) -> T {
	const auto a = ((T(2) - g) / g) * M;
	return sqrt(a + T(2)) / T(2);
}

template <typename T> [[nodiscard]]
auto R2p1(T _2R, T M) -> T {
	const auto a = max(T(1), _2R);
	const auto b = sqrt(pow(a, T(2)) - T(1));
	const auto c = a + b;
	if (M >= T(1)) { return c; }
	else           { return reciprocal(c); }
}

template <typename T> [[nodiscard]]
auto a22R(T v) -> T {
	return cos(v) * T(2);
}

template <class T> [[nodiscard]]
auto approx_tan_rational(T v) -> T {
	auto a = (v * T(-0.0896638)) + T(0.0388452);
	auto b = (v * T(-0.430871)) + T(0.0404318);
	auto c = (a * v) + T(1.00005);
	auto d = (b * v) + T(1);
	return (c * v) / d;
}

template <typename T> [[nodiscard]]
auto compute_wdiv2(T f, float SR) -> T {
	f = max(f, T(0.01));
	f *= std::numbers::pi_v<float> / SR;
	f = min(f, T(1.50845));
	return approx_tan_rational(f);
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

template <typename T>
struct state {
	std::complex<T> fbk;
	std::complex<T> integrator1_fbk;
	std::complex<T> integrator2_fbk;
};

template <typename T>
struct tilt_coeffs {
	T damp;
	T M;
	T recip_M;
	T BW;
	T wMid;
};

template <typename T> struct coeffs {
	T F;
	T recip_Mpow2;
	T wdiv2;
	tilt_coeffs<T> tilt1;
	tilt_coeffs<T> tilt2;
};

template <typename T>
struct zdf_signal {
	std::complex<T> state;
	T signal;
};

template <typename T> [[nodiscard]]
auto make_coeffs(T f, T gain, T slope, T bw, float SR) -> coeffs<T> {
	const auto slope_mod            = std::pow(T(2), slope);
	const auto M                    = blink::math::convert::db_to_linear(gain / T(-2));
	const auto Mpow2                = pow(M, T(2));
	const auto recip_Mpow2          = reciprocal(Mpow2);
	const auto [slope_lo, slope_hi] = bi_clip(slope_mod, T(2));
	const auto slope_hi_half        = slope_hi / T(2);
	const auto g2R_M                = Mpow2 + recip_Mpow2;
	const auto g2R1                 = g2R(slope_lo, g2R_M);
	const auto g2R2                 = g2R(slope_hi_half, g2R_M);
	const auto g2R1_clamp           = min(T(1), g2R1);
	const auto g2R2_clamp           = min(T(1), g2R2);
	const auto R1x                  = g2R1 / g2R1_clamp;
	const auto a1                   = acos(g2R1_clamp);
	const auto a2                   = acos(g2R2_clamp);
	const auto g2R2_R2p1            = R2p1(g2R2, M);
	const auto tilt1_M              = M * g2R2_R2p1;
	const auto tilt2_M              = M / g2R2_R2p1;
	const auto a2div2               = a2 / T(2);
	const auto tilt1_damp           = a22R(a1 + a2div2) * R1x;
	const auto tilt2_damp           = a22R(a1 - a2div2);
	const auto tilt1_wMid           = sqrt(tilt1_M);
	const auto tilt2_wMid           = sqrt(tilt2_M);
	const auto BWpow2               = pow(T(2), bw);
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

template <typename T> [[nodiscard]]
auto svf_fbk_solve(T x, std::complex<T> zdf) -> std::tuple<T, T> {
	const auto a = x * zdf.real();
	const auto b = a + zdf.imag();
	const auto c = T(1) - zdf.real();
	const auto d = T(1.0) / c;
	return std::make_tuple(b * d, d);
}

template <typename T> [[nodiscard]]
auto svf_input_fbk(zdf_signal<T> in) -> zdf_signal<T> {
	const auto [y, _] = svf_fbk_solve(in.signal, in.state * T(-1));
	return {
		.state  = {{T{1}, T{0}}},
		.signal = in.signal + y
	};
}

template <typename T> [[nodiscard]]
auto svf_integrator_input_fbk(zdf_signal<T> in, std::complex<T> fbk) -> std::tuple<zdf_signal<T>, std::complex<T>> {
	fbk *= T(-1);
	const auto [y, G] = svf_fbk_solve(in.signal, fbk);
	const auto out    = zdf_signal{{T{1}, T{0}}, in.signal + y};
	const auto ext    = std::complex<T>{T{G + in.state.real()}, T{(fbk.imag() + in.state.imag()) * G}};
	return std::make_tuple(out, ext);
}

template <typename T> [[nodiscard]]
auto gain(zdf_signal<T> p, T v) -> zdf_signal<T> {
	p.signal    *= v;
	p.state.zdf *= v;
	return p;
}

template <typename T> [[nodiscard]]
auto zdf_output(zdf_signal<T> zdf, std::complex<T> ext) -> zdf_signal<T> {
	zdf.state = {
		zdf.state.real() * ext.real(),
		zdf.state.imag() + (zdf.state.real() * ext.state.img())
	};
	return zdf;
}

template <typename T> [[nodiscard]]
auto lp_to_bp_integrator(zdf_signal<T> in, std::complex<T> fbk, T wdiv2, T BW) -> std::tuple<zdf_signal<T>, std::complex<T>> {
	in                       = gain(in, BW);
	const auto [signal, ext] = svf_integrator_input_fbk(in, fbk);
	signal                   = integrator(signal, wdiv2);
	signal                   = integrator(signal, wdiv2);
	fbk                      = signal.state;
	return std::make_tuple(zdf_output(signal, ext), fbk);
}

template <typename T> [[nodiscard]]
auto svf(band_shelf::state<T> state, T in, T wdiv2, T damp, T BW) -> std::tuple<band_shelf::state<T>, T, T, T> {
	const auto HP           = svf_input_fbk({state.fbk, in});
	const auto [BP, BP_fbk] = lp_to_bp_integrator(HP, state.integrator1_fbk, wdiv2, BW);
	const auto [LP, LP_fbk] = lp_to_bp_integrator(BP, state.integrator2_fbk, wdiv2, BW);
	const auto BPn = gain(BP, damp);
	const auto fbk = LP.state + BPn.state;
	state = {
		.fbk             = LP.state + BPn.state,
		.integrator1_fbk = BP_fbk,
		.integrator2_fbk = LP_fbk,
	};
	return std::make_tuple(state, LP.signal, HP.signal, BPn.signal);
}

template <typename T> [[nodiscard]]
auto tilt(band_shelf::state<T> state, T in, T wdiv2, tilt_coeffs<T> coeffs) -> std::tuple<band_shelf::state<T>, T> {
	const auto [new_state, LP, HP, BPn] = svf(state, in, wdiv2, coeffs.damp, coeffs.BW);
	return std::make_tuple(new_state, (LP / coeffs.M) + (HP * coeffs.M) + BPn);
}

template <typename T> [[nodiscard]]
auto process(band_shelf::state<T> state, T in, const band_shelf::coeffs<T>& coeffs) -> std::tuple<band_shelf::state<T>, T> {
	auto signal = in;
	std::tie(state, signal) = tilt(state, signal, coeffs._wdiv2, coeffs.tilt1);
	std::tie(state, signal) = tilt(state, signal, coeffs._wdiv2, coeffs.tilt2);
	signal *= coeffs._recip_Mpow2;
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
auto tilt(std::complex<float> t, float f, float F, tilt_coeffs<float> coeffs) -> std::complex<float> {
	const auto [LP, BPn, HP] = svf(f, F, coeffs.damp, coeffs.BW);
	const auto a = LP * coeffs.recip_M;
	const auto b = HP * coeffs.M;
	const auto c = a + b;
	const auto d = c + BPn;
	return t * d;
}

[[nodiscard]]
auto compute(const band_shelf::coeffs<float>& coeffs, float f) -> float {
	auto s = std::complex<float>{1.0f, 0.0f};
	s = tilt(s, f, coeffs.F, coeffs.tilt1);
	s = tilt(s, f, coeffs.F, coeffs.tilt2);
	s *= coeffs.recip_Mpow2;
	return std::abs(s);
}

auto compute(const band_shelf::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert (in_omega.size() == out_magnitude.size());
	auto fn_compute = [&coeffs](float omega) -> float { return compute(coeffs, omega); };
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(), fn_compute);
}

} // eq::filters::shelving::band_shelf::tf

namespace eq::filters::shelving::hilo {

template <typename T>
struct tilt_coeffs {
	T damp;
	T M;
	T recip_M;
	T wMid;
};

template <typename T> struct coeffs {
	T F;
	T Mpow2;
	T recip_Mpow2;
    T wdiv2;
	tilt_coeffs<T> tilt1;
	tilt_coeffs<T> tilt2;
	tilt_coeffs<T> tilt3;
	tilt_coeffs<T> tilt4;
};

template <typename T> struct state {};

template <typename T> [[nodiscard]]
auto make_coeffs(T f, T gain, T slope, float SR) -> hilo::coeffs<T> {
	const auto slope_mod              = std::pow(T(2), slope);
	const auto M                      = blink::math::convert::db_to_linear(gain / T(-2));
	const auto Mpow2                  = pow(M, T(2));
	const auto recip_Mpow2            = reciprocal(Mpow2);
	const auto [slope1_lo, slope1_hi] = bi_clip(slope_mod, T(2));
	const auto slope1_hi_half         = slope1_hi / T(2);
	const auto [slope2_lo, slope2_hi] = bi_clip(slope1_hi_half, T(2));
	const auto slope2_hi_half         = slope2_hi / T(2);
	const auto g2R_M                  = Mpow2 + recip_Mpow2;
	const auto g2R1                   = g2R(slope1_lo, g2R_M);
	const auto g2R2                   = g2R(slope2_lo, g2R_M);
	const auto g2R3                   = g2R(slope2_hi_half, g2R_M);
	const auto g2R1_clamp             = min(T(1), g2R1);
	const auto g2R2_clamp             = min(T(1), g2R2);
	const auto g2R3_clamp             = min(T(1), g2R3);
	const auto R1x                    = g2R1 / g2R1_clamp;
	const auto a1                     = acos(g2R1_clamp);
	const auto a2                     = acos(g2R2_clamp);
	const auto a3                     = acos(g2R3_clamp);
	const auto g2R2_R2p1              = R2p1(g2R2, M);
	const auto g2R3_R2p1              = R2p1(g2R3, M);
	const auto M1                     = sqrt(M * g2R3_R2p1);
	const auto tilt34_M               = sqrt(M / g2R3_R2p1);
	const auto tilt1_M                = M * g2R2_R2p1;
	const auto tilt2_M                = M / g2R2_R2p1;
	const auto a2div2                 = a2 / T(2);
	const auto a3div4                 = a3 / T(4);
	const auto tilt1_damp             = a22R(a1 + a2div2) * R1x;
	const auto tilt2_damp             = a22R(a1 - a2div2);
	const auto tilt3_damp             = a22R((a1 + a2div2) - a3div4);
	const auto tilt4_damp             = a22R((a1 - a2div2) - a3div4);
	const auto tilt1_recip_M          = reciprocal(tilt1_M);
	const auto tilt2_recip_M          = reciprocal(tilt2_M);
	const auto tilt34_recip_M         = reciprocal(tilt34_M);
	const auto tilt1_wMid             = sqrt(tilt1_M);
	const auto tilt2_wMid             = sqrt(tilt2_M);
	const auto tilt34_wMid            = sqrt(tilt34_M);
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

template <typename T> [[nodiscard]]
auto process(hilo::state<T> state, T in, const band_shelf::coeffs<T>& coeffs) -> std::tuple<hilo::state<T>, T> {
	auto signal = in;
	std::tie(state, signal) = tilt(state, signal, coeffs._wdiv2, coeffs.tilt1);
	std::tie(state, signal) = tilt(state, signal, coeffs._wdiv2, coeffs.tilt2);
	std::tie(state, signal) = tilt(state, signal, coeffs._wdiv2, coeffs.tilt3);
	std::tie(state, signal) = tilt(state, signal, coeffs._wdiv2, coeffs.tilt4);
	signal *= coeffs._recip_Mpow2;
	return std::make_tuple(state, signal);
}

} // eq::filters::shelving::hilo

namespace eq::filters::shelving::hilo::tf {

[[nodiscard]]
auto svf(float f, float F, float _2R) -> svf_result {
	return shelving::tf::svf({0.0f, f / F}, _2R);
}

[[nodiscard]]
auto tilt(std::complex<float> t, float f, float F, tilt_coeffs<float> coeffs) -> std::complex<float> {
	const auto [LP, BPn, HP] = svf(f, F * coeffs.wMid, coeffs.damp);
	const auto a = LP * coeffs.recip_M;
	const auto b = HP * coeffs.M;
	const auto c = a + b;
	const auto d = c + BPn;
	return t * d;
}

[[nodiscard]]
auto compute(const hilo::coeffs<float>& coeffs, float f) -> std::complex<float> {
	auto s = std::complex<float>{1.0f, 0.0f};
	s = hilo::tf::tilt(s, f, coeffs.F, coeffs.tilt1);
	s = hilo::tf::tilt(s, f, coeffs.F, coeffs.tilt2);
	s = hilo::tf::tilt(s, f, coeffs.F, coeffs.tilt3);
	s = hilo::tf::tilt(s, f, coeffs.F, coeffs.tilt4);
	return s;
}

} // eq::filters::shelving::hilo::tf

namespace eq::filters::shelving::hi_shelf {

template <typename T> using coeffs = hilo::coeffs<T>;
template <typename T> using state  = hilo::state<T>;

template <typename T> [[nodiscard]]
auto make_coeffs(T f, T gain, T slope, float SR) -> coeffs<T> {
	return hilo::make_coeffs(f, gain, slope, SR);
}

template <typename T> [[nodiscard]]
auto process(hi_shelf::state<T> state, T in, const band_shelf::coeffs<T>& coeffs) -> std::tuple<hi_shelf::state<T>, T> {
	auto signal = in;
	std::tie(state, signal) = tilt(state, signal, coeffs._wdiv2, coeffs.tilt1);
	std::tie(state, signal) = tilt(state, signal, coeffs._wdiv2, coeffs.tilt2);
	std::tie(state, signal) = tilt(state, signal, coeffs._wdiv2, coeffs.tilt3);
	std::tie(state, signal) = tilt(state, signal, coeffs._wdiv2, coeffs.tilt4);
	signal *= coeffs._recip_Mpow2;
	return std::make_tuple(state, signal);
}

} // eq::filters::shelving::hi_shelf

namespace eq::filters::shelving::hi_shelf::tf {

[[nodiscard]]
auto compute(const hi_shelf::coeffs<float>& coeffs, float f) -> float {
	return std::abs(hilo::tf::compute(coeffs, f)) * coeffs.Mpow2;
}

auto compute(const hi_shelf::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert (in_omega.size() == out_magnitude.size());
	auto fn_compute = [&coeffs](float omega) -> float { return compute(coeffs, omega); };
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(), fn_compute);
}

} // eq::filters::shelving::hi_shelf::tf

namespace eq::filters::shelving::lo_shelf {

template <typename T> using coeffs = hilo::coeffs<T>;
template <typename T> using state  = hilo::state<T>;

template <typename T> [[nodiscard]]
auto make_coeffs(T f, T gain, T slope, float SR) -> coeffs<T> {
	return hilo::make_coeffs(f, gain, slope, SR);
}

template <typename T> [[nodiscard]]
auto process(lo_shelf::state<T> state, T in, const band_shelf::coeffs<T>& coeffs) -> std::tuple<lo_shelf::state<T>, T> {
	auto signal = in;
	std::tie(state, signal) = tilt(state, signal, coeffs._wdiv2, coeffs.tilt1);
	std::tie(state, signal) = tilt(state, signal, coeffs._wdiv2, coeffs.tilt2);
	std::tie(state, signal) = tilt(state, signal, coeffs._wdiv2, coeffs.tilt3);
	std::tie(state, signal) = tilt(state, signal, coeffs._wdiv2, coeffs.tilt4);
	signal *= coeffs._recip_Mpow2;
	return std::make_tuple(state, signal);
}

} // eq::filters::shelving::lo_shelf

namespace eq::filters::shelving::lo_shelf::tf {

[[nodiscard]]
auto compute(const lo_shelf::coeffs<float>& coeffs, float f) -> float {
	return std::abs(hilo::tf::compute(coeffs, f)) * coeffs.recip_Mpow2;
}

auto compute(const lo_shelf::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert (in_omega.size() == out_magnitude.size());
	auto fn_compute = [&coeffs](float omega) -> float { return compute(coeffs, omega); };
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(), fn_compute);
}

} // eq::filters::shelving::lo_shelf::tf

