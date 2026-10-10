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

auto transfer_j(std::complex<float> v) -> float {
	const auto r = v.real();
	const auto i = v.imag();
	return std::sqrt((r * r) + (i * i));
}

} // eq::filters

namespace eq::filters::band_shelf {

template <typename T> struct coeffs {
	float _F;
	float _2R1;
	float _2R2;
	float _Mpow2_1;
	float _Mpow2_2;
	float _recip_Mpow2;
	float _BW2R;
};

template <typename T> struct state  { T ic1eq, ic2eq; };

template <typename T> [[nodiscard]]
auto g2R(T g, T Mpow2, T recip_Mpow2) -> T {
	const auto a = ((T(2) - g) / g) * (Mpow2 + recip_Mpow2);
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

struct transfer_svf_to_bp_result {
	std::complex<float> LP;
	std::complex<float> BPn;
	std::complex<float> HP;
};

[[nodiscard]]
auto transfer_lp_to_bp(std::complex<float> v, float _2R) -> std::complex<float> {
	const auto a = 1.0f / v;
	const auto b = v + a;
	return b * (1.0f / _2R);
}

[[nodiscard]]
auto transfer_svf_to_bp(float f, float F, float _2R, float BW2R) -> transfer_svf_to_bp_result {
	const auto s     = transfer_lp_to_bp({0.0f, f / F}, _2R);
	const auto spow2 = s * s;
	const auto a     = spow2 + 1.0f;
	const auto b     = s * _2R;
	const auto c     = a + b;
	const auto LP    = 1.0f / c;
	const auto BPn   = (LP * s) * _2R;
	const auto HP    = BPn * s;
	return {
		.LP  = LP,
		.BPn = BPn,
		.HP  = HP
	};
}

[[nodiscard]]
auto transfer_tilt2_to_BP(std::complex<float> t, float f, float F, float _2R, float Mpow2, float recip_Mpow2, float BW2R) -> std::complex<float> {
	const auto [LP, BPn, HP] = transfer_svf_to_bp(f, F, _2R, BW2R);
	const auto a = LP * recip_Mpow2;
	const auto b = HP * Mpow2;
	const auto c = a + b;
	const auto d = c + BPn;
	return t * d;
}

template <typename T> [[nodiscard]]
auto make_coeffs(T f, T gain, T slope, T bw) -> band_shelf::coeffs<T> {
	const auto slope_clamp    = clamp(slope, T(1), T(2));
	const auto M              = blink::math::convert::db_to_linear(gain / T(-2));
	const auto Mpow2          = pow(M, T(2));
	const auto recip_Mpow2    = reciprocal(Mpow2);
	const auto Mpow2_biclip   = bi_clip(slope_clamp, T(2));
	const auto Mpow2_low      = Mpow2_biclip.first;
	const auto Mpow2_high     = Mpow2_biclip.second;
	const auto g2R_low        = g2R(Mpow2_low, Mpow2, recip_Mpow2);
	const auto g2R_high       = g2R(Mpow2_high / T(2), Mpow2, recip_Mpow2);
	const auto g2R_low_clamp  = min(T(1), g2R_low);
	const auto g2R_high_clamp = min(T(1), g2R_high);
	const auto R1x            = g2R_low / g2R_low_clamp;
	const auto a1             = acos(g2R_low_clamp);
	const auto a2             = acos(g2R_high_clamp);
	const auto g2R_high_R2p1  = R2p1(g2R_high, M);
	const auto Mpow2_1        = M * g2R_high_R2p1;
	const auto Mpow2_2        = M / g2R_high_R2p1;
	const auto a2div2         = a2 / T(2);
	const auto _2R1           = a22R(a1 + a2div2) * R1x;
	const auto _2R2           = a22R(a1 - a2div2);
	const auto BWpow2         = pow(bw, T(2));
	const auto wMid           = sqrt(Mpow2);
	const auto BW2R           = (BWpow2 - reciprocal(BWpow2)) * wMid;
	return {
		._F           = f,
		._2R1         = _2R1,
		._2R2         = _2R2,
		._Mpow2_1     = Mpow2_1,
		._Mpow2_2     = Mpow2_2,
		._recip_Mpow2 = recip_Mpow2,
		._BW2R        = BW2R
	};
}

template <typename T> [[nodiscard]]
auto process(band_shelf::state<T> state, T in, const band_shelf::coeffs<T>& coeffs) -> std::tuple<band_shelf::state<T>, T> {
	return {};
}

[[nodiscard]]
auto transfer(const band_shelf::coeffs<float>& coeffs, float f) -> float {
	const auto a = std::complex<float>{1.0f, 0.0f};
	const auto b = transfer_tilt2_to_BP(a, f, coeffs._F, coeffs._2R1, coeffs._Mpow2_1, coeffs._recip_Mpow2, coeffs._BW2R);
	const auto c = transfer_tilt2_to_BP(b, f, coeffs._F, coeffs._2R2, coeffs._Mpow2_2, coeffs._recip_Mpow2, coeffs._BW2R);
	const auto d = c * coeffs._recip_Mpow2;
	return transfer_j(d);
}

auto transfer(const band_shelf::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert (in_omega.size() == out_magnitude.size());
	auto fn_transfer = [&coeffs](float omega) -> float { return transfer(coeffs, omega); };
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(), fn_transfer);
}

} // eq::filters::band_shelf
