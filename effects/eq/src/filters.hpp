#pragma once

#include <cassert>
#include <complex>
#include <cmath>
#include <DSP/MLDSPOps.h>
#include <span>
#include <algorithm>

namespace eq {

[[nodiscard]]
auto tan(ml::DSPVector x) -> ml::DSPVector {
	return ml::sin(x) / ml::cos(x);
}

} // eq

namespace eq::filters::bell {

template <typename T> struct coeffs { T a1, a2, a3, m1; };
template <typename T> struct state  { T ic1eq, ic2eq; };

template <typename T> [[nodiscard]]
auto make_coeffs(T omega, T k, T A) -> bell::coeffs<T> {
	auto kc      = k / A;
	auto piOmega = ml::kPi * omega;
	auto g       = tan(piOmega);
	auto a1      = 1.f / (1.f + g * (g + kc));
	auto a2      = g * a1;
	auto a3      = g * a2;
	auto m1      = kc * (A * A - 1.f);
	return {a1, a2, a3, m1};
}

template <typename T> [[nodiscard]]
auto process(bell::state<T> state, T in, const bell::coeffs<T>& coeffs) -> std::tuple<bell::state<T>, T> {
	auto v3 = in - state.ic2eq;
	auto v1 = coeffs.a1 * state.ic1eq + coeffs.a2 * v3;
	auto v2 = state.ic2eq + coeffs.a2 * state.ic1eq + coeffs.a3 * v3;
	state.ic1eq = 2.f * v1 - state.ic1eq;
	state.ic2eq = 2.f * v2 - state.ic2eq;
	return std::make_tuple(state, in + coeffs.m1 * v1);
}

[[nodiscard]]
auto transfer(const bell::coeffs<float>& coeffs, float omega) -> float {
	static constexpr auto ONE = std::complex<float>{1.f, 0.f};
	auto w           = ml::kPi * 2.f * omega;
	auto z           = std::exp(std::complex<float>(0.f, w));
	auto z_inv       = 1.f / z;
	auto z_inv2      = z_inv * z_inv;
	auto numerator   = ONE + coeffs.m1 * coeffs.a2 * (ONE + coeffs.a2 * z_inv);
	auto denominator = ONE + (coeffs.a2 + coeffs.a3) * z_inv + coeffs.a3 * z_inv2;
	auto H           = numerator / denominator;
	return std::abs(H);
}

auto transfer(const bell::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert (in_omega.size() == out_magnitude.size());
	auto fn_transfer = [&coeffs](float omega) -> float { return transfer(coeffs, omega); };
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(), fn_transfer);
}

} // eq::filters::bell
