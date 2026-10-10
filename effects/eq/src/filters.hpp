#pragma once

#include "adis-filters.h"
#include <cassert>
#include <complex>
#include <DSP/MLDSPOps.h>
#include <span>
#include <algorithm>
#include <numbers>

namespace eq::filters {

struct states {
	std::complex<float> v1;
	std::complex<float> v2;
};

[[nodiscard]]
auto solve(
	std::complex<float> a,
	std::complex<float> b,
	std::complex<float> c,
	std::complex<float> d,
	std::complex<float> r1,
	std::complex<float> r2
) -> std::array<std::complex<float>, 2> {
	const auto determinant = a * d - b * c;
	return {
		(r1 * d - b * r2) / determinant,
		(a * r2 - r1 * c) / determinant
	};
}

[[nodiscard]]
auto shelf_states(float a1, float a2, float a3, float omega) -> states {
    const auto z = std::polar<float>(
        1.0f, 2.0f * std::numbers::pi_v<float> * omega);

    const auto s = solve(
        z + 1.0f - 2.0f * a1,
        2.0f * a2,
        -2.0f * a2,
        z - 1.0f + 2.0f * a3,
        2.0f * a2,
        2.0f * a3);

    const auto scale = (1.0f + z) * 0.5f;
    return {scale * s[0], scale * s[1]};
}

} // eq::filters

namespace eq::filters::bell {

template <typename T> struct coeffs { T a1, a2, a3, m1; };
template <typename T> struct state  { T ic1eq, ic2eq; };

template <typename T> [[nodiscard]]
auto make_coeffs(T omega, T k, T A) -> bell::coeffs<T> {
	const auto kc       = k / A;
	const auto pi_omega = ml::kPi * omega;
	const auto g        = tan(pi_omega);
	const auto a1       = 1.f / (1.f + g * (g + kc));
	const auto a2       = g * a1;
	const auto a3       = g * a2;
	const auto m1       = kc * (A * A - 1.f);
	return {a1, a2, a3, m1};
}

template <typename T> [[nodiscard]]
auto process(bell::state<T> state, T in, const bell::coeffs<T>& coeffs) -> std::tuple<bell::state<T>, T> {
	const auto v3 = in - state.ic2eq;
	const auto v1 = coeffs.a1 * state.ic1eq + coeffs.a2 * v3;
	const auto v2 = state.ic2eq + coeffs.a2 * state.ic1eq + coeffs.a3 * v3;
	state.ic1eq = 2.f * v1 - state.ic1eq;
	state.ic2eq = 2.f * v2 - state.ic2eq;
	return std::make_tuple(state, in + coeffs.m1 * v1);
}

[[nodiscard]]
auto transfer(const bell::coeffs<float>& coeffs, float omega) -> float {
	const auto s = shelf_states(coeffs.a1, coeffs.a2, coeffs.a3, omega);
	return std::abs(1.0f + coeffs.m1 * s.v1);
}

auto transfer(const bell::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert (in_omega.size() == out_magnitude.size());
	auto fn_transfer = [&coeffs](float omega) -> float { return transfer(coeffs, omega); };
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(), fn_transfer);
}

} // eq::filters::bell

namespace eq::filters::pass_hi {

template <typename T> struct coeffs {};
template <typename T> struct state  {};

template <typename T> [[nodiscard]]
auto make_coeffs(T omega, T k) -> pass_hi::coeffs<T> {
	// @TODO:
	return {};
}

template <typename T> [[nodiscard]]
auto process(pass_hi::state<T> state, T in, const pass_hi::coeffs<T>& coeffs) -> std::tuple<pass_hi::state<T>, T> {
	// @TODO:
	return {};
}

[[nodiscard]]
auto transfer(const pass_hi::coeffs<float>& coeffs, float omega) -> float {
	// @TODO:
	return {};
}

auto transfer(const pass_hi::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert(in_omega.size() == out_magnitude.size());
	auto fn_transfer = [&coeffs](float omega) -> float { return transfer(coeffs, omega); };
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(), fn_transfer);
}

} // eq::filters::pass_hi

namespace eq::filters::pass_lo {

template <typename T> struct coeffs {};
template <typename T> struct state  {};

template <typename T> [[nodiscard]]
auto make_coeffs(T omega, T k) -> pass_lo::coeffs<T> {
	// @TODO:
	return {};
}

template <typename T> [[nodiscard]]
auto process(pass_lo::state<T> state, T in, const pass_lo::coeffs<T>& coeffs) -> std::tuple<pass_lo::state<T>, T> {
	// @TODO:
	return {};
}

[[nodiscard]]
auto transfer(const pass_lo::coeffs<float>& coeffs, float omega) -> float {
	// @TODO:
	return {};
}

auto transfer(const pass_lo::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert(in_omega.size() == out_magnitude.size());
	auto fn_transfer = [&coeffs](float omega) -> float { return transfer(coeffs, omega); };
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(), fn_transfer);
}

} // eq::filters::pass_lo

namespace eq::filters::shelf_hi {

template <typename T> struct coeffs {};
template <typename T> struct state  {};

template <typename T> [[nodiscard]]
auto make_coeffs(T omega, T k, T A) -> shelf_hi::coeffs<T> {
	// @TODO:
	return {};
}

template <typename T> [[nodiscard]]
auto process(shelf_hi::state<T> state, T in, const shelf_hi::coeffs<T>& coeffs) -> std::tuple<shelf_hi::state<T>, T> {
	// @TODO:
	return {};
}

[[nodiscard]]
auto transfer(const shelf_hi::coeffs<float>& coeffs, float omega) -> float {
	// @TODO:
	return {};
}

auto transfer(const shelf_hi::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert(in_omega.size() == out_magnitude.size());
	auto fn_transfer = [&coeffs](float omega) -> float { return transfer(coeffs, omega); };
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(), fn_transfer);
}

} // eq::filters::shelf_hi

namespace eq::filters::shelf_lo {

template <typename T> struct coeffs {};
template <typename T> struct state  {};

template <typename T> [[nodiscard]]
auto make_coeffs(T omega, T k, T A) -> shelf_lo::coeffs<T> {
}

template <typename T> [[nodiscard]]
auto process(shelf_lo::state<T> state, T in, const shelf_lo::coeffs<T>& coeffs) -> std::tuple<shelf_lo::state<T>, T> {
	// @TODO:
	return {};
}

[[nodiscard]]
auto transfer(const shelf_lo::coeffs<float>& coeffs, float omega) -> float {
	// @TODO:
	return {};
}

auto transfer(const shelf_lo::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert(in_omega.size() == out_magnitude.size());
	auto fn_transfer = [&coeffs](float omega) -> float { return transfer(coeffs, omega); };
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(), fn_transfer);
}

} // eq::filters::shelf_lo

/*
 *
#include <cmath>
#include <complex>
#include <numbers>

namespace filter_response {

using complex = std::complex<double>;

[[nodiscard]]
auto magnitude(const complex& z) -> double {
	return std::abs(z);
}

[[nodiscard]]
auto response_2(
	const float* A,
	const float* d1,
	const float* d2,
	int count,
	double frequency,
	double sample_rate,
	bool high_pass,
	double gain = 1.0
) -> double {
	const auto q = std::polar(1.0, -2.0 * std::numbers::pi * frequency / sample_rate);
	const auto q2 = q * q;
	auto result = complex{gain, 0.0};

	for (int i = 0; i < count; ++i) {
		const auto numerator = high_pass
			? complex{A[i], 0.0} * (1.0 - q) * (1.0 - q)
			: complex{A[i], 0.0} * (1.0 + q) * (1.0 + q);

		const auto denominator = 1.0 - double(d1[i]) * q - double(d2[i]) * q2;

		result *= numerator / denominator;
	}

	return magnitude(result);
}

[[nodiscard]]
auto response_4(
	const float* A,
	const float* d1,
	const float* d2,
	const float* d3,
	const float* d4,
	int count,
	double frequency,
	double sample_rate,
	bool band_stop,
	double r = 0.0,
	double s = 0.0,
	double gain = 1.0
) -> double {
	const auto q = std::polar(1.0, -2.0 * std::numbers::pi * frequency / sample_rate);
	const auto q2 = q * q;
	const auto q3 = q2 * q;
	const auto q4 = q2 * q2;

	auto result = complex{gain, 0.0};

	for (int i = 0; i < count; ++i) {
		const auto denominator =
			1.0
			- double(d1[i]) * q
			- double(d2[i]) * q2
			- double(d3[i]) * q3
			- double(d4[i]) * q4;

		const auto numerator = band_stop
			? double(A[i]) * (1.0 - r * q + s * q2 - r * q3 + q4)
			: double(A[i]) * (1.0 - 2.0 * q2 + q4);

		result *= numerator / denominator;
	}

	return magnitude(result);
}

// Butterworth

[[nodiscard]]
auto bw_low_pass(const BWLowPass& f, double hz, double fs) -> double {
	return response_2(f.A, f.d1, f.d2, f.n, hz, fs, false);
}

[[nodiscard]]
auto bw_high_pass(const BWHighPass& f, double hz, double fs) -> double {
	return response_2(f.A, f.d1, f.d2, f.n, hz, fs, true);
}

[[nodiscard]]
auto bw_band_pass(const BWBandPass& f, double hz, double fs) -> double {
	return response_4(f.A, f.d1, f.d2, f.d3, f.d4, f.n, hz, fs, false);
}

[[nodiscard]]
auto bw_band_stop(const BWBandStop& f, double hz, double fs) -> double {
	return response_4(f.A, f.d1, f.d2, f.d3, f.d4, f.n, hz, fs,
		true, f.r, f.s);
}

// Chebyshev

[[nodiscard]]
auto che_low_pass(const CHELowPass& f, double hz, double fs) -> double {
	return response_2(f.A, f.d1, f.d2, f.m, hz, fs, false, f.ep);
}

[[nodiscard]]
auto che_high_pass(const CHEHighPass& f, double hz, double fs) -> double {
	return response_2(f.A, f.d1, f.d2, f.m, hz, fs, true, f.ep);
}

[[nodiscard]]
auto che_band_pass(const CHEBandPass& f, double hz, double fs) -> double {
	return response_4(f.A, f.d1, f.d2, f.d3, f.d4, f.m, hz, fs,
		false, 0.0, 0.0, f.ep);
}

[[nodiscard]]
auto che_band_stop(const CHEBandStop& f, double hz, double fs) -> double {
	return response_4(f.A, f.d1, f.d2, f.d3, f.d4, f.m, hz, fs,
		true, f.r, f.s, f.ep);
}

}
*/
