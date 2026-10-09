#pragma once

#include <cassert>
#include <complex>
#include <cmath>
#include <DSP/MLDSPOps.h>
#include <span>
#include <algorithm>

namespace eq::filters::bell {

template <typename T> struct coeffs { T b0, b1, b2; T a1, a2; };
template <typename T> struct state  { T x1 = 0; T x2 = 0; T y1 = 0; T y2 = 0; };

template <typename T> [[nodiscard]]
auto make_coeffs(T omega, T q, T A) -> bell::coeffs<T> {
	const auto w0     = ml::kPi * omega;
	const auto cos_w0 = std::cos(w0);
	const auto sin_w0 = std::sin(w0);
	const auto alpha  = sin_w0 / (T(2) * q);
	const auto a      = A;
	const auto a0     = T(1) + alpha / a;
	const auto b0     = (T(1) + alpha * a) / a0;
	const auto b1     = -T(2) * cos_w0 / a0;
	const auto b2     = (T(1) - alpha * a) / a0;
	const auto a1     = -T(2) * cos_w0 / a0;
	const auto a2     = (T(1) - alpha / a) / a0;
	return {b0, b1, b2, a1, a2};
}

template <typename T> [[nodiscard]]
auto process(bell::state<T> state, T in, const bell::coeffs<T>& coeffs) -> std::tuple<bell::state<T>, T> {
	const auto out =
		coeffs.b0 * in
		+ coeffs.b1 * state.x1
		+ coeffs.b2 * state.x2
		- coeffs.a1 * state.y1
		- coeffs.a2 * state.y2;
	state.x2 = state.x1;
	state.x1 = in;
	state.y2 = state.y1;
	state.y1 = out;
	return {state, out};
}

[[nodiscard]]
auto transfer(const bell::coeffs<float>& coeffs, float omega) -> float {
	const auto w           = ml::kPi * omega;
	const auto z           = std::polar(1.f, -w);
	const auto z2          = z * z;
	const auto numerator   = coeffs.b0 + coeffs.b1 * z + coeffs.b2 * z2;
	const auto denominator = 1.f + coeffs.a1 * z + coeffs.a2 * z2;
	return std::abs(numerator / denominator);
}

auto transfer(const bell::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert(in_omega.size() == out_magnitude.size());
	auto fn_transfer = [&coeffs](float omega) -> float { return transfer(coeffs, omega); };
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(), fn_transfer);
}

} // eq::filters::bell

namespace eq::filters::pass_hi {

template <typename T> struct coeffs {};
template <typename T> struct state  {};

template <typename T> [[nodiscard]]
auto make_coeffs(T omega, T k) -> pass_hi::coeffs<T> {
	// @TODO: implement
	return {};
}

template <typename T> [[nodiscard]]
auto process(pass_hi::state<T> state, T in, const pass_hi::coeffs<T>& coeffs) -> std::tuple<pass_hi::state<T>, T> {
	// @TODO: implement
	return {};
}

[[nodiscard]]
auto transfer(const pass_hi::coeffs<float>& coeffs, float omega) -> float {
	// @TODO: implement
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
	return {};
}

template <typename T> [[nodiscard]]
auto process(pass_lo::state<T> state, T in, const pass_lo::coeffs<T>& coeffs) -> std::tuple<pass_lo::state<T>, T> {
	return {};
}

[[nodiscard]]
auto transfer(const pass_lo::coeffs<float>& coeffs, float omega) -> float {
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
	return {};
}

template <typename T> [[nodiscard]]
auto process(shelf_hi::state<T> state, T in, const shelf_hi::coeffs<T>& coeffs) -> std::tuple<shelf_hi::state<T>, T> {
	return {};
}

[[nodiscard]]
auto transfer(const shelf_hi::coeffs<float>& coeffs, float omega) -> float {
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
	// @TODO: implement
}

template <typename T> [[nodiscard]]
auto process(shelf_lo::state<T> state, T in, const shelf_lo::coeffs<T>& coeffs) -> std::tuple<shelf_lo::state<T>, T> {
	// @TODO: implement
}

[[nodiscard]]
auto transfer(const shelf_lo::coeffs<float>& coeffs, float omega) -> float {
	// @TODO: implement
}

auto transfer(const shelf_lo::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert(in_omega.size() == out_magnitude.size());
	auto fn_transfer = [&coeffs](float omega) -> float { return transfer(coeffs, omega); };
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(), fn_transfer);
}

} // eq::filters::shelf_lo
