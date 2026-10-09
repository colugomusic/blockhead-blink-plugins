#pragma once

#include <cassert>
#include <complex>
#include <cmath>
#include <DSP/MLDSPOps.h>
#include <span>
#include <algorithm>

namespace eq::filters::bell {

template <typename T>
struct coeffs {
};

template <typename T>
struct state {
};

template <typename T>
[[nodiscard]]
auto make_coeffs(T omega, T q, T A) -> bell::coeffs<T> {
	// @TODO:
	return {};
}

template <typename T>
[[nodiscard]]
auto process(bell::state<T> state, T in, const bell::coeffs<T>& coeffs) -> std::tuple<bell::state<T>, T> {
	// @TODO:
	return {};
}

[[nodiscard]]
auto transfer(const bell::coeffs<float>& coeffs, float omega) -> float {
	// @TODO:
	return {};
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
