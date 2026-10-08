#pragma once

#include <cassert>
#include <complex>
#include <cmath>
#include <DSP/MLDSPOps.h>
#include <span>
#include <algorithm>

namespace eq {

constexpr auto OMEGA_MIN = 0.00001f;
constexpr auto OMEGA_MAX = 0.49f;

template <typename T> [[nodiscard]]
auto clamp_omega(T omega) -> T {
  return clamp(omega, T{OMEGA_MIN}, T{OMEGA_MAX});
}

[[nodiscard]]
auto clamp_omega(float omega) -> float {
  return std::clamp(omega, OMEGA_MIN, OMEGA_MAX);
}

[[nodiscard]]
auto tan(ml::DSPVector x) -> ml::DSPVector {
	return ml::sin(x) / ml::cos(x);
}

[[nodiscard]]
auto tan(float x) -> float {
	return std::tan(x);
}

[[nodiscard]]
auto transfer_2pole(
	float a00, float a01,
	float a10, float a11,
	float b0,  float b1,
	float c0,  float c1,
	float d,
	float omega) -> float
{
	const auto z   = std::polar(1.0f, ml::kPi * clamp_omega(omega));
	const auto m00 = z - a00;
	const auto m01 = -a01;
	const auto m10 = -a10;
	const auto m11 = z - a11;
	const auto det = m00 * m11 - m01 * m10;
	const auto s0  = (m11 * b0 - m01 * b1) / det;
	const auto s1  = (-m10 * b0 + m00 * b1) / det;
	const auto h   = d + c0 * s0 + c1 * s1;
	return std::abs(h);
}

} // eq

namespace eq::filters::bell {

template <typename T>
struct coeffs {
	T b0, b1, b2;
	T a1, a2;
};

template <typename T>
struct state {
	T x1 = 0;
	T x2 = 0;
	T y1 = 0;
	T y2 = 0;
};

template <typename T>
[[nodiscard]]
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

template <typename T>
[[nodiscard]]
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

template <typename T> struct coeffs { T g0, g1, g2, gk; };
template <typename T> struct state  { T ic1eq, ic2eq; };

template <typename T> [[nodiscard]]
auto make_coeffs(T omega, T k) -> pass_hi::coeffs<T> {
	const auto pi_omega = ml::kPi * clamp_omega(omega);
	const auto s1       = sin(pi_omega);
	const auto s2       = sin(2.f * pi_omega);
	const auto nrm      = 1.f / (2.f + k * s2);
	const auto cg0      = s2 * nrm;
	const auto cg1      = (-2.0f * s1 * s1 - k * s2) * nrm;
	const auto cg2      = (2.0f * s1 * s1) * nrm;
	return {cg0, cg1, cg2, k};
}

template <typename T> [[nodiscard]]
auto process(pass_hi::state<T> state, T in, const pass_hi::coeffs<T>& coeffs) -> std::tuple<pass_hi::state<T>, T> {
	const auto t0 = in - state.ic2eq;
	const auto t1 = coeffs.g0 * t0 + coeffs.g1 * state.ic1eq;
	const auto t2 = coeffs.g2 * t0 + coeffs.g0 * state.ic1eq;
	const auto v1 = t1 + state.ic1eq;
	const auto v2 = t2 + state.ic2eq;
	state.ic1eq += 2.f * t1;
	state.ic2eq += 2.f * t2;
	return std::make_tuple(state, in - coeffs.gk * v1 - v2);
}

[[nodiscard]]
auto transfer(const pass_hi::coeffs<float>& coeffs, float omega) -> float {
// 	static constexpr auto ONE = std::complex<float>{1.f, 0.f};
// 	const auto w           = ml::kPi * 2.f * omega;
// 	const auto z           = std::exp(std::complex<float>(0.f, w));
// 	const auto z_inv       = 1.f / z;
// 	const auto z_inv2      = z_inv * z_inv;
// 	const auto numerator   = ONE - 2.f * z_inv + z_inv2;
// 	const auto denominator = ONE + (coeffs.g0 + coeffs.gk * coeffs.g1) * z_inv + (coeffs.g0 * coeffs.g1 - coeffs.gk * coeffs.g2 - coeffs.g2) * z_inv2;
// 	const auto H           = numerator / denominator;
// 	return std::abs(H);
// ::
	const auto a00 = 1.f + 2.f * coeffs.g1;
	const auto a01 = -2.f * coeffs.g0;
	const auto a10 = 2.f * coeffs.g0;
	const auto a11 = 1.f - 2.f * coeffs.g2;

	const auto b0 = 2.f * coeffs.g0;
	const auto b1 = 2.f * coeffs.g2;

	const auto c0 =
		-coeffs.g0
		- coeffs.g1 * coeffs.gk
		- coeffs.gk;

	const auto c1 =
		coeffs.g0 * coeffs.gk
		+ coeffs.g2
		- 1.f;

	const auto d =
		1.f
		- coeffs.g0 * coeffs.gk
		- coeffs.g2;

	return eq::transfer_2pole(
		a00, a01, a10, a11,
		b0, b1,
		c0, c1,
		d,
		omega);
}

auto transfer(const pass_hi::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert(in_omega.size() == out_magnitude.size());
	auto fn_transfer = [&coeffs](float omega) -> float { return transfer(coeffs, omega); };
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(), fn_transfer);
}

} // eq::filters::pass_hi

namespace eq::filters::pass_lo {

template <typename T> struct coeffs { T g0, g1, g2; };
template <typename T> struct state  { T ic1eq, ic2eq; };

template <typename T> [[nodiscard]]
auto make_coeffs(T omega, T k) -> pass_lo::coeffs<T> {
	const auto pi_omega = ml::kPi * clamp_omega(omega);
	const auto s1       = sin(pi_omega);
	const auto s2       = sin(2.f * pi_omega);
	const auto nrm      = 1.f / (2.f + k * s2);
	const auto cg0      = s2 * nrm;
	const auto cg1      = (-2.0f * s1 * s1 - k * s2) * nrm;
	const auto cg2      = (2.0f * s1 * s1) * nrm;
	return {cg0, cg1, cg2};
}

template <typename T> [[nodiscard]]
auto process(pass_lo::state<T> state, T in, const pass_lo::coeffs<T>& coeffs) -> std::tuple<pass_lo::state<T>, T> {
	const auto t0 = in - state.ic2eq;
	const auto t1 = coeffs.g0 * t0 + coeffs.g1 * state.ic1eq;
	const auto t2 = coeffs.g2 * t0 + coeffs.g0 * state.ic1eq;
	const auto v2 = t2 + state.ic2eq;
	state.ic1eq += 2.f * t1;
	state.ic2eq += 2.f * t2;
	return v2;
}

[[nodiscard]]
auto transfer(const pass_lo::coeffs<float>& coeffs, float omega) -> float {
	// static constexpr auto ONE = std::complex<float>{1.f, 0.f};
	// const auto w           = ml::kPi * 2.f * omega;
	// const auto z           = std::exp(std::complex<float>(0.f, w));
	// const auto z_inv       = 1.f / z;
	// const auto z_inv2      = z_inv * z_inv;
	// const auto numerator   = coeffs.g0 * (ONE + 2.f * z_inv + z_inv2);
	// const auto denominator = ONE + (coeffs.g0 + coeffs.g1) * z_inv + (coeffs.g0 * coeffs.g1 - coeffs.g2) * z_inv2;
	// const auto H           = numerator / denominator;
	// return std::abs(H);
    const auto a00 = 1.f + 2.f * coeffs.g1;
	const auto a01 = -2.f * coeffs.g0;
	const auto a10 = 2.f * coeffs.g0;
	const auto a11 = 1.f - 2.f * coeffs.g2;

	const auto b0 = 2.f * coeffs.g0;
	const auto b1 = 2.f * coeffs.g2;

	const auto c0 = coeffs.g0;
	const auto c1 = 1.f - coeffs.g2;
	const auto d  = coeffs.g2;

	return eq::transfer_2pole(
		a00, a01, a10, a11,
		b0, b1,
		c0, c1,
		d,
		omega);
}

auto transfer(const pass_lo::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert(in_omega.size() == out_magnitude.size());
	auto fn_transfer = [&coeffs](float omega) -> float { return transfer(coeffs, omega); };
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(), fn_transfer);
}

} // eq::filters::pass_lo

namespace eq::filters::shelf_hi {

template <typename T> struct coeffs { T a1, a2, a3, m0, m1, m2; };
template <typename T> struct state  { T ic1eq, ic2eq; };

template <typename T> [[nodiscard]]
auto make_coeffs(T omega, T k, T A) -> shelf_hi::coeffs<T> {
	const auto pi_omega = ml::kPi * clamp_omega(omega);
	const auto g        = tan(pi_omega) * sqrt(A);
	const auto ca1      = 1.f / (1.f + g * (g + k));
	const auto ca2      = g * ca1;
	const auto ca3      = g * ca2;
	const auto cm0      = A * A;
	const auto cm1      = k * (1.f - A) * A;
	const auto cm2      = 1.f - A * A;
	return {ca1, ca2, ca3, cm0, cm1, cm2};
}

template <typename T> [[nodiscard]]
auto process(shelf_hi::state<T> state, T in, const shelf_hi::coeffs<T>& coeffs) -> std::tuple<shelf_hi::state<T>, T> {
	const auto v3 = in - state.ic2eq;
	const auto v1 = coeffs.a1 * state.ic1eq + coeffs.a2 * v3;
	const auto v2 = state.ic2eq + coeffs.a2 * state.ic1eq + coeffs.a3 * v3;
	state.ic1eq = 2.f * v1 - state.ic1eq;
	state.ic2eq = 2.f * v2 - state.ic2eq;
	return std::make_tuple(state, coeffs.m0 * in + coeffs.m1 * v1 + coeffs.m2 * v2);
}

[[nodiscard]]
auto transfer(const shelf_hi::coeffs<float>& coeffs, float omega) -> float {
	// static constexpr auto ONE = std::complex<float>{1.f, 0.f};
	// const auto w           = ml::kPi * 2.f * omega;
	// const auto z           = std::exp(std::complex<float>(0.f, w));
	// const auto z_inv       = 1.f / z;
	// const auto z_inv2      = z_inv * z_inv;
	// const auto numerator   = coeffs.m0 + coeffs.m1 * coeffs.a2 * (ONE + coeffs.a2 * z_inv) + coeffs.m2 * (ONE + 2.f * coeffs.a2 * z_inv + coeffs.a2 * coeffs.a3 * z_inv2);
	// const auto denominator = ONE + (coeffs.a2 + coeffs.a3) * z_inv + coeffs.a3 * z_inv2;
	// const auto H           = numerator / denominator;
	// return std::abs(H);
	const auto a00 = 2.f * coeffs.a1 - 1.f;
	const auto a01 = -2.f * coeffs.a2;
	const auto a10 = 2.f * coeffs.a2;
	const auto a11 = 1.f - 2.f * coeffs.a3;

	const auto b0 = 2.f * coeffs.a2;
	const auto b1 = 2.f * coeffs.a3;

	const auto c0 =
		coeffs.a1 * coeffs.m1
		+ coeffs.a2 * coeffs.m2;

	const auto c1 =
		-coeffs.a2 * coeffs.m1
		- coeffs.a3 * coeffs.m2
		+ coeffs.m2;

	const auto d =
		coeffs.a2 * coeffs.m1
		+ coeffs.a3 * coeffs.m2
		+ coeffs.m0;

	return eq::transfer_2pole(
		a00, a01, a10, a11,
		b0, b1,
		c0, c1,
		d,
		omega);
}

auto transfer(const shelf_hi::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert(in_omega.size() == out_magnitude.size());
	auto fn_transfer = [&coeffs](float omega) -> float { return transfer(coeffs, omega); };
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(), fn_transfer);
}

} // eq::filters::shelf_hi

namespace eq::filters::shelf_lo {

template <typename T> struct coeffs { T a1, a2, a3, m1, m2; };
template <typename T> struct state  { T ic1eq, ic2eq; };

template <typename T> [[nodiscard]]
auto make_coeffs(T omega, T k, T A) -> shelf_lo::coeffs<T> {
	const auto pi_omega = ml::kPi * clamp_omega(omega);
	const auto g        = tan(pi_omega) / sqrt(A);
	const auto ca1      = 1.f / (1.f + g * (g + k));
	const auto ca2      = g * ca1;
	const auto ca3      = g * ca2;
	const auto cm1      = k * (A - 1.f);
	const auto cm2      = A * A - 1.f;
	return {ca1, ca2, ca3, cm1, cm2};
}

template <typename T> [[nodiscard]]
auto process(shelf_lo::state<T> state, T in, const shelf_lo::coeffs<T>& coeffs) -> std::tuple<shelf_lo::state<T>, T> {
	const auto v3 = in - state.ic2eq;
	const auto v1 = coeffs.a1 * state.ic1eq + coeffs.a2 * v3;
	const auto v2 = state.ic2eq + coeffs.a2 * state.ic1eq + coeffs.a3 * v3;
	state.ic1eq = 2.f * v1 - state.ic1eq;
	state.ic2eq = 2.f * v2 - state.ic2eq;
	return std::make_tuple(state, in + coeffs.m1 * v1 + coeffs.m2 * v2);
}

[[nodiscard]]
auto transfer(const shelf_lo::coeffs<float>& coeffs, float omega) -> float {
	// static constexpr auto ONE = std::complex<float>{1.f, 0.f};
	// const auto w           = ml::kPi * 2.f * omega;
	// const auto z           = std::exp(std::complex<float>(0.f, w));
	// const auto z_inv       = 1.f / z;
	// const auto z_inv2      = z_inv * z_inv;
	// const auto numerator   = ONE + coeffs.m1 * coeffs.a2 * (ONE + coeffs.a2 * z_inv) + coeffs.m2 * (ONE + 2.f * coeffs.a2 * z_inv + coeffs.a2 * coeffs.a3 * z_inv2);
	// const auto denominator = ONE + (coeffs.a2 + coeffs.a3) * z_inv + coeffs.a3 * z_inv2;
	// const auto H           = numerator / denominator;
	// return std::abs(H);
    const auto a00 = 2.f * coeffs.a1 - 1.f;
	const auto a01 = -2.f * coeffs.a2;
	const auto a10 = 2.f * coeffs.a2;
	const auto a11 = 1.f - 2.f * coeffs.a3;

	const auto b0 = 2.f * coeffs.a2;
	const auto b1 = 2.f * coeffs.a3;

	const auto c0 =
		coeffs.a1 * coeffs.m1
		+ coeffs.a2 * coeffs.m2;

	const auto c1 =
		-coeffs.a2 * coeffs.m1
		- coeffs.a3 * coeffs.m2
		+ coeffs.m2;

	const auto d =
		1.f
		+ coeffs.a2 * coeffs.m1
		+ coeffs.a3 * coeffs.m2;

	return eq::transfer_2pole(
		a00, a01, a10, a11,
		b0, b1,
		c0, c1,
		d,
		omega);
}

auto transfer(const shelf_lo::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert(in_omega.size() == out_magnitude.size());
	auto fn_transfer = [&coeffs](float omega) -> float { return transfer(coeffs, omega); };
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(), fn_transfer);
}

} // eq::filters::shelf_lo
