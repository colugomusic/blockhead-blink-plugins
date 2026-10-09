#pragma once

#include <cassert>
#include <complex>
#include <cmath>
#include <DSP/MLDSPOps.h>
#include <span>
#include <algorithm>
#include <array>

namespace eq::filters {

static constexpr int MAX_STAGES = 8;  // supports up to 96 dB/oct (8 × 12 dB)

template <typename T>
struct biquad_coeffs { T b0, b1, b2, a1, a2; };

template <typename T>
struct biquad_state  { T x1 = 0, x2 = 0, y1 = 0, y2 = 0; };

template <typename T> [[nodiscard]]
auto process_biquad(biquad_state<T> s, T in, const biquad_coeffs<T>& c) -> std::tuple<biquad_state<T>, T> {
	const auto out = c.b0*in + c.b1*s.x1 + c.b2*s.x2 - c.a1*s.y1 - c.a2*s.y2;
	s.x2 = s.x1; s.x1 = in;
	s.y2 = s.y1; s.y1 = out;
	return {s, out};
}

[[nodiscard]]
auto transfer_biquad(const biquad_coeffs<float>& c, float omega) -> std::complex<float> {
	const auto z  = std::polar(1.f, -ml::kPi * omega);
	const auto z2 = z * z;
	return (c.b0 + c.b1*z + c.b2*z2) / (1.f + c.a1*z + c.a2*z2);
}

template <typename T>
struct multi_coeffs {
	int n = 1;
	std::array<biquad_coeffs<T>, MAX_STAGES> stages{};
};

template <typename T>
struct multi_state {
	std::array<biquad_state<T>, MAX_STAGES> stages{};
};

template <typename T> [[nodiscard]]
auto process_multi(multi_state<T> s, T in, const multi_coeffs<T>& c) -> std::tuple<multi_state<T>, T> {
	T sig = in;
	for (int i = 0; i < c.n; ++i) {
		auto [ns, out] = process_biquad(s.stages[i], sig, c.stages[i]);
		s.stages[i] = ns;
		sig = out;
	}
	return {s, sig};
}

[[nodiscard]]
auto transfer_multi(const multi_coeffs<float>& c, float omega) -> float {
	auto H = std::complex<float>{1.f, 0.f};
	for (int i = 0; i < c.n; ++i)
		H *= transfer_biquad(c.stages[i], omega);
	return std::abs(H);
}

} // eq::filters

namespace eq::filters::bell {

template <typename T>
struct coeffs {
	int n        = 1;
	T   omega_c  = T(0.5);  // centre frequency (for analytical transfer)
	T   q        = T(0.7);  // Q (for analytical transfer)
	T   A        = T(1);    // sqrt of linear gain, i.e. 10^(dBgain/40)
	std::array<filters::biquad_coeffs<T>, filters::MAX_STAGES> stages{};  // for process()
};

template <typename T> using state = filters::multi_state<T>;

template <typename T> [[nodiscard]]
auto make_coeffs(T omega, T q, T A, int order = 1) -> bell::coeffs<T> {
	auto c    = bell::coeffs<T>{};
	c.n       = std::clamp(order, 1, filters::MAX_STAGES);
	c.omega_c = omega;
	c.q       = q;
	c.A       = A;
	// Design Butterworth-distributed peaking EQ stages for audio processing.
	const auto w0  = T(ml::kPi) * omega;
	const auto cw  = std::cos(w0);
	const auto sw  = std::sin(w0);
	const auto As  = std::pow(A, T(1) / T(c.n));
	for (int i = 0; i < c.n; ++i) {
		const auto q_bw  = T(1) / (T(2) * std::cos(T(ml::kPi) * T(2*i+1) / T(4*c.n)));
		const auto q_eff = q_bw * q / T(0.7071067811865476);
		const auto al    = sw / (T(2) * q_eff);
		const auto a0    = T(1) + al / As;
		c.stages[i] = { (T(1)+al*As)/a0, -T(2)*cw/a0, (T(1)-al*As)/a0, -T(2)*cw/a0, (T(1)-al/As)/a0 };
	}
	return c;
}

template <typename T> [[nodiscard]]
auto process(bell::state<T> s, T in, const bell::coeffs<T>& c) -> std::tuple<bell::state<T>, T> {
	T sig = in;
	for (int i = 0; i < c.n; ++i) {
		auto [ns, out] = filters::process_biquad(s.stages[i], sig, c.stages[i]);
		s.stages[i] = ns;
		sig = out;
	}
	return {s, sig};
}

// Transfer function uses the analytical Nth-order Butterworth bandpass shape:
//   H(ω) = 1 + (A²−1) / (1 + u^{2N}),  u = Q·(f/f0 − f0/f)
// Using the ratio-form u keeps the bell symmetric in log-frequency (octaves),
// and higher order flattens the top while steepening the skirts.
[[nodiscard]]
auto transfer(const bell::coeffs<float>& c, float omega) -> float {
	const float f    = std::clamp(omega, 1.0e-6f, 0.999999f);
	const float f0   = std::clamp(c.omega_c, 1.0e-6f, 0.999999f);
	const float u    = c.q * ((f / f0) - (f0 / f));
	const float u2n  = std::pow(u * u, static_cast<float>(c.n));
	const float g_bp = 1.0f / (1.0f + u2n);
	return 1.0f + (c.A * c.A - 1.0f) * g_bp;
}

auto transfer(const bell::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert(in_omega.size() == out_magnitude.size());
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(),
		[&coeffs](float omega) { return transfer(coeffs, omega); });
}

} // eq::filters::bell

namespace eq::filters::pass_hi {

template <typename T>
struct coeffs {
	int n       = 1;
	T   omega_c = T(0.5);
	T   k       = T(0.7071067811865476);
	std::array<filters::biquad_coeffs<T>, filters::MAX_STAGES> stages{};
};

template <typename T> using state  = filters::multi_state<T>;

template <typename T> [[nodiscard]]
auto make_coeffs(T omega, T k, int order = 1) -> pass_hi::coeffs<T> {
	const auto w0  = ml::kPi * omega;
	const auto cw  = std::cos(w0);
	const auto sw  = std::sin(w0);
	auto c   = pass_hi::coeffs<T>{};
	c.n      = std::clamp(order, 1, filters::MAX_STAGES);
	c.omega_c = omega;
	c.k       = k;
	for (int i = 0; i < c.n; ++i) {
		const auto q_bw = T(1) / (T(2) * std::cos(T(ml::kPi) * T(2*i + 1) / T(4 * c.n)));
		const auto q    = q_bw * k / T(0.7071067811865476);
		const auto al   = sw / (T(2) * q);
		const auto a0   = T(1) + al;
		c.stages[i]     = { (T(1)+cw)/(T(2)*a0), -(T(1)+cw)/a0, (T(1)+cw)/(T(2)*a0), -T(2)*cw/a0, (T(1)-al)/a0 };
	}
	return c;
}

template <typename T> [[nodiscard]]
auto process(pass_hi::state<T> state, T in, const pass_hi::coeffs<T>& coeffs) -> std::tuple<pass_hi::state<T>, T> {
	T sig = in;
	for (int i = 0; i < coeffs.n; ++i) {
		auto [ns, out] = filters::process_biquad(state.stages[i], sig, coeffs.stages[i]);
		state.stages[i] = ns;
		sig = out;
	}
	return {state, sig};
}

[[nodiscard]]
auto transfer(const pass_hi::coeffs<float>& coeffs, float omega) -> float {
	const auto f   = std::clamp(omega, 1.0e-6f, 0.999999f);
	const auto f0  = std::clamp(coeffs.omega_c, 1.0e-6f, 0.999999f);
	const auto k   = std::max(coeffs.k, 0.05f);
	const auto p   = (2.0f * static_cast<float>(coeffs.n)) / k;
	const auto u2n = std::pow(f / f0, p);
	return std::sqrt(u2n / (1.0f + u2n));
}

auto transfer(const pass_hi::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert(in_omega.size() == out_magnitude.size());
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(),
		[&coeffs](float omega) { return transfer(coeffs, omega); });
}

} // eq::filters::pass_hi

namespace eq::filters::pass_lo {

template <typename T>
struct coeffs {
	int n       = 1;
	T   omega_c = T(0.5);
	T   k       = T(0.7071067811865476);
	std::array<filters::biquad_coeffs<T>, filters::MAX_STAGES> stages{};
};

template <typename T> using state  = filters::multi_state<T>;

template <typename T> [[nodiscard]]
auto make_coeffs(T omega, T k, int order = 1) -> pass_lo::coeffs<T> {
	const auto w0  = ml::kPi * omega;
	const auto cw  = std::cos(w0);
	const auto sw  = std::sin(w0);
	auto c   = pass_lo::coeffs<T>{};
	c.n      = std::clamp(order, 1, filters::MAX_STAGES);
	c.omega_c = omega;
	c.k       = k;
	for (int i = 0; i < c.n; ++i) {
		const auto q_bw = T(1) / (T(2) * std::cos(T(ml::kPi) * T(2*i + 1) / T(4 * c.n)));
		const auto q    = q_bw * k / T(0.7071067811865476);
		const auto al   = sw / (T(2) * q);
		const auto a0   = T(1) + al;
		c.stages[i]     = { (T(1)-cw)/(T(2)*a0), (T(1)-cw)/a0, (T(1)-cw)/(T(2)*a0), -T(2)*cw/a0, (T(1)-al)/a0 };
	}
	return c;
}

template <typename T> [[nodiscard]]
auto process(pass_lo::state<T> state, T in, const pass_lo::coeffs<T>& coeffs) -> std::tuple<pass_lo::state<T>, T> {
	T sig = in;
	for (int i = 0; i < coeffs.n; ++i) {
		auto [ns, out] = filters::process_biquad(state.stages[i], sig, coeffs.stages[i]);
		state.stages[i] = ns;
		sig = out;
	}
	return {state, sig};
}

[[nodiscard]]
auto transfer(const pass_lo::coeffs<float>& coeffs, float omega) -> float {
	const auto f   = std::clamp(omega, 1.0e-6f, 0.999999f);
	const auto f0  = std::clamp(coeffs.omega_c, 1.0e-6f, 0.999999f);
	const auto k   = std::max(coeffs.k, 0.05f);
	const auto p   = (2.0f * static_cast<float>(coeffs.n)) / k;
	const auto u2n = std::pow(f / f0, p);
	return 1.0f / std::sqrt(1.0f + u2n);
}

auto transfer(const pass_lo::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert(in_omega.size() == out_magnitude.size());
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(),
		[&coeffs](float omega) { return transfer(coeffs, omega); });
}

} // eq::filters::pass_lo

namespace eq::filters::shelf_hi {

template <typename T> using coeffs = filters::multi_coeffs<T>;
template <typename T> using state  = filters::multi_state<T>;

template <typename T> [[nodiscard]]
auto make_coeffs(T omega, T k, T A, int order = 1) -> shelf_hi::coeffs<T> {
	const auto w0  = ml::kPi * omega;
	const auto cw  = std::cos(w0);
	const auto sw  = std::sin(w0);
	auto c   = shelf_hi::coeffs<T>{};
	c.n      = std::clamp(order, 1, filters::MAX_STAGES);
	const auto As  = std::pow(A, T(1) / T(c.n));
	const auto sAs = std::sqrt(As);
	for (int i = 0; i < c.n; ++i) {
		const auto q_bw = T(1) / (T(2) * std::cos(T(ml::kPi) * T(2*i + 1) / T(4 * c.n)));
		const auto q    = q_bw * k / T(0.7071067811865476);
		const auto al   = sw / (T(2) * q);
		const auto a0   =       (As+T(1)) - (As-T(1))*cw + T(2)*sAs*al;
		c.stages[i].b0  =  As*((As+T(1)) + (As-T(1))*cw + T(2)*sAs*al) / a0;
		c.stages[i].b1  = -T(2)*As*((As-T(1)) + (As+T(1))*cw)          / a0;
		c.stages[i].b2  =  As*((As+T(1)) + (As-T(1))*cw - T(2)*sAs*al) / a0;
		c.stages[i].a1  =  T(2)*((As-T(1)) - (As+T(1))*cw)             / a0;
		c.stages[i].a2  =      ((As+T(1)) - (As-T(1))*cw - T(2)*sAs*al) / a0;
	}
	return c;
}

template <typename T> [[nodiscard]]
auto process(shelf_hi::state<T> state, T in, const shelf_hi::coeffs<T>& coeffs) -> std::tuple<shelf_hi::state<T>, T> {
	return filters::process_multi(state, in, coeffs);
}

[[nodiscard]]
auto transfer(const shelf_hi::coeffs<float>& coeffs, float omega) -> float {
	return filters::transfer_multi(coeffs, omega);
}

auto transfer(const shelf_hi::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert(in_omega.size() == out_magnitude.size());
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(),
		[&coeffs](float omega) { return transfer(coeffs, omega); });
}

} // eq::filters::shelf_hi

namespace eq::filters::shelf_lo {

template <typename T> using coeffs = filters::multi_coeffs<T>;
template <typename T> using state  = filters::multi_state<T>;

template <typename T> [[nodiscard]]
auto make_coeffs(T omega, T k, T A, int order = 1) -> shelf_lo::coeffs<T> {
	const auto w0  = ml::kPi * omega;
	const auto cw  = std::cos(w0);
	const auto sw  = std::sin(w0);
	auto c   = shelf_lo::coeffs<T>{};
	c.n      = std::clamp(order, 1, filters::MAX_STAGES);
	const auto As  = std::pow(A, T(1) / T(c.n));
	const auto sAs = std::sqrt(As);
	for (int i = 0; i < c.n; ++i) {
		const auto q_bw = T(1) / (T(2) * std::cos(T(ml::kPi) * T(2*i + 1) / T(4 * c.n)));
		const auto q    = q_bw * k / T(0.7071067811865476);
		const auto al   = sw / (T(2) * q);
		const auto a0   =       (As+T(1)) + (As-T(1))*cw + T(2)*sAs*al;
		c.stages[i].b0  =  As*((As+T(1)) - (As-T(1))*cw + T(2)*sAs*al) / a0;
		c.stages[i].b1  =  T(2)*As*((As-T(1)) - (As+T(1))*cw)          / a0;
		c.stages[i].b2  =  As*((As+T(1)) - (As-T(1))*cw - T(2)*sAs*al) / a0;
		c.stages[i].a1  = -T(2)*((As-T(1)) + (As+T(1))*cw)             / a0;
		c.stages[i].a2  =      ((As+T(1)) + (As-T(1))*cw - T(2)*sAs*al) / a0;
	}
	return c;
}

template <typename T> [[nodiscard]]
auto process(shelf_lo::state<T> state, T in, const shelf_lo::coeffs<T>& coeffs) -> std::tuple<shelf_lo::state<T>, T> {
	return filters::process_multi(state, in, coeffs);
}

[[nodiscard]]
auto transfer(const shelf_lo::coeffs<float>& coeffs, float omega) -> float {
	return filters::transfer_multi(coeffs, omega);
}

auto transfer(const shelf_lo::coeffs<float>& coeffs, std::span<const float> in_omega, std::span<float> out_magnitude) -> void {
	assert(in_omega.size() == out_magnitude.size());
	std::transform(in_omega.begin(), in_omega.end(), out_magnitude.begin(),
		[&coeffs](float omega) { return transfer(coeffs, omega); });
}

} // eq::filters::shelf_lo
