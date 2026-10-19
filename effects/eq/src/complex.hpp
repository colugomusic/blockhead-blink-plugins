/* @TODO: delete
#pragma once

#include <algorithm>
#include <cmath>
#include <DSP/MLDSPOps.h>

namespace eq {

template <size_t N>
auto acos(ml::DSPVectorArray<N> v) -> ml::DSPVectorArray<N> {
	for (size_t r = 0; r < N; r++) {
		auto& row = v.row(r);
		for (size_t i = 0; i < kFloatsPerDSPVector; i++) {
			row[i] = std::acos(row[i]);
		}
	}
	return v;
}

auto acos(float v) -> float         { return std::acos(v); }
auto cos(float v) -> float          { return std::cos(v); }
auto pow(float v, float n) -> float { return std::pow(v, n); }
auto sqrt(float v) -> float         { return std::sqrt(v); }
auto min(float a, float b) -> float { return std::min(a, b); }
auto max(float a, float b) -> float { return std::max(a, b); }

template <typename T>
struct complex {
	T r = T{0.0};
	T i = T{0.0};
};

template <typename T>
auto operator+(complex<T> a, complex<T> b) -> complex<T> {
	return {
		.r = a.r + b.r,
		.i = a.i + b.i
	};
}

template <typename T>
auto operator*(complex<T> a, complex<T> b) -> complex<T> {
	return {
		.r = (a.r * b.r) - (a.i * b.i),
		.i = (a.r * b.i) + (a.i * b.r)
	};
}

template <typename T>
auto operator+(complex<T> a, T b) -> complex<T> {
	return {
		.r = a.r + b,
		.i = a.i
	};
}

template <typename T>
auto operator*(complex<T> a, T b) -> complex<T> {
	return {
		.r = a.r * b,
		.i = a.i * b
	};
}

template <typename T>
auto operator/(T a, complex<T> b) -> complex<T> {
	const auto denom = (b.r * b.r) + (b.i * b.i);
	return {
		.r = (a * b.r) / denom,
		.i = (-a * b.i) / denom
	};
}

template <typename T>
auto abs(complex<T> a) -> T {
	return sqrt((a.r * a.r) + (a.i * a.i));
}

template <typename T> auto operator*=(complex<T>& a, complex<T> b) -> complex<T>& { a = a * b; return a; }
template <typename T> auto operator*=(complex<T>& a, T b) -> complex<T>&          { a = a * b; return a; }

} // eq
*/
