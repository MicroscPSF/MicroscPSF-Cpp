#pragma once
#include <cmath>

namespace units {

template <typename Float>
struct Meter {
    Float value{};

    Meter(Float) = delete;

    constexpr operator double() const { return value; }
};

template <typename Float>
struct Micrometer {
    Float value{};

    Micrometer(Float) = delete;

    constexpr operator Meter<double>() const { return Meter<double>{value * 1e-6}; }
};

constexpr Meter<double>
operator*(const Meter<double>& a, double factor) {
    return {static_cast<double>(a) * factor};
}

constexpr double
operator/(const Meter<double>& a, const Meter<double>& b) {
    return a.value / b.value;
}

constexpr Micrometer<double>
operator*(double factor, const Micrometer<double>& a) {
    return {a.value * factor};
}

namespace literals {

constexpr ::units::Meter<double>
operator""_m(long double v) {
    return {static_cast<double>(v)};
}

constexpr ::units::Micrometer<double>
operator""_um(long double v) {
    return {static_cast<double>(v)};
}

static_assert(std::abs(Meter<double>(1.5_um) - 1.5e-6) < 1e-6, "Micrometer value out of scale");
static_assert(std::abs((1.5_m).value - 1.5) < 1e-6, "Micrometer value out of scale");

}  // namespace literals
}  // namespace units
