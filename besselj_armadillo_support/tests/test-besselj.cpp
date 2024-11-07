#include <cmath>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEST_CASE("Toolchain has std::cyl_bessel_j") {
    using std::cyl_bessel_j;
    using std::isnan;

    REQUIRE(!isnan(cyl_bessel_j(0.0, 0.0)));
}

TEST_CASE("Besselj spot check") {
    using std::cyl_bessel_j;
    using Catch::Matchers::WithinAbs;
    using Catch::Matchers::WithinRel;

    REQUIRE_THAT(cyl_bessel_j(0.0, 0.0), WithinAbs(1.0, 1e-6));
    REQUIRE_THAT(cyl_bessel_j(1.0, 0.0), WithinAbs(0.0, 1e-6));
}