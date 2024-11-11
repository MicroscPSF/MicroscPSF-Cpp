#include <armadillo>
#include <armadillo_besselj_support.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>

namespace {

arma::vec
arange(double start, double stop, double step) {
    const auto count = static_cast<uint32_t>((stop - start) / step);

    arma::vec values(count);
    for (uint32_t i = 0; i < count; i++) {
        values(i) = (stop - start) * step * i + start;
    }

    return values;
}
}  // namespace

TEST_CASE("Besselj(vector)") {
    using namespace arma;

    const auto A = arange(0.0, 30.0, 0.5);

    REQUIRE(vec(besselj<0>(A)).is_finite());
    REQUIRE(vec(besselj<1>(A)).is_finite());
}

TEST_CASE("Besselj(vector + scalar)") {
    using namespace arma;
    const mat A = arange(0.0, 30.0, 0.5);

    REQUIRE(mat(besselj<0>(A + 1.0)).is_finite());
    REQUIRE(mat(besselj<1>(A + 1.0)).is_finite());

    REQUIRE(mat(besselj<0>(A * 1.0)).is_finite());
    REQUIRE(mat(besselj<1>(A * 1.0)).is_finite());
}

TEST_CASE("Besselj(matrix)") {
    using namespace arma;
    const mat A = arange(0.0, 30.0, 0.5);
    const mat M = A * A.t();

    // Always buffer the intermediate matrix.
    REQUIRE(mat(besselj<0>(M)).is_finite());
    REQUIRE(mat(besselj<1>(M)).is_finite());

    // Deferred execution
    REQUIRE(mat(besselj<0>(A * A.t())).is_finite());
    REQUIRE(mat(besselj<1>(A * A.t())).is_finite());
}
