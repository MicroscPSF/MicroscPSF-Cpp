#include "make_psf.h"

// This must come after header file <armadillo>.
#include <besselj_support.hpp>
namespace {

// Return the range [0, N), excluding N.
arma::vec
iota(double initial_value, uint32_t N) {
    arma::vec v(N);

    std::iota(v.begin(), v.end(), initial_value);
    return v;
}

inline arma::vec
iota(uint32_t N) {
    return iota(0.0, N);
}

arma::vec
linspace(double start, double stop, uint32_t n_samples) {
    return iota(n_samples) * (stop - start) / (n_samples - 1);
}

std::pair<arma::mat, arma::mat>
meshgrid(arma::vec&& x, arma::vec&& y) {
    arma::mat XX(x.n_elem, y.n_elem);
    XX.each_row() = x;

    arma::mat YY(x.n_elem, y.n_elem);
    YY.each_col() = y;
    return {XX, YY};
}

/** Piecewise linear interpolation. */
arma::Cube<double>
cylToRectTransform(const arma::mat& PSF0, const arma::vec& R, double over_sample,
                   microscPSF::pair_t<int32_t> volume) {
    using namespace arma;

    const double x0 = (volume.x - 1) / 2;
    const double y0 = x0;

    const auto [X, Y] = meshgrid(iota(volume.x) - x0, iota(volume.x) - y0);

    const mat rPixel = sqrt(pow(X, 2) + pow(Y, 2));
    const umat index = conv_to<umat>::from(rPixel * over_sample);

    const auto disR = (rPixel - R(index)) * over_sample;

    Cube<double> PSF(volume.x, volume.x, volume.z);
    for (uint32_t zi = 0; zi < uint32_t(volume.z); zi++) {
        const vec h = PSF0.col(zi);
        PSF.slice(zi) = h(index + 1) % disR + h(index) % (1.0 - disR);
    }

    return PSF;
}
}  // namespace

namespace microscPSF {

arma::Cube<double>
makePSF(params_li2017_t params, pair_t<Micron> voxel, pair_t<int32_t> volume,
        precision_li2017_t precision) {
    const double x0 = (volume.x - 1) / 2;
    const double y0 = x0;
    const double z0 = (volume.z - 1) / 2;

    using namespace arma;
    const double max_radius = round(abs(cx_double{volume.x - x0, volume.x - y0})) + 1;

    const vec R = iota(params.sf * max_radius) / params.sf;

    const vec Ti = Meter(params.ti0) + precision.res_axial * (iota(volume.z) - z0);

    const double a = 0.0;
    const auto b = std::min({
        1.0f,                    //
        params.ns / params.NA,   //
        params.ni0 / params.NA,  //
        params.ni / params.NA,   //
        params.ng0 / params.NA,  //
        params.ng / params.NA    //
    });

    const auto L = precision.num_samp;

    const vec Rho = linspace(a, b, L);

    // Approximate function exp(j omega) as  Bessel series
    using ::units::literals::operator""_m;

    const auto NN = precision.num_basis;
    const auto k0 = datum::pi * 2.0 * (1.0_m / Meter(params.lambda));
    const auto r = R * precision.res_lateral;

    const auto A = k0 * params.NA * r;
    const auto Ab = pow(A, 2) * b;

    // Min wavelength
    const auto k00 = datum::pi * 2 / (1.0_m / 545e-9_m);
    const auto factor1 = k0 / k00;

    // Max NA
    const auto factor2 = params.NA / 1.4;

    const auto an = (iota(1.0, NN + 1) * 3 - 2) * factor1 * factor2;

    const auto anRho = an * Rho.t();
    const mat J = besselj<0>(anRho);

    const auto J0A = besselj<0>(Ab);
    const auto J1A = A % besselj<1>(Ab);

    const auto anJ0A = J0A * an.t();

    const auto B0anb = besselj<0>(an * b);
    const auto B1anb = besselj<1>(an * b);

    // bsxfun(@minus, an2', A2);
    const mat domin = pow(an.t(), 2) - pow(A, 2);
    const auto Ele = (anJ0A * B1anb.t() - J1A * B0anb.t()) * b / domin;

    const auto C1 = params.ns * Meter(params.pz);
    const auto C2 = params.ni * (Ti - Meter(params.ti0));
    const auto C3 = params.ng * (Meter(params.tg) - Meter(params.tg0));

    const auto OPDs = C1 * sqrt(1.0 - pow(params.NA * Rho / params.ns, 2));
    const auto OPDi = C2 % sqrt(1.0 - pow(params.NA * Rho / params.ni, 2));
    const auto OPDg = C3 * sqrt(1 - pow(params.NA * Rho / params.ng, 2));

    // bsxfun(plus, OPDi, OPDs + OPDg)
    const auto OPD = OPDi + OPDs + OPDg;

    // Determine the coefficients
    constexpr auto j = cx_double{0.0, 1.0};

    const auto W = k0 * OPD;
    const cx_vec Ffun = cos(W) + j * sin(W);

    // Armadillo does not have solver for complex valued vector.
    const auto Ci = pinv(J) * Ffun;

    // Get PSF for each slice
    const auto ciEle = Ele.t() * Ci;
    const mat PSF0 = real(ciEle % conj(ciEle));

    auto PSF = cylToRectTransform(PSF0, R, params.sf, volume);

    // Normalize the intensity.
    PSF /= accu(PSF);

    return PSF;
}

}  // namespace microscPSF