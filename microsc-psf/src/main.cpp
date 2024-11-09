#include <cassert>

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
    XX.each_row() = x.t();

    arma::mat YY(x.n_elem, y.n_elem);
    YY.each_col() = y;
    return {XX, YY};
}

/** Piecewise linear interpolation. */
arma::Cube<double>
cylToRectTransform(const arma::mat& PSF0, const arma::vec& R, double over_sample,
                   microscPSF::pair_t<int32_t> volume) {
    using namespace arma;
    assert(PSF0.n_cols == uint32_t(volume.z));

    const double x0 = (volume.x - 1) / 2;
    const double y0 = x0;

    const auto [X, Y] = meshgrid(iota(volume.x) - x0, iota(volume.x) - y0);

    const vec rPixel = sqrt(pow(X, 2) + pow(Y, 2)).as_col();
    const uvec index = conv_to<umat>::from(rPixel * over_sample);

#define disR (rPixel - R(index)) * over_sample

    Cube<double> PSF(volume.x, volume.x, volume.z);
    for (uint32_t zi = 0; zi < uint32_t(volume.z); zi++) {
        const mat h = PSF0.col(zi);
        PSF.slice(zi) =
            resize(h(index + 1) % disR + h(index) % (1.0 - disR), SizeMat{X.n_rows, X.n_cols});
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

    const double a = 0.0;
    const float b = std::min({
        1.0f,                    //
        params.ns / params.NA,   //
        params.ni0 / params.NA,  //
        params.ni / params.NA,   //
        params.ng0 / params.NA,  //
        params.ng / params.NA    //
    });

    const auto L = precision.num_samp;

    const rowvec Rho = linspace(a, b, L).t();

    // Approximate function exp(j omega) as  Bessel series
    using ::units::literals::operator""_m;

    const auto k0 = datum::pi * 2.0 * (1.0_m / Meter(params.lambda));

    const vec R = iota(params.sf * max_radius) / params.sf;

    vec A, Ab;
    {
        const auto r = R * precision.res_lateral;
        A = k0 * params.NA * r;
        Ab = pow(A, 2) * b;
    }

    vec an;
    {
        const auto NN = precision.num_basis;

        // Min wavelength
        const double k00 = datum::pi * 2 / (1.0_m / 545e-9_m);
        const double factor1 = k0 / k00;

        // Max NA
        const double factor2 = params.NA / 1.4;

        an = (iota(1.0, NN + 1) * 3 - 2) * factor1 * factor2;
    }

    mat Ele;
    {
        // bsxfun(@minus, an2', A2);
        mat domin(A.n_elem, an.n_elem);
        domin.each_row() = pow(an.t(), 2);
        domin.each_col() -= pow(A, 2);

        Ele = (                                                         //
                  besselj<0>(Ab) * (an.t() % besselj<1>(an * b).t()) -  //
                  (A % besselj<1>(Ab)) * besselj<0>(an * b).t()         //
                  ) *
              b / domin;
    }

    cx_mat Ffun;
    {
#define Ti (Meter(params.ti0) + precision.res_axial*(iota(volume.z) - z0))

        const double C1 = params.ns * Meter(params.pz);
#define C2 (params.ni * (Ti - Meter(params.ti0)))
        const double C3 = params.ng * (Meter(params.tg) - Meter(params.tg0));

#define OPDs (C1 * sqrt(1.0 - pow(params.NA * Rho / params.ns, 2)))
#define OPDi (C2 * sqrt(1.0 - pow(params.NA * Rho / params.ni, 2)))
#define OPDg (C3 * sqrt(1 - pow(params.NA * Rho / params.ng, 2)))

        // bsxfun(plus, OPDi, OPDs + OPDg)
        mat OPD = OPDi;
        OPD.each_row() += OPDs + OPDg;

        // Determine the coefficients
        constexpr auto j = cx_double{0.0, 1.0};

        const auto W = k0 * OPD;
        Ffun = exp(j * W);
    }

    mat PSF0;
    {
        const mat J = besselj<0>(an * Rho);

        // Armadillo does not have solver for complex valued vector.
        #define Ci (pinv(J.t()) * Ffun.t())

        const cx_mat ciEle = Ele * Ci;
        PSF0 = real(ciEle % conj(ciEle));

    }
    std::cout << "PSF0 = " << PSF0.n_rows << ',' << PSF0.n_cols << std::endl;

    auto PSF = cylToRectTransform(PSF0, R, params.sf, volume);

    // Normalize the intensity.
    PSF /= accu(PSF);

    return PSF;
}

}  // namespace microscPSF