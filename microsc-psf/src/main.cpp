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

/** Piecewise linear interpolation. */
arma::Cube<double>
cylToRectTransform(const arma::mat& PSF0, const arma::vec& R, microscPSF::pair_t<int32_t> volume) {
    using namespace arma;
    assert(PSF0.n_cols == uint32_t(volume.z));

    vec rPixel;
    {
        const double x0 = (volume.x - 1) / 2;

        // meshgrid(0:nx - x0, 0:ny - y0)
        mat X(volume.x, volume.x);
        X.each_row() = (iota(volume.x) - x0).t();
#define Y X.t()

        rPixel = sqrt(X % X + Y % Y).as_col();
    }

    Cube<double> PSF(volume.x, volume.x, volume.z);
    for (uint32_t zi = 0; zi < uint32_t(volume.z); zi++) {
        // Memory map the PSF slice to an 1D vector without memory copy.
        vec interpolated{&(PSF(0, 0, zi)), uint32_t(volume.x * volume.x), false, true};

        interp1(R, PSF0.col(zi), rPixel, interpolated, "linear");
        // PSF.slice(zi) = resize(interpolated, SizeMat{volume.x, volume.x});
    }

    return PSF;
}
}  // namespace

namespace microscPSF {

arma::Cube<double>
makePSF(params_li2017_t params, pair_t<Micron> voxel, pair_t<int32_t> volume,
        precision_li2017_t precision) {
    using ::units::literals::operator""_m;

    const double x0 = (volume.x - 1) / 2;
    const double y0 = x0;
    const double z0 = (volume.z - 1) / 2;

    using namespace arma;

    // Max radius is the length of the diagonal of the volume xy plane.
    const double max_radius = round(abs(cx_double{volume.x - x0, volume.x - y0})) + 1;

    const float max_rho = std::min({
        1.0f,                    //
        params.ns / params.NA,   //
        params.ni0 / params.NA,  //
        params.ni / params.NA,   //
        params.ng0 / params.NA,  //
        params.ng / params.NA    //
    });

    // Wavenumber of emitted light.
    const auto k0 = datum::pi * 2.0 * (1.0_m / Meter(params.lambda));

    vec scaling_factor;
    {
        constexpr auto min_wavelength = 436e-9_m;
        constexpr auto max_NA = 1.4;

        // Convert min wavelength to max wavenumber
        const double k00 = datum::pi * 2 / (1.0_m / min_wavelength);

        // Scaling factors for the Fourier-Bessel series expansion
        const auto [factor1, factor2] = std::pair{k0 / k00, params.NA / max_NA};

        scaling_factor = (iota(1.0, precision.num_basis + 1) * 3 - 2) * factor1 * factor2;
    }

    const rowvec Rho = linspace<rowvec>(0.0, max_rho, precision.rho_samples);

    // Calculate phase aberration term
    cx_mat phase;
    {
/** Coverglass z offset in meter */
#define Ti (Meter(params.ti0) + precision.res_axial * (iota(volume.z) - z0))

        const double C1 = params.ns * Meter(params.pz);
#define C2 (params.ni * (Ti - Meter(params.ti0)))
        const double C3 = params.ng * (Meter(params.tg) - Meter(params.tg0));

#define OPDsample (C1 * sqrt(clamp(1.0 - square(params.NA * Rho / params.ns), 0.0, 1.0)))
#define OPDimmersion (C2 * sqrt(clamp(1.0 - square(params.NA * Rho / params.ni), 0.0, 1.0)))
#define OPDglass (C3 * sqrt(clamp(1.0 - square(params.NA * Rho / params.ng), 0.0, 1.0)))

        mat OPD = OPDimmersion;
        OPD.each_row() += OPDsample + OPDglass;

        // Sample the optical phase.
        constexpr auto j = cx_double{0.0, 1.0};
        phase = exp(j * k0 * OPD);
    }

    // Radius coordinates in the PSF RZ plane.
    const vec R = iota(precision.sf * max_radius) / precision.sf;

    // Approximate function exp(j omega) as  Bessel series
    // See equation 5 in Li, Xue, and Blu 2017.
    mat Ele;
    {
        const vec A = k0 * params.NA * (R * precision.res_lateral);

        // bsxfun(@minus, an2', A2);
        mat domin(A.n_elem, scaling_factor.n_elem);
        domin.each_row() = pow(scaling_factor.t(), 2);
        domin.each_col() -= pow(A, 2);

        Ele = (  //
                  besselj<0>(A * max_rho) *
                      (scaling_factor.t() % besselj<1>(scaling_factor * max_rho).t()) -     //
                  (A % besselj<1>(A * max_rho)) * besselj<0>(scaling_factor * max_rho).t()  //
                  ) *
              max_rho / domin;
    }

    // Calculate radial G-L at specified radius and z values (cylindrical
    // coordinates). This models the PSF you would measure by scanning the
    // partical relative to the microscope focus.
    mat PSF0;
    {
        // Define the basis of Bessel functions.
        // Shape: number of basis function by number of rho samples.
        auto&& J = besselj<0>(scaling_factor * Rho);

        // Compute the approximation to the sampled pupil phase by finding the least squares
        // solution to the complex coefficients of the Fourier-Bessel expansion.
        // Shape of C is (number of basis functions by number of z samples).
        // Note the matrix transposes to get the dimensions correct.
        //
        // Note: Armadillo does not have solver for real-valued matrix and complex-valued vector.
        // #define Ci (pinv(J.t()) * phase.t());
        cx_mat&& Ci = solve(conv_to<cx_mat>::from(J.t()), phase.t());

        const cx_mat ciEle = Ele * Ci;
        PSF0 = real(ciEle % conj(ciEle));
    }

    {
        const Mat<uint8_t> PSF0_normalized =
            conv_to<Mat<uint8_t>>::from(PSF0 * 255 / max(max(PSF0)));
        PSF0_normalized.save("psf.pgm", pgm_binary);
    }

    auto PSF = cylToRectTransform(PSF0, R, volume);

    // Normalize the intensity.
    PSF /= accu(PSF);

    return PSF;
}

}  // namespace microscPSF