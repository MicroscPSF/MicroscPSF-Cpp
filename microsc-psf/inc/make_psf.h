#pragma once
#include <armadillo>
#include <cstdint>

#include "units.h"

namespace microscPSF {

using Micron = ::units::Micrometer<double>;
using Meter = ::units::Meter<double>;
using ::units::literals::operator""_um;

struct microscope_params_t {
    Micron ti0 = 150.0_um;     //!< Expected immersion medium thickness
    float ni0 = 1.5;           //!< Expected immersion medium refractive index
    float ni = 1.5;            //!< Measured immersion medium refractive index
    Micron tg0 = 170.0_um;     //!< Expected coverglass thickness
    Micron tg = 170.0_um;      //!< Measured coverglass thickness
    float ng0 = 1.5;           //!< Expected coverglass refractive index
    float ng = 1.5;            //!< Measured coverglass refractive index
    float ns = 1.33;           //!< Sample refractive index
    float NA = 1.4;            //!< Numerical aperture

    Micron pz = 2.0_um;  //!< Particle z distance from the sample-to-coverglass interface
};

/** Solver of linear equations */
enum linear_solver_option_t {
    // Sanderson and Curtin 2020, an adaptive solver for systems of linear equations.
    // https://arma.sourceforge.net/armadillo_solver_2020.pdf
    SandersonAndCurtin2020,

    // Penrose pseudo inverse.
    PenroseInverse
};

struct precision_li2017_t {
    int over_sampling = 2;                //!< Oversample factor
    uint32_t rho_samples = 1000;
    uint32_t num_basis = 100;
    Micron min_wavelength = 0.436_um;
    linear_solver_option_t solver = SandersonAndCurtin2020;
};

template <typename T>
struct pair_t {
    T x{};
    T z{};
};

arma::Cube<double> makePSF(microscope_params_t, pair_t<Micron> voxel, pair_t<int32_t> volume,
                           Micron wavelength = 0.530_um, precision_li2017_t = {});
}  // namespace microscPSF