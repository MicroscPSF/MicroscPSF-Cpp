#pragma once
#include <armadillo>
#include <cstdint>

#include "units.h"

namespace microscPSF {

using Micron = ::units::Micrometer<double>;
using Meter = ::units::Meter<double>;
using namespace ::units::literals;

struct params_li2017_t {
    Micron ti0 = 150.0_um;     //!< Expected immersion medium thickness
    float ni0 = 1.5;           //!< Expected immersion medium refractive index
    float ni = 1.5;            //!< Measured immersion medium refractive index
    Micron tg0 = 170.0_um;     //!< Expected coverglass thickness
    Micron tg = 170.0_um;      //!< Measured coverglass thickness
    float ng0 = 1.5;           //!< Expected coverglass refractive index
    float ng = 1.5;            //!< Measured coverglass refractive index
    float ns = 1.33;           //!< Sample refractive index
    Micron lambda = 0.610_um;  //!< Emission wavelength
    float NA = 1.4;            //!< Numerical aperture

    Micron pz = 2.0_um;  //!< Particle z distance from the sample-to-coverglass interface
};

struct precision_li2017_t {
    int sf = 2;                //!< Oversample factor
    uint32_t rho_samples = 1000;
    Meter res_axial = 0.25_um;
    Meter res_lateral = 0.1_um;
    uint32_t num_basis = 100;
};

template <typename T>
struct pair_t {
    T x{};
    T z{};
};

arma::Cube<double> makePSF(params_li2017_t, pair_t<Micron> voxel, pair_t<int32_t> volume,
                           precision_li2017_t);
}  // namespace microscPSF