#pragma once
#include <armadillo>
#include <cstdint>

#include "units.h"

namespace microscPSF {

using Micron = ::units::Micrometer<double>;
using Meter = ::units::Meter<double>;
using namespace ::units::literals;

struct params_li2017_t {
    Micron ti0 = 150.0_um;
    float ni0 = 1.5;
    float ni = 1.5;
    Micron tg0 = 170.0_um;
    Micron tg = 170.0_um;
    float ng0 = 1.5;
    float ng = 1.5;
    float ns = 1.33;
    Micron lambda = 0.610_um;
    float NA = 1.4;
    int sf = 2;
    int mode = 1;

    Micron pz = 2.0_um;
};

struct precision_li2017_t {
    uint32_t num_samp = 1000;
    Meter res_axial = 250e-9_m;
    Meter res_lateral = 100e-9_m;
    uint8_t num_basis = 100;
};

template <typename T>
struct pair_t {
    T x{};
    T z{};
};

arma::Cube<double> makePSF(params_li2017_t, pair_t<Micron> voxel, pair_t<int32_t> volume,
                           precision_li2017_t);
}  // namespace microscPSF