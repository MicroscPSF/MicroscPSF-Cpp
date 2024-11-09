#include "make_psf.h"

int main() {
    using microscPSF::makePSF;
    using microscPSF::params_li2017_t;
    using microscPSF::precision_li2017_t;
    using arma::hdf5_name;
    namespace hdf5_opts = arma::hdf5_opts;
    using namespace ::units::literals;

    const auto psf = makePSF(params_li2017_t{}, {0.1_um, 0.25_um}, {256, 128}, precision_li2017_t{});
    psf.save(hdf5_name("psf.h5", "psf", hdf5_opts::trans));

    return 0;
}