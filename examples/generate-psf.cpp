#include <armadillo>

#include "make_psf.h"

int main() {
    using microscPSF::makePSF;
    using microscPSF::microscope_params_t;
    using microscPSF::precision_li2017_t;
    using arma::hdf5_name;
    namespace hdf5_opts = arma::hdf5_opts;
    using namespace ::units::literals;

    const auto psf =
        makePSF(microscope_params_t{}, {0.1_um, 0.25_um}, {120, 63}, 0.530_um, precision_li2017_t{});
#ifdef ARMA_USE_HDF5
    psf.save(hdf5_name("psf.h5", "psf", hdf5_opts::trans));
#endif

    {
        using namespace arma;
        mat xy_plane = sqrt(psf.slice(32));
        xy_plane *= 255 / max(xy_plane.as_col());
        xy_plane.save("psf_xy.pgm", pgm_binary);
    }

    {
        using namespace arma;
        mat xz_plane = sqrt(psf.col(60));
        mat(xz_plane.t() * 255 / max(xz_plane.as_col())).save("psf_xz.pgm", pgm_binary);
    }
    return 0;
}