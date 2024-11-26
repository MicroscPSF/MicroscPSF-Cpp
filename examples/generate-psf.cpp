#include <armadillo>

#include "make_psf.h"

int main() {
    using arma::hdf5_name;
    using microsc_psf::makePSF;
    using microsc_psf::microscope_params_t;
    using microsc_psf::precision_li2017_t;
    namespace hdf5_opts = arma::hdf5_opts;
    using namespace ::units::literals;

    microscope_params_t params{};
    params.NA = 1.4;
    params.ti0 = 150.0_um;
    params.ni = 1.5;
    params.ni0 = 1.5;
    params.pz = 2.0_um;

    precision_li2017_t precision{};
    precision.num_basis = 153;
    precision.rho_samples = 1000;

    const auto psf =
        makePSF(params, {0.1_um, 0.25_um}, {256, 128}, 0.610_um, precision);
#ifdef  ARMA_USE_HDF5
    std::cout << "Saving volume to HDF5...\n";
    psf.save(hdf5_name("psf.h5", "psf", hdf5_opts::trans));
    std::cout << R"(Done.
Import to Matlab: h5read("psf.h5", "psf");

Import to Python:
    import h5py
    with h5py.File("psf.h5", "r") as h5file:
        psf = h5file["psf"][()]

)";
#endif

    {
        using namespace arma;
        mat xy_plane = sqrt(psf.slice(32));
        xy_plane *= 255 / max(xy_plane.as_col());
        xy_plane.save("psf_xy.pgm", pgm_binary);
        std::cout << "PSF XY plane saved to 'psf_xy.pgm'.\n";
    }

    {
        using namespace arma;
        mat xz_plane = sqrt(psf.col(60));
        mat(xz_plane.t() * 255 / max(xz_plane.as_col())).save("psf_xz.pgm", pgm_binary);
        std::cout << "PSF XZ plane saved to 'psf_xz.pgm'.\n";
    }
    return 0;
}