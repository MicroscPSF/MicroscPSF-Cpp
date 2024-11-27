#include "linsolver.h"

#include <Eigen/Dense>
#include <cassert>

namespace microsc_psf {

namespace internal {

template <bool transpose_b>
arma::cx_mat
solveWithEigen(const arma::mat& A_buffer, const arma::cx_mat& b_buffer) {
    using arma::cx_double;
    using Eigen::ComputeThinU;
    using Eigen::ComputeThinV;
    using Eigen::Map;
    using Eigen::MatrixX;
    using Eigen::MatrixXd;

    const auto m = A_buffer.n_rows;
    const auto n = A_buffer.n_cols;
    const auto k = (transpose_b ? b_buffer.n_rows : b_buffer.n_cols);

    if constexpr (transpose_b) {
        assert(b_buffer.n_cols == m);
    } else {
        assert(b_buffer.n_rows == m);
    }

    // First, map to Eigen's primary data structure.
    const Map<const MatrixXd> A(A_buffer.memptr(), m, n);
    const Map<const MatrixX<cx_double>> b(b_buffer.memptr(), b_buffer.n_rows, b_buffer.n_cols);

    // Allocate the result buffer
    arma::cx_mat xopt_buffer(n, k);

    // Next, solve with Eigen's SVD solver
    Map<Eigen::MatrixX<cx_double>> xopt(xopt_buffer.memptr(), n, k);
    if constexpr (transpose_b) {
        // Eigen may have an efficient Hermitian transposed solver. Move the transpose step into
        // Eigen.
        xopt = A.bdcSvd(ComputeThinU | ComputeThinV).solve(b.transpose());
    } else {
        xopt = A.bdcSvd(ComputeThinU | ComputeThinV).solve(b);
    }

    return xopt_buffer;
}

template arma::cx_mat solveWithEigen<true>(const arma::mat&, const arma::cx_mat&);
template arma::cx_mat solveWithEigen<false>(const arma::mat&, const arma::cx_mat&);
}  // namespace internal
}  // namespace microsc_psf
