#pragma once

#include <armadillo>

namespace microsc_psf {
namespace internal {

/** Solve min || A * x - b ||^2 for x , using Eigen3 solver.
 *
 * Typical ussage: simulate PSF on machines not having LAPACK support, e.g.
 * Windows on ARM64 CPU.
 *
 * Note: This may seem contrived to have two near-identical C++ libraries
 * serving the linear algebra code, but Armadillo has a nice Matlab-like syntax
 * that simplifies the code review process a lot.
 *
 * @tparam tranpose_b True if the signal b needs a Hermitian transpose ahead of
 * the least squares solver.
 */
template <bool transpose_b>
arma::cx_mat solveWithEigen(const arma::mat& A, const arma::cx_mat& b);
}  // namespace internal
}  // namespace microsc_psf