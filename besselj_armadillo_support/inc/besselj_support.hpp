/** @file
 * A patch to armadillo-code: C++ matrix library to support Bessel function of
 * the first kind.
 */
#pragma once
#include <cmath>
#include <cstdint>

namespace arma {

// Enable OpenMP multithreading support.
template<uint8_t order>
class eop_besselj : public eop_core<eop_besselj<order>>, public eop_use_mp_true {};

/** Enable deferred execution of the bessel function, so that we can write
 * things like:
 *
 * ```c++
 * arma::vec result = besselj0(a * b) * c + d;
 * ```
 *
 * Note(Antony): We don't have a nice way to explicitly instantiate templates
 * over a range of integers. Help wanted.
 *
 * Reference: https://stackoverflow.com/a/67529191
 */
template <>
template <typename eT>
arma_inline eT
eop_core<eop_besselj<0>>::process(const eT val, const eT) {
    return eop_aux::besselj<0>(val);
}

template <>
template <typename eT>
arma_inline eT
eop_core<eop_besselj<1>>::process(const eT val, const eT) {
    return eop_aux::besselj<1>(val);
}

// Implement element-wise bessel function of a vector or matrix.
template <uint8_t order, typename T1>
arma_warn_unused arma_inline
    typename enable_if2<is_arma_type<T1>::value, const eOp<T1, eop_besselj<order>>>::result
    besselj(const T1& A) {
    arma_debug_sigprint();

    static_assert(order <= 1, "Higher order bessel function not yet tested. Help wanted.");
    return eOp<T1, eop_besselj<order>>(A);
}

// TODO(Antony): Implement element-wise bessel function of a tensor, aka 3D matrix.

}  // namespace arma