// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*!\file k2_tree_helper.hpp
 * \brief k2_tree_helper.hpp contains helper functions and definitions for a k^2-tree implementation.
 * \author Francisco Montoto
 */
#ifndef INCLUDED_SDSL_K2_TREE_HELPER
#define INCLUDED_SDSL_K2_TREE_HELPER

#include <cmath>
#include <deque>
#include <stdint.h>
#include <utility>
#include <vector>

#include <sdsl/int_vector.hpp>

#include <ext/alloc_traits.h>

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! Namespace for the k2_tree
namespace k2_tree_ns
{

typedef int_vector<>::size_type idx_type;
typedef int_vector<>::size_type size_type;

template <typename t_bv = bit_vector>
int _build_from_matrix(const std::vector<std::vector<int>> & matrix,
                       const uint8_t k,
                       int n,
                       const int height,
                       int l,
                       int p,
                       int q,
                       std::vector<std::deque<t_bv>> & acc)
{
    unsigned i, j, b_size = pow(k, 2);
    t_bv b(b_size, 0);
    bool is_leaf = (l == height);

    if (is_leaf)
    {
        for (i = 0; i < k; i++)
            for (j = 0; j < k; j++)
                if (p + i < matrix.size() && q + j < matrix.size() && matrix[p + i][q + j] == 1) b[i * k + j] = 1;
    }
    else
    { // Internal node
        for (i = 0; i < k; i++)
            for (j = 0; j < k; j++)
                b[i * k +
                  j] = _build_from_matrix(matrix, k, n / k, height, l + 1, p + i * (n / k), q + j * (n / k), acc);
    }

    // TODO There must be a better way to check if there is a 1 at b.
    for (i = 0; i < b_size; i++)
        if (b[i] == 1) break;
    if (i == b_size) // If there are not 1s at b.
        return 0;

    acc[l].push_back(std::move(b));
    return 1;
}

/*! Get the chunk index ([0, k^2[) of a submatrix point.
 *
 * Gets a point in the global matrix and returns its corresponding chunk
 * in the submatrix specified.
 *
 * \param v Row of the point in the global matrix.
 * \param u Column of the point in the global matrix.
 * \param c_0 Column offset of the submatix in the global matrix.
 * \param r_0 Row offset of the submatrix in the global matrix.
 * \param l size of the chunk at the submatrix.
 * \param k the k parameter from the k^2 tree.
 * \returns the index of the chunk containing the point at the submatrix.
 */
inline uint16_t get_chunk_idx(idx_type v, idx_type u, idx_type c_0, idx_type r_0, size_type l, uint8_t k)
{
    return ((v - r_0) / l) * k + (u - c_0) / l;
}

template <typename t_bv = bit_vector>
void build_template_vector(bit_vector & k_t_, bit_vector & k_l_, t_bv & k_t, t_bv & k_l)
{
    k_t = t_bv(k_t_);
    k_l = t_bv(k_l_);
}

template <>
inline void build_template_vector<bit_vector>(bit_vector & k_t_, bit_vector & k_l_, bit_vector & k_t, bit_vector & k_l)
{
    std::swap(k_t, k_t_);
    std::swap(k_l, k_l_);
}

} // end namespace k2_tree_ns
} // end namespace sdsl

#endif
