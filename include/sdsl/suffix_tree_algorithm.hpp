// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*!\file suffix_tree_algorithm.hpp
 * \brief suffix_tree_algorithm.hpp contains algorithms on CSTs
 * \author Simon Gog
 */
#ifndef INCLUDED_SDSL_SUFFIX_TREE_ALGORITHM
#define INCLUDED_SDSL_SUFFIX_TREE_ALGORITHM

#include <math.h>
#include <set>
#include <stddef.h>
#include <type_traits>
#include <utility>

#include <sdsl/config.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/sdsl_concepts.hpp>

// clang-format off
// Cyclic includes start
#include <sdsl/suffix_array_algorithm.hpp>
// Cyclic includes end
// clang-format on

namespace sdsl
{

//! Forward search for a character c on the path on depth \f$d\f$ to node \f$v\f$.
/*!
 * \param cst      The CST object
 * \param v        The node at the endpoint of the current edge.
 * \param d        The current depth of the path starting with 0.
 * \param c        The character c which should be matched with pathlabel_{root()..v}[d]
 * \param char_pos T[char_pos-d+1..char_pos] matches with the already matched pattern P[0..d-1] or
 *                 d=0 => char_pos=0.
 * \return The number of the matching substrings in T.
 *
 *    \par Time complexity
 *        \f$ \Order{ t_{\Psi} } \f$ or \f$ \Order{t_{cst.child}} \f$
 */
template <class t_cst>
typename t_cst::size_type forward_search(
    t_cst const & cst,
    typename t_cst::node_type & v,
    const typename t_cst::size_type d,
    const typename t_cst::char_type c,
    typename t_cst::size_type & char_pos,
    SDSL_UNUSED
    typename std::enable_if<std::is_same<cst_tag, typename t_cst::index_category>::value, cst_tag>::type x = cst_tag())
{
    auto cc = cst.csa.char2comp[c]; // check if c occurs in the text of the csa
    if (cc == 0 and cc != c)        //   "    " "    "     "  "    "   "  "   "
        return 0;
    typename t_cst::size_type depth_node = cst.depth(v);
    if (d < depth_node)
    { // in an edge, no  branching
        char_pos = cst.csa.psi[char_pos];
        if (char_pos < cst.csa.C[cc] or char_pos >= cst.csa.C[cc + 1])
            return 0;
        return cst.size(v);
    }
    else if (d == depth_node)
    { // at a node,  branching
        v = cst.child(v, c, char_pos);
        if (v == cst.root())
            return 0;
        else
            return cst.size(v);
    }
    else
    {
        return 0;
    }
}

//! Forward search for a pattern pat on the path on depth \f$d\f$ to node \f$v\f$.
/*!
 *    \param cst      The compressed suffix tree.
 *    \param v        The node at the endpoint of the current edge.
 *    \param d        The current depth of the path. 0 = first character on each edge of the root node.
 *    \param pat      The character c which should be matched at the path on depth \f$d\f$ to node \f$v\f$.
 *    \param char_pos One position in the text, which corresponds to the text that is already matched. If v=cst.root()
 * and d=0 => char_pos=0.
 *
 *    \par Time complexity
 *        \f$ \Order{ t_{\Psi} } \f$ or \f$ \Order{t_{cst.child}} \f$
 */
template <class t_cst, class t_pat_iter>
typename t_cst::size_type forward_search(
    t_cst const & cst,
    typename t_cst::node_type & v,
    typename t_cst::size_type d,
    t_pat_iter begin,
    t_pat_iter end,
    typename t_cst::size_type & char_pos,
    SDSL_UNUSED
    typename std::enable_if<std::is_same<cst_tag, typename t_cst::index_category>::value, cst_tag>::type x = cst_tag())
{
    if (begin == end)
        return cst.size(v);
    typename t_cst::size_type size = 0;
    t_pat_iter it = begin;
    while (it != end and (size = forward_search(cst, v, d, *it, char_pos)))
    {
        ++d;
        ++it;
    }
    return size;
}

//! Counts the number of occurrences of a pattern in a CST.
/*!
 * \tparam t_cst      CST type.
 * \tparam t_pat_iter Pattern iterator type.
 *
 * \param cst   The CST object.
 * \param begin Iterator to the begin of the pattern (inclusive).
 * \param end   Iterator to the end of the pattern (exclusive).
 * \return The number of occurrences of the pattern in the CSA.
 *
 * \par Time complexity
 *        \f$ \Order{ t_{backward\_search} } \f$
 */
template <class t_cst, class t_pat_iter>
typename t_cst::size_type count(t_cst const & cst, t_pat_iter begin, t_pat_iter end, cst_tag)
{
    return count(cst.csa, begin, end);
}

//! Calculates all occurrences of a pattern pat in a CST.
/*!
 * \tparam t_cst      CST type.
 * \tparam t_pat_iter Pattern iterator type.
 * \tparam t_rac      Resizeable random access container.
 *
 * \param cst   The CST object.
 * \param begin Iterator to the begin of the pattern (inclusive).
 * \param end   Iterator to the end of the pattern (exclusive).
 * \return A vector containing the occurrences of the pattern  in the CST.
 *
 * \par Time complexity
 *        \f$ \Order{ t_{backward\_search} + z \cdot t_{SA} } \f$, where \f$z\f$ is the number of
 *         occurrences of pattern in the CST.
 */
template <class t_cst, class t_pat_iter, class t_rac = int_vector<64>>
t_rac locate(t_cst const & cst,
             t_pat_iter begin,
             t_pat_iter end,
             SDSL_UNUSED
             typename std::enable_if<std::is_same<cst_tag, typename t_cst::index_category>::value, cst_tag>::type x =
                 cst_tag())
{
    return locate(cst.csa, begin, end);
}

//! Calculate the concatenation of edge labels from the root to the node v of a CST.
/*!
 * \tparam t_cst       CST type.
 * \tparam t_text_iter Random access iterator type.
 *
 * \param cst   The CST object.
 * \param v     The node where the concatenation of the edge label ends.
 * \param text  Random access iterator pointing to the start of an container, which can hold at least (end-begin+1)
 * character. \returns The length of the extracted edge label. \pre text has to be initialized with enough memory (\f$
 * cst.depth(v)+1\f$ bytes) to hold the extracted text.
 */
template <class t_cst, class t_text_iter>
typename t_cst::size_type extract(
    t_cst const & cst,
    const typename t_cst::node_type & v,
    t_text_iter text,
    SDSL_UNUSED
    typename std::enable_if<std::is_same<cst_tag, typename t_cst::index_category>::value, cst_tag>::type x = cst_tag())
{
    if (v == cst.root())
    {
        text[0] = 0;
        return 0;
    }
    // first get the suffix array entry of the leftmost leaf in the subtree rooted at v
    typename t_cst::size_type begin = cst.csa[cst.lb(v)];
    // then call the extract method on the compressed suffix array
    extract(cst.csa, begin, begin + cst.depth(v) - 1, text);
}

//! Calculate the concatenation of edge labels from the root to the node v of of c CST.
/*!
 * \tparam t_rac Random access container which should hold the result.
 * \tparam t_cst CSA type.
 *
 * \param cst The CST object.
 * \return A t_rac object holding the extracted edge label.
 * \return The string of the concatenated edge labels from the root to the node v.
 */
template <class t_cst>
typename t_cst::csa_type::string_type extract(
    t_cst const & cst,
    const typename t_cst::node_type & v,
    SDSL_UNUSED
    typename std::enable_if<std::is_same<cst_tag, typename t_cst::index_category>::value, cst_tag>::type x = cst_tag())
{
    typedef typename t_cst::csa_type::string_type t_rac;
    if (v == cst.root())
    {
        return t_rac{};
    }
    // first get the suffix array entry of the leftmost leaf in the subtree rooted at v
    typename t_cst::size_type begin = cst.csa[cst.lb(v)];
    // then call the extract method on the compressed suffix array
    return extract(cst.csa, begin, begin + cst.depth(v) - 1);
}

//! Calculate the zeroth order entropy of the text that follows a certain substring s
/*!
 * \param v     A suffix tree node v. The label of the path from the root to v is s.
 * \param cst   The suffix tree of v.
 * \return      The zeroth order entropy of the concatenation of all characters that follow
 * s in the original text.
 */
template <class t_cst>
double H0(const typename t_cst::node_type & v, t_cst const & cst)
{
    if (cst.is_leaf(v))
    {
        return 0;
    }
    else
    {
        double h0 = 0;
        auto n = cst.size(v);
        for (auto const & child : cst.children(v))
        {
            double p = ((double)cst.size(child)) / n;
            h0 -= p * log2(p);
        }
        return h0;
    }
}

//! Calculate the k-th order entropy of a text
/*!
 * \param cst The suffix tree.
 * \param k   Parameter k for which H_k should be calculated.
 * \return    H_k and the number of contexts.
 */
template <class t_cst>
std::pair<double, size_t> Hk(t_cst const & cst, typename t_cst::size_type k)
{
    double hk = 0;
    size_t context = 0;
    std::set<typename t_cst::size_type> leafs_with_d_smaller_k;
    for (typename t_cst::size_type d = 1; d < k; ++d)
    {
        leafs_with_d_smaller_k.insert(cst.csa.isa[cst.csa.size() - d]);
    }
    for (typename t_cst::const_iterator it = cst.begin(), end = cst.end(); it != end; ++it)
    {
        if (it.visit() == 1)
        {
            if (!cst.is_leaf(*it))
            {
                typename t_cst::size_type d = cst.depth(*it);
                if (d >= k)
                {
                    if (d == k)
                    {
                        hk += cst.size(*it) * H0(*it, cst);
                    }
                    ++context;
                    it.skip_subtree();
                }
            }
            else
            {
                // if d of leaf is >= k, add context
                if (leafs_with_d_smaller_k.find(cst.lb(*it)) == leafs_with_d_smaller_k.end())
                {
                    ++context;
                }
            }
        }
    }
    hk /= cst.size();
    return {hk, context};
}

} // namespace sdsl
#endif
