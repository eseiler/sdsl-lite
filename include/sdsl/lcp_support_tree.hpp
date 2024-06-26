// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
#ifndef INCLUDED_SDSL_SUPPORT_LCP_TREE
#define INCLUDED_SDSL_SUPPORT_LCP_TREE

#include <iostream>
#include <string>

#include <sdsl/bits.hpp>
#include <sdsl/cereal.hpp>
#include <sdsl/config.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/io.hpp>
#include <sdsl/iterators.hpp>
#include <sdsl/lcp_wt.hpp>
#include <sdsl/ram_fs.hpp>
#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/sorted_multi_stack_support.hpp>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

namespace sdsl
{

inline void construct_first_child_lcp(int_vector_buffer<> & lcp_buf, int_vector<> & fc_lcp)
{
    typedef int_vector_size_type size_type;
    size_type n = lcp_buf.size();
    if (n == 0)
    { // if n == 0 we are done
        fc_lcp = int_vector<>(0);
        return;
    }
    fc_lcp = int_vector<>(n, 0, bits::hi(n) + 1);
    size_type fc_cnt = 0; // first child counter
    sorted_multi_stack_support vec_stack(n);
    size_type y;
    for (size_type i = 0, x; i < n; ++i)
    {
        x = lcp_buf[i];
        while (!vec_stack.empty() and x < vec_stack.top())
        {
            y = vec_stack.top();
            if (vec_stack.pop())
            {
                fc_lcp[fc_cnt++] = y;
            }
        }
        vec_stack.push(x);
    }

    while (!vec_stack.empty())
    {
        y = vec_stack.top();
        if (vec_stack.pop())
        {
            fc_lcp[fc_cnt++] = y;
        }
    }
    if (fc_cnt < fc_lcp.size())
    {
        fc_lcp.resize(fc_cnt);
        fc_lcp.shrink_to_fit();
    }
}

/*! This class composes a virtual LCP array from a LCP arrays which is in suffix array order
 * (e.g. lcp_byte or lcp_bitcompressed) and a CST.
 *    The time consumption of the []-operator depends on:
 *    - The time consumption of the tlcp_idx function of the CST
 *    - The access time to the suffix array ordered LCP array
 *
 * \tparam t_lcp Type of the underlying LCP array. Must be an suffix array ordered one.
 * \tparam t_cst Type of the underlying CST.
 */
template <class t_lcp, class t_cst>
class _lcp_support_tree
{
public:
    typedef typename t_lcp::value_type value_type;
    typedef random_access_const_iterator<_lcp_support_tree> const_iterator;
    typedef const_iterator iterator;
    typedef const value_type const_reference;
    typedef const_reference reference;
    typedef const_reference * pointer;
    typedef const pointer const_pointer;
    typedef typename t_lcp::size_type size_type;
    typedef typename t_lcp::difference_type difference_type;

    typedef lcp_tree_compressed_tag lcp_category;

    enum
    {
        fast_access = 0,
        text_order = t_lcp::text_order,
        sa_order = t_lcp::sa_order
    };

    template <class CST>
    struct type
    {
        typedef _lcp_support_tree lcp_type;
    };

private:
    t_cst const * m_cst;
    t_lcp m_lcp;

public:
    //! Default constructor
    _lcp_support_tree() = default;

    // Destructor
    ~_lcp_support_tree() = default;

    //! Copy/Move constructor
    _lcp_support_tree(_lcp_support_tree const &) = default;
    _lcp_support_tree(_lcp_support_tree &&) = default;
    _lcp_support_tree & operator=(_lcp_support_tree const &) = default;
    _lcp_support_tree & operator=(_lcp_support_tree &&) = default;

    //! Constructor
    /*!
     *  \param config  Cache configuration.
     *  \param cst     A pointer to the CST.
     */
    _lcp_support_tree(cache_config & config, t_cst const * cst = nullptr)
    {
        m_cst = cst;
        std::string fc_lcp_key = "fc_lcp_" + util::to_string(util::id());
        std::string tmp_file = cache_file_name(fc_lcp_key, config);
        {
            int_vector<0> temp_lcp;
            int_vector_buffer<> lcp_buf(cache_file_name(conf::KEY_LCP, config));
            construct_first_child_lcp(lcp_buf, temp_lcp);
            // TODO: store LCP values directly
            store_to_file(temp_lcp, tmp_file);
        }
        {
            {
                m_lcp = t_lcp(config, fc_lcp_key); // works for lcp_kurtz, lcp_wt and lcp_bitcompressed
            }
        }
        sdsl::remove(tmp_file);
    }

    size_type size() const
    {
        return m_cst->size();
    }

    void set_cst(t_cst const * cst)
    {
        m_cst = cst;
    }

    static size_type max_size()
    {
        return t_lcp::max_size();
    }

    size_type empty() const
    {
        return m_lcp.empty();
    }

    //! Returns a const_iterator to the first element.
    const_iterator begin() const
    {
        return const_iterator(this, 0);
    }

    //! Returns a const_iterator to the element after the last element.
    const_iterator end() const
    {
        return const_iterator(this, size());
    }

    //! []-operator
    /*!\param i Index of the value. \f$ i \in [0..size()-1]\f$.
     * \par Time complexity
     *     \f$ \Order{t_{find\_close} + t_{rank}} \f$
     */
    inline value_type operator[](size_type i) const
    {
        return m_lcp[m_cst->tlcp_idx(i)];
    }

    //! Serialize to a stream.
    size_type serialize(std::ostream & out, structure_tree_node * v = nullptr, std::string name = "") const
    {
        size_type written_bytes = 0;
        structure_tree_node * child = structure_tree::add_child(v, name, util::class_name(*this));
        written_bytes += m_lcp.serialize(out, child, "lcp");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    //! Load from a stream.
    void load(std::istream & in, t_cst const * cst = nullptr)
    {
        m_lcp.load(in); // works for lcp_byte and lcp_bitcompressed
        m_cst = cst;
    }

    template <typename archive_t>
    void CEREAL_SAVE_FUNCTION_NAME(archive_t & ar) const
    {
        ar(CEREAL_NVP(m_lcp));
    }

    template <typename archive_t>
    void CEREAL_LOAD_FUNCTION_NAME(archive_t & ar)
    {
        ar(CEREAL_NVP(m_lcp));
    }

    //! Equality operator.
    bool operator==(_lcp_support_tree const & other) const noexcept
    {
        return (m_lcp == other.m_lcp);
    }

    //! Inequality operator.
    bool operator!=(_lcp_support_tree const & other) const noexcept
    {
        return !(*this == other);
    }
};

//! Helper class which provides _lcp_support_tree the context of a CST.
template <class t_lcp = lcp_wt<>>
struct lcp_support_tree
{
    template <class t_cst>
    using type = _lcp_support_tree<t_lcp, t_cst>;
};

} // namespace sdsl
#endif
