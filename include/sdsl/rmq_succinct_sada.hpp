// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*!\file rmq_succinct_sada.hpp
 * \brief rmq_succinct_sada.hpp contains the class rmq_succinct_sada which supports range minimum or range maximum
 * queries on a random access container in constant time and \f$4 n+o(n) bits\f$ space. \author Simon Gog
 */
#ifndef INCLUDED_SDSL_RMQ_SUCCINCT_SADA
#define INCLUDED_SDSL_RMQ_SUCCINCT_SADA

#include <assert.h>
#include <iosfwd>
#include <stack>
#include <stdint.h>
#include <string>

#include <sdsl/bp_support_sada.hpp>
#include <sdsl/cereal.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support_v5.hpp>
#include <sdsl/rmq_succinct_sct.hpp>
#include <sdsl/select_support_mcl.hpp>
#include <sdsl/select_support_scan.hpp>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

//! Namespace for the succinct data structure library.
namespace sdsl
{

template <bool t_min = true,
          class t_bp_support = bp_support_sada<256, 32, rank_support_v5<1>, select_support_scan<1>>,
          class t_rank_10 = rank_support_v<10, 2>,
          class t_select_10 = select_support_mcl<10, 2>>
class rmq_succinct_sada;

template <class t_bp_support = bp_support_sada<256, 32, rank_support_v5<>, select_support_scan<1>>,
          class t_rank_10 = rank_support_v<10, 2>,
          class t_select_10 = select_support_mcl<10, 2>>
struct range_maximum_support_sada
{
    typedef rmq_succinct_sada<false, t_bp_support, t_rank_10, t_select_10> type;
};

//! A class to support range minimum or range maximum queries on a random access container.
/*!
 *  \tparam t_min        Specifies whether the data structure should answer range min/max queries (mimumum=true)
 *  \tparam t_bp_support Type of Support structure for the BPS-DFS
 *  \tparam t_rank_10    Type of rank structure for it pattern `10`.
 *  \tparam t_select_10  Type of select structure for bit pattern `10`.
 * \par Time complexity
 *        \f$ \Order{1} \f$ for the range min/max queries if the BPS  supports constant time operations.
 * \par Space complexity:
 *        \f$ 4n+o(n) \f$ bits for the data structure ( \f$ n=size() \f$ ).
 *
 * \par Reference
 *      TODO: Kunihiko Sadakane
 *
 */
template <bool t_min, class t_bp_support, class t_rank_10, class t_select_10>
class rmq_succinct_sada
{
    bit_vector m_ect_bp;           //!< A bit vector which contains the BP-DFS
    t_bp_support m_ect_bp_support; //!< Support structure for the BP-DFS
    t_rank_10 m_ect_bp_rank10;     //!< A rank support for bit pattern `10`
    t_select_10 m_ect_bp_select10; //!< A select support for bit pattern `10`

public:
    typedef typename bit_vector::size_type size_type;
    typedef typename bit_vector::size_type value_type;

    typedef t_bp_support bp_support_type;
    typedef t_rank_10 rank_support10_type;
    typedef t_select_10 select_support10_type;

    bit_vector const & ect_bp = m_ect_bp;
    bp_support_type const & ect_bp_support = m_ect_bp_support;
    rank_support10_type const & ect_bp_rank10 = m_ect_bp_rank10;
    select_support10_type const & ect_bp_select10 = m_ect_bp_select10;

private:
    typedef rmq_succinct_sct<t_min> rmq_construct_helper_type;

    // helper class for the construction
    struct state
    {
        size_type l, r; // left and right interval
        size_type m;    // index of the rmq
        uint8_t visit;  // 1==first, 2==second, 3==third visit
        state(size_type fl = 0, size_type fr = 0, size_type fm = 0, uint8_t fvisit = 0) :
            l(fl),
            r(fr),
            m(fm),
            visit(fvisit)
        {}
    };

    //! Constructor
    /*!\tparam t_rac A random access container.
     *  \param  v     Ponter to container object.
     */
    template <class t_rac>
    void construct_bp_of_extended_cartesian_tree(t_rac const * v, rmq_construct_helper_type const & rmq_helper)
    {
        m_ect_bp.resize(4 * v->size());
        if (v->size() > 0)
        {
            size_type bp_cnt = 0;
            size_type l = 0, r = v->size() - 1;
            std::stack<state> state_stack;
            state_stack.push(state(l, r, rmq_helper(l, r), 1));
            while (!state_stack.empty())
            {
                state s = state_stack.top();
                state_stack.pop();
                if (1 == s.visit)
                {
                    m_ect_bp[bp_cnt++] = 1; // write beginning of inner node
                    state_stack.push(state(s.l, s.r, s.m, 2));
                    if (s.m > s.l)
                    {
                        state_stack.push(state(s.l, s.m - 1, rmq_helper(s.l, s.m - 1), 1));
                    }
                }
                else if (2 == s.visit)
                {
                    m_ect_bp[bp_cnt++] = 1; // write leaf
                    m_ect_bp[bp_cnt++] = 0;
                    state_stack.push(state(s.l, s.r, s.m, 3));
                    if (s.m < s.r)
                    {
                        state_stack.push(state(s.m + 1, s.r, rmq_helper(s.m + 1, s.r), 1));
                    }
                }
                else if (3 == s.visit)
                {
                    m_ect_bp[bp_cnt++] = 0; // write end of inner node
                }
            }
            assert(bp_cnt == 4 * v->size());
        }
    }

public:
    //! Default Constructor
    rmq_succinct_sada()
    {}

    //! Constructor
    template <class t_rac>
    rmq_succinct_sada(t_rac const * v = nullptr)
    {
        if (v != nullptr)
        {
            rmq_construct_helper_type rmq_helper(v);
            m_ect_bp.resize(4 * v->size());
            construct_bp_of_extended_cartesian_tree(v, rmq_helper);
            m_ect_bp_support = bp_support_type(&m_ect_bp);
            util::init_support(m_ect_bp_rank10, &m_ect_bp);
            util::init_support(m_ect_bp_select10, &m_ect_bp);
        }
    }

    //! Copy constructor
    rmq_succinct_sada(rmq_succinct_sada const & rm) :
        m_ect_bp(rm.m_ect_bp),
        m_ect_bp_support(rm.m_ect_bp_support),
        m_ect_bp_rank10(rm.m_ect_bp_rank10),
        m_ect_bp_select10(rm.m_ect_bp_select10)
    {
        m_ect_bp_support.set_vector(&m_ect_bp);
        m_ect_bp_rank10.set_vector(&m_ect_bp);
        m_ect_bp_select10.set_vector(&m_ect_bp);
    }

    //! Move constructor
    rmq_succinct_sada(rmq_succinct_sada && rm)
    {
        *this = std::move(rm);
    }

    //! Destructor
    ~rmq_succinct_sada()
    {}

    rmq_succinct_sada & operator=(rmq_succinct_sada const & rm)
    {
        if (this != &rm)
        {
            rmq_succinct_sada tmp(rm);
            *this = std::move(tmp);
        }
        return *this;
    }

    rmq_succinct_sada & operator=(rmq_succinct_sada && rm)
    {
        if (this != &rm)
        {
            m_ect_bp = std::move(rm.m_ect_bp);
            m_ect_bp_support = std::move(rm.m_ect_bp_support);
            m_ect_bp_support.set_vector(&m_ect_bp);
            m_ect_bp_rank10 = std::move(rm.m_ect_bp_rank10);
            m_ect_bp_rank10.set_vector(&m_ect_bp);
            m_ect_bp_select10 = std::move(rm.m_ect_bp_select10);
            m_ect_bp_select10.set_vector(&m_ect_bp);
        }
        return *this;
    }

    //! Range minimum/maximum query for the supported random access container v.
    /*!
     * \param l Leftmost position of the interval \f$[\ell..r]\f$.
     * \param r Rightmost position of the interval \f$[\ell..r]\f$.
     * \return The minimal index i with \f$\ell \leq i \leq r\f$ for which \f$ v[i] \f$ is minimal/maximal.
     * \pre
     *   - r < size()
     *   - \f$ \ell \leq r \f$
     * \par Time complexity
     *      \f$ \Order{1} \f$
     */
    size_type operator()(const size_type l, const size_type r) const
    {
        assert(l <= r);
        assert(r < size());
        if (l == r)
            return l;
        size_type x = m_ect_bp_select10(l + 1);
        size_type y = m_ect_bp_select10(r + 1);
        size_type z = m_ect_bp_support.rmq(x, y);
        size_type f = z + 1 - 2 * (m_ect_bp[z]);
        return m_ect_bp_rank10(f - 1);
    }

    size_type size() const
    {
        return m_ect_bp.size() / 4;
    }

    size_type serialize(std::ostream & out, structure_tree_node * v = nullptr, std::string name = "") const
    {
        structure_tree_node * child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += m_ect_bp.serialize(out, child, "ect_bp");
        written_bytes += m_ect_bp_support.serialize(out, child, "ect_bp_support");
        written_bytes += m_ect_bp_rank10.serialize(out, child, "ect_bp_rank10");
        written_bytes += m_ect_bp_select10.serialize(out, child, "ect_bp_select10");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void load(std::istream & in)
    {
        m_ect_bp.load(in);
        m_ect_bp_support.load(in, &m_ect_bp);
        m_ect_bp_rank10.load(in, &m_ect_bp);
        m_ect_bp_select10.load(in, &m_ect_bp);
    }

    template <typename archive_t>
    void CEREAL_SAVE_FUNCTION_NAME(archive_t & ar) const
    {
        ar(CEREAL_NVP(m_ect_bp));
        ar(CEREAL_NVP(m_ect_bp_support));
        ar(CEREAL_NVP(m_ect_bp_rank10));
        ar(CEREAL_NVP(m_ect_bp_select10));
    }

    template <typename archive_t>
    void CEREAL_LOAD_FUNCTION_NAME(archive_t & ar)
    {
        ar(CEREAL_NVP(m_ect_bp));
        ar(CEREAL_NVP(m_ect_bp_support));
        m_ect_bp_support.set_vector(&m_ect_bp);
        ar(CEREAL_NVP(m_ect_bp_rank10));
        m_ect_bp_rank10.set_vector(&m_ect_bp);
        ar(CEREAL_NVP(m_ect_bp_select10));
        m_ect_bp_select10.set_vector(&m_ect_bp);
    }

    //! Equality operator.
    bool operator==(rmq_succinct_sada const & other) const noexcept
    {
        return (m_ect_bp == other.m_ect_bp) && (m_ect_bp_support == other.m_ect_bp_support)
            && (m_ect_bp_rank10 == other.m_ect_bp_rank10) && (m_ect_bp_select10 == other.m_ect_bp_select10);
    }

    //! Inequality operator.
    bool operator!=(rmq_succinct_sada const & other) const noexcept
    {
        return !(*this == other);
    }
};

} // end namespace sdsl
#endif
