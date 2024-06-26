// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*!\file wt_gmr.hpp
 * \brief wt_gmr.hpp contains a specialized class to support select, rank
 * and access on inputs over a large alphabet.
 * \author Alexander Diehm, Timo Beller, Simon Gog
 */
#ifndef INCLUDED_SDSL_WT_GMR
#define INCLUDED_SDSL_WT_GMR

#include <algorithm>
#include <assert.h>
#include <cstdint>
#include <iosfwd>
#include <iterator>
#include <string>
#include <type_traits>
#include <utility>

#include <sdsl/bits.hpp>
#include <sdsl/cereal.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/io.hpp>
#include <sdsl/iterators.hpp>
#include <sdsl/ram_fs.hpp>
#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! Class inv_multi_perm_support adds access to the inverse of permutations.
/*!
 * \tparam t_s    Sampling parameter of the inverse permutation.
 * \tparam t_rac  Type of the random access container used for storing the permutation.
 * \tparam t_bv   Type of the bitvector used to indicate back-pointers.
 * \tparam t_rank Type of rank_support to rank the indicator bitvector.
 *
 * This support class adds access to the inverse of permutations in at
 * most \(t_s\) steps.
 *
 * \par References
 *      [1] J. Munro, R. Raman, V. Raman, S. Rao: ,,Succinct representation
 *          of permutations'', Proceedings of ICALP 2003
 */
template <uint64_t t_s = 32,
          class t_rac = int_vector<>,
          class t_bv = bit_vector,
          class t_rank = typename t_bv::rank_1_type>
class inv_multi_perm_support
{
public:
    typedef t_rac iv_type;
    typedef typename iv_type::size_type size_type;
    typedef typename iv_type::value_type value_type;
    typedef typename iv_type::difference_type difference_type;
    typedef t_bv bit_vector_type;
    typedef t_rank rank_type;
    typedef random_access_const_iterator<inv_multi_perm_support> const_iterator;

private:
    iv_type const * m_perm = nullptr; // pointer to supported permutation
    uint64_t m_chunksize;             // size of one permutation
    int_vector<> m_back_pointer;      // back pointers
    bit_vector_type m_marked;         // back pointer marking
    rank_type m_marked_rank;          // rank support for back pointer marking

public:
    //! Default constructor
    inv_multi_perm_support(){};

    //! Constructor
    inv_multi_perm_support(iv_type const * perm, int_vector<> & iv, uint64_t chunksize) :
        m_perm(perm),
        m_chunksize(chunksize)
    {
        bit_vector marked(iv.size(), 0);
        bit_vector done(m_chunksize, 0);

        size_type max_back_pointer = 0;
        for (size_type i = 0, off = 0; i < iv.size(); ++i)
        {
            if (i == off + chunksize)
            {
                off = i;
                util::set_to_value(done, 0);
            }
            if (!done[i - off])
            {
                done[i - off] = 1;
                size_type back_pointer = i, j = i, j_new = 0;
                uint64_t steps = 0, all_steps = 0;
                while ((j_new = (iv[j] + off)) != i)
                {
                    j = j_new;
                    done[j - off] = 1;
                    ++steps;
                    ++all_steps;
                    if (t_s == steps)
                    {
                        max_back_pointer = std::max(max_back_pointer, back_pointer - off);
                        marked[j] = 1;
                        steps = 0;
                        back_pointer = j;
                    }
                }
                if (all_steps > t_s)
                {
                    marked[i] = 1;
                    max_back_pointer = std::max(max_back_pointer, back_pointer - off);
                }
            }
        }

        m_marked = t_bv(std::move(marked));
        util::init_support(m_marked_rank, &m_marked);

        util::set_to_value(done, 0);
        size_type n_bp = m_marked_rank(iv.size());
        m_back_pointer = int_vector<>(n_bp, 0, bits::hi(max_back_pointer) + 1);

        for (size_type i = 0, off = 0; i < iv.size(); ++i)
        {
            if (i == off + chunksize)
            {
                off = i;
                util::set_to_value(done, 0);
            }
            if (!done[i - off])
            {
                done[i - off] = 1;
                size_type back_pointer = i, j = i, j_new = 0;
                uint64_t steps = 0, all_steps = 0;
                while ((j_new = (iv[j] + off)) != i)
                {
                    j = j_new;
                    done[j - off] = 1;
                    ++steps;
                    ++all_steps;
                    if (t_s == steps)
                    {
                        m_back_pointer[m_marked_rank(j)] = back_pointer - off;
                        steps = 0;
                        back_pointer = j;
                    }
                }
                if (all_steps > t_s)
                {
                    m_back_pointer[m_marked_rank(i)] = back_pointer - off;
                }
            }
        }
    }

    //! Copy constructor
    inv_multi_perm_support(inv_multi_perm_support const & p) :
        m_perm(p.m_perm),
        m_chunksize(p.m_chunksize),
        m_back_pointer(p.m_back_pointer),
        m_marked(p.m_marked),
        m_marked_rank(p.m_marked_rank)
    {
        m_marked_rank.set_vector(&m_marked);
    }

    //! Move constructor
    inv_multi_perm_support(inv_multi_perm_support && p)
    {
        *this = std::move(p);
    }

    //! Assignment operation
    inv_multi_perm_support & operator=(inv_multi_perm_support const & p)
    {
        if (this != &p)
        {
            m_perm = p.m_perm;
            m_chunksize = p.m_chunksize;
            m_back_pointer = p.m_back_pointer;
            m_marked = p.m_marked;
            m_marked_rank = p.m_marked_rank;
            m_marked_rank.set_vector(&m_marked);
        }
        return *this;
    }

    //! Assignment move operation
    inv_multi_perm_support & operator=(inv_multi_perm_support && p)
    {
        if (this != &p)
        {
            m_perm = std::move(p.m_perm);
            m_chunksize = std::move(p.m_chunksize);
            m_back_pointer = std::move(p.m_back_pointer);
            m_marked = std::move(p.m_marked);
            m_marked_rank = std::move(p.m_marked_rank);
            m_marked_rank.set_vector(&m_marked);
        }
        return *this;
    }

    //! Returns the size of the original vector.
    size_type size() const
    {
        return nullptr == m_perm ? 0 : m_perm->size();
    }

    //! Returns whether the original vector contains no data.
    bool empty() const
    {
        return size() == 0;
    }

    //! Access operator
    /*
     *  \par Time complexity
     *       \f$ \Order{t_s} \f$
     */
    value_type operator[](size_type i) const
    {
        size_type off = (i / m_chunksize) * m_chunksize;
        size_type j = i, j_new = 0;
        while ((j_new = ((*m_perm)[j]) + off) != i)
        {
            if (m_marked[j])
            {
                j = m_back_pointer[m_marked_rank(j)] + off;
                while ((j_new = ((*m_perm)[j]) + off) != i)
                    j = j_new;
            }
            else
            {
                j = j_new;
            }
        }
        return j;
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

    void set_vector(iv_type const * v)
    {
        m_perm = v;
    }

    //! Serialize into stream
    size_type serialize(std::ostream & out, structure_tree_node * v = nullptr, std::string name = "") const
    {
        structure_tree_node * child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += write_member(m_chunksize, out, child, "chunksize");
        written_bytes += m_back_pointer.serialize(out, child, "back_pointer");
        written_bytes += m_marked.serialize(out, child, "marked");
        written_bytes += m_marked_rank.serialize(out, child, "marked_rank");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    //! Load sampling from disk
    void load(std::istream & in, iv_type const * v = nullptr)
    {
        set_vector(v);
        read_member(m_chunksize, in);
        m_back_pointer.load(in);
        m_marked.load(in);
        m_marked_rank.load(in, &m_marked);
    }

    //!\brief Serialise (save) via cereal
    template <typename archive_t>
    void CEREAL_SAVE_FUNCTION_NAME(archive_t & ar) const
    {
        ar(CEREAL_NVP(m_chunksize));
        ar(CEREAL_NVP(m_back_pointer));
        ar(CEREAL_NVP(m_marked));
        ar(CEREAL_NVP(m_marked_rank));
    }

    //!\brief Load via cereal
    template <typename archive_t>
    void CEREAL_LOAD_FUNCTION_NAME(archive_t & ar)
    {
        ar(CEREAL_NVP(m_chunksize));
        ar(CEREAL_NVP(m_back_pointer));
        ar(CEREAL_NVP(m_marked));
        ar(CEREAL_NVP(m_marked_rank));
        m_marked_rank.set_vector(&m_marked);
    }

    //! Equality operator.
    bool operator==(inv_multi_perm_support const & other) const noexcept
    {
        return (m_chunksize == other.m_chunksize) && (m_back_pointer == other.m_back_pointer)
            && (m_marked == other.m_marked) && (m_marked_rank == other.m_marked_rank);
    }

    //! Inequality operator.
    bool operator!=(inv_multi_perm_support const & other) const noexcept
    {
        return !(*this == other);
    }
};

template <class t_rac>
void _transform_to_compressed(int_vector<> & iv,
                              typename std::enable_if<!(std::is_same<t_rac, int_vector<>>::value), t_rac>::type & rac,
                              const std::string filename)
{
    std::string tmp_file_name = tmp_file(filename, "_compress_int_vector");
    store_to_file(iv, tmp_file_name);
    util::clear(iv);
    int_vector_buffer<> buf(tmp_file_name, std::ios::in, 1024 * 1024, iv.width());
    rac = t_rac(buf);
    buf.close(true); // delete tmp_file
}

template <class t_rac>
void _transform_to_compressed(int_vector<> & iv,
                              typename std::enable_if<std::is_same<t_rac, int_vector<>>::value, t_rac>::type & rac,
                              const std::string)
{
    rac = std::move(iv);
}

//! A wavelet tree class for integer sequences.
/*!
 *  \tparam t_rac         Type of the random access container used for E.
 *  \tparam t_bitvector   Type of the bitvector used for storing B.
 *  \tparam t_select      Type of the support structure for select on pattern `1`.
 *  \tparam t_select_zero Type of the support structure for select on pattern `0`.
 *
 * This is an implementation of the first proposal in the SODA paper of Golynski et. al.
 * which support fast rank and select, but not fast access.
 *
 * \par References
 * [1] A. Golynski, J. Munro and S. Rao:
 *     ,,Rank/select operations on large alphabets: a tool for text indexing''
 *     Proceedings of SODA 2006.
 *
 *   @ingroup wt
 */
template <class t_rac = int_vector<>,
          class t_bitvector = bit_vector,
          class t_select = typename t_bitvector::select_1_type,
          class t_select_zero = typename t_bitvector::select_0_type>
class wt_gmr_rs
{
public:
    typedef int_vector<>::size_type size_type;
    typedef int_vector<>::value_type value_type;
    typedef typename t_bitvector::difference_type difference_type;
    typedef random_access_const_iterator<wt_gmr_rs> const_iterator;
    typedef const_iterator iterator;
    typedef wt_tag index_category;
    typedef int_alphabet_tag alphabet_category;
    enum
    {
        lex_ordered = 0
    };

private:
    t_bitvector m_bv_blocks;
    t_rac m_e;
    t_select m_bv_blocks_select1;
    t_select_zero m_bv_blocks_select0;
    uint64_t m_size;           // input length
    uint64_t m_block_size = 0; // size of the blocks
    uint64_t m_blocks;         // blocks per character
    uint64_t m_sigma = 0;

public:
    size_type const & sigma = m_sigma;

    //! Default constructor
    wt_gmr_rs() = default;

    //! Construct the wavelet tree from a sequence defined by two interators
    /*!
     * \param begin   Iterator to the start of the input.
     * \param end     Iterator one past the end of the input.
     * \param tmp_dir Temporary directory for intermediate results.
     *    \par Time complexity
     *        \f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
     *        I.e. we need \Order{n\log n} if rac is a permutation of 0..n-1.
     */
    template <typename t_it>
    wt_gmr_rs(t_it begin, t_it end, std::string tmp_dir = ram_file_name("")) : m_size(std::distance(begin, end))
    {
        // Determine max. symbol
        for (auto it = begin; it != end; ++it)
        {
            value_type value = *it;
            if (m_block_size < value)
                m_block_size = value;
        }
        ++m_block_size;

        // Create and fill m_bv_blocks
        m_blocks = (m_size + m_block_size - 1) / m_block_size;
        bit_vector b(m_size + m_block_size * m_blocks + 1, 0);
        int_vector<> symbols(m_block_size, 0, bits::hi(m_size) + 1);
        {
            int_vector<> tmp(m_block_size * m_blocks, 0, bits::hi(m_block_size) + 1);
            uint64_t j = 0, offset = 0;
            for (auto it = begin; it != end; ++it, ++j)
            {
                if (j == m_block_size)
                {
                    ++offset;
                    j = 0;
                }
                ++tmp[(*it) * m_blocks + offset];
            }

            for (uint64_t i = 0; i < symbols.size(); ++i)
            {
                for (uint64_t j = m_blocks * i; j < (i + 1) * m_blocks; ++j)
                {
                    symbols[i] += tmp[j];
                }
            }

            for (uint64_t i = 0, l = 1; i < tmp.size(); ++i, ++l)
            {
                for (uint64_t j = 0; j < tmp[i]; ++j)
                    b[l++] = 1;
            }

            // calc m_sigma
            bool write = true;
            uint64_t blocks = 0;
            for (uint64_t i = 1; i < b.size(); ++i)
            {
                if (blocks == m_blocks)
                {
                    blocks = 0;
                    write = true;
                }
                if (b[i])
                {
                    if (write)
                    {
                        ++m_sigma;
                        write = false;
                    }
                }
                else
                    ++blocks;
            }

            m_bv_blocks = t_bitvector(std::move(b));
        }

        // Create and fill e
        int_vector<> positions(m_size, 0, bits::hi(m_block_size) + 1);
        for (uint64_t i = 0, tmp = 0, sum = 0; i < m_block_size; ++i)
        {
            tmp = symbols[i];
            symbols[i] = sum;
            sum += tmp;
        }
        for (auto it = begin; it != end;)
        {
            for (uint64_t j = 0; j < m_block_size and it != end; ++it, ++j)
            {
                positions[symbols[*it]++] = j;
            }
        }
        _transform_to_compressed<t_rac>(positions, m_e, tmp_dir);

        util::init_support(m_bv_blocks_select0, &m_bv_blocks);
        util::init_support(m_bv_blocks_select1, &m_bv_blocks);
    }

    //! Copy constructor
    wt_gmr_rs(wt_gmr_rs const & wt) :
        m_bv_blocks(wt.m_bv_blocks),
        m_e(wt.m_e),
        m_bv_blocks_select1(wt.m_bv_blocks_select1),
        m_bv_blocks_select0(wt.m_bv_blocks_select0),
        m_size(wt.m_size),
        m_block_size(wt.m_block_size),
        m_blocks(wt.m_blocks),
        m_sigma(wt.m_sigma)
    {
        m_bv_blocks_select1.set_vector(&m_bv_blocks);
        m_bv_blocks_select0.set_vector(&m_bv_blocks);
    }

    //! Move copy constructor
    wt_gmr_rs(wt_gmr_rs && wt) :
        m_bv_blocks(std::move(wt.m_bv_blocks)),
        m_e(std::move(wt.m_e)),
        m_bv_blocks_select1(std::move(wt.m_bv_blocks_select1)),
        m_bv_blocks_select0(std::move(wt.m_bv_blocks_select0)),
        m_size(wt.m_size),
        m_block_size(wt.m_block_size),
        m_blocks(wt.m_blocks),
        m_sigma(wt.m_sigma)
    {
        m_bv_blocks_select1.set_vector(&m_bv_blocks);
        m_bv_blocks_select0.set_vector(&m_bv_blocks);
    }

    //! Assignment operator
    wt_gmr_rs & operator=(wt_gmr_rs const & wt)
    {
        wt_gmr_rs tmp(wt);
        *this = std::move(tmp);
        return *this;
    }

    //! Move assignment operator
    wt_gmr_rs & operator=(wt_gmr_rs && wt)
    {
        m_bv_blocks = std::move(wt.m_bv_blocks);
        m_e = std::move(wt.m_e);
        m_bv_blocks_select1 = std::move(wt.m_bv_blocks_select1);
        m_bv_blocks_select1.set_vector(&m_bv_blocks);
        m_bv_blocks_select0 = std::move(wt.m_bv_blocks_select0);
        m_bv_blocks_select0.set_vector(&m_bv_blocks);
        m_size = wt.m_size;
        m_block_size = wt.m_block_size;
        m_blocks = wt.m_blocks;
        m_sigma = wt.m_sigma;
        return *this;
    }

    //! Returns the size of the original vector.
    size_type size() const
    {
        return m_size;
    }

    //! Returns whether the wavelet tree contains no data.
    bool empty() const
    {
        return m_size == 0;
    }

    //! Recovers the i-th symbol of the original vector.
    /*!\param i The index of the symbol in the original vector.
     *  \returns The i-th symbol of the original vector.
     *  \par Time complexity
     *       \f$ \Order{|\Sigma|} \f$
     *  \par Precondition
     *       \f$ i < size() \f$
     */
    value_type operator[](size_type i) const
    {
        assert(i < m_size);
        size_type block = i / m_block_size + 1, val = i % m_block_size, search_begin, search_end, j;
        while (true)
        {
            j = m_bv_blocks_select0(block) + 1;
            search_begin = j - block;
            if (m_bv_blocks[j])
            {
                search_end = m_bv_blocks_select0(block + 1) - (block);
                if (search_end - search_begin < 50)
                { // After a short test, this seems to be a good threshold
                    while (search_begin < search_end and m_e[search_begin] <= val)
                    {
                        if (m_e[search_begin] == val)
                        {
                            return (block - 1) / m_blocks;
                        }
                        ++search_begin;
                    }
                }
                else
                {
                    if (std::binary_search(m_e.begin() + search_begin, m_e.begin() + search_end, val))
                    {
                        return (block - 1) / m_blocks;
                    }
                }
            }
            block += m_blocks;
        }
    }

    //! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
    /*!
     *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
     *  \param c The symbol to count the occurrences in the prefix.
     *    \returns The number of occurrences of symbol c in the prefix [0..i-1] of the supported vector.
     *  \par Time complexity
     *       \f$ \Order{\log |\Sigma|} \f$
     *  \par Precondition
     *       \f$ i \leq size() \f$
     */
    size_type rank(size_type i, value_type c) const
    {
        if (0 == i or c > m_block_size - 1)
        {
            return 0;
        }

        size_type offset = 0;
        size_type ones_before_cblock = m_bv_blocks_select0(c * m_blocks + 1) - c * m_blocks;

        auto begin = m_e.begin() + m_bv_blocks_select0(c * m_blocks + (i - 1) / m_block_size + 1)
                   - (c * m_blocks + (i - 1) / m_block_size + 1) + 1;
        auto end = m_e.begin() + m_bv_blocks_select0(c * m_blocks + (i - 1) / m_block_size + 2)
                 - (c * m_blocks + (i - 1) / m_block_size + 1);

        size_type val = (i - 1) % m_block_size;
        if (end - begin < 50)
        { // After a short test, this seems to be a good threshold
            offset = std::find_if(begin,
                                  end,
                                  [&val](auto const && x)
                                  {
                                      return x > val;
                                  })
                   - begin;
        }
        else
        {
            offset = std::lower_bound(begin, end, val + 1) - begin;
        }
        return (begin - m_e.begin()) + offset - ones_before_cblock;
    }

    //! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
    /*!
     *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
     *  \param c The symbol to count the occurrences in the prefix.
     *    \returns The number of occurrences of symbol c in the prefix [0..i-1] of the supported vector.
     *  \par Time complexity
     *       \f$ \Order{|\Sigma|} \f$
     *  \par Precondition
     *       \f$ i \leq size() \f$
     */
    std::pair<size_type, value_type> inverse_select(size_type i) const
    {
        assert(i < m_size);
        size_type block = i / m_block_size + 1, val = i % m_block_size, offset = 0, search_begin, search_end, j;
        while (true)
        {
            j = m_bv_blocks_select0(block) + 1;
            search_begin = j - block;
            if (m_bv_blocks[j])
            {
                search_end = m_bv_blocks_select0(block + 1) - (block);
                offset = 0;
                if (search_end - search_begin < 50)
                { // After a short test, this seems to be a good threshold
                    while (search_begin < search_end and m_e[search_begin] <= val)
                    {
                        if (m_e[search_begin] == val)
                        {
                            value_type c = (block - 1) / m_blocks;
                            size_type ones_before_cblock = m_bv_blocks_select0(c * m_blocks + 1) - (c * m_blocks);
                            size_type r = search_begin - ones_before_cblock;
                            return std::make_pair(r, c);
                        }
                        ++search_begin;
                    }
                }
                else
                {
                    offset = std::lower_bound(m_e.begin() + search_begin, m_e.begin() + search_end, val) - m_e.begin();
                    if (offset < search_end)
                    {
                        if (m_e[offset] == val)
                        {
                            value_type c = (block - 1) / m_blocks;
                            size_type ones_before_cblock = m_bv_blocks_select0(c * m_blocks + 1) - (c * m_blocks);
                            size_type r = offset - ones_before_cblock;
                            return std::make_pair(r, c);
                        }
                    }
                }
            }
            block += m_blocks;
        }
    }

    //! Calculates the i-th occurrence of the symbol c in the supported vector.
    /*!
     *  \param i The i-th occurrence.
     *  \param c The symbol c.
     *  \par Time complexity
     *       \f$ \Order{1} \f$
     *  \par Precondition
     *       \f$ 1 \leq i \leq rank(size(), c) \f$
     */
    size_type select(size_type i, value_type c) const
    {
        size_type k = m_bv_blocks_select0(c * m_blocks + 1) - (c * m_blocks) + i;
        return (m_bv_blocks_select1(k) - k) * m_block_size + m_e[k - 1] - c * m_blocks * m_block_size;
    }

    //! Serializes the data structure into the given ostream
    size_type serialize(std::ostream & out, structure_tree_node * v = nullptr, std::string name = "") const
    {
        structure_tree_node * child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += write_member(m_size, out, child, "size");
        written_bytes += write_member(m_block_size, out, child, "block_size");
        written_bytes += write_member(m_blocks, out, child, "blocks");
        written_bytes += write_member(m_sigma, out, child, "sigma");
        written_bytes += m_e.serialize(out, child, "E");
        written_bytes += m_bv_blocks.serialize(out, child, "bv_blocks");
        written_bytes += m_bv_blocks_select0.serialize(out, child, "bv_blocks_select0");
        written_bytes += m_bv_blocks_select1.serialize(out, child, "bv_blocks_select1");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    //! Loads the data structure from the given istream.
    void load(std::istream & in)
    {
        read_member(m_size, in);
        read_member(m_block_size, in);
        read_member(m_blocks, in);
        read_member(m_sigma, in);
        m_e.load(in);
        m_bv_blocks.load(in);
        m_bv_blocks_select0.load(in, &m_bv_blocks);
        m_bv_blocks_select1.load(in, &m_bv_blocks);
    }

    //!\brief Serialise (save) via cereal
    template <typename archive_t>
    void CEREAL_SAVE_FUNCTION_NAME(archive_t & ar) const
    {
        ar(CEREAL_NVP(m_size));
        ar(CEREAL_NVP(m_block_size));
        ar(CEREAL_NVP(m_blocks));
        ar(CEREAL_NVP(m_sigma));
        ar(CEREAL_NVP(m_e));
        ar(CEREAL_NVP(m_bv_blocks));
        ar(CEREAL_NVP(m_bv_blocks_select0));
        ar(CEREAL_NVP(m_bv_blocks_select1));
    }

    //!\brief Load via cereal
    template <typename archive_t>
    void CEREAL_LOAD_FUNCTION_NAME(archive_t & ar)
    {
        ar(CEREAL_NVP(m_size));
        ar(CEREAL_NVP(m_block_size));
        ar(CEREAL_NVP(m_blocks));
        ar(CEREAL_NVP(m_sigma));
        ar(CEREAL_NVP(m_e));
        ar(CEREAL_NVP(m_bv_blocks));
        ar(CEREAL_NVP(m_bv_blocks_select0));
        m_bv_blocks_select0.set_vector(&m_bv_blocks);
        ar(CEREAL_NVP(m_bv_blocks_select1));
        m_bv_blocks_select1.set_vector(&m_bv_blocks);
    }

    iterator begin()
    {
        return {this, 0};
    };
    const_iterator end()
    {
        return {this, size()};
    };
    iterator begin() const
    {
        return {this, 0};
    };
    const_iterator end() const
    {
        return {this, size()};
    };

    //! Equality operator.
    bool operator==(wt_gmr_rs const & other) const noexcept
    {
        return (m_size == other.m_size) && (m_block_size == other.m_block_size) && (m_blocks == other.m_blocks)
            && (m_sigma == other.m_sigma) && (m_e == other.m_e) && (m_bv_blocks == other.m_bv_blocks)
            && (m_bv_blocks_select0 == other.m_bv_blocks_select0) && (m_bv_blocks_select1 == other.m_bv_blocks_select1);
    }

    //! Inequality operator.
    bool operator!=(wt_gmr_rs const & other) const noexcept
    {
        return !(*this == other);
    }
};

//! A wavelet tree class for integer sequences.
/*!
 *  \tparam t_rac         Type of the random access container used for storing the permutation.
 *  \tparam t_inv_support Type of the support structure for inverse permutation
 *  \tparam t_bitvector   Type of the bitvector used for storing B and X.
 *  \tparam t_select      Type of the support structure for select on pattern `1`.
 *  \tparam t_select_zero Type of the support structure for select on pattern `0`.
 *
 * This is an implementation of the second proposal in the SODA paper of Golynski et. al.
 * which supports fast access, inverse select, rank, and select.
 *
 * \par References
 * [1] A. Golynski, J. Munro and S. Rao:
 *     ,,Rank/select operations on large alphabets: a tool for text indexing''
 *     Proceedings of SODA 2006.
 *
 *   @ingroup wt
 */
template <class t_rac = int_vector<>,
          class t_inverse_support = inv_multi_perm_support<32, t_rac>,
          class t_bitvector = bit_vector,
          class t_select = typename t_bitvector::select_1_type,
          class t_select_zero = typename t_bitvector::select_0_type>
class wt_gmr
{
public:
    typedef typename t_rac::size_type size_type;
    typedef typename t_rac::value_type value_type;
    typedef typename t_bitvector::difference_type difference_type;
    typedef random_access_const_iterator<wt_gmr> const_iterator;
    typedef const_iterator iterator;
    typedef wt_tag index_category;
    typedef int_alphabet_tag alphabet_category;
    enum
    {
        lex_ordered = 0
    };

private:
    t_bitvector m_bv_blocks; // 0 indicates end of block. Corresponds to B in the paper.
    t_bitvector m_bv_chunks; // 0 indicates end of symbol in chunk. Corresponds to X in the paper.

    t_rac m_perm;            // Contains permutation of each chunk. Corresponds to \f$ \pi \f$ in the paper.
    t_inverse_support m_ips; // Support for inverse permutation

    t_select m_bv_blocks_select1, m_bv_chunks_select1;
    t_select_zero m_bv_blocks_select0, m_bv_chunks_select0;

    uint64_t m_size;           // input length
    uint64_t m_max_symbol = 0; // maximum character + 1
    uint64_t m_chunks;         // number of chunks
    uint64_t m_chunksize;
    uint64_t m_sigma = 0;

public:
    size_type const & sigma = m_sigma;

    //! Default constructor
    wt_gmr() = default;

    //! Construct the wavelet tree from a sequence defined by two interators
    /*!
     * \param begin   Iterator to the start of the input.
     * \param end     Iterator one past the end of the input.
     * \param tmp_dir Temporary directory for intermediate results.
     *    \par Time complexity
     *        \f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
     *        I.e. we need \Order{n\log n} if rac is a permutation of 0..n-1.
     */
    template <typename t_it>
    wt_gmr(t_it begin, t_it end, std::string tmp_dir = ram_file_name("")) : m_size(std::distance(begin, end))
    {
        // Determine max. symbol
        for (auto it = begin; it != end; ++it)
        {
            value_type value = *it;
            if (m_max_symbol < value)
                m_max_symbol = value;
        }
        ++m_max_symbol;
        m_chunksize = (1ULL << (bits::hi(m_max_symbol - 1) + 1)); // In some cases this is better than m_max_smbol
        m_chunks = (m_size + m_chunksize - 1) / m_chunksize;

        // calc m_bv_blocks
        {
            bit_vector b(m_size + m_max_symbol * m_chunks + 1, 0);
            int_vector<> tmp(m_max_symbol * m_chunks, 0, bits::hi(m_max_symbol - 1) + 2);

            uint64_t offset = 0, j = 0;
            for (auto it = begin; it != end; ++it, ++j)
            {
                if (j == m_chunksize)
                {
                    ++offset;
                    j = 0;
                }
                ++tmp[(*it) * m_chunks + offset];
            }

            for (uint64_t i = 0, l = 1; i < tmp.size(); ++i, ++l)
                for (uint64_t j = 0; j < tmp[i]; ++j)
                    b[l++] = 1;

            // calc m_sigma
            bool write = true;
            uint64_t blocks = 0;
            for (uint64_t i = 1; i < b.size(); ++i)
            {
                if (blocks == m_chunks)
                {
                    blocks = 0;
                    write = true;
                }
                if (b[i])
                {
                    if (write)
                    {
                        ++m_sigma;
                        write = false;
                    }
                }
                else
                    ++blocks;
            }

            m_bv_blocks = t_bitvector(std::move(b));
        }

        // Calc perm and bv_chunks
        {
            uint64_t x_pos = 0;
            bit_vector x(m_size + m_chunks * m_max_symbol + 1, 0);

            // fill perm and m_bv_chunks for every chunk
            int_vector<> perm(m_size, 0, bits::hi(m_max_symbol - 1) + 1);
            for (uint64_t i = 0; i < m_chunks; ++i)
            {
                int_vector<> symbols(m_max_symbol, 0, bits::hi(m_max_symbol - 1) + 2);

                // calc symbols
                for (uint64_t j = i * m_chunksize; j < (i + 1) * m_chunksize and j < m_size; ++j)
                {
                    ++symbols[*(begin + j)];
                }
                // calc m_bv_chunks
                for (uint64_t j = 0; j < m_max_symbol; ++j, ++x_pos)
                    for (uint64_t k = 0; k < symbols[j]; ++k)
                        x[++x_pos] = 1;

                // calc symbols prefix sum
                for (uint64_t j = 0, tmp = 0, sum = 0; j < m_max_symbol; ++j)
                {
                    tmp = symbols[j];
                    symbols[j] = sum;
                    sum += tmp;
                }
                // calc perm
                for (uint64_t j = i * m_chunksize, k = 0; j < (i + 1) * m_chunksize and j < m_size; ++j, ++k)
                {
                    perm[i * m_chunksize + (symbols[*(begin + j)]++)] = k;
                }
            }
            m_bv_chunks = t_bitvector(std::move(x));
            m_ips = t_inverse_support(&m_perm, perm, m_chunksize);
            _transform_to_compressed<t_rac>(perm, m_perm, tmp_dir);
            m_ips.set_vector(&m_perm);
        }
        util::init_support(m_bv_chunks_select1, &m_bv_chunks);
        util::init_support(m_bv_chunks_select0, &m_bv_chunks);
        util::init_support(m_bv_blocks_select1, &m_bv_blocks);
        util::init_support(m_bv_blocks_select0, &m_bv_blocks);
    }

    //! Copy constructor
    wt_gmr(wt_gmr const & wt) :
        m_bv_blocks(wt.m_bv_blocks),
        m_bv_chunks(wt.m_bv_chunks),
        m_perm(wt.m_perm),
        m_ips(wt.m_ips),
        m_bv_blocks_select1(wt.m_bv_blocks_select1),
        m_bv_chunks_select1(wt.m_bv_chunks_select1),
        m_bv_blocks_select0(wt.m_bv_blocks_select0),
        m_bv_chunks_select0(wt.m_bv_chunks_select0),
        m_size(wt.m_size),
        m_max_symbol(wt.m_max_symbol),
        m_chunks(wt.m_chunks),
        m_chunksize(wt.m_chunksize),
        m_sigma(wt.m_sigma)
    {
        m_ips.set_vector(&m_perm);
        m_bv_blocks_select1.set_vector(&m_bv_blocks);
        m_bv_chunks_select1.set_vector(&m_bv_chunks);
        m_bv_blocks_select0.set_vector(&m_bv_blocks);
        m_bv_chunks_select0.set_vector(&m_bv_chunks);
    }

    //! Move constructor
    wt_gmr(wt_gmr && wt) :
        m_bv_blocks(std::move(wt.m_bv_blocks)),
        m_bv_chunks(std::move(wt.m_bv_chunks)),
        m_perm(std::move(wt.m_perm)),
        m_ips(std::move(wt.m_ips)),
        m_bv_blocks_select1(std::move(wt.m_bv_blocks_select1)),
        m_bv_chunks_select1(std::move(wt.m_bv_chunks_select1)),
        m_bv_blocks_select0(std::move(wt.m_bv_blocks_select0)),
        m_bv_chunks_select0(std::move(wt.m_bv_chunks_select0)),
        m_size(wt.m_size),
        m_max_symbol(wt.m_max_symbol),
        m_chunks(wt.m_chunks),
        m_chunksize(wt.m_chunksize),
        m_sigma(wt.m_sigma)
    {
        m_ips.set_vector(&m_perm);
        m_bv_blocks_select1.set_vector(&m_bv_blocks);
        m_bv_chunks_select1.set_vector(&m_bv_chunks);
        m_bv_blocks_select0.set_vector(&m_bv_blocks);
        m_bv_chunks_select0.set_vector(&m_bv_chunks);
    }

    //! Assignment operator
    wt_gmr & operator=(wt_gmr const & wt)
    {
        wt_gmr tmp(wt);
        *this = std::move(tmp);
        return *this;
    }

    //! Move assignment operator
    wt_gmr & operator=(wt_gmr && wt)
    {
        m_bv_blocks = std::move(wt.m_bv_blocks);
        m_bv_chunks = std::move(wt.m_bv_chunks);
        m_perm = std::move(wt.m_perm);
        m_ips = std::move(wt.m_ips);
        m_ips.set_vector(&m_perm);
        m_bv_blocks_select1 = std::move(wt.m_bv_blocks_select1);
        m_bv_blocks_select1.set_vector(&m_bv_blocks);
        m_bv_chunks_select1 = std::move(wt.m_bv_chunks_select1);
        m_bv_chunks_select1.set_vector(&m_bv_chunks);
        m_bv_blocks_select0 = std::move(wt.m_bv_blocks_select0);
        m_bv_blocks_select0.set_vector(&m_bv_blocks);
        m_bv_chunks_select0 = std::move(wt.m_bv_chunks_select0);
        m_bv_chunks_select0.set_vector(&m_bv_chunks);
        m_size = wt.m_size;
        m_max_symbol = wt.m_max_symbol;
        m_chunks = wt.m_chunks;
        m_chunksize = wt.m_chunksize;
        m_sigma = wt.m_sigma;
        return *this;
    }

    //! Returns the size of the original vector.
    size_type size() const
    {
        return m_size;
    }

    //! Returns whether the wavelet tree contains no data.
    bool empty() const
    {
        return m_size == 0;
    }

    //! Recovers the i-th symbol of the original vector.
    /*!\param i The index of the symbol in the original vector.
     *  \returns The i-th symbol of the original vector.
     *  \par Time complexity
     *       \f$ \Order{1} + 1 Access to the inverse permutation \f$
     *  \par Precondition
     *       \f$ i < size() \f$
     */
    value_type operator[](size_type i) const
    {
        assert(i < size());
        uint64_t chunk = i / m_chunksize;
        uint64_t x = m_ips[i];
        return m_bv_chunks_select1(x + 1) - x - (chunk * m_max_symbol) - 1;
    }

    //! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
    /*!
     *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
     *  \param c The symbol to count the occurrences in the prefix.
     *    \returns The number of occurrences of symbol c in the prefix [0..i-1] of the supported vector.
     *  \par Time complexity
     *       \f$ \Order{\log |\Sigma|} \f$
     *  \par Precondition
     *       \f$ i \leq size() \f$
     */
    size_type rank(size_type i, value_type c) const
    {
        assert(i <= size());

        if (0 == i or c > m_max_symbol - 1)
        {
            return 0;
        }

        uint64_t chunk = (i - 1) / m_chunksize;
        uint64_t ones_before_c = m_bv_blocks_select0(c * m_chunks + 1) - (c * m_chunks + 1) + 1;
        uint64_t c_ones_before_chunk =
            m_bv_blocks_select0(c * m_chunks + chunk + 1) - (c * m_chunks + chunk + 1) + 1 - ones_before_c;

        uint64_t c_ones_in_chunk = 0;
        auto begin =
            m_perm.begin() + m_bv_chunks_select0(chunk * m_max_symbol + 1 + c) - (chunk * m_max_symbol + 1 + c) + 1;
        auto end =
            m_perm.begin() + m_bv_chunks_select0(chunk * m_max_symbol + 2 + c) - (chunk * m_max_symbol + 2 + c) + 1;

        size_type val = (i - 1) % m_chunksize;
        if (end - begin < 50)
        { // After a short test, this seems to be a good threshold
            c_ones_in_chunk = std::find_if(begin,
                                           end,
                                           [&val](auto const && x)
                                           {
                                               return x > val;
                                           })
                            - begin;
        }
        else
        {
            c_ones_in_chunk = std::lower_bound(begin, end, val + 1) - begin;
        }
        return c_ones_before_chunk + c_ones_in_chunk;
    }

    //! Calculates how many occurrences of symbol input[i] are in the prefix [0..i-1] of the original input.
    /*!
     *  \param i The index of the symbol.
     *  \return  Pair (rank(input[i],i), input[i])
     *  \par Time complexity
     *       \f$ \Order{1} + One access to the inverse permutation \f$
     *  \par Precondition
     *       \f$ i < size() \f$
     */
    std::pair<size_type, value_type> inverse_select(size_type i) const
    {
        assert(i < size());
        uint64_t chunk = i / m_chunksize;
        uint64_t x = m_ips[i];
        uint64_t tmp = m_bv_chunks_select1(x + 1);
        uint64_t c = tmp - x - (chunk * m_max_symbol) - 1;

        uint64_t ones_before_c = m_bv_blocks_select0(c * m_chunks + 1) - (c * m_chunks + 1) + 1;
        uint64_t c_before_chunk =
            m_bv_blocks_select0(c * m_chunks + chunk + 1) - (c * m_chunks + chunk + 1) + 1 - ones_before_c;
        uint64_t c_in_chunk = tmp - m_bv_chunks_select0(c + 1 + chunk * m_max_symbol) - 1;
        return std::make_pair(c_before_chunk + c_in_chunk, c);
    }

    //! Calculates the i-th occurrence of the symbol c in the supported vector.
    /*!
     *  \param i The i-th occurrence.
     *  \param c The symbol c.
     *  \par Time complexity
     *       \f$ \Order{1} \f$
     *  \par Precondition
     *       \f$ 1 \leq i \leq rank(size(), c) \f$
     */
    size_type select(size_type i, value_type c) const
    {
        assert(1 <= i and i <= rank(size(), c));

        uint64_t ones_before_c = m_bv_blocks_select0(c * m_chunks + 1) - (c * m_chunks);
        uint64_t chunk = m_bv_blocks_select1(ones_before_c + i) - ones_before_c - (c * m_chunks + 1) - i + 1;
        uint64_t c_ones_before_chunk =
            m_bv_blocks_select0(c * m_chunks + chunk + 1) - (c * m_chunks + chunk) - ones_before_c;
        uint64_t pi_pos = m_bv_chunks_select0(chunk * m_max_symbol + c + 1) + (i - c_ones_before_chunk)
                        - chunk * m_max_symbol - c - 1;

        return m_perm[pi_pos] + chunk * m_chunksize;
    }

    //! Serializes the data structure into the given ostream
    size_type serialize(std::ostream & out, structure_tree_node * v = nullptr, std::string name = "") const
    {
        structure_tree_node * child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += write_member(m_size, out, child, "size");
        written_bytes += write_member(m_max_symbol, out, child, "max_symbol");
        written_bytes += write_member(m_chunks, out, child, "chunks");
        written_bytes += write_member(m_chunksize, out, child, "chunksize");
        written_bytes += write_member(m_sigma, out, child, "sigma");
        written_bytes += m_bv_blocks.serialize(out, child, "bv_blocks");
        written_bytes += m_bv_blocks_select0.serialize(out, child, "bv_blocks_select0");
        written_bytes += m_bv_blocks_select1.serialize(out, child, "bv_blocks_select1");
        written_bytes += m_bv_chunks.serialize(out, child, "bv_chunks");
        written_bytes += m_bv_chunks_select0.serialize(out, child, "bv_chunks_select0");
        written_bytes += m_bv_chunks_select1.serialize(out, child, "bv_chunks_select1");
        written_bytes += m_perm.serialize(out, child, "permutation");
        written_bytes += m_ips.serialize(out, child, "inverse_permutation_support");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    //! Loads the data structure from the given istream.
    void load(std::istream & in)
    {
        read_member(m_size, in);
        read_member(m_max_symbol, in);
        read_member(m_chunks, in);
        read_member(m_chunksize, in);
        read_member(m_sigma, in);
        m_bv_blocks.load(in);
        m_bv_blocks_select0.load(in, &m_bv_blocks);
        m_bv_blocks_select1.load(in, &m_bv_blocks);
        m_bv_chunks.load(in);
        m_bv_chunks_select0.load(in, &m_bv_chunks);
        m_bv_chunks_select1.load(in, &m_bv_chunks);
        m_perm.load(in);
        m_ips.load(in, &m_perm);
    }

    //!\brief Serialise (save) via cereal
    template <typename archive_t>
    void CEREAL_SAVE_FUNCTION_NAME(archive_t & ar) const
    {
        ar(CEREAL_NVP(m_size));
        ar(CEREAL_NVP(m_max_symbol));
        ar(CEREAL_NVP(m_chunks));
        ar(CEREAL_NVP(m_chunksize));
        ar(CEREAL_NVP(m_sigma));
        ar(CEREAL_NVP(m_bv_blocks));
        ar(CEREAL_NVP(m_bv_blocks_select0));
        ar(CEREAL_NVP(m_bv_blocks_select1));
        ar(CEREAL_NVP(m_bv_chunks));
        ar(CEREAL_NVP(m_bv_chunks_select0));
        ar(CEREAL_NVP(m_bv_chunks_select1));
        ar(CEREAL_NVP(m_perm));
        ar(CEREAL_NVP(m_ips));
    }

    //!\brief Load via cereal
    template <typename archive_t>
    void CEREAL_LOAD_FUNCTION_NAME(archive_t & ar)
    {
        ar(CEREAL_NVP(m_size));
        ar(CEREAL_NVP(m_max_symbol));
        ar(CEREAL_NVP(m_chunks));
        ar(CEREAL_NVP(m_chunksize));
        ar(CEREAL_NVP(m_sigma));
        ar(CEREAL_NVP(m_bv_blocks));
        ar(CEREAL_NVP(m_bv_blocks_select0));
        m_bv_blocks_select0.set_vector(&m_bv_blocks);
        ar(CEREAL_NVP(m_bv_blocks_select1));
        m_bv_blocks_select1.set_vector(&m_bv_blocks);
        ar(CEREAL_NVP(m_bv_chunks));
        ar(CEREAL_NVP(m_bv_chunks_select0));
        m_bv_chunks_select0.set_vector(&m_bv_chunks);
        ar(CEREAL_NVP(m_bv_chunks_select1));
        m_bv_chunks_select1.set_vector(&m_bv_chunks);
        ar(CEREAL_NVP(m_perm));
        ar(CEREAL_NVP(m_ips));
        m_ips.set_vector(&m_perm);
    }

    iterator begin()
    {
        return {this, 0};
    };
    const_iterator end()
    {
        return {this, size()};
    };
    iterator begin() const
    {
        return {this, 0};
    };
    const_iterator end() const
    {
        return {this, size()};
    };

    //! Equality operator.
    bool operator==(wt_gmr const & other) const noexcept
    {
        return (m_size == other.m_size) && (m_max_symbol == other.m_max_symbol) && (m_chunks == other.m_chunks)
            && (m_chunksize == other.m_chunksize) && (m_sigma == other.m_sigma) && (m_bv_blocks == other.m_bv_blocks)
            && (m_bv_blocks_select0 == other.m_bv_blocks_select0) && (m_bv_blocks_select1 == other.m_bv_blocks_select1)
            && (m_bv_chunks == other.m_bv_chunks) && (m_bv_chunks_select0 == other.m_bv_chunks_select0)
            && (m_bv_chunks_select1 == other.m_bv_chunks_select1) && (m_perm == other.m_perm) && (m_ips == other.m_ips);
    }

    //! Inequality operator.
    bool operator!=(wt_gmr const & other) const noexcept
    {
        return !(*this == other);
    }
};
} // namespace sdsl

#endif
