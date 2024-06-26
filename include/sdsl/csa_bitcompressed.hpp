// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*!\file csa_bitcompressed.hpp
 * \brief csa_bitcompressed.hpp contains a bitcompressed suffix array.
 * \author Simon Gog
 */
#ifndef INCLUDED_SDSL_CSA_UNCOMPRESSED
#define INCLUDED_SDSL_CSA_UNCOMPRESSED

#include <iostream>
#include <stddef.h>
#include <stdint.h>
#include <string>

#include <sdsl/cereal.hpp>
#include <sdsl/config.hpp>
#include <sdsl/csa_alphabet_strategy.hpp>
#include <sdsl/csa_sampling_strategy.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/io.hpp>
#include <sdsl/iterators.hpp>
#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/structure_tree.hpp>
#include <sdsl/suffix_array_helper.hpp>
#include <sdsl/util.hpp>

namespace sdsl
{

//! A class for the uncompressed suffix array (SA).
/*!
 * This class stores the information of the suffix array and the inverse suffix array in uncompressed form.
 * In contrast to this class, classes like sdsl::csa_sada, and sdsl::csa_wt store
 * the suffix array and inverse suffix array data in compressed form.
 *
 * The interface of this class is exactly the same as for the compressed indexes. This is the reason
 * why it is in the group of compressed suffix arrays.
 *
 *  \tparam t_alphabet_strat  Policy for alphabet representation.
 *
 * \par Space complexity
 *        \f$ 2n\cdot \log n\f$ bits, where \f$n\f$ equals the \f$size()\f$ of the suffix array.
 * \sa sdsl::csa_sada, sdsl::csa_wt
 * @ingroup csa
 */
template <class t_alphabet_strat = byte_alphabet>
class csa_bitcompressed
{
    friend class bwt_of_csa_psi<csa_bitcompressed>;

public:
    typedef uint64_t value_type;                                            // STL Container requirement
    typedef random_access_const_iterator<csa_bitcompressed> const_iterator; // STL Container requirement
    typedef const_iterator iterator;                                        // STL Container requirement
    typedef const value_type const_reference;
    typedef const_reference reference;
    typedef const_reference * pointer;
    typedef const pointer const_pointer;
    typedef int_vector<>::size_type size_type; // STL Container requirement
    typedef size_type csa_size_type;
    typedef ptrdiff_t difference_type; // STL Container requirement
    typedef traverse_csa_saisa<csa_bitcompressed, true> psi_type;
    typedef traverse_csa_saisa<csa_bitcompressed, false> lf_type;
    typedef bwt_of_csa_psi<csa_bitcompressed> bwt_type;
    typedef text_of_csa<csa_bitcompressed> text_type;
    typedef first_row_of_csa<csa_bitcompressed> first_row_type;
    typedef _sa_order_sampling<csa_bitcompressed, 0> sa_sample_type;
    typedef _isa_sampling<csa_bitcompressed, 0> isa_sample_type;
    typedef isa_sample_type isa_type;
    typedef t_alphabet_strat alphabet_type;
    typedef typename alphabet_type::char_type char_type; // Note: This is the char type of the CSA not the WT!
    typedef typename alphabet_type::comp_char_type comp_char_type;
    typedef typename alphabet_type::string_type string_type;
    typedef typename alphabet_type::alphabet_category alphabet_category;
    typedef csa_bitcompressed csa_type;

    typedef csa_tag index_category;
    typedef psi_tag extract_category;

    enum
    {
        sa_sample_dens = 1,
        isa_sample_dens = 1
    };

private:
    sa_sample_type m_sa;   // vector for suffix array values
    isa_sample_type m_isa; // vector for inverse suffix array values
    alphabet_type m_alphabet;

public:
    const typename alphabet_type::char2comp_type & char2comp = m_alphabet.char2comp;
    const typename alphabet_type::comp2char_type & comp2char = m_alphabet.comp2char;
    const typename alphabet_type::C_type & C = m_alphabet.C;
    const typename alphabet_type::sigma_type & sigma = m_alphabet.sigma;
    const psi_type psi = psi_type(*this);
    const lf_type lf = lf_type(*this);
    const bwt_type bwt = bwt_type(*this);
    const bwt_type L = bwt_type(*this);
    isa_type const & isa = m_isa;
    const first_row_type F = first_row_type(*this);
    const text_type text = text_type(*this);
    sa_sample_type const & sa_sample = m_sa;
    isa_sample_type const & isa_sample = m_isa;

    //! Default constructor
    csa_bitcompressed()
    {}
    //! Copy constructor
    csa_bitcompressed(csa_bitcompressed const & csa) : m_sa(csa.m_sa), m_isa(csa.m_isa), m_alphabet(csa.m_alphabet)
    {}

    //! Move constructor
    csa_bitcompressed(csa_bitcompressed && csa)
    {
        *this = std::move(csa);
    }

    //! Constructor
    csa_bitcompressed(cache_config & config)
    {
        std::string text_file = cache_file_name(key_text<alphabet_type::int_width>(), config);
        int_vector_buffer<alphabet_type::int_width> text_buf(text_file);
        int_vector_buffer<> sa_buf(cache_file_name(conf::KEY_SA, config));
        size_type n = text_buf.size();

        m_alphabet = alphabet_type(text_buf, n);
        m_sa = sa_sample_type(config);
        m_isa = isa_sample_type(config);
    }

    //! Number of elements in the instance.
    /*! Required for the Container Concept of the STL.
     *  \sa max_size, empty
     */
    size_type size() const
    {
        return m_sa.size();
    }

    //! Returns the largest size that csa_bitcompressed can ever have.
    /*! Required for the Container Concept of the STL.
     *  \sa size
     */
    static size_type max_size()
    {
        return int_vector<>::max_size();
    }

    //! Returns if the data structure is empty.
    /*! Required for the Container Concept of the STL.
     * \sa size
     */
    bool empty() const
    {
        return m_sa.empty();
    }

    //! Returns a const_iterator to the first element.
    /*! Required for the STL Container Concept.
     *  \sa end
     */
    const_iterator begin() const
    {
        return const_iterator(this, 0);
    }

    //! Returns a const_iterator to the element after the last element.
    /*! Required for the STL Container Concept.
     *  \sa begin.
     */
    const_iterator end() const
    {
        return const_iterator(this, size());
    }

    //! []-operator
    /*!\param i Index of the value. \f$ i \in [0..size()-1]\f$.
     *
     * Required for the STL Random Access Container Concept.
     */
    inline value_type operator[](size_type i) const
    {
        return m_sa[i];
    }

    //! Assignment Operator.
    /*!
     *    Required for the Assignable Concept of the STL.
     */
    csa_bitcompressed & operator=(csa_bitcompressed const & csa)
    {
        if (this != &csa)
        {
            csa_bitcompressed tmp(csa);
            *this = std::move(tmp);
        }
        return *this;
    }

    //! Assignment Move Operator.
    /*!
     *    Required for the Assignable Concept of the STL.
     */
    csa_bitcompressed & operator=(csa_bitcompressed && csa)
    {
        if (this != &csa)
        {
            m_sa = std::move(csa.m_sa);
            m_isa = std::move(csa.m_isa);
            m_alphabet = std::move(csa.m_alphabet);
        }
        return *this;
    }

    //! Equality operator.
    bool operator==(csa_bitcompressed const & other) const noexcept
    {
        return (m_sa == other.m_sa) && (m_isa == other.m_isa) && (m_alphabet == other.m_alphabet);
    }

    //! Inequality operator.
    bool operator!=(csa_bitcompressed const & other) const noexcept
    {
        return !(*this == other);
    }

    //! Serialize to a stream.
    /*!\param out Output stream to write the data structure.
     *  \return The number of written bytes.
     */
    size_type serialize(std::ostream & out, structure_tree_node * v = nullptr, std::string name = "") const
    {
        structure_tree_node * child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += m_sa.serialize(out, child, "m_sa");
        written_bytes += m_isa.serialize(out, child, "m_isa");
        written_bytes += m_alphabet.serialize(out, child, "m_alphabet");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void load(std::istream & in)
    {
        m_sa.load(in);
        m_isa.load(in);
        m_alphabet.load(in);
    }

    template <typename archive_t>
    void CEREAL_SAVE_FUNCTION_NAME(archive_t & ar) const
    {
        ar(CEREAL_NVP(m_sa));
        ar(CEREAL_NVP(m_isa));
        ar(CEREAL_NVP(m_alphabet));
    }

    template <typename archive_t>
    void CEREAL_LOAD_FUNCTION_NAME(archive_t & ar)
    {
        ar(CEREAL_NVP(m_sa));
        ar(CEREAL_NVP(m_isa));
        ar(CEREAL_NVP(m_alphabet));
    }

    size_type get_sample_dens() const
    {
        return 1;
    }

private:
    // Calculates how many symbols c are in the prefix [0..i-1] of the BWT of the original text.
    /*
     *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in [0..size()]\f$.
     *  \param c The symbol to count the occurrences in the prefix.
     *    \returns The number of occurrences of symbol c in the prefix [0..i-1] of the BWT.
     *  \par Time complexity
     *        \f$ \Order{\log n} \f$
     */
    size_type rank_bwt(size_type i, const char_type c) const
    {
        // TODO: special case if c == BWT[i-1] we can use LF to get a constant time answer
        comp_char_type cc = char2comp[c];
        if (cc == 0 and c != 0) // character is not in the text => return 0
            return 0;
        // binary search the interval [C[cc]..C[cc+1]-1] for the result
        size_type lower_b = C[cc],
                  upper_b = C[((size_type)1) + cc]; // lower_b inclusive, upper_b exclusive
        while (lower_b + 1 < upper_b)
        {
            size_type mid = (lower_b + upper_b) / 2;
            if (psi[mid] >= i)
                upper_b = mid;
            else
                lower_b = mid;
        }
        if (lower_b > C[cc])
            return lower_b - C[cc] + 1;
        else
        {                            // lower_b == m_C[cc]
            return psi[lower_b] < i; // 1 if psi[lower_b]<i, 0 otherwise
        }
    }

    // Calculates the i-th occurrence of symbol c in the BWT of the original text.
    /*
     *  \param i The i-th occurrence. \f$i\in [1..rank(size(),c)]\f$.
     *  \param c Character c.
     *    \returns The i-th occurrence of c in the BWT or size() if c does
     *           not occur t times in BWT>
     *  \par Time complexity
     *        \f$ \Order{t_{\Psi}} \f$
     */
    size_type select_bwt(size_type i, const char_type c) const
    {
        comp_char_type cc = char2comp[c];
        if (cc == 0 and c != 0) // character is not in the text => return size()
            return size();
        if (C[cc] + i - 1 < C[((size_type)1) + cc])
        {
            return psi[C[cc] + i - 1];
        }
        return size();
    }
};

} // end namespace sdsl
#endif
