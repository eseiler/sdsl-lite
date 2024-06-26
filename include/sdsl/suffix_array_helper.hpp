// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*!\file suffix_array_helper.hpp
 * \brief suffix_array_helper.hpp contains some helper classes for CSTs
 * \author Simon Gog
 */
#ifndef INCLUDED_SDSL_SUFFIX_ARRAY_HELPER
#define INCLUDED_SDSL_SUFFIX_ARRAY_HELPER

#include <cassert>

#include <sdsl/iterators.hpp>
#include <sdsl/sdsl_concepts.hpp>

namespace sdsl
{

//! Get the symbol at position i in the first row of the sorted suffixes of CSA
/*
 * \param i   Position in the first row.
 * \param csa CSA
 * \par Time complexity
 *    \f$ \Order{\log \sigma} \f$
 *  TODO: add hinted binary search? Two way binary search?
 */
template <typename t_csa>
typename t_csa::char_type first_row_symbol(const typename t_csa::size_type i, t_csa const & csa)
{
    assert(i < csa.size());
    if (csa.sigma < 16)
    { //<- if sigma is small search linear
        typename t_csa::size_type res = 1;
        while (res < csa.sigma and csa.C[res] <= i)
            ++res;
        return csa.comp2char[res - 1];
    }
    else
    {
        // binary search the character with C
        typename t_csa::size_type upper_c = csa.sigma,
                                  lower_c = 0; // lower_c inclusive, upper_c exclusive
        typename t_csa::size_type res = 0;
        do
        {
            res = (upper_c + lower_c) / 2;
            if (i < csa.C[res])
            {
                upper_c = res;
            }
            else if (i >= csa.C[res + 1])
            {
                lower_c = res + 1;
            }
        }
        while (i < csa.C[res] or i >= csa.C[res + 1]); // i is not in the interval
        return csa.comp2char[res];
    }
}

// psi[] trait
template <typename t_csa, bool t_direction>
struct traverse_csa_psi_trait
{
    typedef typename t_csa::value_type value_type;
    typedef typename t_csa::size_type size_type;
    static value_type access(t_csa const & csa, size_type i)
    {
        return csa.psi[i];
    }
};

// lf[] trait
template <typename t_csa>
struct traverse_csa_psi_trait<t_csa, false>
{
    typedef typename t_csa::value_type value_type;
    typedef typename t_csa::size_type size_type;
    static value_type access(t_csa const & csa, size_type i)
    {
        // TODO: in case of a very sparse sampling of SA it may be faster to
        //  use \sigma binary searches on PSI function to determine the
        // LF values.
        return csa.isa[(csa[i] + csa.size() - 1) % csa.size()];
    }
};

template <typename t_csa, bool t_direction>
class traverse_csa_psi
{
public:
    typedef typename t_csa::value_type value_type;
    typedef typename t_csa::size_type size_type;
    typedef typename t_csa::difference_type difference_type;
    typedef random_access_const_iterator<traverse_csa_psi> const_iterator;
    typedef csa_member_tag category;
    typedef int_alphabet_tag alphabet_category;

private:
    t_csa const & m_csa;

public:
    //! Constructor
    traverse_csa_psi(t_csa const & csa_psi) : m_csa(csa_psi)
    {}
    //! Copy constructor
    traverse_csa_psi(traverse_csa_psi const & tcsa) : m_csa(tcsa.m_csa)
    {}

    //! Calculate the \f$\Psi\f$ or \f$LF\f$ value at position i.
    /*!\param i The index for which the \f$\Psi\f$ or \f$LF\f$ value should be calculated, \f$i\in [0..size()-1]\f$.
     */
    value_type operator[](size_type i) const
    {
        assert(i < size());
        return traverse_csa_psi_trait<t_csa, t_direction>::access(m_csa, i);
    }

    //! Returns the size of the \f$\Psi\f$ function.
    size_type size() const
    {
        return m_csa.size();
    }

    //! Returns if the \f$\Psi\f$ function is empty.
    size_type empty() const
    {
        return m_csa.empty();
    }

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
};

// psi[] trait
template <typename t_csa, bool t_direction>
struct traverse_csa_saisa_trait
{
    typedef typename t_csa::value_type value_type;
    typedef typename t_csa::size_type size_type;
    static value_type access(t_csa const & csa, size_type i)
    {
        // \f$\Psi[i] = SA^{-1}[SA[i]+1 \mod n]\f$, where \f$n\f$ is the length of the suffix array SA
        return csa.isa[(csa[i] + 1) % csa.size()];
    }
};

// lf[] trait
template <typename t_csa>
struct traverse_csa_saisa_trait<t_csa, false>
{
    typedef typename t_csa::value_type value_type;
    typedef typename t_csa::size_type size_type;
    static value_type access(t_csa const & csa, size_type i)
    {
        // TODO: in case of a very sparse sampling of SA it may be faster to
        //  use \sigma binary searches on PSI function to determine the
        // LF values.
        return csa.isa[(csa[i] + csa.size() - 1) % csa.size()];
    }
};

//! A helper class for the \f$\Psi\f$ function for (compressed) suffix arrays which provide also the inverse suffix
//! array values (like sdsl::csa_bitcompressed).
template <typename t_csa, bool t_direction>
class traverse_csa_saisa
{
public:
    typedef typename t_csa::value_type value_type;
    typedef typename t_csa::size_type size_type;
    typedef typename t_csa::difference_type difference_type;
    typedef random_access_const_iterator<traverse_csa_saisa> const_iterator; // STL Container requirement
    typedef csa_member_tag category;
    typedef int_alphabet_tag alphabet_category;

private:
    t_csa const & m_csa;

public:
    //! Constructor
    traverse_csa_saisa(t_csa const & csa) : m_csa(csa)
    {}

    // Copy constructor
    traverse_csa_saisa(traverse_csa_saisa const & tcsa) : m_csa(tcsa.m_csa)
    {}

    //! Calculate the \f$\Psi\f$ value at position i.
    /*!	\param i The index for which the \f$\Psi\f$ value should be calculated, \f$i\in [0..size()-1]\f$.
     *	\par Time complexity
     *		\f$ \Order{\saaccess+\isaaccess} \f$
     */
    value_type operator[](size_type i) const
    {
        assert(i < size());
        return traverse_csa_saisa_trait<t_csa, t_direction>::access(m_csa, i);
    }

    //! Returns the size of the \f$\Psi\f$ function.
    size_type size() const
    {
        return m_csa.size();
    }

    //! Returns if the \f$\Psi\f$ function is empty.
    size_type empty() const
    {
        return m_csa, empty();
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
};

//! A wrapper for the bwt of a compressed suffix array that is based on the \f$\psi\f$ function.
template <typename t_csa>
class bwt_of_csa_psi
{
public:
    typedef typename t_csa::char_type value_type;
    typedef typename t_csa::size_type size_type;
    typedef typename t_csa::char_type char_type;
    typedef typename t_csa::difference_type difference_type;
    typedef random_access_const_iterator<bwt_of_csa_psi> const_iterator;
    typedef csa_member_tag category;
    typedef typename t_csa::alphabet_category alphabet_category;

private:
    t_csa const & m_csa; //<- pointer to the (compressed) suffix array that is based on the \f$\Psi\f$ function.

public:
    //! Constructor
    bwt_of_csa_psi(t_csa const & csa) : m_csa(csa)
    {}

    //! Calculate the Burrows Wheeler Transform (BWT) at position i.
    /*!	\param i The index for which the BWT value should be calculated, \f$i\in [0..size()-1]\f$.
     *	\par Time complexity
     *		\f$ \Order{\log |\Sigma|} \f$
     */
    value_type operator[](size_type i) const
    {
        assert(i < size());
        size_type pos = m_csa.lf[i];
        return first_row_symbol(pos, m_csa);
    }

    //! Calculates how many symbols c are in the prefix [0..i-1]
    /*!
     *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in [0..size()]\f$.
     *  \param c The symbol to count the occurrences in the prefix.
     *    \returns The number of occurrences of symbol c in the prefix [0..i-1].
     *  \par Time complexity
     *        \f$ \Order{\log n t_{\Psi}} \f$
     */
    size_type rank(size_type i, const char_type c) const
    {
        return m_csa.rank_bwt(i, c);
    }

    //! Calculates the position of the i-th c.
    /*!
     *  \param i The i-th occurrence. \f$i\in [1..rank(size(),c)]\f$.
     *  \param c Symbol c.
     *    \returns The position of the i-th c or size() if c does occur less then i times.
     *  \par Time complexity
     *        \f$ \Order{t_{\Psi}} \f$
     */
    size_type select(size_type i, const char_type c) const
    {
        return m_csa.select_bwt(i, c);
    }

    //! Returns the size of the \f$\Psi\f$ function.
    size_type size() const
    {
        return m_csa.size();
    }

    //! Returns if the bwt is empty.
    size_type empty() const
    {
        return m_csa.empty();
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
};

// psi[] trait
template <typename t_csa, bool t_direction>
struct traverse_csa_wt_traits
{
    typedef typename t_csa::value_type value_type;
    typedef typename t_csa::char_type char_type;
    typedef typename t_csa::size_type size_type;
    static value_type access(t_csa const & csa, size_type i)
    {
        char_type c = csa.F[i];
        return csa.wavelet_tree.select(i - csa.C[csa.char2comp[c]] + 1, c);
    }
};

// lf[] trait
template <typename t_csa>
struct traverse_csa_wt_traits<t_csa, false>
{
    typedef typename t_csa::value_type value_type;
    typedef typename t_csa::char_type char_type;
    typedef typename t_csa::size_type size_type;
    static value_type access(t_csa const & csa, size_type i)
    {
        typename t_csa::char_type c;
        auto rc = csa.wavelet_tree.inverse_select(i);
        size_type j = rc.first;
        c = rc.second;
        return csa.C[csa.char2comp[c]] + j;
    }
};

//! A wrapper class for the \f$\Psi\f$ and LF function for (compressed) suffix arrays that are based on a wavelet tree
//! (like sdsl::csa_wt).
template <typename t_csa, bool t_direction>
class traverse_csa_wt
{
public:
    typedef typename t_csa::value_type value_type;
    typedef typename t_csa::size_type size_type;
    typedef typename t_csa::char_type char_type;
    typedef typename t_csa::difference_type difference_type;
    typedef random_access_const_iterator<traverse_csa_wt> const_iterator;
    typedef csa_member_tag category;
    typedef int_alphabet_tag alphabet_category;

private:
    t_csa const & m_csa; //<- pointer to the (compressed) suffix array that is based on a wavelet tree
    traverse_csa_wt(){}; // disable default constructor

public:
    //! Constructor
    traverse_csa_wt(t_csa const & csa_wt) : m_csa(csa_wt)
    {}
    //! Calculate the \f$\Psi\f$ value at position i.
    /*!	\param i The index for which the \f$\Psi\f$ value should be calculated, \f$i\in [0..size()-1]\f$.
     *	\par Time complexity
     *		\f$ \Order{\log |\Sigma|} \f$
     */
    value_type operator[](size_type i) const
    {
        assert(i < m_csa.size());
        return traverse_csa_wt_traits<t_csa, t_direction>::access(m_csa, i);
    }

    //! Returns the size of the \f$\Psi\f$ function.
    size_type size() const
    {
        return m_csa.size();
    }
    //! Returns if the \f$\Psi\f$ function is empty.
    size_type empty() const
    {
        return m_csa.empty();
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
};

template <typename t_csa>
class bwt_of_csa_wt
{
public:
    typedef const typename t_csa::char_type value_type;
    typedef typename t_csa::size_type size_type;
    typedef typename t_csa::char_type char_type;
    typedef typename t_csa::difference_type difference_type;
    typedef random_access_const_iterator<bwt_of_csa_wt> const_iterator;
    typedef csa_member_tag category;
    typedef typename t_csa::alphabet_category alphabet_category;

private:
    t_csa const & m_csa; //<- pointer to the (compressed) suffix array that is based on a wavelet tree
    bwt_of_csa_wt(){};   // disable default constructor

public:
    //! Constructor
    bwt_of_csa_wt(t_csa const & csa_wt) : m_csa(csa_wt)
    {}
    //! Calculate the Burrows Wheeler Transform (BWT) at position i.
    /*!	\param i The index for which the \f$\Psi\f$ value should be calculated, \f$i\in [0..size()-1]\f$.
     *	\par Time complexity
     *		\f$ \Order{\log |\Sigma|} \f$
     */
    value_type operator[](size_type i) const
    {
        assert(i < size());
        return m_csa.wavelet_tree[i];
    }
    //! Returns the size of the BWT function.
    size_type size() const
    {
        return m_csa.size();
    }

    //! Calculates how many symbols c are in the prefix [0..i-1].
    /*!
     *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in [0..size()]\f$.
     *  \param c The symbol to count the occurrences in the prefix.
     *    \returns The number of occurrences of symbol c in the prefix [0..i-1].
     *  \par Time complexity
     *        \f$ \Order{\log |\Sigma|} \f$
     */
    size_type rank(size_type i, const char_type c) const
    {
        return m_csa.rank_bwt(i, c);
    }

    //! Calculates the position of the i-th c.
    /*!
     *  \param i The i-th occurrence. \f$i\in [1..rank(size(),c)]\f$.
     *  \param c Symbol c.
     *    \returns The position of the i-th c or size() if c does occur less then i times.
     *  \par Time complexity
     *        \f$ \Order{t_{\Psi}} \f$
     */
    size_type select(size_type i, const char_type c) const
    {
        return m_csa.select(i, c);
    }

    //! Returns if the BWT function is empty.
    size_type empty() const
    {
        return m_csa.empty();
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
};

template <typename t_csa>
class isa_of_csa_wt
{
public:
    typedef typename t_csa::value_type value_type;
    typedef typename t_csa::size_type size_type;
    typedef typename t_csa::difference_type difference_type;
    typedef random_access_const_iterator<isa_of_csa_wt> const_iterator;
    typedef csa_member_tag category;
    typedef int_alphabet_tag alphabet_category;

private:
    t_csa const & m_csa; //<- pointer to the (compressed) suffix array that is based on a wavelet tree
    isa_of_csa_wt(){};   // disable default constructor

public:
    //! Constructor
    isa_of_csa_wt(t_csa const & csa_wt) : m_csa(csa_wt)
    {}

    //! Access operator to ISA.
    /*!\param i Index \f$i\in [0..size()-1]\f$.
     */
    value_type operator[](size_type i) const
    {
        assert(i < size());
        auto sample = m_csa.isa_sample.sample_qeq(i);
        value_type result = std::get<0>(sample);
        if (std::get<1>(sample) < i)
        {
            i = std::get<1>(sample) + m_csa.size() - i;
        }
        else
        {
            i = std::get<1>(sample) - i;
        }
        while (i--)
        {
            result = m_csa.lf[result];
        }
        return result;
    }

    //! Returns the size of the CSA.
    size_type size() const
    {
        return m_csa.size();
    }
    //! Returns if the CSA is empty.
    size_type empty() const
    {
        return m_csa.empty();
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
};

template <typename t_csa>
class isa_of_csa_psi
{
public:
    typedef typename t_csa::value_type value_type;
    typedef typename t_csa::size_type size_type;
    typedef typename t_csa::difference_type difference_type;
    typedef random_access_const_iterator<isa_of_csa_psi> const_iterator;
    typedef csa_member_tag category;
    typedef int_alphabet_tag alphabet_category;

private:
    t_csa const & m_csa; //<- pointer to the (compressed) suffix array that is based on a wavelet tree
    isa_of_csa_psi(){};  // disable default constructor

public:
    //! Constructor
    isa_of_csa_psi(t_csa const & csa_wt) : m_csa(csa_wt)
    {}

    //! Access operator to ISA.
    /*!\param i Index \f$i\in [0..size()-1]\f$.
     */
    value_type operator[](size_type i) const
    {
        assert(i < size());
        // get the rightmost sampled isa value to the left of i
        auto sample = m_csa.isa_sample.sample_leq(i);
        value_type result = std::get<0>(sample);
        i = i - std::get<1>(sample);
        while (i--)
        {
            result = m_csa.psi[result];
        }
        return result;
    }
    //! Returns the size of the CSA.
    size_type size() const
    {
        return m_csa.size();
    }
    //! Returns if the CSA is empty.
    size_type empty() const
    {
        return m_csa.empty();
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
};

template <typename t_csa>
class first_row_of_csa
{
public:
    typedef const typename t_csa::char_type value_type;
    typedef typename t_csa::size_type size_type;
    typedef typename t_csa::difference_type difference_type;
    typedef random_access_const_iterator<first_row_of_csa> const_iterator;
    typedef csa_member_tag category;
    typedef typename t_csa::alphabet_category alphabet_category;

private:
    t_csa const & m_csa;

public:
    //! Constructor
    first_row_of_csa(t_csa const & csa) : m_csa(csa)
    {}
    //! Calculate F[i]
    /*!\param i The index for which the \f$F\f$ value should be calculated, \f$i\in [0..size()-1]\f$.
     *  \par Time complexity
     *      \f$ \Order{\log |\Sigma|} \f$
     */
    value_type operator[](size_type i) const
    {
        assert(i < size());
        return first_row_symbol(i, m_csa);
    }
    //! Returns the size of the F column.
    size_type size() const
    {
        return m_csa.size();
    }
    //! Returns if the F column is empty.
    size_type empty() const
    {
        return m_csa.empty();
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
};

template <typename t_csa>
class text_of_csa
{
public:
    typedef typename t_csa::char_type value_type;
    typedef typename t_csa::size_type size_type;
    typedef typename t_csa::difference_type difference_type;
    typedef random_access_const_iterator<text_of_csa> const_iterator;
    typedef csa_member_tag category;
    typedef typename t_csa::alphabet_category alphabet_category;

private:
    t_csa const & m_csa;
    text_of_csa()
    {}

public:
    //! Constructor
    text_of_csa(t_csa const & csa) : m_csa(csa)
    {}

    //! Character at index \f$i\f$ of the original text.
    /*!	\param i Text position , \f$i\in [0..size()-1]\f$.
     *	\par Time complexity
     *		\f$ t_{ISA} \log\sigma \f$
     */
    value_type operator[](size_type i) const
    {
        assert(i < size());
        return first_row_symbol(m_csa.isa[i], m_csa);
    }

    //! Returns the size of the original text.
    size_type size() const
    {
        return m_csa.size();
    }

    //! Returns if text text has size 0.
    size_type empty() const
    {
        return m_csa.empty();
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
};
} // namespace sdsl

#endif
