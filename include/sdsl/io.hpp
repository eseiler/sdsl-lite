// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*!\file io.hpp
 * \brief io.hpp contains some methods for reading/writing sdsl structures.
 * \author Simon Gog
 */
#ifndef INCLUDED_SDSL_IO
#define INCLUDED_SDSL_IO

#include <algorithm>
#include <cctype>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <unordered_map>
#include <utility>
#include <vector>

#include <sdsl/bits.hpp>
#include <sdsl/config.hpp>
#include <sdsl/platform.hpp>
#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/sfstream.hpp>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

namespace sdsl
{
template <uint8_t = 0u>
class int_vector;

int remove(std::string const &);

template <typename T>
void load_vector(std::vector<T> &, std::istream &);

template <typename T>
uint64_t
serialize_vector(std::vector<T> const &, std::ostream &, sdsl::structure_tree_node * v = nullptr, std::string = "");

// has_serialize<X>::value is true if class X has
// implement method serialize
// Adapted solution from jrok's proposal:
// http://stackoverflow.com/questions/87372/check-if-a-class-has-a-member-function-of-a-given-signature
template <typename X>
struct has_serialize
{
    template <typename T>
    static constexpr auto check(T *) ->
        typename std::is_same<decltype(std::declval<T>().serialize(std::declval<std::ostream &>(),
                                                                   std::declval<structure_tree_node *>(),
                                                                   std::declval<std::string>())),
                              typename T::size_type>::type
    {
        return std::true_type();
    }
    template <typename>
    static constexpr std::false_type check(...)
    {
        return std::false_type();
    }
    typedef decltype(check<X>(nullptr)) type;
    static constexpr bool value = type::value;
};

// has_load<X>::value is true if class X has
// implement method load
template <typename X>
struct has_load
{
    template <typename T>
    static constexpr auto check(T *) ->
        typename std::is_same<decltype(std::declval<T>().load(std::declval<std::istream &>())), void>::type
    {
        return std::true_type();
    }
    template <typename>
    static constexpr std::false_type check(...)
    {
        return std::false_type();
    }
    typedef decltype(check<X>(nullptr)) type;
    static constexpr bool value = type::value;
};

// Writes primitive-typed variable t to stream out
template <typename T>
size_t write_member(T const & t, std::ostream & out, sdsl::structure_tree_node * v = nullptr, std::string name = "")
{
    sdsl::structure_tree_node * child = sdsl::structure_tree::add_child(v, name, util::class_name(t));
    out.write((char *)&t, sizeof(t));
    size_t written_bytes = sizeof(t);
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

// Specialization for std::string
template <>
inline size_t
write_member<std::string>(std::string const & t, std::ostream & out, sdsl::structure_tree_node * v, std::string name)
{
    structure_tree_node * child = structure_tree::add_child(v, name, util::class_name(t));
    size_t written_bytes = 0;
    written_bytes += write_member(t.size(), out, child, "length");
    out.write(t.c_str(), t.size());
    written_bytes += t.size();
    structure_tree::add_size(v, written_bytes);
    return written_bytes;
}

// Writes primitive-typed variable t to stream out
template <typename T>
void read_member(T & t, std::istream & in)
{
    in.read((char *)&t, sizeof(t));
}

// Specialization for std::string
template <>
inline void read_member<std::string>(std::string & t, std::istream & in)
{
    std::string::size_type size;
    read_member(size, in);
    char * buf = new char[size];
    in.read(buf, size);
    std::string temp(buf, size);
    delete[] buf;
    t = std::move(temp);
}

template <typename X>
typename std::enable_if<has_serialize<X>::value, typename X::size_type>::type
serialize(X const & x, std::ostream & out, structure_tree_node * v = nullptr, std::string name = "")
{
    return x.serialize(out, v, name);
}

template <typename X>
typename std::enable_if<std::is_standard_layout<X>::value && std::is_trivial<X>::value, uint64_t>::type
serialize(X const & x, std::ostream & out, structure_tree_node * v = nullptr, std::string name = "")
{
    return write_member(x, out, v, name);
}

template <typename X>
uint64_t
serialize(std::vector<X> const & x, std::ostream & out, structure_tree_node * v = nullptr, std::string name = "")
{

    return serialize(x.size(), out, v, name) + serialize_vector(x, out, v, name);
}

template <typename X>
typename std::enable_if<has_load<X>::value, void>::type load(X & x, std::istream & in)
{
    x.load(in);
}

template <typename X>
typename std::enable_if<std::is_standard_layout<X>::value && std::is_trivial<X>::value, void>::type
load(X & x, std::istream & in)
{
    read_member(x, in);
}

template <typename X>
void load(std::vector<X> & x, std::istream & in)
{
    typename std::vector<X>::size_type size;
    load(size, in);
    x.resize(size);
    load_vector(x, in);
}

//! Load sdsl-object v from a file.
/*!
 * \param v sdsl-
 * \param file Name of the serialized file.
 */
template <typename T>
bool load_from_file(T & v, std::string const & file);

//! Load an int_vector from a plain array of `num_bytes`-byte integers with X in \{0, 1,2,4,8\} from disk.
// TODO: Remove ENDIAN dependency.
template <typename t_int_vec>
bool load_vector_from_file(t_int_vec & v, std::string const & file, uint8_t num_bytes = 1, uint8_t max_int_width = 64)
{
    if ((uint8_t)0 == num_bytes)
    { // if byte size is variable read int_vector<0> from file
        return load_from_file(v, file);
    }
    else if (num_bytes == 'd')
    {
        uint64_t x = 0, max_x = 0;
        isfstream in(file, std::ios::in | std::ios::binary);
        if (!in)
        {
            return false;
        }
        else
        {
            std::vector<uint64_t> tmp;
            while (in >> x)
            {
                tmp.push_back(x);
                max_x = std::max(x, max_x);
            }
            v.width(bits::hi(max_x) + 1);
            v.resize(tmp.size());
            for (size_t i = 0; i < tmp.size(); ++i)
            {
                v[i] = tmp[i];
            }
            return true;
        }
    }
    else
    {
        off_t file_size = util::file_size(file);
        if (file_size == 0)
        {
            v.resize(0);
            return true;
        }
        if (file_size % num_bytes != 0)
        {
            throw std::logic_error("file size " + util::to_string(file_size) + " of \"" + file
                                   + "\" is not a multiple of " + util::to_string(num_bytes));
            return false;
        }
        isfstream in(file, std::ios::in | std::ios::binary);
        if (in)
        {
            v.width(std::min((int)8 * num_bytes, (int)max_int_width));
            v.resize(file_size / num_bytes);
            if (8 == t_int_vec::fixed_int_width and 1 == num_bytes)
            { // if int_vector<8> is created from byte alphabet file
                in.read((char *)v.data(), file_size);
            }
            else
            {
                size_t idx = 0;
                const size_t block_size = conf::SDSL_BLOCK_SIZE * num_bytes;
                std::vector<uint8_t> buf(block_size);
                // TODO: check for larger alphabets with num_bytes*8 = v::fixed_int_width

                uint64_t x = 0; // value
                uint8_t cur_byte = 0;
                do
                {
                    in.read((char *)buf.data(), block_size);
                    size_t read = in.gcount();
                    uint8_t * begin = buf.data();
                    uint8_t * end = begin + read;
                    while (begin < end)
                    {
                        x |= ((uint64_t)(*begin)) << (cur_byte * 8);
                        ++cur_byte;
                        if (cur_byte == num_bytes)
                        {
                            v[idx++] = x;
                            cur_byte = 0;
                            x = 0ULL;
                        }
                        ++begin;
                    }
                }
                while (idx < v.size());
                in.close();
            }
            return true;
        }
        else
        {
            return false;
        }
    }
}

//! Store a data structure to a file.
/*! The data structure has to provide a serialize function.
 *  \param v Data structure to store.
 *  \param file Name of the file where to store the data structure.
 *  \param Return if the data structure was stored successfully
 */
template <typename T>
bool store_to_file(T const & v, std::string const & file);

//! Specialization of store_to_file for a char array
inline bool store_to_file(char const * v, std::string const & file)
{
    osfstream out(file, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out)
    {
        if (util::verbose)
        {
            std::cerr << "ERROR: store_to_file(const char *v, const std::string&)" << std::endl;
            return false;
        }
    }
    uint64_t n = strlen((char const *)v);
    out.write(v, n);
    out.close();
    return true;
}

//! Specialization of store_to_file for int_vector
template <uint8_t t_width>
bool store_to_file(int_vector<t_width> const & v, std::string const & file);

//! Store an int_vector as plain int_type array to disk
template <typename int_type, typename t_int_vec>
bool store_to_plain_array(t_int_vec & v, std::string const & file)
{
    osfstream out(file, std::ios::out | std::ios::binary);
    if (out)
    {
        for (typename t_int_vec::size_type i = 0; i < v.size(); ++i)
        {
            int_type x = v[i];
            out.write((char *)&x, sizeof(int_type));
        }
        return true;
    }
    else
    {
        return false;
    }
}

template <typename T>
size_t
serialize_empty_object(std::ostream &, structure_tree_node * v = nullptr, std::string name = "", T const * t = nullptr)
{
    structure_tree_node * child = structure_tree::add_child(v, name, util::class_name(*t));
    size_t written_bytes = 0;
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

//! Get the size of a data structure in bytes.
/*!
 *  \param v A reference to the data structure for which the size in bytes should be calculated.
 */
template <typename T>
typename T::size_type size_in_bytes(T const & t);

//! Get the size of a data structure in mega bytes (MiB).
/*!
 *  \param t A reference to the data structure for which the size in bytes should be calculated.
 */
template <typename T>
double size_in_mega_bytes(T const & t);

struct nullstream : std::ostream
{
    struct nullbuf : std::streambuf
    {
        int overflow(int c)
        {
            return traits_type::not_eof(c);
        }
        int xputc(int)
        {
            return 0;
        }
        std::streamsize xsputn(char const *, std::streamsize n)
        {
            return n;
        }
        int sync()
        {
            return 0;
        }
    } m_sbuf;
    nullstream() : std::ios(&m_sbuf), std::ostream(&m_sbuf), m_sbuf()
    {}
};

//! Serialize each element of an std::vector
/*!
 * \param vec The vector which should be serialized.
 * \param out Output stream to which should be written.
 * \param v   Structure tree node. Note: If all elements have the same
 *            structure, then it is tried to combine all elements (i.e.
 *            make one node w with size set to the cumulative sum of all
 *           sizes of the children)
 */
template <typename T>
uint64_t
serialize_vector(std::vector<T> const & vec, std::ostream & out, sdsl::structure_tree_node * v, std::string name)
{
    if (vec.size() > 0)
    {
        sdsl::structure_tree_node * child =
            sdsl::structure_tree::add_child(v, name, "std::vector<" + util::class_name(vec[0]) + ">");
        size_t written_bytes = 0;
        for (auto const & x : vec)
        {
            written_bytes += serialize(x, out, child, "[]");
        }
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }
    else
    {
        return 0;
    }
}

//! Load all elements of a vector from a input stream
/*!\param vec  Vector whose elements should be loaded.
 *  \param in   Input stream.
 *  \par Note
 *   The vector has to be resized prior the loading
 *   of its elements.
 */
template <typename T>
void load_vector(std::vector<T> & vec, std::istream & in)
{
    for (typename std::vector<T>::size_type i = 0; i < vec.size(); ++i)
    {
        load(vec[i], in);
    }
}

template <format_type F, typename X>
void write_structure(X const & x, std::ostream & out)
{
    std::unique_ptr<structure_tree_node> st_node(new structure_tree_node("name", "type"));
    nullstream ns;
    serialize(x, ns, st_node.get(), "");
    if (st_node.get()->children.size() > 0)
    {
        for (auto const & child : st_node.get()->children)
        {
            sdsl::write_structure_tree<F>(child.second.get(), out);
        }
    }
}

template <format_type F, typename X>
void write_structure(X const & x, std::string file)
{
    std::ofstream out(file);
    write_structure<F>(x, out);
}

template <format_type F, typename... Xs>
void write_structure(std::ostream & out, Xs... xs)
{
    typedef std::unique_ptr<structure_tree_node> up_stn_type;
    up_stn_type st_node(new structure_tree_node("name", "type"));
    _write_structure(st_node, xs...);
    sdsl::write_structure_tree<F>(st_node.get(), out);
}

template <typename X, typename... Xs>
void _write_structure(std::unique_ptr<structure_tree_node> & st_node, X x, Xs... xs)
{
    nullstream ns;
    serialize(x, ns, st_node.get(), "");
    _write_structure(st_node, xs...);
}

inline void _write_structure(std::unique_ptr<structure_tree_node> &)
{}

//! Internal function used by csXprintf
inline uint64_t _parse_number(std::string::const_iterator & c, std::string::const_iterator const & end)
{
    std::string::const_iterator s = c;
    while (c != end and isdigit(*c))
        ++c;
    if (c > s)
    {
        return std::stoull(std::string(s, c));
    }
    else
    {
        return 0;
    }
}

//! Internal function used by csXprintf
template <typename t_csa>
t_csa const & _idx_csa(t_csa const & t, csa_tag)
{
    return t;
}

//! Internal function used by csXprintf
template <typename t_cst>
const typename t_cst::csa_type & _idx_csa(t_cst const & t, cst_tag)
{
    return t.csa;
}

//! Internal function used by csXprintf
template <typename t_csa>
std::string _idx_lcp_val(t_csa const &, uint64_t, uint64_t, csa_tag)
{
    return "";
}

//! Internal function used by csXprintf
template <typename t_cst>
std::string _idx_lcp_val(t_cst const & t, uint64_t i, uint64_t w, cst_tag)
{
    return util::to_string(t.lcp[i], w);
}

template <typename t_csx, typename t_alph = typename t_csx::alphabet_category>
struct default_sentinel
{
    static char const value = '$';
};

template <typename t_csx>
struct default_sentinel<t_csx, byte_alphabet_tag>
{
    static char const value = '$';
};

template <typename t_csx>
struct default_sentinel<t_csx, int_alphabet_tag>
{
    static char const value = '0';
};

//! Prints members of CSAs and CSTs
/*! This is a printf like method to write members of CSAs and CSTs into an outstream.
 * \tparam t_idx   Type of the index. Class should be of concept csa_tag or cst_tag.
 * \param out      Output stream.
 * \param format   Format string. See explanation below.
 * \param idx      CSA or CST object.
 * \param sentinel Character which should replace the 0-symbol in BWT/ TEXT.
 *
 * \par Format string
 *   Each line of the output will be formatted according to the format string.
 *   All content, except tokens which start with `%` will be copied. Tokens
 *   which start with `%` will be replaced as follows (let w be a positive
 *    number. setw(w) is used to format single numbers):
 *
 *      Token      |  Replacement | Comment
 *      -----------------------------------------------------------------------
 *       %[w]I     | Row index i.                           |
 *       %[w]S     | SA[i]                                  |
 *       %[w]s     | ISA[i]                                 |
 *       %[w]P     | PSI[i]                                 |
 *       %[w]p     | LF[i]                                  |
 *       %[w]L     | LCP[i]                                 | only for CSTs
 *       %[w]B     | BWT[i]                                 |
 *       %[w[:W]]T | Print min(idx.size(),w) chars of each  |
 *                 | suffix, each char formatted by setw(W).|
 *       %%        | %                                      |
 */
template <typename t_idx>
void csXprintf(std::ostream & out,
               std::string const & format,
               t_idx const & idx,
               char sentinel = default_sentinel<t_idx>::value)
{
    typename t_idx::index_category cat;
    const typename t_idx::csa_type & csa = _idx_csa(idx, cat);
    std::vector<std::string> res(csa.size());
    bool truncate = false;
    for (std::string::const_iterator c = format.begin(), s = c; c != format.end(); s = c)
    {
        while (c != format.end() and *c != '%')
            ++c; // string before the next `%`
        if (c > s)
        { // copy format string part
            std::vector<std::string> to_copy(csa.size(), std::string(s, c));
            transform(res.begin(), res.end(), to_copy.begin(), res.begin(), std::plus<std::string>());
        }
        if (c == format.end())
            break;
        ++c;                                         // skip `%`
        uint64_t w = _parse_number(c, format.end()); // element width
        if (c == format.end())
            break;
        uint64_t W = 0; // character width
        if (':' == *c)
        {
            ++c;
            W = _parse_number(c, format.end());
        }
        if (c == format.end())
            break;
        for (uint64_t i = 0; i < csa.size(); ++i)
        {
            switch (*c)
            {
            case 'I':
                res[i] += util::to_string(i, w);
                break;
            case 'S':
                res[i] += util::to_string(csa[i], w);
                break;
            case 's':
                res[i] += util::to_string(csa.isa[i], w);
                break;
            case 'P':
                res[i] += util::to_string(csa.psi[i], w);
                break;
            case 'p':
                res[i] += util::to_string(csa.lf[i], w);
                break;
            case 'L':
                res[i] += _idx_lcp_val(idx, i, w, cat);
                break;
            case 'B':
                if (0 == csa.bwt[i])
                {
                    res[i] += util::to_string(sentinel, w);
                }
                else
                {
                    res[i] += util::to_string(csa.bwt[i], w);
                }
                break;
            case 'U':
                truncate = true;
                SDSL_FALLTHROUGH
            case 'T':
                for (uint64_t k = 0; (w > 0 and k < w) or (0 == w and k < csa.size()); ++k)
                {
                    if (0 == csa.text[(csa[i] + k) % csa.size()])
                    {
                        res[i] += util::to_string(sentinel, W);
                        if (truncate)
                        {
                            truncate = false;
                            break;
                        }
                    }
                    else
                    {
                        res[i] += util::to_string(csa.text[(csa[i] + k) % csa.size()], W);
                    }
                }
                break;
            case 'u':
                truncate = true;
                SDSL_FALLTHROUGH
            case 't':
                for (uint64_t k = 0; (w > 0 and k < w) or (0 == w and k < csa.size()); ++k)
                {
                    if (0 == csa.text[(i + k) % csa.size()])
                    {
                        res[i] += util::to_string(sentinel, W);
                        if (truncate)
                        {
                            truncate = false;
                            break;
                        }
                    }
                    else
                    {
                        res[i] += util::to_string(csa.text[(i + k) % csa.size()], W);
                    }
                }
                break;
            case '%':
                res[i] += "%";
                break;
            }
        }
        ++c;
    }
    for (size_t i = 0; i < res.size(); ++i)
        out << res[i] << std::endl;
}

//! Returns the file name of the resource.
/*!
 * \param  key        Resource key.
 * \param  config    Cache configuration.
 * \return The file name of the resource.
 */
inline std::string cache_file_name(std::string const & key, cache_config const & config)
{
    if (config.file_map.count(key) != 0)
    {
        return config.file_map.at(key);
    }
    return config.dir + "/" + key + "_" + config.id + ".sdsl";
}

//! Returns the file name of the resource.
/*!
 * \param  key        Resource key.
 * \param  config    Cache configuration.
 * \return The file name of the resource.
 */
template <typename T>
std::string cache_file_name(std::string const & key, cache_config const & config)
{
    return cache_file_name(key + "_" + util::class_to_hash(T()), config);
}

//! Register the existing resource specified by the key to the cache
/*!
 *  \param key        Resource key.
 *  \param config    Cache configuration.
 *
 *  Note: If the resource does not exist under the given key,
 *  it will be not added to the cache configuration.
 */
inline void register_cache_file(std::string const & key, cache_config & config)
{
    std::string file_name = cache_file_name(key, config);
    isfstream in(file_name);
    if (in)
    { // if file exists, register it.
        config.file_map[key] = file_name;
    }
}

//! Checks if the resource specified by the key exists in the cache.
/*!
 * \param key    Resource key.
 * \param config Cache configuration.
 * \return True, if the file exists, false otherwise.
 */
inline bool cache_file_exists(std::string const & key, cache_config const & config)
{
    std::string file_name = cache_file_name(key, config);
    isfstream in(file_name);
    if (in)
    {
        in.close();
        return true;
    }
    return false;
}

//! Checks if the resource specified by the key and type exists in the cache.
/*!
 * \tparam T     Type.
 * \param key    Resource key.
 * \param config Cache configuration.
 * \return True, if the file exists, false otherwise.
 */
template <typename T>
bool cache_file_exists(std::string const & key, cache_config const & config)
{
    return cache_file_exists(key + "_" + util::class_to_hash(T()), config);
}

//! Returns a name for a temporary file. I.e. the name was not used before.
inline std::string tmp_file(cache_config const & config, std::string name_part = "")
{
    return config.dir + "/" + util::to_string(util::pid()) + "_" + util::to_string(util::id()) + name_part + ".sdsl";
}

//! Returns a name for a temporary file. I.e. the name was not used before.
inline std::string tmp_file(std::string const & filename, std::string name_part = "")
{
    return util::dirname(filename) + "/" + util::to_string(util::pid()) + "_" + util::to_string(util::id()) + name_part
         + ".sdsl";
}

template <typename T>
bool load_from_cache(T & v, std::string const & key, cache_config const & config, bool add_type_hash = false)
{
    std::string file;
    if (add_type_hash)
    {
        file = cache_file_name<T>(key, config);
    }
    else
    {
        file = cache_file_name(key, config);
    }
    if (load_from_file(v, file))
    {
        if (util::verbose)
        {
            std::cerr << "Load `" << file << std::endl;
        }
        return true;
    }
    else
    {
        std::cerr << "WARNING: Could not load file '";
        std::cerr << file << "'" << std::endl;
        return false;
    }
}

//! Stores the object v as a resource in the cache.
/*!
 *  \param
 */
template <typename T>
bool store_to_cache(T const & v, std::string const & key, cache_config & config, bool add_type_hash = false)
{
    std::string file;
    if (add_type_hash)
    {
        file = cache_file_name<T>(key, config);
    }
    else
    {
        file = cache_file_name(key, config);
    }
    if (store_to_file(v, file))
    {
        config.file_map[std::string(key)] = file;
        return true;
    }
    else
    {
        std::cerr << "WARNING: store_to_cache: could not store file `" << file << "`" << std::endl;
        return false;
    }
}

template <typename T>
bool remove_from_cache(std::string const & key, cache_config & config, bool add_type_hash = false)
{
    std::string file;
    if (add_type_hash)
    {
        file = cache_file_name<T>(key, config);
    }
    else
    {
        file = cache_file_name(key, config);
    }
    config.file_map.erase(key);
    if (sdsl::remove(file) == 0)
    {
        return true;
    }
    else
    {
        std::cerr << "WARNING: delete_from_cache: could not delete file `" << file << "`" << std::endl;
        return false;
    }
}

//==================== Template functions ====================

template <typename T>
typename T::size_type size_in_bytes(T const & t)
{
    nullstream ns;
    return serialize(t, ns);
}

template <typename T>
double size_in_mega_bytes(T const & t)
{
    return size_in_bytes(t) / (1024.0 * 1024.0);
}

template <typename T>
void add_hash(T const & t, std::ostream & out)
{
    uint64_t hash_value = util::hashvalue_of_classname(t);
    write_member(hash_value, out);
}

template <typename T>
bool store_to_file(T const & t, std::string const & file)
{
    osfstream out(file, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out)
    {
        if (util::verbose)
        {
            std::cerr << "ERROR: store_to_file not successful for: `" << file << "`" << std::endl;
        }
        return false;
    }
    serialize(t, out);
    out.close();
    if (util::verbose)
    {
        std::cerr << "INFO: store_to_file: `" << file << "`" << std::endl;
    }
    return true;
}

template <typename T>
bool store_to_checked_file(T const & t, std::string const & file)
{
    std::string checkfile = file + "_check";
    osfstream out(checkfile, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out)
    {
        if (util::verbose)
        {
            std::cerr << "ERROR: store_to_checked_file not successful for: `" << checkfile << "`" << std::endl;
        }
        return false;
    }
    add_hash(t, out);
    out.close();
    return store_to_file(t, file);
}

inline bool store_to_checked_file(char const * v, std::string const & file)
{
    std::string checkfile = file + "_check";
    osfstream out(checkfile, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out)
    {
        if (util::verbose)
        {
            std::cerr << "ERROR: store_to_checked_file(const char *v, const std::string&)" << std::endl;
            return false;
        }
    }
    add_hash(v, out);
    out.close();
    return store_to_file(v, file);
}

inline bool store_to_file(std::string const & v, std::string const & file)
{
    osfstream out(file, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out)
    {
        if (util::verbose)
        {
            std::cerr << "ERROR: store_to_file(const std::string& v, const std::string&)" << std::endl;
            return false;
        }
    }
    out.write(v.data(), v.size());
    out.close();
    return true;
}

template <uint8_t t_width>
bool store_to_file(int_vector<t_width> const & v, std::string const & file)
{
    osfstream out(file, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out)
    {
        std::cerr << "ERROR: util::store_to_file:: Could not open file `" << file << "`" << std::endl;
        return false;
    }
    else
    {
        if (util::verbose)
        {
            std::cerr << "INFO: store_to_file: `" << file << "`" << std::endl;
        }
    }
    v.serialize(out, nullptr, "");
    out.close();
    return true;
}

template <uint8_t t_width>
bool store_to_checked_file(int_vector<t_width> const & v, std::string const & file)
{
    std::string checkfile = file + "_check";
    osfstream out(checkfile, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out)
    {
        std::cerr << "ERROR: util::store_to_checked_file: Could not open check file `" << checkfile << "`" << std::endl;
        return false;
    }
    else
    {
        if (util::verbose)
        {
            std::cerr << "INFO: store_to_checked_file: `" << checkfile << "`" << std::endl;
        }
    }
    add_hash(v, out);
    out.close();
    return store_to_file(v, file);
}

template <typename T>
bool load_from_file(T & v, std::string const & file)
{
    isfstream in(file, std::ios::binary | std::ios::in);
    if (!in)
    {
        if (util::verbose)
        {
            std::cerr << "Could not load file `" << file << "`" << std::endl;
        }
        return false;
    }
    load(v, in);
    in.close();
    if (util::verbose)
    {
        std::cerr << "Load file `" << file << "`" << std::endl;
    }
    return true;
}

template <typename T>
bool load_from_checked_file(T & v, std::string const & file)
{
    isfstream in(file + "_check", std::ios::binary | std::ios::in);
    if (!in)
    {
        if (util::verbose)
        {
            std::cerr << "Could not load check file `" << file << "_check`" << std::endl;
        }
        return false;
    }
    uint64_t hash_value;
    read_member(hash_value, in);
    if (hash_value != util::hashvalue_of_classname(v))
    {
        if (util::verbose)
        {
            std::cerr << "File `" << file << "` is not an instance of the class `"
                      << sdsl::util::demangle2(typeid(T).name()) << "`" << std::endl;
        }
        return false;
    }
    return load_from_file(v, file);
}

template <typename t_iv>
inline typename std::enable_if<std::is_same<typename t_iv::index_category, iv_tag>::value
                                   or std::is_same<typename t_iv::index_category, csa_tag>::value
                                   or std::is_same<typename t_iv::index_category, lcp_tag>::value,
                               std::ostream &>::type
operator<<(std::ostream & os, t_iv const & v)
{
    for (auto it = v.begin(), end = v.end(); it != end; ++it)
    {
        os << *it;
        if (it + 1 != end)
            os << " ";
    }
    return os;
}

template <typename t_iv>
inline typename std::enable_if<std::is_same<typename t_iv::index_category, wt_tag>::value, std::ostream &>::type
operator<<(std::ostream & os, t_iv const & v)
{
    for (auto it = v.begin(), end = v.end(); it != end; ++it)
    {
        os << *it;
        if (it + 1 != end and std::is_same<typename t_iv::alphabet_category, int_alphabet_tag>::value)
            os << " ";
    }
    return os;
}

template <typename t_int>
inline typename std::enable_if<std::is_integral<t_int>::value, std::ostream &>::type
operator<<(std::ostream & os, std::vector<t_int> const & v)
{
    for (auto it = v.begin(), end = v.end(); it != end; ++it)
    {
        os << *it;
        if (it + 1 != end)
            os << " ";
    }
    return os;
}

template <typename t_iv>
inline typename std::enable_if<std::is_same<typename t_iv::category, csa_member_tag>::value, std::ostream &>::type
operator<<(std::ostream & os, t_iv const & v)
{
    for (auto it = v.begin(), end = v.end(); it != end; ++it)
    {
        os << *it;
        if (it + 1 != end and std::is_same<typename t_iv::alphabet_category, int_alphabet_tag>::value)
            os << " ";
    }
    return os;
}
} // namespace sdsl
#endif
