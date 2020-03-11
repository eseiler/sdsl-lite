#include "sdsl/wavelet_trees.hpp"
#include "common.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <algorithm> // for std::min
#include <random>

namespace
{

using namespace sdsl;
using namespace std;

typedef int_vector<>::size_type size_type;

string temp_file;
string temp_dir;
int_vector<8> text;

template<class T>
class wt_byte_epr_test : public ::testing::Test { };

using testing::Types;

typedef Types<
    wt_epr<4>,
    wt_epr<5>,
    wt_epr<6>,
    wt_epr<7>,
    wt_epr<12>,
    wt_epr<13>,
    wt_epr<17>
> Implementations;

TYPED_TEST_CASE(wt_byte_epr_test, Implementations);

TYPED_TEST(wt_byte_epr_test, create_and_store)
{
    static_assert(sdsl::util::is_regular<TypeParam>::value, "Type is not regular");
    std::mt19937_64 gen(12345); // random but fixed text
    std::uniform_int_distribution<> dis(0, TypeParam::fixed_alphabet_size-1); // no 0 allowed

    text.resize(500000);

    for (auto it = text.begin(); it != text.end(); ++it)
        *it = dis(gen);
    //
    // for (auto && x : text)
    //     std::cout <<  (int)x << ',';
    // std::cout << '\n';

    TypeParam wt(text.begin(), text.end());

    ASSERT_TRUE(store_to_file(wt, temp_file));
}

//! Test sigma
TYPED_TEST(wt_byte_epr_test, sigma)
{
    TypeParam wt;
    ASSERT_TRUE(load_from_file(wt, temp_file));
    ASSERT_EQ(text.size(), wt.size());
    bit_vector occur(256, 0);
    uint16_t sigma = 0;
    for (size_type j=0; j<text.size(); ++j) {
        if (!occur[(unsigned char)text[j]]) {
            occur[(unsigned char)text[j]] = 1;
            ++sigma;
        }
    }
    ASSERT_EQ(sigma, wt.sigma);
}

template<class t_wt>
void compare_wt(const int_vector<8>& text, const t_wt& wt)
{
    ASSERT_EQ(text.size(), wt.size());
    for (size_type j=0; j<text.size(); ++j) {
        ASSERT_EQ((typename t_wt::value_type)text[j], wt[j])<<" j="<<j;
    }
}

//! Test Access method, Copy-construtor, Move-constructor, Copy-assign and Move-assign
TYPED_TEST(wt_byte_epr_test, access_copy_move_and_swap)
{
    TypeParam wt;
    ASSERT_TRUE(load_from_file(wt, temp_file));
    compare_wt(text, wt);

    // Copy-constructor
    TypeParam wt2(wt);
    compare_wt(text, wt2);

    // Move-constructor
    TypeParam wt3(std::move(wt2));
    compare_wt(text, wt3);

    // Copy-Assign
    TypeParam wt4;
    wt4 = wt3;
    compare_wt(text, wt4);

    // Move-Assign
    TypeParam wt5;
    wt5 = std::move(wt4);
    compare_wt(text, wt5);
}

//! Test rank methods
TYPED_TEST(wt_byte_epr_test, rank)
{
    TypeParam wt;
    ASSERT_TRUE(load_from_file(wt, temp_file));

    ASSERT_EQ(text.size(), wt.size());

    // for (auto && x : text)
    //     std::cout <<  (int)x << ',';
    // std::cout << '\n';

    // Test rank(i, c) for each character c and position i
    std::vector<size_t> cnt_prefix_rank(wt.sigma, 0);
    // 234 332 143 211 345
    // auto x = wt.rank(15, 4);
    // (void) x;
    for (size_t i = 0; i < text.size() + 1; ++i)
    {
        for (size_t c = 0; c < wt.sigma; ++c)
        {
            if (i > 0 && text[i - 1] <= c) // 12 - 12
                ++cnt_prefix_rank[c];

            // auto const rank = rb(i, c);
            // std::cout << "i=" << i << " c=" << c << '\n';
            if (c > 0)
                ASSERT_EQ(cnt_prefix_rank[c] - cnt_prefix_rank[c - 1], wt.rank(i, c)) << "i=" << i << " c=" << c << "/" << wt.sigma - 1;
            else
                ASSERT_EQ(cnt_prefix_rank[c], wt.rank(i, c)) << "i=" << i << " c=" << c << "/" << wt.sigma - 1;

            if (c > 0)
            {
                auto lex = wt.lex_smaller_count(i, c);
                ASSERT_EQ(cnt_prefix_rank[c - 1], std::get<1>(lex)) << "i=" << i << " c=" << c << "/" << wt.sigma - 1;
            }
        }
    }
}

#if SDSL_HAS_CEREAL
template <typename in_archive_t, typename out_archive_t, typename TypeParam>
void do_serialisation(TypeParam const & l)
{
	{
		std::ofstream os{temp_file, std::ios::binary};
		out_archive_t oarchive{os};
		oarchive(l);
	}

	{
		TypeParam in_l{};
		std::ifstream is{temp_file, std::ios::binary};
		in_archive_t iarchive{is};
		iarchive(in_l);
		EXPECT_EQ(l, in_l);
	}
}

TYPED_TEST(wt_byte_epr_test, cereal)
{
	if (temp_dir != "@/")
	{
		TypeParam wt;
	        ASSERT_TRUE(load_from_file(wt, temp_file));

		do_serialisation<cereal::BinaryInputArchive,         cereal::BinaryOutputArchive>        (wt);
		do_serialisation<cereal::PortableBinaryInputArchive, cereal::PortableBinaryOutputArchive>(wt);
		do_serialisation<cereal::JSONInputArchive,           cereal::JSONOutputArchive>          (wt);
		do_serialisation<cereal::XMLInputArchive,            cereal::XMLOutputArchive>           (wt);
	}
}
#endif // SDSL_HAS_CEREAL

TYPED_TEST(wt_byte_epr_test, delete_)
{
    sdsl::remove(temp_file);
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 2) {
        // LCOV_EXCL_START
        std::cout << "Usage: " << argv[0] << " tmp_dir" << std::endl;
        return 1;
        // LCOV_EXCL_STOP
    }
    temp_dir = argv[1];
    temp_file = temp_dir + "/wt_epr";

    return RUN_ALL_TESTS();
}
