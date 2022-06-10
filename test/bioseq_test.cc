#include "lbcb/biotypes//bioseq.h"

#include <iostream>
#include <string>
#include <string_view>

#include "../include/lbcb/biotypes/bioseq.h"
#include "catch2/catch_test_macros.hpp"

namespace lbcb {}

namespace lbcb::detail {
TEST_CASE("biotypes/bioseq") {
  SECTION("compress (basic data_0)") {
    std::string_view basic_data_0{"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"};
    std::uint64_t expected_data = 0b00000000000000000000000000000000;
    auto res = Compress(basic_data_0);
    std::cout << basic_data_0 << std::endl;
    std::cout << expected_data << std::endl;
    std::cout << res[0] << std::endl;
    REQUIRE(res[0] == expected_data);
  }

  SECTION("compress (basic data_1)") {
    std::string_view basic_data_1{"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"};
    std::uint64_t expected_data =
        0b0101010101010101010101010101010101010101010101010101010101010101;
    auto res1 = Compress(basic_data_1);
    std::cout << basic_data_1 << std::endl;
    std::cout << expected_data << std::endl;
    std::cout << res1[0] << std::endl;
    REQUIRE(res1[0] == expected_data);
  }
  SECTION("compress (basic data_2)") {
    std::string_view basic_data_2{"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"};
    std::uint64_t expected_data =
        0b1010101010101010101010101010101010101010101010101010101010101010;
    auto res2 = Compress(basic_data_2);
    std::cout << basic_data_2 << std::endl;
    std::cout << expected_data << std::endl;
    std::cout << res2[0] << std::endl;
    REQUIRE(res2[0] == expected_data);
  }
  SECTION("compress (basic data_3)") {
    std::string_view basic_data_3{"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"};
    std::uint64_t expected_data =
        0b1111111111111111111111111111111111111111111111111111111111111111;
    auto res3 = Compress(basic_data_3);
    std::cout << basic_data_3 << std::endl;
    std::cout << expected_data << std::endl;
    std::cout << res3[0] << std::endl;
    REQUIRE(res3[0] == expected_data);
  }
  SECTION("compress (basic data_4)") {
    std::string_view basic_data_4{"ATCTTGTATTCCTTAATTTTTTTTTTTTACGT"};
    std::uint64_t expected_data =
        0b0011011111101100111101011111000011111111111111111111111100011011;
    auto res4 = Compress(basic_data_4);
    std::cout << basic_data_4 << std::endl;
    std::cout << expected_data << std::endl;
    std::cout << res4[0] << std::endl;
    REQUIRE(res4[0] == expected_data);
  }
  SECTION("compress (basic data_5)") {
    std::string_view basic_data_5{"ATCTTGTATTCCTTAATTTTTTTTTTTTACGT"};
    std::uint64_t expected_data =
        0b0011011111101100111101011111000011111111111111111111111100011011;
    auto res5 = Compress(basic_data_5);
    std::cout << basic_data_5 << std::endl;
    std::cout << expected_data << std::endl;
    std::cout << res5[0] << std::endl;
    REQUIRE(res5[0] == expected_data);
  }
  SECTION("compress (basic data_6)") {
    std::string_view basic_data_6{
        "ATCTTGTATTCCTTAATTTTTTTTTTTTACGTATCTTGTATTCCTTAATTTTTTTTTTTTACGT"};
    std::uint64_t expected_data =
        0b0011011111101100111101011111000011111111111111111111111100011011;
    auto res6 = Compress(basic_data_6);
    std::cout << basic_data_6 << std::endl;
    std::cout << expected_data << std::endl;
    std::cout << res6[0] << std::endl;
    REQUIRE(res6[0] == expected_data);
    REQUIRE(res6[1] == expected_data);
  }
}
}  // namespace lbcb::detail
