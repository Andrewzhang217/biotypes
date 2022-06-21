#include "lbcb/biotypes//bioseq.h"

#include <algorithm>
#include <iterator>
#include <string_view>

#include "catch2/catch_test_macros.hpp"

namespace lbcb::testing {
static constexpr std::array<std::string_view, 32> data_array{
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG",
    "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
    "ATCTTGTATTCCTTAATTTTTTTTTTTTACGT",
    "ATCTTGTATTCCTTAATTTTTTTTTTTTACGTATCTTGTATTCCTTAATTTTTTTTTTTTACGT",
    "ATCTTGTATTCCTTAATTTTTTTTTTTTACGTATCTTGTATTCCTTAATTTTTTTTTTTTACGTGG",
    "ATCTTGTATTCCTTAATTTTTTTTTTTTACGTACG",
    "ATCTTGTATTCCTTAATTTTTTTTTTTTACGT"};

static constexpr std::array<std::string_view, 32> quality_array{
    "7?CIGJB:D:-F7LA:GI9FDHBIJ7,GHGJB",
    "7?CIGJB:D:-F7LA:GI9FDHBIJ7,GHGJBKHNI7IN,EML8IFIA7HN7J6,L6686LCJE?"
    "JKA6G7AK6GK5C6@6IK+++?5+=<;227*6054"};

static constexpr std::array<std::uint64_t, 32> bits_array{
    0b0000000000000000000000000000000000000000000000000000000000000000,
    0b0101010101010101010101010101010101010101010101010101010101010101,
    0b1010101010101010101010101010101010101010101010101010101010101010,
    0b1111111111111111111111111111111111111111111111111111111111111111,
    0b0011011111101100111101011111000011111111111111111111111100011011,
};
static constexpr std::uint64_t trailing(
    0b0000000000000000000000000000000000000000000000000000000000001010);

static Sequence sequence0("sequence0", data_array[0]);
static Sequence sequence4("sequence4", data_array[4]);

TEST_CASE("biotypes/bioseq") {
  SECTION("Compress (basic data)") {
    auto res0 = detail::Compress(data_array[0]);
    REQUIRE(res0[0] == bits_array[0]);
    auto res1 = detail::Compress(data_array[1]);
    REQUIRE(res1[0] == bits_array[1]);
    auto res2 = detail::Compress(data_array[2]);
    REQUIRE(res2[0] == bits_array[2]);
    auto res3 = detail::Compress(data_array[3]);
    REQUIRE(res3[0] == bits_array[3]);
    auto res4 = detail::Compress(data_array[4]);
    REQUIRE(res4[0] == bits_array[4]);
    auto res5 = detail::Compress(data_array[5]);
    REQUIRE(res5[0] == bits_array[4]);
    REQUIRE(res5[1] == bits_array[4]);
    auto res6 = detail::Compress(data_array[6]);
    REQUIRE(res6[0] == bits_array[4]);
    REQUIRE(res6[1] == bits_array[4]);
    REQUIRE(res6[2] == trailing);
  }
  SECTION("API: atValue") {
    REQUIRE(sequence0.atValue(1) == 'A');
    Sequence sequence4("sequence4", data_array[4]);
    REQUIRE(sequence4.atValue(5) == 'G');
    REQUIRE(sequence4.atValue(15) == 'A');
    REQUIRE(sequence4.atValue(31) == 'T');
    Sequence sequence6("sequence6", data_array[6]);
    REQUIRE(sequence6.atValue(32) == 'A');
  }
  SECTION("API: atQuality") {}
  SECTION("API: atBase") {}
  SECTION("API: begin") {
    REQUIRE(sequence0.begin() == Sequence::iterator{sequence0, 0});
  }
  SECTION("API: end") {
    REQUIRE(sequence0.end() ==
            Sequence::iterator{sequence0, data_array[0].size()});
  }
  SECTION("iterator") {
    auto iterator = sequence4.begin();
    iterator += 2;
    REQUIRE((*iterator).value == 'C');
    iterator -= 1;
    REQUIRE((*iterator).value == 'T');
    iterator = iterator + 6;
    REQUIRE((*iterator).value == 'A');
    iterator = iterator - 2;
    REQUIRE((*iterator).value == 'G');
    iterator++;
    REQUIRE(iterator->value == 'T');
    int counter = 0;
    for (auto it = sequence4.begin(); it != sequence4.end(); ++it) {
      REQUIRE((*it).value == data_array[4][counter]);
      REQUIRE(it->value == data_array[4][counter++]);
    }
    auto it = sequence4.begin();
    REQUIRE(it[0].value == 'A');
    REQUIRE(it[1].value == 'T');
    REQUIRE(it[2].value == 'C');

    auto buff = std::string(sequence4.size(), '0');
    std::transform(sequence4.begin(), sequence4.end(), buff.begin(),
                   [](Base const base) -> char { return base.value; });
    REQUIRE(buff == data_array[4]);
  }
}
};  // namespace lbcb::testing
