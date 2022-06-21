#include "lbcb/biotypes//bioseq.h"

#include <algorithm>
#include <iterator>
#include <string_view>
#include <vector>

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
  SECTION("API: AtValue") {
    REQUIRE(sequence0.AtValue(1) == 'A');
    Sequence sequence4("sequence4", data_array[4]);
    REQUIRE(sequence4.AtValue(5) == 'G');
    REQUIRE(sequence4.AtValue(15) == 'A');
    REQUIRE(sequence4.AtValue(31) == 'T');
    Sequence sequence6("sequence6", data_array[6]);
    REQUIRE(sequence6.AtValue(32) == 'A');
  }
  SECTION("API: AtQuality") {}
  SECTION("API: AtBase") {}
  SECTION("API: Begin") {
    REQUIRE(sequence0.Begin() == Sequence::Iterator{sequence0, 0});
  }
  SECTION("API: End") {
    REQUIRE(sequence0.End() ==
            Sequence::Iterator{sequence0, data_array[0].size()});
  }
  SECTION("Iterator") {
    auto iterator = sequence4.Begin();
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
    for (auto it = sequence4.Begin(); it != sequence4.End(); ++it) {
      REQUIRE((*it).value == data_array[4][counter]);
      REQUIRE(it->value == data_array[4][counter++]);
    }
    std::vector<Base> from_vector;
    counter = 0;
    for (auto it = sequence4.Begin(); it != sequence4.End(); ++it) {
      REQUIRE(it->value == data_array[4][counter++]);
      from_vector.emplace_back(*it);
    }
    std::vector<Base> to_vector;
    std::copy(from_vector.begin(), from_vector.end(),
              std::back_inserter(to_vector));
    counter = 0;
    for (auto& it : to_vector) {
      REQUIRE(it.value == data_array[4][counter++]);
    }
  }
}
};  // namespace lbcb::testing
