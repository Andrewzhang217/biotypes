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
    "7?CIGJB:D:-F7LA:GI9FDHBIJ7,GHGJBKHNI7IN,EML8IFIA7HN7J6,L6686LCJE"
    "JKA6G7AK6GK5C6@6IK+++?5+=<;227*6054"};

static const Sequence sequence0("sequence0", data_array[0]);
static const Sequence sequence1("sequence1", data_array[1]);
static const Sequence sequence2("sequence2", data_array[2]);
static const Sequence sequence3("sequence3", data_array[3]);
static const Sequence sequence4("sequence4", data_array[4]);
static const Sequence sequence_with_quality("sequenceQ", data_array[4],
                                            quality_array[0]);
static const Sequence sequence_with_quality1("sequenceQ1", data_array[5],
                                             quality_array[1]);
static const Sequence* ptr0{&sequence0};
static const Sequence* ptr1(&sequence1);
static const Sequence* ptr2(&sequence2);
static const Sequence* ptr3(&sequence3);
static const Sequence* ptr4(&sequence4);

TEST_CASE("biotypes/bioseq::Sequence") {
  SECTION("API: atValue") {
    REQUIRE(sequence0.atValue(1) == 'A');
    Sequence sequence4("sequence4", data_array[4]);
    REQUIRE(sequence4.atValue(5) == 'G');
    REQUIRE(sequence4.atValue(15) == 'A');
    REQUIRE(sequence4.atValue(31) == 'T');
    Sequence sequence6("sequence6", data_array[6]);
    REQUIRE(sequence6.atValue(32) == 'A');
  }
  SECTION("API: atQuality") {
    REQUIRE(sequence_with_quality.atQuality(0) == 'A');
    REQUIRE(sequence_with_quality.atQuality(1) == 'A');
    REQUIRE(sequence_with_quality.atQuality(31) == 'A');
    REQUIRE(sequence_with_quality1.atQuality(32) == 'B');
    REQUIRE(sequence_with_quality1.atQuality(33) == 'B');
    REQUIRE(sequence_with_quality1.atQuality(63) == 'B');
  }
  SECTION("API: atBase") {
    REQUIRE(sequence_with_quality.atBase(5) == Base{'G', 'A'});
    REQUIRE(sequence3.atBase(19) == Base{'T', 127});
  }
  SECTION("API: begin") { REQUIRE(sequence4.begin()->value == 'A'); }
  SECTION("API: end") {
    REQUIRE((*sequence4.end()).value == Sequence::kSentinel.value);
    REQUIRE((*sequence4.end()).phred33 == Sequence::kSentinel.phred33);
  }
  SECTION("Sequence copy") {
    SECTION("copy construction") {
      auto const seq0_cpy(sequence0);
      REQUIRE(seq0_cpy.size() == sequence0.size());
      for (auto it = 0; it < seq0_cpy.size(); ++it) {
        REQUIRE(seq0_cpy.atValue(it) == sequence0.atValue(it));
      }
    }

    SECTION("copy assignment") {
      auto const seq0_cpy = sequence0;
      REQUIRE(seq0_cpy.size() == sequence0.size());
      for (auto it = 0; it < seq0_cpy.size(); ++it) {
        REQUIRE(seq0_cpy.atValue(it) == sequence0.atValue(it));
      }
    }

    SECTION("move construction") {
      auto const seq1_cpy(std::move(sequence1));

      CHECK(sequence1.name() == "sequence1");
      CHECK(sequence1.size() == 32);
      for (auto i = 0; i < 32; ++i) {
        CHECK(sequence1.atValue(i) == data_array[1][i]);
      }

      CHECK(seq1_cpy.name() == sequence1.name());
      CHECK(seq1_cpy.size() == sequence1.size());
      for (auto i = 0; i < 32; ++i) {
        CHECK(seq1_cpy.atValue(i) == data_array[1][i]);
      }
    }

    SECTION("move assignment") {
      auto const seq4_cpy = std::move(sequence4);

      CHECK(sequence4.name() == "sequence4");
      CHECK(sequence4.size() == 32);
      for (auto i = 0; i < 32; ++i) {
        CHECK(sequence4.atValue(i) == data_array[4][i]);
      }

      CHECK(seq4_cpy.name() == sequence4.name());
      CHECK(seq4_cpy.size() == sequence4.size());
      for (auto i = 0; i < 32; ++i) {
        CHECK(seq4_cpy.atValue(i) == data_array[4][i]);
      }
    }
  }
}
TEST_CASE("biotypes/bioseq::iterators") {
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
    iterator--;
    REQUIRE(iterator->value == 'G');
    auto it0 = sequence0.begin();
    it0 = sequence1.begin();
    REQUIRE((*it0).value == 'C');
    auto it1 = sequence1.begin();
    auto it4{it1};
    REQUIRE((*it4).value == 'C');
    auto it2 = sequence2.begin();
    auto it3 = std::move(it2);
    REQUIRE((*it3).value == 'G');
    auto temp{std::move(it3)};
    REQUIRE((*temp).value == 'G');
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
