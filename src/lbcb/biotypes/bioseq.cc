#include "lbcb/biotypes/bioseq.h"

#include <iostream>
#include <vector>

namespace lbcb {

Sequence::Sequence(std::string_view name, std::string_view data)
    : name_((name)), data_(data) {
  // TODO: magic
}
Sequence::Sequence(std::string_view name, std::string_view data,
                   std::string_view quality)
    : name_(name), data_(data), quality_(quality) {}
}  // namespace lbcb

namespace lbcb::detail {
std::vector<std::uint64_t> Compress(std::string_view data) {
  std::uint64_t block;
  std::vector<std::uint64_t> compressed_data;
  compressed_data.reserve(data.size() / 32 + 1);
  int counter = 0;
  for (char it : data) {
    counter++;
    block <<= 2;
    switch (it) {
      case 'A':
        block += 0b00;
        break;
      case 'C':
        block += 0b01;
        break;
      case 'G':
        block += 0b10;
        break;
      case 'T':
        block += 0b11;
        break;
      default:
        assert(false);
        std::cout << "Wrong input: not a nucleotide" << std::endl;
        return compressed_data;
    }

    if (counter == 32) {
      compressed_data.emplace_back(block);
      counter = 0;
      block = 0;
    }
  }
  // input size is not the multiple of 32
  //  block <<= (32 - counter);
  //  compressed_data.emplace_back(block);
  return compressed_data;
}
}  // namespace lbcb::detail
