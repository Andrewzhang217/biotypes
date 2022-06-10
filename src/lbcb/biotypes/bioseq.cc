#include "lbcb/biotypes/bioseq.h"

#include <vector>

namespace lbcb {

Sequence::Sequence(std::string_view name, std::string_view data)
    : name_((name)), data_(data) {
  // TODO: magic
}
Sequence::Sequence(std::string_view name, std::string_view data,
                   std::string_view quality)
    : name_(name), data_(data), quality_(quality) {}

std::vector<std::uint64_t> Sequence::Compress(std::string& data) {
  std::uint64_t block;
  std::vector<std::uint64_t> compressed_data;
  int counter = 0;
  for (char& it : data) {
    counter++;
    if (counter % 32 == 1) {
      compressed_data.emplace_back(block);
      block = 0;
    }
    switch (it) {
      case 'A':
        block += 0;
        block <<= 2;
      case 'C':
        block += 1;
        block <<= 2;
      case 'G':
        block += 2;
        block <<= 2;
      case 'T':
        block += 3;
        block <<= 2;
      default:
        assert(false);
    }
  }
  return compressed_data;
}

}  // namespace lbcb
