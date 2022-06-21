#include "lbcb/biotypes/bioseq.h"

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace lbcb {
Sequence::Sequence(std::string_view name, std::string_view data)
    : name_((name)),
      compressed_data_(detail::Compress(data)),
      size_(data.size()) {}
Sequence::Sequence(std::string_view name, std::string_view data,
                   std::string_view quality)
    : name_(name),
      compressed_data_(detail::Compress(data)),
      compressed_quality_(detail::Compress(quality)),
      size_(data.size()) {}
Base Sequence::AtBase(std::size_t pos) const {
  assert(pos < size_);
  return {AtValue(pos), AtQuality(pos)};
}
char Sequence::AtValue(std::size_t pos) const {
  assert(pos < size_);
  std::uint64_t block = compressed_data_[pos >> 5];
  block <<= 2 * (pos & 31);
  return detail::kNucleotideDecoder[block >> 62];
}
char Sequence::AtQuality(std::size_t pos) const { return {}; }
Sequence::Iterator Sequence::Begin() { return {*this, 0}; }
Sequence::Iterator Sequence::End() { return {*this, size_}; }

}  // namespace lbcb

namespace lbcb::detail {
std::vector<std::uint64_t> Compress(std::string_view src) {
  std::vector<std::uint64_t> dst;
  dst.resize(std::ceil(static_cast<double>(src.size()) / 32));
  std::uint64_t active_block = 0U;
  auto index = 0U;

  for (auto i = 0U; i < src.size(); ++i) {
    if (i > 0 && (i & 31) == 0) {
      dst[index] = std::exchange(active_block, 0);
      index++;
    }
    active_block = (active_block << 2U) | kBaseEncodings[src[i]];
  }

  dst[index] = active_block;
  return dst;
}
std::string Decompress(const std::vector<std::uint64_t>& dst) {
  std::string src;
  for (unsigned long long it : dst) {
    std::string curr;
    int index = 0;
    while (index < 32) {
      std::uint64_t temp = (it & (3ULL << 62)) >> 62;
      curr += kNucleotideDecoder[temp];
      index++;
      it <<= 2;
    }
    src += curr;
    curr = "";
  }
  return src;
}

}  // namespace lbcb::detail
