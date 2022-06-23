#include "lbcb/biotypes/bioseq.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>

namespace lbcb::detail {

/* clang-format off */
static constexpr auto kBaseEncodings = std::array<std::uint8_t, 256> {
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255,   0, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255,   0,   1 ,  1,   0, 255, 255,   2,
    3, 255, 255,   2, 255,   1,   0, 255,
    255, 255,   0,   1,   3,   3,   2,   0,
    255,   3, 255, 255, 255, 255, 255, 255,
    255,   0,   1,   1,   0, 255, 255,   2,
    3, 255, 255,   2, 255,   1,   0, 255,
    255, 255,   0,   1,   3,   3,   2,   0,
    255,   3, 255, 255, 255, 255, 255, 255
};
/* clang-format on */
/* Mapping : - 0
             A 0 a 0
             B 1 b 1
             C 1 c 1
             D 0 d 0
             G 2 g 2
             H 3 h 3
             K 2 k 2
             M 1 m 1
             N 0 n 0
             R 0 r 0
             S 1 s 1
             T 3 t 3
             U 3 u 3
             V 2 v 2
             W 0 w 0
             Y 3 y 3
 */
static constexpr char kNucleotideDecoder[] = {'A', 'C', 'G', 'T'};

std::vector<std::uint64_t> CompressData(std::string_view src) {
  std::vector<std::uint64_t> dst;
  dst.resize(std::ceil(static_cast<double>(src.size()) / 32));
  std::uint64_t active_block = 0U;
  auto index = 0U;

  for (auto i = 0U; i < src.size(); ++i) {
    if (i > 0 && (i & 31) == 0) {
      dst[index++] = std::exchange(active_block, 0);
    }
    active_block = (active_block << 2U) | kBaseEncodings[src[i]];
  }

  dst[index] = active_block;
  return dst;
}
std::vector<char> CompressQuality(std::string_view src) {
  std::vector<char> dst;
  dst.resize(std::ceil(static_cast<double>(src.size()) / 32));
  int index = 0;
  int sum = 0;

  for (auto i = 0; i < src.size(); ++i) {
    if (i > 0 && (i & 31) == 0) {
      dst[index] = static_cast<char>(sum / 32);
      index++;
      sum = 0;
    }
    sum += src[i];
  }
  if ((src.size() & 31) == 0) {
    dst[index] = static_cast<char>(sum / 32);
  } else {
    dst[index] = static_cast<char>(sum / (src.size() & 31));
  }
  return dst;
}
std::string DecompressData(const std::vector<std::uint64_t>& dst) {
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

namespace lbcb {
bool Base::operator==(const Base& other) const {
  return value == other.value && phred33 == other.phred33;
}
Sequence::iterator::iterator(const lbcb::Sequence* sequence,
                             std::size_t start_index)
    : sequence_(sequence), pos_(start_index) {
  if (start_index < sequence->size_) {
    base_ = sequence->atBase(start_index);
  } else {
    base_ = Sequence::kSentinel;
  }
}
const Base& Sequence::iterator::operator*() noexcept { return base_; }
Sequence::iterator::pointer Sequence::iterator::operator->() const noexcept {
  return &base_;
}
Sequence::iterator::value_type Sequence::iterator::operator[](
    std::size_t pos) noexcept {
  return sequence_->atBase(pos);
}
Sequence::iterator& Sequence::iterator::operator+=(std::size_t diff) noexcept {
  pos_ += diff;
  if (pos_ >= 0 && pos_ < sequence_->size()) {
    base_ = sequence_->atBase(pos_);
  }
  return *this;
}
Sequence::iterator& Sequence::iterator::operator-=(std::size_t diff) noexcept {
  pos_ -= diff;
  if (pos_ >= 0 && pos_ < sequence_->size()) {
    base_ = sequence_->atBase(pos_);
  }
  return *this;
}
Sequence::iterator& Sequence::iterator::operator++() noexcept {
  pos_++;
  if (pos_ >= 0 && pos_ < sequence_->size()) {
    base_ = sequence_->atBase(pos_);
  }
  return *this;
}
Sequence::iterator& Sequence::iterator::operator--() noexcept {
  pos_--;
  if (pos_ >= 0 && pos_ < sequence_->size()) {
    base_ = sequence_->atBase(pos_);
  }
  return *this;
}
const Sequence::iterator Sequence::iterator::operator++(int) noexcept {
  auto tmp = *this;
  ++(*this);
  return tmp;
}
const Sequence::iterator Sequence::iterator::operator--(int) noexcept {
  auto tmp = *this;
  --(*this);
  return tmp;
}
Sequence::iterator& Sequence::iterator::operator+(
    std::size_t diff) const noexcept {
  auto it = new iterator(sequence_, pos_ + diff);
  return *it;
}
Sequence::iterator& Sequence::iterator::operator-(
    std::size_t diff) const noexcept {
  auto it = new iterator(sequence_, pos_ - diff);
  return *it;
}
bool Sequence::iterator::operator==(
    const Sequence::iterator& rhs) const noexcept {
  return sequence_ == rhs.sequence_ && pos_ == rhs.pos_;
}
bool Sequence::iterator::operator!=(
    const Sequence::iterator& rhs) const noexcept {
  return !(*this == rhs);
}
bool Sequence::iterator::operator>(
    const Sequence::iterator& rhs) const noexcept {
  return sequence_ == rhs.sequence_ && pos_ > rhs.pos_;
}
bool Sequence::iterator::operator<(
    const Sequence::iterator& rhs) const noexcept {
  return sequence_ == rhs.sequence_ && pos_ <= rhs.pos_;
}
bool Sequence::iterator::operator>=(
    const Sequence::iterator& rhs) const noexcept {
  return sequence_ == rhs.sequence_ && pos_ >= rhs.pos_;
}
bool Sequence::iterator::operator<=(
    const Sequence::iterator& rhs) const noexcept {
  return sequence_ == rhs.sequence_ && pos_ <= rhs.pos_;
}
Sequence::Sequence(std::string_view name, std::string_view data)
    : name_((name)),
      compressed_data_(detail::CompressData(data)),
      compressed_quality_(),
      size_(data.size()) {}
Sequence::Sequence(std::string_view name, std::string_view data,
                   std::string_view quality)
    : name_(name),
      compressed_data_(detail::CompressData(data)),
      compressed_quality_(detail::CompressQuality(quality)),
      size_(data.size()) {}
Base Sequence::atBase(std::size_t pos) const {
  assert(pos >= 0 && pos < size_);
  if (compressed_quality_.empty())
    return {atValue(pos), static_cast<char>(127)};
  return {atValue(pos), atQuality(pos)};
}
char Sequence::atValue(std::size_t pos) const {
  assert(pos >= 0 && pos < size_);
  std::uint64_t block = compressed_data_[pos >> 5];
  block <<= 2 * (pos & 31);
  return detail::kNucleotideDecoder[block >> 62];
}
char Sequence::atQuality(std::size_t pos) const {
  assert(pos >= 0 && pos < size_);
  assert(!compressed_quality_.empty());
  return compressed_quality_[pos / 32];
}
Sequence::iterator Sequence::begin() const { return {this, 0}; }
Sequence::iterator Sequence::end() const { return {this, size_}; }
std::string Sequence::name() const noexcept { return name_; }
std::size_t Sequence::size() const noexcept { return size_; }

}  // namespace lbcb
