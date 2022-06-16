#ifndef LBCB_BIOTYPES_BIOSEQ_H_
#define LBCB_BIOTYPES_BIOSEQ_H_

#include <array>
#include <string>
#include <string_view>
#include <vector>

namespace lbcb {

struct Base {
  char value;
  char phred33;
};

template <typename T>
class Iterator {
 public:
  Iterator() : pos_(0) {}
  explicit Iterator(std::size_t pos) : pos_(pos) {}

  std::size_t GetPos() { return pos_; };
  Iterator& operator+=(std::size_t diff) {
    pos_ += diff;
    return *this;
  }
  Iterator& operator-=(std::size_t diff) {
    pos_ -= diff;
    return *this;
  }
  Iterator& operator++() {
    pos_++;
    return *this;
  }
  Iterator& operator--() {
    pos_--;
    return *this;
  }
  Iterator operator+(std::size_t diff) const { return Iterator(pos_ + diff); }
  Iterator operator-(std::size_t diff) const { return Iterator(pos_ - diff); }

  bool operator==(const Iterator& rhs) const { return pos_ == rhs.pos_; }
  bool operator!=(const Iterator& rhs) const { return pos_ != rhs.pos_; }
  bool operator>(const Iterator& rhs) const { return pos_ > rhs.pos_; }
  bool operator<(const Iterator& rhs) const { return pos_ <= rhs.pos_; }
  bool operator>=(const Iterator& rhs) const { return pos_ >= rhs.pos_; }
  bool operator<=(const Iterator& rhs) const { return pos_ <= rhs.pos_; }

 private:
  std::size_t pos_;
};

class Sequence {
 public:
  Sequence(std::string_view name, std::string_view data);
  Sequence(std::string_view name, std::string_view data,
           std::string_view quality);

  [[nodiscard]] Base AtBase(std::size_t pos);
  [[nodiscard]] char AtValue(std::size_t pos);
  [[nodiscard]] char AtQuality(std::size_t pos);
  [[nodiscard]] Iterator<Base> Begin();
  [[nodiscard]] Iterator<Base> End();

 private:
  std::string name_;
  std::vector<std::uint64_t> compressed_data_;
  std::vector<std::uint64_t> compressed_quality_;
  std::size_t size_;
  Iterator<Base> iterator_;
};

}  // namespace lbcb

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
inline char constexpr kNucleotideDecoder[] = {'A', 'C', 'G', 'T'};

std::vector<std::uint64_t> Compress(std::string_view data);
std::string Decompress(const std::vector<std::uint64_t>& compressed_data);

}  // namespace lbcb::detail

#endif /* LBCB_BIOTYPES_BIOSEQ_H_  */
