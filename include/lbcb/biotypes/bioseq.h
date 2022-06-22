#ifndef LBCB_BIOTYPES_BIOSEQ_H_
#define LBCB_BIOTYPES_BIOSEQ_H_

#include <array>
#include <iostream>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace lbcb {

struct Base {
  char value;
  char phred33;
};

class Sequence {
 public:
  class iterator {
   public:
    using value_type = Base;
    using reference = const Base&;
    using pointer = const Base*;

    iterator(Sequence* sequence, std::size_t start_index)
        : sequence_(sequence), pos_(start_index) {
      if (start_index < sequence->size_) {
        base_ = sequence->atBase(start_index);
      } else {
        base_ = kSentinel;
      }
    }
    reference operator*() noexcept { return base_; }
    pointer operator->() const noexcept { return &base_; }
    value_type operator[](std::size_t pos) noexcept {
      return sequence_->atBase(pos);
    }
    iterator& operator=(const iterator& other) = default;
    iterator& operator+=(std::size_t diff) noexcept {
      pos_ += diff;
      if (pos_ >= 0 && pos_ < sequence_->size()) {
        base_ = sequence_->atBase(pos_);
      }
      return *this;
    }
    iterator& operator-=(std::size_t diff) noexcept {
      pos_ -= diff;
      if (pos_ >= 0 && pos_ < sequence_->size()) {
        base_ = sequence_->atBase(pos_);
      }
      return *this;
    }
    iterator& operator++() noexcept {
      pos_++;
      if (pos_ >= 0 && pos_ < sequence_->size()) {
        base_ = sequence_->atBase(pos_);
      }
      return *this;
    }
    iterator& operator--() noexcept {
      pos_--;
      if (pos_ >= 0 && pos_ < sequence_->size()) {
        base_ = sequence_->atBase(pos_);
      }
      return *this;
    }
    const iterator operator++(int) noexcept {
      auto tmp = *this;
      ++(*this);
      return tmp;
    }
    const iterator operator--(int) noexcept {
      auto tmp = *this;
      --(*this);
      return tmp;
    }
    iterator& operator+(std::size_t diff) const noexcept {
      auto it = new iterator(sequence_, pos_ + diff);
      return *it;
    }
    iterator& operator-(std::size_t diff) const noexcept {
      auto it = new iterator(sequence_, pos_ - diff);
      return *it;
    }
    bool operator==(const iterator& rhs) const noexcept {
      return std::addressof(*sequence_) == std::addressof(*rhs.sequence_) &&
             pos_ == rhs.pos_;
    }
    bool operator!=(const iterator& rhs) const noexcept {
      return !(*this == rhs);
    }
    bool operator>(const iterator& rhs) const noexcept {
      return std::addressof(*sequence_) == std::addressof(*rhs.sequence_) &&
             pos_ > rhs.pos_;
    }
    bool operator<(const iterator& rhs) const noexcept {
      return std::addressof(*sequence_) == std::addressof(*rhs.sequence_) &&
             pos_ <= rhs.pos_;
    }
    bool operator>=(const iterator& rhs) const noexcept {
      return std::addressof(*sequence_) == std::addressof(*rhs.sequence_) &&
             pos_ >= rhs.pos_;
    }
    bool operator<=(const iterator& rhs) const noexcept {
      return std::addressof(*sequence_) == std::addressof(*rhs.sequence_) &&
             pos_ <= rhs.pos_;
    }

   private:
    Base base_;
    Sequence* sequence_;
    std::size_t pos_;
  };
  Sequence(std::string_view name, std::string_view data);
  Sequence(std::string_view name, std::string_view data,
           std::string_view quality);

  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] Base atBase(std::size_t pos) const;
  [[nodiscard]] char atValue(std::size_t pos) const;
  [[nodiscard]] char atQuality(std::size_t pos) const;
  [[nodiscard]] Sequence::iterator begin();
  [[nodiscard]] Sequence::iterator end();

 private:
  std::string name_;
  std::vector<std::uint64_t> compressed_data_;
  std::vector<std::uint64_t> compressed_quality_;
  std::size_t size_;
  static constexpr Base kSentinel{static_cast<char>(255),
                                  static_cast<char>(255)};
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
