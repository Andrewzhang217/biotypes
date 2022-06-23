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

    iterator(const Sequence* sequence, std::size_t start_index);
    iterator(const iterator& other) = default;
    iterator(iterator&&) = default;
    iterator& operator=(const iterator& other) = default;
    iterator& operator=(iterator&& other) = default;
    reference operator*() noexcept;
    pointer operator->() const noexcept;
    value_type operator[](std::size_t pos) noexcept;
    iterator& operator+=(std::size_t diff) noexcept;
    iterator& operator-=(std::size_t diff) noexcept;
    iterator& operator++() noexcept;
    iterator& operator--() noexcept;
    const iterator operator++(int) noexcept;
    const iterator operator--(int) noexcept;
    iterator& operator+(std::size_t diff) const noexcept;
    iterator& operator-(std::size_t diff) const noexcept;
    bool operator==(const iterator& rhs) const noexcept;
    bool operator!=(const iterator& rhs) const noexcept;
    bool operator>(const iterator& rhs) const noexcept;
    bool operator<(const iterator& rhs) const noexcept;
    bool operator>=(const iterator& rhs) const noexcept;
    bool operator<=(const iterator& rhs) const noexcept;

   private:
    Base base_;
    const Sequence* sequence_;
    std::size_t pos_;
  };
  Sequence(std::string_view name, std::string_view data);
  Sequence(std::string_view name, std::string_view data,
           std::string_view quality);

  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] Base atBase(std::size_t pos) const;
  [[nodiscard]] char atValue(std::size_t pos) const;
  [[nodiscard]] char atQuality(std::size_t pos) const;
  [[nodiscard]] Sequence::iterator begin() const;
  [[nodiscard]] Sequence::iterator end() const;

 private:
  static constexpr Base kSentinel{static_cast<char>(255),
                                  static_cast<char>(255)};
  std::string name_;
  std::vector<std::uint64_t> compressed_data_;
  std::vector<std::uint64_t> compressed_quality_;
  std::size_t size_;
};
}  // namespace lbcb

#endif /* LBCB_BIOTYPES_BIOSEQ_H_  */
