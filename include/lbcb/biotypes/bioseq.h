#ifndef LBCB_BIOTYPES_BIOSEQ_H_
#define LBCB_BIOTYPES_BIOSEQ_H_

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
  Sequence(std::string_view name, std::string_view data);
  Sequence(std::string_view name, std::string_view data,
           std::string_view quality);

  // TODO: make iterator with value_type = Base

 private:
  std::string name_;
  std::string data_;
  std::string quality_;
};

}  // namespace lbcb

namespace lbcb::detail {
std::vector<std::uint64_t> Compress(std::string_view data);
}

#endif /* LBCB_BIOTYPES_BIOSEQ_H_  */
