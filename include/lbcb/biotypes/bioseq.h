#ifndef LBCB_BIOTYPES_BIOSEQ_H_
#define LBCB_BIOTYPES_BIOSEQ_H_

#include <string_view>

namespace lbcb::biotypes {

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
};

}  // namespace lbcb::biotypes

#endif /* LBCB_BIOTYPES_BIOSEQ_H_  */
