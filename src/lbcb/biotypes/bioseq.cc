#include "lbcb/biotypes/bioseq.h"

namespace lbcb {

Sequence::Sequence(std::string_view name, std::string_view data)
    : name_((name)), data_(data) {
  // TODO: magic
}
Sequence::Sequence(std::string_view name, std::string_view data,
                   std::string_view quality)
    : name_(name), data_(data), quality_(quality) {}

}  // namespace lbcb
