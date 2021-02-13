#include "f4_hash.hpp"

ostream & operator << (ostream & os, const F4_Hash & ht) {
  for (size_t i = 0; i < ht.MAXIMUM; ++i) {
    if (ht.table[i].size() != 0) {
      os << i << ": ";
      for (auto iter = ht.table[i].begin(); iter != ht.table[i].end(); ++iter)
        os << '(' << *(iter->first)
            << ", {" << iter->second.location << ", " << iter->second.matrix_id
            << "} ) ";
      os << endl;
    }
  }
  return os;
}