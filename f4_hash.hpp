#ifndef __F4_HASH_
#define __F4_HASH_

#include <random>
using std::default_random_engine;
using std::uniform_int_distribution;

#include <chrono>

#include <list>
using std::list;

#include <utility>
using std::pair;

#include "system_constants.hpp"

#include "monomial.hpp"

class F4_Hash {

  protected:

    NVAR_TYPE n;
    WT_TYPE * weights;
    static const size_t MAXIMUM = 1 << 18;
    list<pair<const Monomial *, size_t> > table[MAXIMUM];

  public:

    F4_Hash(NVAR_TYPE num_vars) {
      n = num_vars;
      weights = new WT_TYPE[n];
      default_random_engine generator(
          std::chrono::system_clock::now().time_since_epoch().count()
      );
      uniform_int_distribution<int> rand(0, MAXIMUM - 1);
      for (NVAR_TYPE i = 0; i < n; ++i)
        weights[i] = rand(generator);
    }

    ~F4_Hash() { delete [] weights; }

    /**
      @brief determines the index of @p t in the lookup table
      @param t @c Monomial
    */
    size_t get_index(const Monomial & t) {
      DEG_TYPE index = 0;
      for (NVAR_TYPE i = 0; i < n; ++i)
        index += (t[i] * weights[i]) % MAXIMUM;
      return index % MAXIMUM;
    }

    /**
      @brief determines the index of \f$tu\f$ in the lookup table
      @param t @c Monomial
      @param u @c Monomial
    */
    size_t get_index(const Monomial & t, const Monomial & u) {
      DEG_TYPE index = 0;
      for (NVAR_TYPE i = 0; i < n; ++i)
        index += ((t[i] + u[i]) * weights[i]) % MAXIMUM;
      return index % MAXIMUM;
    }

    /**
      @brief adds a monomial to the lookup table
      @details Use only when creating the hash table.
        This does not check for uniqueness; it simply adds monomials!
        So make sure you are adding a monomial that isn't already there.
      @param t @c Monomial for future lookup
      @param location column for coefficients of @p t
    */
    void add_monomial(const Monomial * t, const size_t location) {
      table[get_index(*t)].emplace_back(t, location);
    }

    /**
      @brief looks up the location of \f$tu\f$ in the array
      @details This will not add \f$tu\f$ if it isn't already there, so
        call it only after the table has been completely set up!
      @param t @c Monomial
      @param t @c Monomial
    */
    size_t lookup_product(
        const Monomial & t, const Monomial & u
    ) {
      auto i = get_index(t, u);
      auto curr = table[i].begin();
      while (curr != table[i].end() and not curr->first->like_multiple(t, u))
        ++curr;
      return curr->second;
    }

    friend ostream & operator << (ostream &, const F4_Hash &);

};

ostream & operator << (ostream & os, const F4_Hash & ht) {
  for (size_t i = 0; i < ht.MAXIMUM; ++i) {
    if (ht.table[i].size() != 0) {
      os << i << ": ";
      for (auto iter = ht.table[i].begin(); iter != ht.table[i].end(); ++iter)
        os << '(' << *(iter->first) << ',' << iter->second << ") ";
      os << endl;
    }
  }
}

#endif