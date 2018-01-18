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

/**
  @ingroup GBComputation
  @brief hash table class, based on discussions with Christian Eder
  @details For a homogeneous system it is advisable to recreate this for each
    matrix. For inhomogeneous systems it may be advantageous to reuse the table
    as one proceeds through the system.
*/
class F4_Hash {

  protected:

    NVAR_TYPE n; /**< @brief number of variables in the monomials stored */
    WT_TYPE * weights; /**< @brief randomized weighting for each exponent */
    static const size_t MAXIMUM = 1 << 18; /**< @brief number of entries in table */
    list<pair<const Monomial *, size_t> > * table; /**< @brief the actual table */

  public:

    /**
      @brief allocates the table and sets up randomized hash function
      @param num_vars number of variables in the monomials this table will check
    */
    F4_Hash(NVAR_TYPE num_vars) {
      table = new list<pair<const Monomial *, size_t> > [MAXIMUM];
      n = num_vars;
      weights = new WT_TYPE[n];
      default_random_engine generator(
          std::chrono::system_clock::now().time_since_epoch().count()
      );
      uniform_int_distribution<int> rand(0, MAXIMUM - 1);
      for (NVAR_TYPE i = 0; i < n; ++i)
        weights[i] = rand(generator);
    }

    /**
      @brief releases the table and the randomized hash function
    */
    ~F4_Hash() {
      delete [] weights;
      delete [] table;
    }

    /**
      @brief determines the index of @p t in the lookup table
      @param t @c Monomial
      @return which table entry contains the list that contains @p t
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
      @return which table entry contains the list that contains \f$tu\f$
    */
    size_t get_index(const Monomial & t, const Monomial & u) {
      get_index(t, u.log());
    }

    /**
      @brief determines the index of \f$tu\f$ in the lookup table
      @param t @c Monomial
      @param u @c exponent vector
      @return which table entry contains the list that contains \f$tu\f$
    */
    size_t get_index(const Monomial & t, const EXP_TYPE * u) {
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
      @brief looks up the location of \f$tu\f$ in the array; that is,
        the second argument to @c add_monomial() when \f$tu\f$ was the first
      @details This will not add \f$tu\f$ if it isn't already there, so
        call it only after the table has been completely set up!
      @param t @c Monomial
      @param u @c Monomial
      @return location of \f$tu\f$ in the array
    */
    size_t lookup_product(const Monomial & t, const Monomial & u) {
      auto i = get_index(t, u);
      auto curr = table[i].begin();
      while (curr != table[i].end() and not curr->first->like_multiple(t, u))
        ++curr;
      return curr->second;
    }

    /**
      @brief looks up the location of \f$tu\f$ in the array; that is,
        the second argument to @c add_monomial() when \f$tu\f$ was the first
      @details This will not add \f$tu\f$ if it isn't already there, so
        call it only after the table has been completely set up!
      @param t @c Monomial
      @param u @c exponent vector
      @return location of \f$tu\f$ in the array
    */
    size_t lookup_product(const Monomial & t, const EXP_TYPE * u) {
      auto i = get_index(t, u);
      auto curr = table[i].begin();
      while (curr != table[i].end() and not curr->first->like_multiple(u, t))
        ++curr;
      return curr->second;
    }

    friend ostream & operator << (ostream &, const F4_Hash &);

};

/**
  @brief prints the non-empty locations as a list of tuples (monomial, location)
*/
ostream & operator << (ostream & os, const F4_Hash & ht);

#endif