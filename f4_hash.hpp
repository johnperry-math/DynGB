#ifndef __F4_HASH_
#define __F4_HASH_

#include <random>
using std::default_random_engine;
using std::uniform_int_distribution;

#include <chrono>

#include <list>
using std::list;
#include <vector>
using std::vector;

#include <utility>
using std::pair;

#include "system_constants.hpp"

#include "monomial.hpp"
#include "polynomial_hashed.hpp"

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
    vector< vector<pair<const Monomial *, size_t> > > table; /**< @brief the actual table */

  public:

    /** @brief maximum length of a list in the table; currently for info only */
    size_t max_size = 0;
    /**
      @brief allocates the table and sets up randomized hash function
      @param num_vars number of variables in the monomials this table will check
    */
    explicit F4_Hash(NVAR_TYPE num_vars) : table(MAXIMUM) {
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
    }

    /**
      @brief determines the index of @p t in the lookup table
      @param t @c Monomial
      @return which table entry contains the list that contains @p t
    */
    size_t get_index(const Monomial & t) {
      DEG_TYPE index = t.weighted_degree(weights);
      return index % MAXIMUM;
    }

    /**
      @brief determines the index of \f$tu\f$ in the lookup table
      @param t @c Monomial
      @param u @c Monomial
      @return which table entry contains the list that contains \f$tu\f$
    */
    size_t get_index(const Monomial & t, const Monomial & u) {
      DEG_TYPE index = t.weighted_degree(weights) + u.weighted_degree(weights);
      return index % MAXIMUM;
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
        index += ((t[i] + u[i]) * weights[i]);
      return index % MAXIMUM;
    }

    /**
      @brief modifies the location of the specified monomial in the lookup table
      @details Use only in conjunction with monomials allocated with
        @c add_monomial(const Monomial *).
      @warning Invoking this on a monomial whose location hasn't changed or
        been reassigned will lead to errors and possibly termination.
        Furthermore, there are no guards for whether the monomial actually
        appears in the table, so if the monomial hasn't been added there will
        be errors and probably termination.
      @param t @c Monomial that has already been added to the table
      @param location column for coefficients of @c t
    */
    void update_location(const Monomial * t, size_t location) {
      auto & list = table[get_index(*t)];
      auto curr = list.begin();
      while (not t->is_like(*(curr->first))) ++curr;
      curr->second = location;
    }

    /**
      @brief indicates whether a monomial has been added to the table
      @return @c true if and only if @p t has been added to the table
      @param t @c Monomial whose presence we'd like to determine
    */
    bool contains(const Monomial * t) {
      auto & list = table[get_index(*t)];
      auto curr = list.begin();
      while (curr != list.end() and not curr->first->is_like(*t)) ++curr;
      return curr != list.end();
    }

    /**
      @brief indicates whether a product has been added to the table
      @return @c true if and only if \f$ tu \f$ has been added to the table
      @param t monomial whose product we're checking
      @param u monomial whose product we're checking
    */
    bool contains_product(const Monomial & t, const Monomial & u) {
      auto & list = table[get_index(t, u)];
      auto curr = list.begin();
      while (curr != list.end() and not curr->first->like_multiple(t, u)) ++curr;
      return curr != list.end();
    }

    /**
      @brief adds a monomial without an index to the lookup table
      @details Use only when creating the hash table.
        This does not check for uniqueness; it simply adds monomials!
        So make sure you are adding a monomial that isn't already there.
        For this version of the command, you must subsequently call
        @c update_location() so that the table has proper references.
      @warning The location of this monomial is <b>invalid</b>.
        Essentially, this function merely indicates that such an object
        needs to be tracked, and will eventually be updated for the table.
      @param t @c Monomial for future lookup
    */
    void add_monomial(const Monomial * t) {
      auto & list = table[get_index(*t)];
      list.emplace_back(t, 0);
      if (list.size() > max_size) max_size = list.size();
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
      auto & list = table[get_index(*t)];
      list.emplace_back(t, location);
      if (list.size() > max_size) max_size = list.size();
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
      auto & list = table[get_index(t, u)];
      auto curr = list.begin();
      while (curr != list.end() and not curr->first->like_multiple(t, u))
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
      auto & list = table[get_index(t, u)];
      auto curr = list.begin();
      while (curr != list.end() and not curr->first->like_multiple(u, t))
        ++curr;
      return curr->second;
    }

    /**
      @brief indicates which location is associated with @p t
      @param t @c Monomial
      @return location of @p t in the array
    */
    unsigned operator[](const Monomial & t) {
      auto & list = table[get_index(t)];
      auto curr = list.begin();
      while (not curr->first->is_like(t)) ++curr;
      return curr->second;
    }

    /** @brief list the monomials stored at this monomial's index */
    void vomit(const Monomial & t) {
      auto & list = table[get_index(t)];
      cout << "bucket " << get_index(t) << endl;
      for (auto storage : list) {
        cout << "\t( " << *storage.first << " , " << storage.second << " )\n";
      }
    }

    friend ostream & operator << (ostream &, const F4_Hash &);

};

/**
  @brief prints the non-empty locations as a list of tuples (monomial, location)
*/
ostream & operator << (ostream & os, const F4_Hash & ht);

#endif
