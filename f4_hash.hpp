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

extern double get_index_time;

/**
  information for the monomial: location in the table, when last used
*/
struct usage_information {
  size_t location = 0;
  size_t matrix_location = 0;
  unsigned matrix_id = -1;
};

const usage_information NOT_FOUND;

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
    /**
      @brief number of entries in table
    */
    static const size_t MAXIMUM = 1 << 18;
    /**
      @brief the actual table
    */
    vector< vector<pair<const Monomial *, usage_information> > > table;
    /**
      @brief signature for the table, to help with caching computations
    */
    WT_TYPE signature;

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
      uniform_int_distribution<int> rand(0, MAXIMUM / n - 1 );
      for (NVAR_TYPE i = 0; i < n; ++i) {
        weights[i] = rand(generator);
      }
      signature = ( rand(generator) + rand(generator) ) % MAXIMUM;
    }

    /**
      @brief releases the table and the randomized hash function
    */
    ~F4_Hash() {
      delete [] weights;
      unsigned unused = 0, max_length = 0;
      for ( auto i = 0; i < table.size(); ++i ) {
        if ( table[i].size() == 0 )
          ++unused;
        else if ( table[i].size() > max_length )
          max_length = table[i].size();
      }
      cout << "unused " << unused << " of " << table.size() << endl;
      cout << "maximum length is " << max_length << endl;
    }

    /**
      @brief determines the index of @p t in the lookup table
      @param t @c Monomial
      @return which table entry contains the list that contains @p t
    */
    size_t get_index(const Monomial & t) const {
      DEG_TYPE index = t.cached_weighted_degree(weights);
      return index % MAXIMUM;
    }

    /**
      @brief determines the index of \f$tu\f$ in the lookup table
      @param t @c Monomial
      @param u @c Monomial
      @return which table entry contains the list that contains \f$tu\f$
    */
    size_t get_index(const Monomial & t, const Monomial & u) const {
      DEG_TYPE index = t.cached_weighted_degree(weights) + u.cached_weighted_degree(weights);
      return index % MAXIMUM;
    }

    /**
      @brief determines the index of \f$tu\f$ in the lookup table
      @param t @c Monomial
      @param u @c exponent vector
      @return which table entry contains the list that contains \f$tu\f$
    */
    //size_t get_index(const Monomial & t, const EXP_TYPE * u) {
    //  DEG_TYPE index = 0;
    //  for (NVAR_TYPE i = 0; i < n; ++i)
    //    index += ((t[i] + u[i]) * weights[i]);
    //  return index % MAXIMUM;
    //}

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
      @param matrix which matrix the location applies to
      @param location column for coefficients of @c t
    */
    void update_location(const Monomial * t, unsigned matrix, size_t location) {
      auto & list = table[get_index(*t)];
      auto curr = list.begin();
      while (not t->is_like(*(curr->first))) ++curr;
      curr->second = { curr->second.location, location, matrix };
    }

    /**
      @brief indicates whether a monomial has been added to the table
      @return @c true if and only if @p t has been added to the table
      @param t @c Monomial whose presence we'd like to determine
      @param matrix which matrix we want @c t to be in
    */
    //usage_information contains_in_matrix(const Monomial * t, unsigned matrix) {
    //  auto & list = table[get_index(*t)];
    //  auto curr = list.begin();
    //  while (curr != list.end() and not curr->first->is_like(*t)) ++curr;
    //  if (curr == list.end()) return NOT_FOUND;
    //  else return curr->second;
    //}

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
      @param matrix which matrix we want @c tu to be in
    */
    usage_information contains_product_in_matrix(
        const Monomial & t, const Monomial & u, unsigned matrix
    ) {
      //std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
      auto idx = get_index(t,u);
      //std::chrono::high_resolution_clock::time_point stop = std::chrono::high_resolution_clock::now();
      //get_index_time += std::chrono::duration_cast<std::chrono::duration<double> >(stop - start).count();
      auto & list = table[idx];
      auto curr = list.begin();
      while (curr != list.end() and not curr->first->like_multiple(t, u)) ++curr;
      if (curr == list.end()) return NOT_FOUND;
      else return curr->second;
    }

    bool contains_product(const Monomial & t, const Monomial & u) {
      auto & list = table[get_index(t, u)];
      auto curr = list.begin();
      //std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
      while (curr != list.end() and not ( curr->first->like_multiple(t, u)) ) ++curr;
      //std::chrono::high_resolution_clock::time_point stop = std::chrono::high_resolution_clock::now();
      //get_index_time += std::chrono::duration_cast<std::chrono::duration<double> >(stop - start).count();
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
      @param matrix which matrix the monomial is for
    */
    void add_monomial_for_matrix(
        const Monomial * t, size_t location, unsigned matrix = 0
    ) {
      auto & list = table[get_index(*t)];
      auto curr = list.begin();
      while (curr != list.end() and not curr->first->is_like(*t)) ++curr;
      if (curr != list.end()) {
        curr->second.location = location;
        curr->second.matrix_id = matrix;
      } else {
        usage_information new_value { location, 0, matrix };
        list.emplace_back(t, new_value );
      }
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
      usage_information new_value { location, 0, 0 };
      list.emplace_back(t, new_value);
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
      return curr->second.matrix_location;
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
    //size_t lookup_product(const Monomial & t, const EXP_TYPE * u) {
    //  auto & list = table[get_index(t, u)];
    //  auto curr = list.begin();
    //  while (curr != list.end() and not curr->first->like_multiple(u, t))
    //    ++curr;
    //  return curr->second.matrix_location;
    //}

    /**
      @brief indicates which location is associated with @p t
      @param t @c Monomial
      @return location of @p t in the array
    */
    unsigned operator[](const Monomial & t) const {
      auto & list = table[get_index(t)];
      auto curr = list.begin();
      while (not curr->first->is_like(t)) ++curr;
      return curr->second.location;
    }

    unsigned matrix_location(const Monomial & t) const {
      auto & list = table[get_index(t)];
      auto curr = list.begin();
      while (not curr->first->is_like(t)) ++curr;
      return curr->second.matrix_location;
    }

    /** @brief list the monomials stored at this monomial's index */
    void vomit(const Monomial & t) {
      auto & list = table[get_index(t)];
      cout << "bucket " << get_index(t) << endl;
      for (auto storage : list) {
        cout << "\t( " << *storage.first << " , { " << storage.second.location
              << ", " << storage.second.matrix_id << " } )\n";
      }
    }

    friend ostream & operator << (ostream &, const F4_Hash &);

    unsigned num_unused() {
      unsigned result = 0;
      for (auto entry : table) {
        if (entry.size() == 0) ++result;
      }
      return result;
    }
    
    unsigned num_total() { return table.size(); }
    
    unsigned load() {
      unsigned result = 0;
      for (auto entry : table) {
        result += entry.size();
      }
      return result / table.size();
    }

};

/**
  @brief prints the non-empty locations as a list of tuples (monomial, location)
*/
ostream & operator << (ostream & os, const F4_Hash & ht);

#endif
