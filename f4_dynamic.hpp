#ifndef __F4_REDUCTION_HPP__
#define __F4_REDUCTION_HPP__

/*****************************************************************************\
* This file is part of DynGB.                                                 *
*                                                                             *
* DynGB is free software: you can redistribute it and/or modify               *
* it under the terms of the GNU General Public License as published by        *
* the Free Software Foundation, either version 2 of the License, or           *
* (at your option) any later version.                                         *
*                                                                             *
* DynGB is distributed in the hope that it will be useful,                    *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with DynGB. If not, see <http://www.gnu.org/licenses/>.               *
\*****************************************************************************/

#include <set>
using std::set;
#include <list>
using std::list;
#include <vector>
using std::vector;
#include <cstdlib>
#include <utility>
using std::pair;
#include <functional>
using std::function;
#include <mutex>
using std::mutex;

#include "system_constants.hpp"

#include "fields.hpp"
#include "monomial.hpp"
#include "monomial_ordering.hpp"
#include "polynomial.hpp"
#include "critical_pair.hpp"
#include "f4_hash.hpp"

#include "lp_solver.hpp"
using LP_Solvers::LP_Solver;
#include "dynamic_engine.hpp"
using Dynamic_Engine::PP_With_Ideal;
using Dynamic_Engine::Dynamic_Heuristic;

/**
  @brief enumeration of styles of row analysis for selecting a new ordering
*/
enum class Analysis {
  row_sequential, /**< analyze only one row (topmost unprocessed) */
  whole_matrix    /**< analyze every row of the matrix */
};

/**
  @ingroup GBComputation
  @brief equivalent to @c buchberger(), but for Faug&egrave;re&rsquo;s F4 algorithm
  @param F generators of a polynomial ideal
  @param hash_table an @c F4_Hash to store monomials computed for the basis
  @param static_algorithm whether the algorithm should run traditionally or
      using dynamic techniques
  @param max_refinements maximum number of refinements per matrix;
      default (0) means no maximum
  @return a Gr&ouml;bner basis of the ideal generated by @f$ F @f$
    with respect to the ordering already assigned to its polynomials
*/
list<Abstract_Polynomial *> f4_control(
    const list<Abstract_Polynomial *> &F,
    F4_Hash & hash_table,
    const bool static_algorithm = true,
    const unsigned max_refinements = 0,
    const Analysis style = Analysis::row_sequential
);

extern unsigned location_of_monomial_index(
    vector< pair< unsigned, COEF_TYPE > > & row, unsigned i
);

/**
  @brief used to compare monomials for STL containers such as @c map
*/
struct MonCmp {
  /**
    @brief returns @c true iff the monomial pointed to by @p t is smaller than
      the monomial pointed to by @p u
    @param t monomial we're comparing
    @param u monomial we're comparing
    @return @c true iff the monomial pointed to by @p t is smaller than
      the monomial pointed to by @p u
  */
  bool operator()(
        const Monomial * t,
        const Monomial * u
  ) const {
    return *t < *u;
  }
};

/**
  @brief Implementation of Faug&egrave;re&rsquo;s F4 algorithm.
  @ingroup GBComputation
  @details Currently computes a Gr&ouml;bner basis by selecting several
    s-polynomials of lowest lcm degree. Data is stored in a semi-sparse matrix
    format, with each row a contiguous array of entries (@c A).
    Each row of @c A begins at the position indicated
    by the corresponding @c offset,
    and its first non-zero entry appears at the position indicated by
    the corresponding @c head.
    So the leading coefficient of the polynomial in row @c k appears in
    <c>A[head[k]]</c> and the leading monomial appears in
    <c>M[offset[k]+head[k]]</c>.
    Starting from <c>head[k]</c>, row @c k is actually dense.
    While this is not sparse, it does succeed in saving more space
    than one might expect.
*/
class F4_Reduction_Data {
public:
  /** @name Construction */
  ///@{
  /**
    @brief encapsulation of one step of the F4 algorithm for the polynomials
        indicated by @p P and @p B
    @param curr_ord current monomial ordering
    @param P list of critical pairs that will create a new basis; matrix will
        have this many rows
    @param B list of polynomials currently in the basis
    @param method heuristic the system should use in choosing leading monomials
  */
  F4_Reduction_Data(
      WGrevlex * curr_ord,
      const list<Critical_Pair_Dynamic *> & P,
      const list<Abstract_Polynomial *> & B,
      Dynamic_Heuristic method
  );
  /**
    @brief adds monomials of @f$ ug @f$ to @c M_builder
    @param curr_ord current monomial ordering
    @param g polynomial, such as an s-polynomial generator or a reducer
    @param u monomial to multiply to @p g to add monomials
    @param new_row if @c true, adds the leading monomial; otherwise, adds only
        trailing monomials
  */
  void add_monomials(
      const WGrevlex * curr_ord,
      const Abstract_Polynomial *g,
      const Monomial & u,
      bool new_row = false
  );
  /**
    @brief creates the matrix
    @details called internally by constructor, which has already performed some
        setup; inadvisable to use elsewhere
    @param P list of critical pairs that will generate the matrix
  */
  void initialize_many(const list<Critical_Pair_Dynamic *> & P);
  ///@}
  /** @name Destruction */
  ///@{
  /**
    @brief releases space for the matrix
  */
  ~F4_Reduction_Data();
  ///@}
  /** @name Basic properties */
  ///@{
  /**
    @brief returns @c true iff all the entries are 0
  */
  bool is_zero();
  /** @brief returns the number of nonzero entries on the indicated row */
  unsigned number_of_nonzero_entries(unsigned i) { return nonzero_entries[i]; }
  /** @brief basic properties */
  unsigned number_of_rows() const { return A.size(); }
  /** @brief set of monomials in row @p i that are not zero (absolute index) */
  void monomials_in_row(unsigned i, list<int> &) const;
  /** @brief current ordering of monomials */
  const WGrevlex * current_ordering() const { return mord; }
  /**
    @brief returns the number of compatible monomials detected in the given row
  */
  unsigned number_of_compatibles(unsigned row) {
    return compatible_pps[row].size();
  }
  /** @brief find the index of the monomial of greatest weight in this row */
  unsigned head_monomial_index(unsigned i, bool static_algorithm = false) {
    const auto & Ai = A[i];
    unsigned result = Ai[0].first;
    if (not static_algorithm) { 
      auto & t = *M[result];
      for (unsigned k = 1; k < Ai.size(); ++k) {
        if (mord->first_smaller(t, *M[Ai[k].first])) {
          result = Ai[k].first;
          t = *M[result];
        }
      }
    }
    return result;
  }
  ///@}
  /** @name Conversion */
  ///@{
  /**
    @brief converts indicated row to a Polynomial_Hashed, returns result
  */
  Polynomial_Hashed * finalize(unsigned, vector< Monomial * > &, F4_Hash &);
  ///@}
  /** @name Modification */
  ///@{
  ///@}
  /** @name Computation */
  ///@{
  /** @brief reduces polynomials */
  void reduce_by_old();
  /** @brief reduces polynomials */
  void reduce_by_new(unsigned i, unsigned lhead_i, const set<unsigned> &);
  /**
    @ingroup GBComputation
    @author John Perry
    @date 2019
    @brief Create constraints for a candidate LPP.
    @param pp_I pair of PP with the ideal it would have.
    @param monomials_for_comparison indices of monomials used to generate constraints with LPP
    @param result the new constraints
  */
  void constraints_for_new_pp(
    const PP_With_Ideal &I,
    const set<int> &monomials_for_comparison,
    vector<Constraint> &result
  );
  /**
    @ingroup GBComputation
    @author John Perry
    @date 2019
    @brief identifies potential leading monomials and sorts them by preference
    @param my_row which row of the matrix we are working on
    @param currentLPP current leading monomial of this polynomial
    @param allPP_indices indices of the support of a polynomial in need of a new choice of LPP
  */
  void create_and_sort_ideals(
      int my_row,
      const list<Monomial> & CurrentLPPs,
      const Dense_Univariate_Integer_Polynomial * current_hilbert_numerator,
      const list<Abstract_Polynomial *> & CurrentPolys,
      const list<Critical_Pair_Dynamic *> & crit_pairs,
      const Ray & w,
      Dynamic_Heuristic method
  );
  /**
    @ingroup GBComputation
    @author John Perry
    @date 2019
    @brief refines the ordering, if possible
    @param my_row the row of the matrix that we are refining
    @param currSkel the current skeleton, corresponding to the choices `CurrentLPPs`
    @param CurrentPolys list of current polynomials in the basis (needed to verify consistency)
  */
  pair< bool, LP_Solver * > refine(
      unsigned my_row,
      LP_Solver * currSkel,
      const list<Abstract_Polynomial *> & CurrentPolys
  );
  /**
    @ingroup GBComputation
    @author John Perry
    @date 2019
    @brief determines the weights of each monomial, according to the given skeleton
    @param skel current skeleton
  */
  void recache_weights(LP_Solver * skel);
  /**
    @ingroup GBComputation
    @author John Perry
    @date 2019
    @brief cleans up various records after we compute a new ordering
    @param my_row the row of the matrix that we are refining
    @param w the old ordering
    @param current_hilbert_numerator the basis' current Hilbert numerator
    @param currSkel the current skeleton
    @param winner the @c PP_With_Ideal with information relating to the selected monomial
    @param ordering_changed whether the ordering has changed (written, not read)
  */
  /**
    @ingroup GBComputation
    @author John Perry
    @date 2019
    @brief prints the weights cached for each monomial, according to most recent
      skeleton used
  */
  void print_cached_weights();
  void reassign(
      unsigned my_row,
      const Ray & w,
      Dense_Univariate_Integer_Polynomial * & current_hilbert_numerator,
      LP_Solver * currSkel,
      const PP_With_Ideal & winner,
      bool & ordering_changed
  );
  /**
    @brief selects leading monomials for one remaining nonzero rows
    @details This functions basically fills the role that select_monomial()
      fills in the dynamic Buchberger algorithm.
      Unlike the other select_dynamic(), however,
      it chooses only one row at a time.
    @param unprocessed set of unprocessed rows;
      the result will from from this set
    @param T list of current leading monomials
    @param G list of current basis polynomials
    @param P list of current critical pairs
    @param curr_ord current monomial ordering
    @param skel the current LP_Solver that defines the term ordering
    @param row_option which approach to take to selecting a row
    @return the row with the recommended new ordering
  */
  unsigned select_dynamic_single(
      set<unsigned> & unprocessed,
      list<Monomial> & T,
      const list<Abstract_Polynomial *> G,
      const list<Critical_Pair_Dynamic *> & P,
      WGrevlex * curr_ord,
      LP_Solver * & skel,
      const Analysis & style = Analysis::row_sequential
  );
  /**
    @brief eliminates duplicates of rows:
      later rows that are identical to earlier rows will be eliminated
    @param in_use rows that are currently still in use by the matrix
  */
  void simplify_identical_rows(set<unsigned> & in_use);
  /** @brief make row monic at indicated index */
  void normalize(unsigned i, unsigned lhead_i) {
    auto & Ai = A[i];
    auto mod = Rx.ground_field().modulus();
    COEF_TYPE a = Rx.ground_field().inverse(Ai[location_of_monomial_index(Ai, lhead_i)].second);
    for (auto k = 0; k < Ai.size(); ++k) {
      Ai[k].second *= a; Ai[k].second %= mod;
    }
  }
  ///@}
  /** @name I/O */
  ///@{
  /** @brief lists the reducers selected for each column, in order */
  void list_reducers();
  /**
    @brief prints the matrix
    @param show_data whether to show the monomials that correspond to each column
  */
  void print_matrix(bool show_data=false);
  /**
    @brief prints indicated row of the matrix
    @details If the second argument is true, the row is printed as a polynomial.
  */
  void print_row(unsigned, bool=false);
  /**
    @brief prints the pairs of monomials and reducers that are being built
  */
  void print_builder();
  ///@}
protected:
  /**
    @brief creates rows of the matrix indexed by the specified pairs,
      starting at the specified row
  */
  void initialize_some_rows(const list<Critical_Pair_Dynamic *> &, unsigned);
  /** @brief builds a reducer for the specified column */
  void build_reducer(unsigned);
  /**
    @brief reduces the specified set of rows, suitable for multithreading
    @param rows the rows to reduce
    @param buffer space to expand each row in @p rows and then reduce
    @param next @c expand needs this to track monomial locations
  */
  void reduce_my_rows(
      const vector<int> &rows,
      vector <COEF_TYPE> &buffer, vector <unsigned> & prev, vector <unsigned> & next
  );
  /**
    @brief reduces the specified set of rows by the specified row,
      suitable for multithreading
    @param i row of the row to use when reducing
    @param lhead_i the location of the (new) leading monomial of @p i
    @param buffer space to expand each row for reduction
    @param next @c expand needs this to track monomial locations
    @param to_reduce which rows of the matrix this thread will reduce
  */
  void reduce_my_new_rows(
      unsigned i,
      unsigned lhead_i,
      vector< COEF_TYPE > & buffer,
      vector< unsigned > & prev, vector< unsigned > & next,
      const set<unsigned> & to_reduce,
      unsigned mod
  );
  /** @brief number of columns in the polynomial */
  unsigned num_cols;
  /** @brief number of rows in the matrix */
  unsigned num_rows;
  /** @brief monomials for each matrix */
  vector< Monomial * > M;
  /** @brief hash table of the monomials in @c M */
  F4_Hash M_table;
  /** @brief locks on the reducers */
  vector<mutex> red_lock;
  /**
    @brief coefficient data in sparse representation; each vector entry 
        indicates position and coefficient
  */
  vector< vector < pair < unsigned, COEF_TYPE > > > A;
  /** @brief whether the row was modified during reduction */
  vector<bool> dirty;
  /** @brief index of the preferred head term of this row (absolute index) */
  vector<unsigned> pref_head;
  /** @brief number of nonzero entries of each expanded row of A */
  vector<unsigned> nonzero_entries;
  /** @brief compatible pp's for each row (absolute index) */
  vector< list<int> > compatible_pps;
  /** @brief cache of monomial weights under current ordering */
  vector< vector<unsigned> > pp_weights;
  /** @brief potential ideals for the given row */
  vector< list<PP_With_Ideal> > potential_ideals;
  /** @brief storage of monomials and reducers while preprocessing */
  map<Monomial *, Abstract_Polynomial *, MonCmp> M_builder;
  /** @brief finalized list of indices of reducers for the corresponding monomials of @c f */
  vector<Abstract_Polynomial *> R;
  /** @brief mutexes for building rows of @c R */
  vector<mutex> red_mutex;
  /** @brief reducers actually generated */
  vector<vector<pair<unsigned, COEF_TYPE> > > R_built;
  /** @brief current basis of ideal */
  const list<Abstract_Polynomial *> & G;
  /** @brief how the monomials are ordered */
  WGrevlex * mord;
  /** @brief polynomial ring */
  Polynomial_Ring & Rx;
  /** @brief strategy data for each polynomial */
  vector<Poly_Sugar_Data *> strategies;
  /**
    @brief each entry contains a set of monomials that could be a leading
      monomial of the corresponding row of the matrix
  */
  vector<list<PP_With_Ideal> > I;
  /** @brief heuristic the system should use in choosing leading monomials */
  Dynamic_Heuristic heur;
  /** @brief function corresponding to @c heur */
  function <bool(PP_With_Ideal &, PP_With_Ideal &)> heuristic_judges_smaller;
  friend void compatible_pp(
    int my_row,
    F4_Reduction_Data & F4,
    const LP_Solver * skel,
    bool & stop,
    vector<bool> & completed
  );

};

#endif