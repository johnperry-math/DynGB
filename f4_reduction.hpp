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

#include "system_constants.hpp"

#include "fields.hpp"
#include "monomial.hpp"
#include "monomial_ordering.hpp"
#include "polynomial.hpp"
#include "critical_pair.hpp"
#include "f4_hash.hpp"

#include <set>
#include <list>
#include <vector>
#include <cstdlib>
#include <utility>

using std::set;
using std::list;
using std::pair;
using std::vector;

/**
  @ingroup GBComputation
  @brief equivalent to @c buchberger(), but for Faug&egrave;re&rsquo;s F4 algorithm
  @param F generators of a polynomial ideal
  @return a Gr&ouml;bner basis of the ideal generated by @f$ F @f$
    with respect to the ordering already assigned to its polynomials
*/
list<Constant_Polynomial *> f4_control(const list<Abstract_Polynomial *> &F);

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
    @param P list of critical pairs that will create a new basis; matrix will
        have this many rows
    @param B list of polynomials currently in the basis
  */
  F4_Reduction_Data(
      const list<Critical_Pair_Basic *> & P,
      const list<Abstract_Polynomial *> & B
  );
  /**
    @brief adds monomials of @f$ ug @f$ to @c M_build
    @details
      If @c g is an s-polynomial generator, set the flag to @c true,
      @c ti to the beginning of @c M_build, @c ri to the beginning of @c R_build,
      and @c u to the first multiple in the critical pair.
      Setting the flag to @c true signals that @c ti and @c ri are not at the
      correct insertion point, so the function needs to find this point.
      In this case, @c ti and @c ri will change.

      If @c g is not a generator, then @c ti and @c ri should point
      to appropriate locations for the monomial divisible by
      the leading term of @c g. The algorithm advances both a copy of @c ti
      and an iterator on @c g to insert new monomials.
      In this case, the algorithm modifies neither @c ti nor @c ri.
    @param ti iterator through @c M_build; see details for more information
    @param ri iterator through @c R_build; see details for more information
    @param g polynomial, such as an s-polynomial generator or a reducer
    @param u monomial to multiply to @p g to add monomials
    @param new_row whether @p g is adding a new row to the basis
  */
  void add_monomials(
      list<Monomial *>::iterator & ti,
      list<Abstract_Polynomial *>::iterator & ri,
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
  void initialize_many(const list<Critical_Pair_Basic *> & P);
  ///@}
  /** @name Destruction */
  ///@{
  /**
    @brief releases space for the matrix and deletes any strategies not already
    set to @c nullptr
  */
  ~F4_Reduction_Data();
  ///@}
  /** @name Basic properties */
  ///@{
  /**
    @brief returns @c true iff all the entries are 0
  */
  bool is_zero();
  /**
    @brief returns the strategies currently in use
  */
  vector<Poly_Sugar_Data *> get_strategies() { return strategies; }
  /** @brief basic properties */
  unsigned number_of_rows() { return A.size(); }
  ///@}
  /** @name Conversion */
  ///@{
  /**
    @brief converts @c this to a vector of Constant_Polynomial
      and returns the result
  */
  vector<Constant_Polynomial *> finalize();
  ///@}
  /** @name Modification */
  ///@{
  /** @brief clears the strategy; do this if you have saved it elsewhere */
  void clear_strategy(unsigned i) { strategies[i] = nullptr; }
  ///@}
  /** @name Computation */
  ///@{
  /** @brief reduces polynomials */
  void reduce_by_old();
  /** @brief reduces polynomials */
  void reduce_by_new();
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
  ///@}
protected:
  /*void check_consistency() {
    for (unsigned i = 0; i < num_rows; ++i) {
      const auto & Ai = A[i];
      unsigned n = 0;
      for (auto a : Ai)
        if (a != 0) ++n;
      if (n != nonzero_entries[i]) cout << "row " << i << " inconsistent\n";
    }
  }*/
  /**
    @brief creates rows of the matrix indexed by the specified pairs,
      starting at the specified row
  */
  void initialize_some_rows(const list<Critical_Pair_Basic *> &, unsigned);
  /** @brief reduces the specified set of rows, suitable for multithreading */
  void reduce_my_rows(const set<unsigned> &);
  /**
    @brief reduces the specified set of rows by the specified row,
      suitable for multithreading
  */
  void reduce_my_new_rows(unsigned, const Prime_Field &, const set<unsigned> &);
  /** @brief number of columns in the polynomial */
  unsigned num_cols;
  /** @brief number of rows in the matrix */
  unsigned num_rows;
  /** @brief monomials for each matrix */
  vector<Monomial *> M;
  /** @brief hash table of the monomials in @c M */
  F4_Hash M_table;
  /**
    @brief coefficient data in sparse representation; each vector entry is a
      subrow of the dense matrix
  */
  vector<vector<COEF_TYPE> > A;
  /** @brief index of the starting point of this row in the dense matrix*/
  vector<unsigned> offset;
  /**
    @brief index of the first non-zero entry of this row in the dense matrix
      (counted from offset)
  */
  vector<unsigned> head;
  /** @brief monomials while building */
  list<Monomial *> M_build;
  /** @brief indices of reducers for the corresponding elements of @c M */
  list<Abstract_Polynomial *> R_build;
  /** @brief finalized list of indices of reducers for the corresponding monomials of @c f */
  vector<Abstract_Polynomial *> R;
  /** @brief reducers actually generated */
  vector<vector<pair<unsigned, COEF_TYPE> > > R_built;
  /** @brief count of threads that have actually read a generated actually built */
  vector<unsigned> num_readers;
  /** @brief current basis of ideal */
  const list<Abstract_Polynomial *> & G;
  /** @brief how the monomials are ordered */
  const Monomial_Ordering * mord;
  /** @brief number of nonzero entries of each row of A */
  vector<unsigned> nonzero_entries;
  /** @brief polynomial ring */
  Polynomial_Ring & Rx;
  /** @brief strategy data for each polynomial */
  vector<Poly_Sugar_Data *> strategies;
};

#endif