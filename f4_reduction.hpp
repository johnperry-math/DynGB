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
* Foobar is distributed in the hope that it will be useful,                   *
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

#include <set>
#include <list>
#include <vector>
#include <cstdlib>

using std::set;
using std::list;
using std::vector;

/**
  @ingroup GBComputation
  @brief equivalent to @c buchberger(), but for Faug&egrave;re&rsquo;s F4 algorithm
*/
list<Constant_Polynomial *> f4_control(
    const list<Abstract_Polynomial *> &F,
    int method,
    unsigned strategy,
    WT_TYPE * strategy_weights
);

/**
  @brief Implementation of Faug&egrave;re&rsquo;s F4 algorithm.
  @ingroup GBComputation
  @details Currently computes a Gr&ouml;bner basis by selecting one s-polynomial
    at a time. Eventually needs to select several at a time, then reduce them.
*/
class F4_Reduction_Data {
public:
  /** @name Construction */
  ///@{
  F4_Reduction_Data(
      Critical_Pair_Basic & p,
      list<Abstract_Polynomial *> & B
  );
  /**
    @brief adds appropriate monomial multiples of @c g to @c M_build
    @details  If @c g is a generator, set @c ti to @c nullptr and
      @c u to the multiple necessary for the s-polynomial. This should be
      defined already in the critical pair.
      <b>Do not leave @c u as @c nullptr in this case!</b>

      If @c g is not a generator, then @c ti should point to the monomial
      divisible by the leading term of @c g. The algorithm advances both @c ti
      and a new iterator through @c g ahead one monomial and begins comparing.

      The iterator @c ri is an iterator through @c R_build.
  */
  void add_monomials(
      list<Monomial *>::iterator * ti,
      list<Abstract_Polynomial *>::iterator & ri,
      const Abstract_Polynomial *g,
      const Monomial * u = nullptr
  );
  /**
    @brief copy temporary to permanent reducers; free list of monomials;
      insert @c p's coefficients, aligned with multiples of @c t
  */
  void initialize(Abstract_Polynomial * p, const Monomial & t);
  ///@}
  /** @name Destruction */
  ///@{
  ~F4_Reduction_Data();
  ///@}
  /** @name Basic properties */
  ///@{
  /**
    @brief returns @c true iff all the entries are 0
  */
  bool is_zero();
  /**
    @brief returns the strategy currently in use
  */
  Poly_Strategy_Data * get_strategy() { return strategy; }
  ///@}
  /** @name Conversion */
  ///@{
  /** @brief converts @c this to a Constant_Polynomial and returns the result */
  Constant_Polynomial * finalize();
  ///@}
  /** @name Modification */
  ///@{
  /** @brief clears the strategy; do this if you have saved it elsewhere */
  void clear_strategy() { strategy = nullptr; }
  ///@}
  /** @name Computation */
  ///@{
  /** @brief advances @c head; useful after a computation */
  inline void advance_head() {
    while (head < len and A[head].is_zero()) ++head;
  }
  /**
    @brief reduces polynomial
    @warning There @b must be at least one monomial in the polynomial;
      we do not verify this initially
  */
  void reduce();
  ///@}
  /** @name I/O */
  ///@{
  void list_reducers();
  ///@}
protected:
  /** @brief head of the polynomial (location of leading term) */
  unsigned head;
  /** @brief length of the polynomial */
  unsigned len;
  /** @brief polynomial data */
  Monomial ** M;
  /** @brief polynomial data */
  Prime_Field_Element * A;
  /** @brief monomials of @c f */
  list<Monomial *> * M_build;
  /** @brief indices of reducers for the corresponding elements of @c M */
  list<Abstract_Polynomial *> * R_build;
  /** @brief finalized list of indices of reducers for the corresponding monomials of @c f */
  vector<Abstract_Polynomial *> R;
  /** @brief current basis of ideal */
  list<Abstract_Polynomial *> & G;
  /** @brief how the monomials are ordered */
  const Monomial_Ordering * mord;
  /** @brief number of nonzero entries of A */
  unsigned nonzero_entries;
  /** @brief polynomial ring */
  Polynomial_Ring * Rx;
  /** @brief strategy in use */
  Poly_Strategy_Data * strategy;
};

#endif