#ifndef __ALGORITHM_BUCHBERGER_BASIC_HPP_
#define __ALGORITHM_BUCHBERGER_BASIC_HPP_

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

#include <set>
#include <list>
#include <iostream>
#include <iterator>

using std::list;
using std::cout; using std::endl;
using std::next;
using std::set;
using std::iterator; using std::forward_iterator_tag;

#include "system_constants.hpp"

#include "fields.hpp"
#include "monomial.hpp"
#include "polynomial.hpp"
#include "critical_pair.hpp"
#include "normal_strategy.hpp"
#include "polynomial_array.hpp"
#include "polynomial_geobucket.hpp"
#include "polynomial_linked_list.hpp"
#include "polynomial_double_buffered.hpp"

#include "sugar_strategy.hpp"
#include "weighted_sugar_strategy.hpp"

/**
  @defgroup GBComputation Gr&ouml;bner basis computation
  @brief classes related directly to Gr&ouml;bner basis computation
*/

/**
  @brief Checks whether @p p is in danger of forming a Buchberger triple
  with some pair listed in @p C. We need this to avoid deleting useful new pairs.
  @ingroup GBComputation
  @param p a critical pair of some sort
  @param C a list of critical pairs
  @return @c true iff @p p does not form a Buchberger triple with some pair in @p C
*/
template <typename T>
bool no_triplet(const T *p, const list<T *>C) {
  bool result = true;
  for (T * c : C) {
    const Monomial & u = c->lcm();
    bool degrees_smaller = true;
    for (NVAR_TYPE i = 0; degrees_smaller and i < u.num_vars(); ++i)
      degrees_smaller = degrees_smaller and u.degree(i) <= p->lcm_degree(i);
    result = result and !degrees_smaller;
    if (not result) break;
  }
  return result;
}

template bool no_triplet<Critical_Pair_Basic>(
    const Critical_Pair_Basic *, const list<Critical_Pair_Basic *>
);

/**
  @brief Checks if the lcm of @p t and @p u is like the lcm stored in @p p.
  @param t a monomial
  @param u a monomial
  @param p a critical pair
  @ingroup GBComputation
  @return \c true iff \f$ lcm(t,u) \f$ is @p p&rsquo;s lcm
*/
bool lcm_alike(
    const Monomial & t,
    const Monomial & u,
    const Critical_Pair_Basic * p
);

/**
  @brief Remove redundant polynomials from <c>G</c>.
  @ingroup GBComputation
  @details (A polynomial <i>g</i> in <i>G</i> is redundant when we can find <i>h</i>
  in <i>G</i> whose leading monomial divides <i>g</i>&rsquo;s leading monomial.)
  @param G list of generators of a polynomial ideal, of which some may be redundant
  @return a list of polynomials without the redundant elements
*/
list<Abstract_Polynomial *> reduce_basis(list<Abstract_Polynomial *>G);

/**
  @brief Implementation of Gebauer-Moeller algorithm.
  Based on description in Becker and Weispfenning (1993).
  @ingroup GBComputation
  @param P list of critical pairs
  @param G current basis
  @param r polynomial to add to basis (and to generate new pairs)
  @param strategy how to sort pairs
*/
void gm_update(
    list<Critical_Pair_Basic *> & P,
    list<Abstract_Polynomial *> & G,
    Abstract_Polynomial * r,
    unsigned strategy
);

/**
  @brief A brief report on the number of critical pairs. If <c>verbose</c> is true,
  also lists them.
  @ingroup GBComputation
  @param P a list of critical pairs
  @param verbose whether to list the pairs as well as report their number
*/
template <typename T>
void report_critical_pairs(const list<T *>P, bool verbose = false) {
  cout << P.size() << " critical pairs remaining\n";
  if (verbose)
    for (T * p : P)
      cout << '\t' << *p << endl;
}

template void report_critical_pairs<Critical_Pair_Basic>(
    const list<Critical_Pair_Basic *>, bool
);

/**
  @brief used to sort polynomials by leading monomial
  @ingroup GBComputation
*/
struct smaller_lm {
  /**
    @brief returns \c true iff \p f&rsquo;s leading monomial is
      smaller than \c g&rsquo;s
    @param f a polynomial of some sort
    @param g a polynomial of some sort
    @return  \c true iff \p f&rsquo;s leading monomial is
      smaller than \c g&rsquo;s
  */
  bool operator()(Abstract_Polynomial *f, Abstract_Polynomial *g) {
    return f->leading_monomial() < g->leading_monomial();
  }
};

/**
  @brief checks that \p G is a Gr&ouml;bner basis by verifying each s-polynomial
    reduces to zero
  @ingroup GBComputation
  @param G list of generators of an ideal
  @param strategy how to select critical pairs 
*/
void check_correctness(list<Constant_Polynomial *>G, StrategyFlags strategy);

/**
  @brief prints the number of polynomials in the basis
  @param verbose print each polynomial&rsquo;s leading term
  @param very_verbose print each polynomial
*/
void report_basis(
    list<Abstract_Polynomial *> G,
    bool verbose=false,
    bool very_verbose=false
);

/**
  @brief gives a summary of information in \c p, with additional information
      depending on \c strategy
*/
void report_front_pair(Critical_Pair_Basic *p, unsigned strategy);

/**
  @brief Implementation of Buchberger&rsquo;s algorithm.
  @ingroup GBComputation
  @param F generators of the ideal whose Gr&ouml;bner basis you&rsquo;d like to
    compute
  @param rep which polynomial representation to use for the s-polynomials
      (default is GEOBUCKETS)
  @param strategy which strategy to use when selecting a critical pair
      (default is SUGAR_STRATEGY)
  @param strategy_weights if using a weighted sugar strategy, place an array
      of weights here
  @return list of polynomials in a Gr&ouml;bner basis of \c F
*/
list<Constant_Polynomial *> buchberger(
    const list<Abstract_Polynomial *> &F,
    SPolyCreationFlags rep = GEOBUCKETS,
    StrategyFlags strategy = SUGAR_STRATEGY,
    WT_TYPE * strategy_weights = nullptr
);

/**
  @brief Applies the strategy to find the &ldquo;smallest&rdquo; critical pair.
  @ingroup GBComputation

  @details Rather than sort each time, which typically shuffles
  about a lot of pairs that will later be deleted,
  we simply pass through the list once, find the smallest one,
  and move it to the front.
  @param P a list of critical pairs
*/
template <typename T>
void sort_pairs_by_strategy(list<T *> & P) {
  // identify lowest lcm
  typename list<T *>::iterator minkey_i = P.begin();
  const Pair_Strategy_Data * minkey_data = (*minkey_i)->pair_key();
  for (typename list<T *>::iterator pi = next(minkey_i);
       pi != P.end();
       ++pi
      )
  {
    const Pair_Strategy_Data * pairkey = (*pi)->pair_key();
    if (*pairkey < *minkey_data) {
      minkey_i = pi;
      minkey_data = pairkey;
    }
  }
  // move it to front
  T * minkey = *minkey_i;
  P.erase(minkey_i);
  P.push_front(minkey);
}

template void sort_pairs_by_strategy<Critical_Pair_Basic>(
    list<Critical_Pair_Basic *> &
);

#endif