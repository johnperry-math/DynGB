#ifndef __REDUCTION_SUPPORT_HPP
#define __REDUCTION_SUPPORT_HPP

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
#include "polynomial.hpp"
#include "monomial.hpp"
#include "fields.hpp"

/**
  @brief Find a polynomial in the basis @p G that can reduce @p r.
  @param r the polynomial we want to reduce
  @param G a set of potential reductors, usually the basis of an ideal
  @return nullptr if no such polynomial exists.
  @ingroup GBComputation
*/
template <typename T>
Abstract_Polynomial * find_reducer(Abstract_Polynomial * r, const T & G) {
  Abstract_Polynomial * result = nullptr;
  const Monomial & t = r->leading_monomial();
  bool found = false;
  // loop through G until we find a result, if one exists
  for (Abstract_Polynomial * g : G)
  {
    if (
        g->leading_monomial() | t
        and (r->strategy() == nullptr or r->strategy()->valid_reduction(*r, *g)))
    {
      result = g;
      found = true;
      break;
    }
  }
  return result;
}

/**
  @ingroup GBComputation
  @brief reduce the polynomial @c **sp by @c *g
  @details This is only a top reduction, though it is repetitive.
*/
void top_reduce(Mutable_Polynomial *s, Abstract_Polynomial * g, int comm_id = 0);

/**
  @brief Reduce the polynomial r over the basis G.
  @ingroup GBComputation
  @param sp the s-polynomial we want to reduce
  @param G the generators of an ideal
  @param comm_id for MPI multiprocessing

  @details This is a complete reduction, not just the head,
  so the value of *sp may change.
*/
template <typename T>
void reduce_over_basis(Mutable_Polynomial **sp, const T & G, int comm_id=0) {
  Abstract_Polynomial * g; // used to loop through G
  Mutable_Polynomial * s = *sp; // s-poly
  Mutable_Polynomial * r = s->zero_polynomial(); // remainder / residue
  bool verbose = false;
  bool very_verbose = false;
  // continue reducing until s is zero
  while (!s->is_zero()) {
    if (verbose) cout << comm_id << " reducing " << s->leading_monomial() << "?\n"; 
    if (very_verbose) s->println();
    if ((g = find_reducer<T>(s, G))) {
      if (verbose) cout << comm_id << " yes! by " << g->leading_monomial() << endl;
      if (very_verbose) { cout << "yes! by "; g->println(); }
      top_reduce(s, g);
      if (very_verbose) {
        cout << "\tresults in ";
        if (s->is_zero())
          cout << 0 << endl;
        else
          s->println();
        if (not s->is_zero())
          cout << "\tresults in " << s->leading_monomial() << endl;
      }
    } else {
      if (verbose or very_verbose) cout << "no!\n";
      // no reducer; move leading term to residue and continue
      Abstract_Polynomial * t = s->detach_head();
      // t should already be smaller than what's already in r, so no danger
      r->add_last(t->leading_coefficient(), t->leading_monomial());
      delete t;
    }
  }
  //cout << comm_id << " wrapping up with " << *r << "\n";
  // move the strategy from s to r
  r->set_strategy(s->strategy());
  s->set_strategy(nullptr);
  delete s;
  *sp = r;
}

#endif