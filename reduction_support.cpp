#ifndef __REDUCTION_SUPPORT_CPP_
#define __REDUCTION_SUPPORT_CPP_

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

#include "reduction_support.hpp"

/**
  @ingroup GBComputation
  @brief reduce the polynomial @c **sp by @c *g
  @details This is only a top reduction, though it is repetitive, as in:
      it continues to reduce @c **sp by @c g until it can be done no more.
  @param s a Mutable_Polynomial, the one you want to reduce
  @param g a polynomial, the one you want to use in reduction
  @param comm_id pass this for verbose messages during multiprocessing/MPI
*/
void top_reduce(Mutable_Polynomial *s, Abstract_Polynomial * g, int comm_id) {
  bool verbose = false;
  bool very_verbose = false;
  Monomial & u = g->leading_monomial();
  while ((not s->is_zero()) and g->leading_monomial() | (s->leading_monomial())) {
    Monomial t(s->leading_monomial());
    if (verbose) cout << "top-reducing " << s->leading_coefficient() << ' ' << t << " by " << g->leading_monomial() << endl;
    if (very_verbose) cout << "top-reducing\n\t" << *s << "\nby\n\t" << *g << endl;
    t /= g->leading_monomial();
    Prime_Field_Element a = s->leading_coefficient();
    Prime_Field_Element b = g->leading_coefficient();
    a /= b;
    // s = s - atg
    //if (verbose) cout << comm_id << " entering pre-reduction tasks\n";
    if (s->strategy() != nullptr)
      s->strategy()->pre_reduction_tasks(t, *g);
    //if (verbose) cout << comm_id << " finished pre-reduction tasks\n";
    s->add_polynomial_multiple(a, t, *g, true);
  }
}

#endif