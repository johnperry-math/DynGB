#ifndef __CYCLIC_N_HPP_
#define __CYCLIC_N_HPP_

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
#include <iostream>

using std::set;
using std::cout; using std::endl;

#include "system_constants.hpp"

#include "fields.hpp"
#include "monomial.hpp"
#include "polynomial.hpp"
#include "polynomial_ring.hpp"
#include "monomial_ordering.hpp"
#include "algorithm_buchberger_basic.hpp"

extern Monomial_Ordering * generic_grevlex_ptr;

/**
  @ingroup polygroup
  @brief generates the Cyclic-@f$ n @f$ system
  @return a set of generators of the Cyclic-@f$ n @f$ system, as pointers to
      Constant_Polynomial
  @param n number of variables
  @param F ground field
  @param homog whether to homogenize the system
      (affects only the last polynomial)
  @param mord a Monomial_Ordering
  @details Generates the Cyclic-@f$ n @f$ system, @f[
      x_0 + x_1 + ... + x_n,\\
      x_0 x_1 + x_1 x_2 + x_2 x_3 + ... + x_n x_1,\\
      x_0 x_1 x_2 + x_1 x_2 x_3 + x_2 x_3 x_4 + x_n x_1 x_2,\\
      \ldots\\
      x_0 x_1 \ldots x_n - 1.@f]
  Use @f$n > 2@f$.
  If <c>homog</c> is <c>true</c>, the last monomial is @f$h^5@f$ instead of 1.
*/
list<Abstract_Polynomial *> cyclic_n(
    NVAR_TYPE n, Prime_Field & F, bool homog,
    Monomial_Ordering * mord = generic_grevlex_ptr
) {
  list <Abstract_Polynomial *> result;
  // set up coefficients and monomials
  Polynomial_Ring * R = (homog) ? new Polynomial_Ring(n+1, F)
                                : new Polynomial_Ring(n, F);
  NVAR_TYPE max_n = (homog) ? n + 1 : n;
  Prime_Field_Element * A = static_cast<Prime_Field_Element *>(
      malloc(sizeof(Prime_Field_Element) * n)
  );
  Monomial * M = static_cast<Monomial *>(calloc(n, sizeof(Monomial)));
  for (NVAR_TYPE i = 0; i < n; ++i) {
    M[i].common_initialization();
    M[i].initialize_exponents(max_n);
    M[i].set_monomial_ordering(mord);
  }
  // ith polynomial for i = 1, ... n-1
  for (NVAR_TYPE i = 0; i < n - 1; ++i) {
    // jth monomial...
    for (NVAR_TYPE j = 0; j < n; ++j)
    {
      A[j] = F.unity();
      // clear exponents first...
      for (NVAR_TYPE k = 0; k < max_n; ++k)
        M[j].set_exponent(k,0);
      // set relevant exponents to 1
      for (NVAR_TYPE k = j; k < i + j + 1; ++k)
      {
        NVAR_TYPE l = (k >= n) ? k - n : k;
        M[j].set_exponent(l, 1);
      }
    }
    Constant_Polynomial * f = new Constant_Polynomial(n, *R, M, A);
    f->sort_by_order();
    result.push_back(f);
  }
  // last polynomial has a different structure so we can't run it in the loop
  for (NVAR_TYPE i = 0; i < max_n; ++i) {
    if (!homog or i < max_n - 1) {
      M[0].set_exponent(i, 1);
      M[1].set_exponent(i, 0);
    }
    else {
      M[0].set_exponent(i, 0);
      M[1].set_exponent(i, n);
    }
  }
  M[0].set_monomial_ordering(mord);
  M[0].set_monomial_ordering(mord);
  A[0] = F.unity();
  A[1] = -A[0];
  result.push_back(new Constant_Polynomial(2, *R, M, A));
  for (DEG_TYPE i = 0; i < n; ++i) M[i].deinitialize();
  free(M);
  free(A);
  return result;
}

#endif