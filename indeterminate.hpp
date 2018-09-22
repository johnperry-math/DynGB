#ifndef __INDETERMINATE_HPP_
#define __INDETERMINATE_HPP_

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

#include <iostream>

#include "system_constants.hpp"

#include "monomial.hpp"
#include "polynomial_ring.hpp"

class Monomial;
class Polynomial_Ring;

/**
  @class Indeterminate
  @author John Perry
  @ingroup polygroup
  @date 2016
  @brief Implementation of indeterminates, for easier building of polynomials.

  @details The main purpose of this class is to help make it easier to build
  polynomials. It is not especially useful otherwise, and is highly inefficient.
  However, it makes the following possible:
  @code
  Prime_Field F = Prime_Field(43);
  string var_names [] = { "x", "y" };
  P = Polynomial_Ring(2, F, var_names);
  Indeterminate x(R, 0);
  Indeterminate y(R, 1);
  Monomial x3y = (x^3) * (x*y);
  x3y.printlncout();
  @endcode
  &hellip;and the result should be @f$x^3y@f$. Alternately, you could do:
  @code
  Indeterminate * X = P.indeterminates();
  Monomial x3y = (X[0]^3) * (X[0]*X[1]);
  x3y.printlncout();
  free(X);
  @endcode
  &hellip;with the same result. Just be careful in the second case to destroy
  the evidence.
  @example test_monomials.cpp

  @warning Keep in mind that the constructor&rsquo;s index should correspond
    to a valid number; that is, <c>xi < P.number_of_variables()</c>!
*/
class Indeterminate {

public:

  /** @name Construction */
  ///@{

  /** @brief @c this will correspond to the <c>xi</c>th indeterminate of @c P. */
  Indeterminate(Polynomial_Ring & P, NVAR_TYPE xi) : R(&P), i(xi) { }

  /** @brief copy constructor */
  Indeterminate(const Indeterminate & other) : R(other.R), i(other.i) { }

  Indeterminate & operator =(const Indeterminate & other);

  ///@}

  /** @name Basic properties */
  ///@{

  /** @brief the Polynomial_Ring @c this lives in */
  Polynomial_Ring & base_ring() const { return *R; }

  /** @brief which variable in base_ring() @c this is*/
  NVAR_TYPE index_in_ring() const { return i; }

  ///@}

  /** @name Computation */
  ///@{

  /** @brief returns @c this to the <c>a</c>th power */
  Monomial operator ^(EXP_TYPE a);

  /** @brief returns the product of @c this and @p y */
  Monomial operator *(Indeterminate y);

  /** @brief returns the product of @c this and @p y */
  Monomial operator *(Monomial t);

  ///@}

  /** @name I/O */
  ///@{

  /** @brief prints @c this with the appropriate name */
  friend ostream & operator << (ostream &, Indeterminate &);

  ///@}

protected:

  /** @brief the ring @c this lives in */
  Polynomial_Ring * R;
  /** @brief which indeterminate in @c R @c this is */
  NVAR_TYPE i;
};

#endif