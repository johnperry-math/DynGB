#include "indeterminate.hpp"

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

Indeterminate & Indeterminate::operator=(const Indeterminate & other) {
  R = other.R; i = other.i;
  return *this;
}

Monomial Indeterminate::operator ^(EXP_TYPE a) {
  Monomial t(R->number_of_variables());
  t.set_exponent(i, a);
  return t;
}

Monomial Indeterminate::operator *(Indeterminate y) {
  Monomial u(R->number_of_variables());
  if (i == y.i)
    u.set_exponent(i, 2);
  else {
    u.set_exponent(i,1);
    u.set_exponent(y.i, 1);
  }
  return u;
}

Monomial Indeterminate::operator *(Monomial t) {
  return t * (*this);
}

ostream & operator << (ostream & os, Indeterminate & x) {
  os << x.R->name(x.i);
  return os;
}
