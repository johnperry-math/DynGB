#ifndef __MONOMIAL_ORDERING_CPP_
#define __MONOMIAL_ORDERING_CPP_

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

#include "monomial_ordering.hpp"

/*
  set_data, first_larger_or_equal, first_larger_or_equal_than_multiple
  defined in monomial.cpp to avoid a compiler error
*/

Monomial_Ordering::~Monomial_Ordering() { }

bool Monomial_Ordering::first_smaller_or_equal(
    const Monomial & t, const Monomial & u
) const {
  return not first_larger(t, u);
}

bool Monomial_Ordering::first_smaller_than_multiple(
    const Monomial &t, const Monomial & u, const Monomial & v
) const {
  return not first_larger_or_equal_than_multiple(t, u, v);
}

bool Monomial_Ordering::first_smaller_or_equal_than_multiple(
    const Monomial & t, const Monomial & u, const Monomial & v
) const {
  return not first_larger_than_multiple(t, u, v);
}

#endif