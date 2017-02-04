#ifndef __STRATEGIES_CPP_
#define __STRATEGIES_CPP_

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

#include "strategies.hpp"

bool Poly_Strategy_Data::operator==(const Poly_Strategy_Data &other) const {
  return equivalent(other);
}

bool Poly_Strategy_Data::operator >(const Poly_Strategy_Data &other) const {
  return first_larger(other);
}

bool Poly_Strategy_Data::operator >=(const Poly_Strategy_Data &other) const {
  return *this == other or *this > other;
}

bool Poly_Strategy_Data::operator <(const Poly_Strategy_Data &other) const {
  return not (*this >= other);
}

bool Poly_Strategy_Data::operator <=(const Poly_Strategy_Data &other) const {
  return not first_larger(other);
}

bool Poly_Strategy_Data::valid_reduction(
    const Abstract_Polynomial & r, const Abstract_Polynomial & g
) const { return true; }

ostream & operator <<(ostream &os, const Poly_Strategy_Data &sd) {
  os << "generic Poly_Strategy_Data";
  return os;
}

Pair_Strategy_Data::~Pair_Strategy_Data() { };

bool Pair_Strategy_Data::operator ==(const Pair_Strategy_Data & sd) const {
  return equivalent(sd);
}

bool Pair_Strategy_Data::operator >(const Pair_Strategy_Data & sd) const {
  return first_larger(sd);
}

bool Pair_Strategy_Data::operator >=(const Pair_Strategy_Data & sd) const {
  return *this > sd or *this == sd;
}

bool Pair_Strategy_Data::operator <(const Pair_Strategy_Data & sd) const {
  return not (*this >= sd);
}

bool Pair_Strategy_Data::operator <=(const Pair_Strategy_Data & sd) const {
  return not first_larger(sd);
}

void Pair_Strategy_Data::pre_spolynomial_tasks() const { }

#endif