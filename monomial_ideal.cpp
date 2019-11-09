#ifndef __MONOMIAL_IDEAL_CPP_
#define __MONOMIAL_IDEAL_CPP_

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

#include "monomial_ideal.hpp"

ostream & operator << (ostream & os, const Monomial_Ideal & I) {
  os << "[ ";
  for (
        list<Monomial>::const_iterator ti = I.gens.begin();
        ti != I.gens.end();
        ++ti
  ) {
    os << *ti << ',';
  }
  os << " ] (";
  if (I.hNum == nullptr) os << I.hNum << ", ";
  else os << *(I.hNum) << ", ";
  if (I.hRedNum == nullptr) os << I.hRedNum << ", ";
  else os << *(I.hRedNum) << ", ";
  if (I.hPol == nullptr) os << I.hPol << ")";
  else os << *(I.hPol) << ")";
  return os;
}

bool first_smaller(const Monomial & t, const Monomial & u) {
  return t < u;
}

bool standard_degree(const Monomial & t, const Monomial & u) {
  return t.total_degree() < u.total_degree();
}

list<Monomial> colon_ideal_without_ideals(
    const list<Monomial> & U, const Monomial & t
) {
  list<Monomial> V;
  for (const Monomial & u : U) {
    Monomial tnew = u.colon(t);
    V.push_back(tnew);
  }
  V.sort(standard_degree);
  for (auto vi = V.begin(); vi != V.end(); ++vi) {
    auto ui { vi };
    for (++ui; ui != V.end(); /* */) {
      if (ui->divisible_by(*vi)) {
        auto wi { ui };
        ++wi;
        V.erase(ui);
        ui = wi;
      } else
        ++ui;
    }
  }
  return V;
}

#endif