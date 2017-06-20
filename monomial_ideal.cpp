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

list<Monomial> colon_ideal_without_ideals(
    const list<Monomial> & U, const Monomial & t
) {
  list<Monomial> V;
  for (const Monomial & u : U) {
    Monomial tnew = u.colon(t);
    // first prune from V monomials that tnew divides
    bool redundant = false;
    for (
          list<Monomial>::const_iterator vi = V.begin();
          vi != V.end();
          ++vi
    ) {
      if ((not redundant) and tnew.divisible_by(*vi))
        redundant = true; 
      while (
              vi != V.end() and (not tnew.is_like(*vi))
                  and vi->divisible_by(tnew)
      ) {
        list<Monomial>::const_iterator wi = vi; ++wi;
        V.erase(vi);
        vi = wi;
      }
    }
    // check that tnew not divisible by anything in V
    if (not redundant)
      V.push_back(tnew);
  }
  return V;
}

#endif