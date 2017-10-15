#ifndef __BETTI_CPP_
#define __BETTI_CPP_

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
using std::set;
#include <utility>
using std::pair;
#include <algorithm>
using std::any_of;

#include "betti.hpp"

map<DEG_TYPE, unsigned long> incremental_betti(
  const list<Monomial> & T, const WT_TYPE * grading
) {
  map<DEG_TYPE, unsigned long> result;
  const Monomial & u = T.back();
  auto Tstop = T.end();
  --Tstop;
  list<pair<Monomial, Monomial>> S;
  for (auto ti = T.begin(); ti != Tstop; ++ti) {
    S.emplace_back(*ti, u);
    S.back().first.set_monomial_ordering(u.monomial_ordering());
    for (auto ui = T.begin(); ui != ti; ++ui) {
      S.emplace_back(*ti, *ui);
      S.back().first.set_monomial_ordering(u.monomial_ordering());
    }
  }
  list<pair<Monomial, Monomial>> S2;
  for (auto p : S) {
    bool found = false;
    if (
        not found and u | p.first.lcm(p.second) and
        any_of(S2.begin(), S2.end(), [p,u](pair<Monomial, Monomial> q){return (q.first == u or q.second == u) and (q.first == p.first or q.second==p.first);}) and
        any_of(S2.begin(), S2.end(), [p,u](pair<Monomial, Monomial> q){return (q.first == u or q.second == u) and (q.first == p.second or q.second==p.second);})
    ) {
      found = true;
    }
    for (auto t : T) {
      if (found) break;
      if (
          t | p.first.lcm(p.second) and
          any_of(S2.begin(), S2.end(), [p,t](pair<Monomial, Monomial> q){return (q.first == t or q.second == t) and (q.first == p.second or q.second==p.second);}) and
          any_of(S2.begin(), S2.end(), [p,t](pair<Monomial, Monomial> q){return (q.first == t or q.second == t) and (q.first == p.second or q.second==p.second);})
      ) {
        found = true;
      }
    }
    if (not found) S2.emplace_back(p);
  }
  for (auto p : S2) {
    DEG_TYPE d = (grading == nullptr)
        ? p.first.lcm(p.second).total_degree()
        : p.first.lcm(p.second).weighted_degree(grading);
    if (result.find(d) != result.end())
      ++result[d];
    else
      result[d] = 1;
  }
  return result;
}

#endif