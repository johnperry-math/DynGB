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
using std::pair; using std::tuple; using std::get;
#include <algorithm>
using std::any_of; using std::min; using std::max; using std::sort;

#include "betti.hpp"

WT_TYPE lcm_degree(const Monomial &t, const Monomial &u, const WT_TYPE * grading = nullptr) {
  auto a = t.packed_log(), b = u.packed_log();
  NVAR_TYPE i = 0, j = 0;
  WT_TYPE result = 0;
  if (grading != nullptr) {
    for (/* */; i < t.packed_size() and j < u.packed_size(); /* */) {
      if (a[i] == b[j]) {
        result += max(a[i+1], b[j+1]) * grading[a[i]];
        i += 2; j += 2;
      } else if (a[i] < b[j]) {
        result += a[i+1] * grading[a[i]];
        i += 2;
      } else {
        result += b[j+1] * grading[b[j]];
        j += 2;
      }
    }
    for (/* */; i < t.packed_size(); i += 2)
      result += a[i+1] * grading[a[i]];
    for (/* */; j < u.packed_size(); j += 2)
      result += b[j+1] * grading[b[j]];
  } else {
    for (/* */; i < t.packed_size() and j < u.packed_size(); /* */) {
      if (a[i] == b[j]) {
        result += max(a[i+1], b[j+1]);
        i += 2; j += 2;
      } else if (a[i] < b[j]) {
        result += a[i+1];
        i += 2;
      } else {
        result += b[j+1];
        j += 2;
      }
    }
    for (/* */; i < t.packed_size(); i += 2)
      result += a[i+1];
    for (/* */; j < u.packed_size(); j += 2)
      result += b[j+1];
  }
  return result;
}

struct Comparer {
  const vector<const Monomial *> * T;
  bool operator()(const pair<int, int> &p, const pair<int, int> &q) {
    return lcm_degree(*((*T)[p.first]), *((*T)[p.second])) <
        lcm_degree(*((*T)[q.first]), *((*T)[q.second])) ;
  }
};

map<DEG_TYPE, unsigned long> incremental_betti(
  const list<Monomial> & T,
  const Monomial & u,
  set< pair<const Monomial *, const Monomial *> > & R,
  const WT_TYPE * grading
) {
  map<DEG_TYPE, unsigned long> result;
  // cancel old triples
  auto ri = R.begin();
  while (ri != R.end()) {
    auto & r = *ri;
    auto rlcm{r.first->lcm(*r.second)};
    if (rlcm == rlcm.lcm(u) and rlcm != r.first->lcm(u) and rlcm != r.second->lcm(u)) {
      auto tmp {ri};
      ++ri;
      R.erase(tmp);
    } else ++ri;
  }
  // cancel new pairs subject to proper divisibility of lcms
  set<tuple<const Monomial *, const Monomial *, const Monomial> > A;
  for (auto & t : T) { A.emplace(&t, &u, t.lcm(u)); }
  auto ai = A.begin();
  while (ai != A.end()) {
    if (any_of(
          A.begin(), A.end(),
          [&ai](auto bi){
              return *ai != bi and get<2>(bi) != get<2>(*ai) and get<2>(bi) | get<2>(*ai);
          } 
    )) {
      auto tmp{ai}; ++ai;
      A.erase(tmp);
    } else ++ai;
  }
  // retain only one pair at each lcm
  ai = A.begin();
  while (ai != A.end()) {
    auto ai2 = ai; ++ai2;
    auto bi { ai2 };
    while (bi != A.end()) {
      if (get<2>(*ai) == get<2>(*bi)) {
        if (get<0>(*bi)->is_coprime(*get<1>(*bi))) {
          A.erase(ai);
          ai = bi;
          ++bi;
        } else {
          auto tmp { bi };
          ++tmp;
          if (ai2 == bi) ai2 = tmp;
          A.erase(bi);
          bi = tmp;
        }
      } else ++bi;
    }
    if (get<0>(*ai)->is_coprime(*get<1>(*ai)))
      A.erase(ai);
    else
      R.emplace(get<0>(*ai), get<1>(*ai));
    ai = ai2;
  }
  for (auto & r : R) {
    //cout << *r.first << ", " << *r.second << "; ";
    DEG_TYPE d = lcm_degree(*(r.first), *(r.second), grading);
    if (result.find(d) != result.end())
      ++result[d];
    else
      result[d] = 1;
  }
  return result;
}


map<DEG_TYPE, unsigned long> full_betti(
  const list<Monomial> & U, const WT_TYPE * grading
) {
  map<DEG_TYPE, unsigned long> result;
  vector<const Monomial *> T(U.size());
  int i = 0;
  for (auto & u : U) { T[i] = &u; ++i; }
  // first add all pairs
  size_t T_size = T.size();
  vector< pair<int, int> > P(T.size() * (T.size() - 1) / 2);
  unsigned k = 0;
  for (int i = 0; i < T.size(); ++i) {
    for (int j = i + 1; j < T.size(); ++j) {
      P[k].first = i; P[k].second = j;
      ++k;
    }
  }
  // now we remove relatively prime pairs,
  // as well as those that they divide
  vector< pair<int, int> > Q(P.size());
  k = 0;
  auto pi = P.begin();
  while (pi != P.end()) {
    if (pi->first != -1) {
      const Monomial & t = *T[pi->first];
      if (not t.is_coprime(*T[pi->second])) {
        Q[k] = *pi;
        ++k;
      } else {
        //cout << "dropping " << *T[pi->first] << ", " << *T[pi->second] << endl;
        const Monomial & u = *T[pi->second];
        auto tu { t.lcm(u) };
        auto pi2 = pi; ++pi2;
        while (pi2 != P.end()) {
          if (pi2->first != -1 and (not tu.divides_lcm(*T[pi2->first], *T[pi2->second])))
            ++pi2;
          else {
            pi2->first = -1;
            ++pi2;
          }
        }
      }
    }
    ++pi;
  }
  // we are now left only with non-relatively prime pairs
  // remove pairs that will be computed after others whose lcm divides theirs
  set< pair<int, int> > R;
  pair<int, int> first, second;
  Comparer compare;
  compare.T = &T;
  sort(Q.begin(), Q.begin() + k, compare);
  for (auto qi = Q.begin(); qi < Q.begin() + k; ++qi) {
    auto & q = *qi;
    bool any_divides = false;
    for (int i = 0; (not any_divides) and i < T.size(); ++i) {
      first.first = min(i, q.first); first.second = max(i, q.first);
      second.first = min(i, q.second); second.second = max(i, q.second);
      if (
          T[i]->divides_lcm(*T[q.first], *T[q.second]) and
          R.find(first) != R.end() and R.find(second) != R.end()
      ) {
        any_divides = true;
        //cout << "erasing " << *T[q.first] << ", " << *T[q.second] << " because of " << *T[i] << endl;
      }
    }
    if (not any_divides)
      R.emplace(q);
  }
  for (auto & r : R) {
    //cout << *T[r.first] << ", " << *T[r.second] << "; ";
    DEG_TYPE d = lcm_degree(*T[r.first], *T[r.second], grading);
    if (result.find(d) != result.end())
      ++result[d];
    else
      result[d] = 1;
  }
  //cout << result.size() << "; ";
  return result;
}

#endif