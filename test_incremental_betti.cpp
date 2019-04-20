#include "betti.hpp"
#include "monomial.hpp"

#include <set>
using std::set;

#include <utility>
using std::pair;

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

void full_test() {
  Monomial w       { 1, 0, 0, 0 };

  Monomial x2      { 0, 2, 0, 0 };

  Monomial x1x22   { 0, 1, 2, 0 };
  Monomial x1x32   { 0, 1, 0, 2 };
  Monomial x33     { 0, 0, 0, 3 };

  Monomial x1x22x3 { 0, 1, 2, 1 };

  Monomial x1x24   { 0, 1, 4, 0 };
  Monomial x22x33  { 0, 0, 2, 3 };
  Monomial x1      { 0, 1, 0, 0 };

  Monomial x24x32  { 0, 0, 4, 2 };

  WT_TYPE wts[4] = { 6, 4, 3, 2 };
  list<Monomial> T { w, x2 };
  cout << w << ", " << x2 << ", " << x1x32 << ": ";
  for (auto p : full_betti(T, x1x32))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << "Should be (4,1)\n";
  cout << endl;

  cout << w << ", " << x2 << ", " << x1x22 << " (6 4 3 2) " << " : ";
  for (auto p : full_betti(T, x1x22, wts))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << "Should be (14,1)\n";
  cout << w << ", " << x2 << ", " << x1x32 << " (6 4 3 2) " << " : ";
  for (auto p : full_betti(T, x1x32, wts))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << "Should be (12,1)\n";
  cout << w << ", " << x2 << ", " << x33 << " (6 4 3 2) " << " : ";
  for (auto p : full_betti(T, x33, wts))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << "Should show nothing\n";
  cout << endl;

  WT_TYPE wts2[4] = { 4, 3, 1, 2 };
  list<Monomial> U { w, x2, x1x32, x1x22x3 };
  cout << w << ", " << x2 << ", " << x1x32 << ", " << x1x22x3 << ", " << x1x24 << " : ";
  for (auto p : full_betti(U, x1x24))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << w << ", " << x2 << ", " << x1x32 << ", " << x1x22x3 << ", " << x22x33 << " : ";
  for (auto p : full_betti(U, x22x33))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << w << ", " << x2 << ", " << x1x32 << ", " << x1x22x3 << ", " << x1 << " : ";
  for (auto p : full_betti(U, x1))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << w << ", " << x2 << ", " << x1x32 << ", " << x1x22x3 << ", " << x1x24 << " (4 3 1 2) " << " : ";
  for (auto p : full_betti(U, x1x24, wts2))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << w << ", " << x2 << ", " << x1x32 << ", " << x1x22x3 << ", " << x22x33 << " (4 3 1 2) " << " : ";
  for (auto p : full_betti(U, x22x33, wts2))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << w << ", " << x2 << ", " << x1x32 << ", " << x1x22x3 << ", " << x1 << " (4 3 1 2) " << " : ";
  for (auto p : full_betti(U, x1, wts2))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << endl;

  list<Monomial> V { w, x2, x1x32, x1x22x3, x22x33, x1x24 };
  cout << w << ", " << x2 << ", " << x1x32 << ", " << x1x22x3 << ", " << x22x33 << ", " << x1x24 << ", " << x24x32 << " : ";
  for (auto p : full_betti(V, x24x32))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << "Should be (4,1), (5,2), (6,3), (7,3)\n";
}

void incremental_test() {
  Monomial w       { 1, 0, 0, 0 };

  Monomial x2      { 0, 2, 0, 0 };

  Monomial x1x22   { 0, 1, 2, 0 };
  Monomial x1x32   { 0, 1, 0, 2 };
  Monomial x33     { 0, 0, 0, 3 };

  Monomial x1x22x3 { 0, 1, 2, 1 };

  Monomial x1x24   { 0, 1, 4, 0 };
  Monomial x22x33  { 0, 0, 2, 3 };
  Monomial x1      { 0, 1, 0, 0 };

  Monomial x24x32  { 0, 0, 4, 2 };

  WT_TYPE wts[4] = { 6, 4, 3, 2 };
  list<Monomial> T { w, x2 };
  cout << w << ", " << x2 << ", " << x1x32 << ": ";
  set<pair<const Monomial *, const Monomial *> > P;
  for (auto p : incremental_betti(T, x1x32, P))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << "Should be (4,1)\n";
  cout << endl;

  P.clear();
  cout << w << ", " << x2 << ", " << x1x22 << " (6 4 3 2) " << " : ";
  for (auto p : incremental_betti(T, x1x22, P, wts))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << "Should be (14,1)\n";
  P.clear();
  cout << w << ", " << x2 << ", " << x1x32 << " (6 4 3 2) " << " : ";
  for (auto p : incremental_betti(T, x1x32, P, wts))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  P.clear();
  cout << "Should be (12,1)\n";
  cout << w << ", " << x2 << ", " << x33 << " (6 4 3 2) " << " : ";
  for (auto p : incremental_betti(T, x33, P, wts))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << "Should show nothing\n";
  cout << endl;

  /*WT_TYPE wts2[4] = { 4, 3, 1, 2 };
  list<Monomial> U { w, x2, x1x32, x1x22x3 };
  cout << w << ", " << x2 << ", " << x1x32 << ", " << x1x22x3 << ", " << x1x24 << " : ";
  for (auto p : full_betti(U, x1x24))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << w << ", " << x2 << ", " << x1x32 << ", " << x1x22x3 << ", " << x22x33 << " : ";
  for (auto p : full_betti(U, x22x33))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << w << ", " << x2 << ", " << x1x32 << ", " << x1x22x3 << ", " << x1 << " : ";
  for (auto p : full_betti(U, x1))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << w << ", " << x2 << ", " << x1x32 << ", " << x1x22x3 << ", " << x1x24 << " (4 3 1 2) " << " : ";
  for (auto p : full_betti(U, x1x24, wts2))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << w << ", " << x2 << ", " << x1x32 << ", " << x1x22x3 << ", " << x22x33 << " (4 3 1 2) " << " : ";
  for (auto p : full_betti(U, x22x33, wts2))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << w << ", " << x2 << ", " << x1x32 << ", " << x1x22x3 << ", " << x1 << " (4 3 1 2) " << " : ";
  for (auto p : full_betti(U, x1, wts2))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << endl;*/

  P.clear();
  P.emplace(&x2, &x1x32);   // 1,2
  P.emplace(&x2, &x1x22x3); // 1,3
  P.emplace(&x2, &x1x24);   // 1,5
  P.emplace(&x1x32, &x1x22x3); // 2,3
  P.emplace(&x1x32, &x22x33);  // 2,4
  P.emplace(&x1x32, &x1x24);   // 2,5
  P.emplace(&x1x22x3, &x1x24); // 3,5
  list<Monomial> V { w, x2, x1x32, x1x22x3, x22x33, x1x24 };
  cout << w << ", " << x2 << ", " << x1x32 << ", " << x1x22x3 << ", " << x22x33 << ", " << x1x24 << ", " << x24x32 << " : ";
  for (auto p : incremental_betti(V, x24x32, P))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << "Should be (4,1), (5,2), (6,3), (7,3)\n"; 
}

int main() {
  full_test();
  cout << endl << endl;
  incremental_test();
}