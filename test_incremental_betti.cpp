#include "betti.hpp"
#include "monomial.hpp"

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

int main() {
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
  for (auto p : incremental_betti(T, x1x32))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << endl;

  for (auto p : incremental_betti(T, x1x22, wts))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  for (auto p : incremental_betti(T, x1x32, wts))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  for (auto p : incremental_betti(T, x33, wts))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << endl;

  WT_TYPE wts2[4] = { 4, 3, 1, 2 };
  list<Monomial> U { w, x2, x1x32, x1x22x3 };
  for (auto p : incremental_betti(U, x1x24))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  for (auto p : incremental_betti(U, x22x33))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  for (auto p : incremental_betti(U, x1))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  for (auto p : incremental_betti(U, x1x24, wts2))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  for (auto p : incremental_betti(U, x22x33, wts2))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  for (auto p : incremental_betti(U, x1, wts2))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
  cout << endl;

  list<Monomial> V { w, x2, x1x32, x1x22x3, x22x33, x1x24 };
  for (auto p : incremental_betti(V, x24x32))
    cout << p.first << ',' << p.second << "; ";
  cout << endl;
}