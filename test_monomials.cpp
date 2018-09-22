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
using std::cout; using std::endl;

#include "fields.hpp"
#include "monomial.hpp"

#include "fields.hpp"
#include "polynomial_ring.hpp"
#include "indeterminate.hpp"

int main()
{
  DEG_TYPE a [] = {1, 0, 5, 2, 0};
  Monomial s { 1, 0, 5, 2, 0 };
  Monomial t { 1, 0, 3, 2, 0 };
  cout << "s = " << s << endl;
  cout << "t = " << t << endl;

  Monomial u { 3, 1, 0, 0, 2 };
  cout << "u = " << u << endl << endl;

  cout << "\nMasks:\n";
  cout << s.mask().to_string() << endl << t.mask() << endl << u.mask() << endl << endl;

  cout << "Does s divide t? " << t.divisible_by(s) << endl;
  cout << "Does t divide s? " << s.divisible_by(t) << endl;
  cout << "Does t divide u? " << u.divisible_by(t) << endl;
  cout << "Does u divide t? " << t.divisible_by(u) << endl;
  cout << "t:u = " << t.colon(u) << endl;
  cout << "u:t = " << u.colon(t) << endl;
  cout << "lcm(t,u) = " << t.lcm(u) << endl;
  cout << "lcm(u,t) = " << u.lcm(t) << endl;
  cout << "gcd(t,u) = " << t.gcd(u) << endl;
  cout << "gcd(u,t) = " << u.gcd(t) << endl;

  Monomial w(t);
  w *= u;
  cout << "tu = " << w << endl;

  Monomial v { 2, 2, 0, 1, 0 };
  cout << endl << "v = " << v << endl;

  Monomial tv(t.lcm(v));
  tv /= t;
  cout << "lcm(t,v)/t = " << tv << endl;
  Monomial larger { 0, 0, 2, 0, 1, 0 };
  Monomial smaller { 0, 1, 0, 1, 1, 0 };
  cout << smaller << " < " << larger << " ? " << (smaller < larger) << endl;
  cout << larger << " < " << smaller << " ? " << (larger < smaller) << endl;

  Prime_Field F(43);
  string var_names [] = { "x", "y" };
  Polynomial_Ring R = Polynomial_Ring(2, F, var_names);
  Indeterminate * X = R.indeterminates();
  Monomial x2 = X[0]^2;
  x2.print(true, cout, var_names);
  Monomial xy = X[0]*X[1];
  xy.print(true, cout, var_names);
  Monomial x3y_first = x2 * xy;
  Monomial x3y_second = (X[0]^2) * xy;
  x3y_first.print(false, cout, var_names); cout << ',';
  x3y_second.print(true, cout, var_names);
  free(X);
}