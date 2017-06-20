#include <set>
#include <iostream>

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

using std::set;
using std::cout; using std::endl;

#include "algorithm_buchberger_basic.hpp"
#include "fields.hpp"
#include "monomial.hpp"
#include "polynomial.hpp"

int main(int argc, char *argv[]) {
  // set up the field
  Prime_Field F43 = Prime_Field(43);
  Polynomial_Ring R(4, F43);
  Prime_Field_Element a = F43.unity();
  // set up our polynomials
  // first poly: linear
  unsigned e1 [] {1,0,0,0}; Monomial t1(4, e1);
  unsigned e2 [] {0,1,0,0}; Monomial t2(4, e2);
  unsigned e3 [] {0,0,1,0}; Monomial t3(4, e3);
  unsigned e4 [] {0,0,0,1}; Monomial t4(4, e4);
  Monomial M1 [] { t1, t2, t3, t4 };
  Prime_Field_Element C1 [] { a, a, a, a };
  Constant_Polynomial f1(4, R, M1, C1);
  // second poly: quadratic
  unsigned e5 [] {1,1,0,0}; Monomial t5(4, e5);
  unsigned e6 [] {0,1,1,0}; Monomial t6(4, e6);
  unsigned e7 [] {0,0,1,1}; Monomial t7(4, e7);
  unsigned e8 [] {1,0,0,1}; Monomial t8(4, e8);
  Monomial M2 [] { t5, t6, t7, t8 };
  Prime_Field_Element C2 [] { a, a, a, a };
  Constant_Polynomial f2(4, R, M2, C2);
  f2.sort_by_order();
  // third poly: cubic
  unsigned e9  [] {1,1,1,0}; Monomial t9 (4, e9);
  unsigned e10 [] {0,1,1,1}; Monomial t10(4, e10);
  unsigned e11 [] {1,0,1,1}; Monomial t11(4, e11);
  unsigned e12 [] {1,1,0,1}; Monomial t12(4, e12);
  Monomial M3 [] { t9, t10, t11, t12 };
  Prime_Field_Element C3 [] { a, a, a, a };
  Constant_Polynomial f3(4, R, M3, C3);
  f3.sort_by_order();
  // fourth poly: quartic
  unsigned e13 [] {1,1,1,1}; Monomial t13(4, e13);
  Monomial t14(4);
  Monomial M4 [] { t13, t14 };
  Prime_Field_Element C4 [] { a, -a };
  Constant_Polynomial f4(2, R, M4, C4);
  f4.sort_by_order();
  // message
  cout << "Computing a Groebner basis for " << f1 << ", " << f2
       << ", " << f3 << ", " << f4 << endl;
  // compute basis
  set<Abstract_Polynomial *> F;
  F.insert(&f1); F.insert(&f2); F.insert(&f3); F.insert(&f4);
  list<Constant_Polynomial *> G = buchberger(F);
  cout << "Basis:\n";
  for (list<Constant_Polynomial *>::iterator g = G.begin(); g != G.end(); ++g)
    cout << '\t' << *(*g) << endl;
  cout << "bye\n";
}
