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

#include "system_constants.hpp"

#include "monomial.hpp"
#include "polynomial.hpp"
#include "algorithm_buchberger_basic.hpp"
#include "algorithm_buchberger_dynamic.hpp"

int main() {
  Prime_Field FF = Prime_Field(32003);
  string X [9] = {"t", "x", "y", "z", "a", "b", "c", "d", "e"} ;
  Polynomial_Ring R(9, FF, X );
  Prime_Field_Element a = FF.unity();
  // set up our polynomials
  // first poly
  Monomial t11 { 0, 32, 0, 32, 0, 0, 0, 0, 0 };
  Monomial t12 { 0, 0, 82, 0, 1, 0, 0, 0, 0 };
  Monomial M1 [] { t11, t12 };
  Prime_Field_Element C1 [] { a, -a };
  Constant_Polynomial f1(2, R, M1, C1);
  f1.sort_by_order();
  // second poly
  Monomial t21 { 0, 45, 0, 0, 0, 0, 0, 0, 0 };
  Monomial t22 { 0, 0, 13, 21, 0, 1, 0, 0, 0 };
  Monomial M2 [] { t21, t22 };
  Prime_Field_Element C2 [] { a, -a };
  Constant_Polynomial f2(2, R, M2, C2);
  f2.sort_by_order();
  // third poly
  Monomial t31 { 0, 41, 0, 0, 0, 0, 1, 0, 0 };
  Monomial t32 { 0, 0, 33, 12, 0, 0, 0, 0, 0 };
  Monomial M3 [] { t31, t32 };
  Prime_Field_Element C3 [] { a, -a };
  Constant_Polynomial f3(2, R, M3, C3);
  f3.sort_by_order();
  // fourth poly
  Monomial t41 { 0, 22, 0, 0, 0, 0, 0, 0, 0 };
  Monomial t42 { 0, 0, 33, 12, 0, 0, 0, 1, 0 };
  Monomial M4 [] { t41, t42 };
  Prime_Field_Element C4 [] { a, -a };
  Constant_Polynomial f4(2, R, M4, C4);
  f4.sort_by_order();
  // fifth poly
  Monomial t51 { 0, 5, 17, 22, 0, 0, 0, 0, 1 };
  Monomial t52 { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  Monomial M5 [] { t51, t52 };
  Prime_Field_Element C5 [] { a, -a };
  Constant_Polynomial f5(2, R, M5, C5);
  f5.sort_by_order();
  // sixth poly
  Monomial t61 { 1, 1, 1, 1, 0, 0, 0, 0, 0 };
  Monomial t62 { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  Monomial M6 [] { t61, t62 };
  Prime_Field_Element C6 [] { a, -a };
  Constant_Polynomial f6(2, R, M6, C6);
  f6.sort_by_order();
  // message
  cout << "Trying an initial analysis on"
       << ":\n\t" << f1
       << ",\n\t" << f2
       << ",\n\t" << f3
       << ",\n\t" << f4
       << ",\n\t" << f5
       << ",\n\t" << f6
       << endl;
  // compute basis
  list<Abstract_Polynomial *> F;
  F.push_back(&f1); F.push_back(&f2); F.push_back(&f3);
  F.push_back(&f4); F.push_back(&f5); F.push_back(&f6);
  Monomial_Ordering * mord;
  Skeleton * skel = new Skeleton(9);
  initial_analysis(F, &mord, skel);
  cout << "Chose leading monomials:\n";
  for (auto f : F) {
    f->set_monomial_ordering(mord);
    cout << f->leading_monomial() << ", ";
  }
  cout << endl;
}