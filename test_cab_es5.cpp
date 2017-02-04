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

#include <set>
#include <cstdlib>
#include <cstring>
#include <iostream>

using std::set;
using std::cout; using std::endl;

#include "system_constants.hpp"

#include "fields.hpp"
#include "monomial.hpp"
#include "polynomial.hpp"

#include "dynamic_engine.hpp"

#include "algorithm_buchberger_basic.hpp"
#include "algorithm_buchberger_dynamic.hpp"

int main(int argc, char *argv[]) {
  if (argc != 3 or (strcmp(argv[2],"stat") and strcmp(argv[2],"dyn"))) {
    cout << "need to know method (usually 2) and then if dynamic (stat or dyn)\n";
    return 1;
  }
  // obtain method -- don't screw it up b/c we don't check it
  SPolyCreationFlags method = (SPolyCreationFlags )atoi(argv[1]);
  bool static_algorithm = true;
  if (!strcmp(argv[2],"dyn")) static_algorithm = false;
  // set up the field
  Prime_Field FF = Prime_Field(32003);
  string X [9] = {"a", "b", "c", "d", "r", "s", "t", "x", "y"} ;
  Polynomial_Ring R(9, FF, X );
  Prime_Field_Element a = FF.unity();
  DEG_TYPE e0 [] { 0, 0, 0, 0, 0, 0 };
  // set up our polynomials
  // first poly
  Monomial bx { 0, 1, 0, 0, 0, 0, 0, 1, 0 };
  Monomial ay { 1, 0, 0, 0, 0, 0, 0, 0, 1 };
  Monomial M1 [] { bx, ay };
  Prime_Field_Element C1 [] { a, -a };
  Constant_Polynomial f1(2, R, M1, C1);
  f1.sort_by_order();
  // second poly
  Monomial dx { 0, 0, 0, 1, 0, 0, 0, 1, 0 };
  Monomial d1 { 0, 0, 0, 1, 0, 0, 0, 0, 0 };
  Monomial cy { 0, 0, 1, 0, 0, 0, 0, 0, 1 };
  Monomial y1 { 0, 0, 0, 0, 0, 0, 0, 0, 1 };
  Monomial M2 [] { dx, d1, cy, y1 };
  Prime_Field_Element C2 [] { a, -a, -a, a };
  Constant_Polynomial f2(4, R, M2, C2);
  f2.sort_by_order();
  // third poly
  Monomial b2 { 0, 2, 0, 0, 0, 0, 0, 0, 0 };
  Monomial a2 { 2, 0, 0, 0, 0, 0, 0, 0, 0 };
  Monomial r2 { 0, 0, 0, 0, 2, 0, 0, 0, 0 };
  Monomial M3 [] { b2, a2, r2 };
  Prime_Field_Element C3 [] { a, a, -a };
  Constant_Polynomial f3(3, R, M3, C3);
  f3.sort_by_order();
  // fourth poly
  Monomial c2 { 0, 0, 2, 0, 0, 0, 0, 0, 0 };
  Monomial c1 { 0, 0, 1, 0, 0, 0, 0, 0, 0 };
  Monomial one { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  Monomial d2 { 0, 0, 0, 2, 0, 0, 0, 0, 0 };
  Monomial s2 { 0, 0, 0, 0, 0, 2, 0, 0, 0 };
  Monomial M4 [] { c2, c1, one, d2, s2 };
  Prime_Field_Element C4 [] { a, -a*2, a, a, -a };
  Constant_Polynomial f4(5, R, M4, C4);
  f4.sort_by_order();
  // fifth poly
  Monomial ac { 1, 0, 1, 0, 0, 0, 0, 0, 0 };
  Monomial bd { 0, 1, 0, 1, 0, 0, 0, 0, 0 };
  Monomial t2 { 0, 0, 0, 0, 0, 0, 2, 0, 0 };
  Monomial M5 [] { a2, ac, c2, b2, bd, d2, t2 };
  Prime_Field_Element C5 [] { a, -a*2, a, a, -a*2, a, -a };
  Constant_Polynomial f5(7, R, M5, C5);
  f5.sort_by_order();
  // message
  cout << "Computing a Groebner basis for\n\t" << f1
       << ",\n\t" << f2
       << ",\n\t" << f3
       << ",\n\t" << f4
       << ",\n\t" << f5
       << endl;
  // compute basis
  list<Abstract_Polynomial *> F;
  F.push_back(&f1); F.push_back(&f2); F.push_back(&f3);
  F.push_back(&f4); F.push_back(&f5);
  list<Constant_Polynomial *> G;
  if (static_algorithm) G = buchberger(F, method, SUGAR_STRATEGY);
  else G = buchberger_dynamic(
      F, method, SUGAR_STRATEGY, nullptr, DynamicHeuristic::ORD_HILBERT_THEN_DEG
  );
  cout << "Basis:\n";
  for (Constant_Polynomial * g : G) {
    cout << '\t';
    g->leading_monomial().print(true, cout, R.name_list());
    delete g;
  }
  cout << "bye\n";
}
