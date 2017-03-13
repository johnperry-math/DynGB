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
  string X [9] = { "t", "u", "v", "w", "x", "y", "z" } ;
  Polynomial_Ring R(7, FF, X );
  Prime_Field_Element a = FF.unity();
  Monomial one { 0, 0, 0, 0, 0, 0, 0 };
  // set up our polynomials
  // first poly
  Monomial x2y1z4 { 0, 0, 0, 0, 2, 1, 4 };
  Monomial t1 { 1, 0, 0, 0, 0, 0, 0 };
  Monomial M1 [] { x2y1z4, t1 };
  Prime_Field_Element C1 [] { -a, a };
  Constant_Polynomial f1(2, R, M1, C1);
  f1.sort_by_order();
  // second poly
  Monomial x5y7 { 0, 0, 0, 0, 5, 7, 0 };
  Monomial u1z2 { 0, 1, 0, 0, 0, 0, 2 };
  Monomial M2 [] { x5y7, u1z2 };
  Prime_Field_Element C2 [] { -a, a };
  Constant_Polynomial f2(2, R, M2, C2);
  f2.sort_by_order();
  // third poly
  Monomial v1x3z1 { 0, 0, 1, 0, 3, 0, 1 };
  Monomial y2 { 0, 0, 0, 0, 0, 2, 0 };
  Monomial M3 [] { v1x3z1, y2 };
  Prime_Field_Element C3 [] { a, -a };
  Constant_Polynomial f3(2, R, M3, C3);
  f3.sort_by_order();
  // fourth poly
  Monomial w1z5 { 0, 0, 0, 1, 0, 0, 5 };
  Monomial x1y3 { 0, 0, 0, 0, 1, 3, 0 };
  Monomial M4 [] { w1z5, x1y3 };
  Prime_Field_Element C4 [] { a, -a };
  Constant_Polynomial f4(2, R, M4, C4);
  f4.sort_by_order();
  // message
  cout << "Computing a Groebner basis for\n\t" << f1
       << ",\n\t" << f2
       << ",\n\t" << f3
       << ",\n\t" << f4
       << endl;
  // compute basis
  list<Abstract_Polynomial *> F;
  F.push_back(&f1); F.push_back(&f2); F.push_back(&f3);
  F.push_back(&f4);
  list<Constant_Polynomial *> G;
  if (static_algorithm) G = buchberger(F, method, StrategyFlags::SUGAR_STRATEGY);
  else G = buchberger_dynamic(
      F, method, StrategyFlags::SUGAR_STRATEGY, nullptr,
      DynamicHeuristic::ORD_HILBERT_THEN_DEG
  );
  cout << "Basis:\n";
  for (Constant_Polynomial * g : G) {
    cout << '\t';
    g->leading_monomial().print(true, cout, R.name_list());
    delete g;
  }
  cout << "bye\n";
}
