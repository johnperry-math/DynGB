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
  string X [6] = {"B", "P", "S", "T", "W", "Z"} ;
  Polynomial_Ring R(6, FF, X );
  Prime_Field_Element a = FF.unity();
  DEG_TYPE e0 [] { 0, 0, 0, 0, 0, 0 };
  // set up our polynomials
  // first poly
  Monomial P { 0, 1, 0, 0, 0, 0 };
  Monomial S { 0, 0, 1, 0, 0, 0 };
  Monomial B { 1, 0, 0, 0, 0, 0 };
  Monomial one { 0, 0, 0 ,0, 0, 0 };
  Monomial M1 [] { P, S, B, one };
  Prime_Field_Element C1 [] { a*45, a*35, -a*165, -a*36 };
  Constant_Polynomial f1(4, R, M1, C1);
  f1.sort_by_order();
  // second poly
  Monomial Z { 0, 0, 0, 0, 0, 1 } ;
  Monomial T { 0, 0, 0, 1, 0, 0 };
  Monomial M2 [] { P, Z, T, S };
  Prime_Field_Element C2 [] { a*35, a*40, a*25, -a*27 };
  Constant_Polynomial f2(4, R, M2, C2);
  f2.sort_by_order();
  // third poly
  Monomial W  { 0, 0, 0, 0, 1, 0 };
  Monomial PS { 0, 1, 1, 0, 0, 0 };
  Monomial B2 { 2, 0, 0, 0, 0, 0 };
  Monomial M3 [] { W, PS, Z, T, B2 };
  Prime_Field_Element C3 [] { a*15, a*25, a*30, -a*18, -a*165 };
  Constant_Polynomial f3(5, R, M3, C3);
  f3.sort_by_order();
  // fourth poly
  Monomial PT { 0, 1, 0, 1, 0, 0 };
  Monomial ZS { 0, 0, 1, 0, 0, 1 };
  Monomial M4 [] { W, PT, ZS };
  Prime_Field_Element C4 [] { -a*9, a*15, a*20 };
  Constant_Polynomial f4(3, R, M4, C4);
  f4.sort_by_order();
  // fifth poly
  Monomial WP { 0, 1, 0, 0, 1, 0 };
  Monomial ZT { 0, 0, 0, 1, 0, 1 };
  Monomial B3 { 3, 0, 0, 0, 0, 0 };
  Monomial M5 [] { WP, ZT, B3 };
  Prime_Field_Element C5 [] { a, a*2, -a*11 };
  Constant_Polynomial f5(3, R, M5, C5);
  f5.sort_by_order();
  // sixth poly
  Monomial SB { 1, 0, 1, 0, 0, 0 };
  Monomial M6 [] { W, SB, B2 };
  Prime_Field_Element C6 [] { a*99, -a*11, a*3 };
  Constant_Polynomial f6(3, R, M6, C6);
  f6.sort_by_order();
  // message
  cout << "Computing a Groebner basis for\n\t" << f1
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
  list<Constant_Polynomial *> G;
  if (static_algorithm) G = buchberger(F, method, StrategyFlags::SUGAR_STRATEGY);
  else G = buchberger_dynamic(
      F, method, StrategyFlags::SUGAR_STRATEGY,
      nullptr, Dynamic_Heuristic::ORD_HILBERT_THEN_DEG
  );
  cout << "Basis:\n";
  for (Constant_Polynomial * g : G) {
    cout << '\t';
    g->leading_monomial().print(true, cout, R.name_list());
    delete g;
  }
  cout << "bye\n";
}
