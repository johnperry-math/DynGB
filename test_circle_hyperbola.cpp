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
  Prime_Field_Element a = F43.unity();
  Prime_Field_Element b(-4, &F43);
  Prime_Field_Element c(-1, &F43);
  // set up our polynomials
  // first poly: circle
  unsigned e1 [] {2,0}; Monomial t1(2, e1); // x^2
  unsigned e2 [] {0,2}; Monomial t2(2, e2); // y^2
  unsigned e3 [] {0,0}; Monomial t3(2, e3); // 1
  Monomial M1 [] { t1, t2, t3 };
  Prime_Field_Element C1 [] { a, a, b };
  Constant_Polynomial f1(3, M1, C1);
  // second poly: hyperbola
  unsigned e4 [] {1,1}; Monomial t4(2, e4); // xy
  unsigned e5 [] {0,0}; Monomial t5(2, e5); // 1
  Monomial M2 [] { t4, t5 };
  Prime_Field_Element C2 [] { a, -a };
  Constant_Polynomial f2(2, M2, C2);
  // message
  cout << "Computing a Groebner basis for " << f1 << ", " << f2 << endl;
  // compute basis
  set<Abstract_Polynomial *> F;
  F.insert(&f1); F.insert(&f2);
  set<Constant_Polynomial *> G = buchberger(F);
  cout << "Basis:\n";
  for (auto g = G.begin(); g != G.end(); ++g)
    cout << '\t' << *(*g) << endl;
  cout << "bye\n";
}
