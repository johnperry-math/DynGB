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

#include "hilbert_functions.hpp"

// prototypes
void test_greuel_pfister();
void test_splitting_case();

int main() {

  test_splitting_case();
  test_greuel_pfister();

  cout << endl;  
}

void test_splitting_case() {

  Prime_Field F { 43 };
  auto R = Polynomial_Ring(10, F);
  Indeterminate a(R, 0), b(R, 1), c(R, 2), d(R, 3), e(R, 4), f(R, 5), g(R, 6),
    h(R, 7), i(R, 8), j(R, 9);

  list<Monomial> B = { a*b, c*d*e, f*g, h*i*j, j^2 };

  cout << "Hilbert numerator for { ";
  for ( auto & t : B ) { cout << t << " , "; }
  cout << " } :\n";
  auto hn =  hilbert_numerator_bigatti(B);
  cout << *hn << endl;
  delete hn;

}

void test_greuel_pfister() {

  Prime_Field F{ 32003 };
  auto R = Polynomial_Ring(5, F);
  Indeterminate t(R, 0), x(R, 1), y(R, 2), z(R, 3), w(R, 4);

  list<Monomial> B = {
    (x^14)*(w^12),
    (x^12)*z*(w^11),
    (t^2)*(x^3)*(y^7)*(w^11),
    t*(x^5)*(y^5)*(w^11),
    (t^2)*(x^4)*(y^4)*(w^11),
    (x^7)*(z^2)*(w^10),
    (x^5)*(y^2)*z*(w^10),
    (x^4)*(y^3)*z*(w^10),
    (x^5)*(z^3)*(w^10),
    (x^6)*(y^2)*z*(w^8),
    (t^3)*x*(y^4)*(w^9),
    (x^4)*y*(z^3)*(w^9),
    (t^4)*x*(y^2)*(w^10),
    (x^3)*(z^4)*(w^10),
    (t^3)*x*(y^5)*(w^7),
    (t^5)*(x^3)*(w^8),
    (t^3)*(x^4)*(w^9),
    (t^2)*(x^5)*(w^9),
    t*(x^6)*(w^9),
    (x^2)*y*(z^4)*(w^9),
    (t^7)*(y^3)*(w^5),
    (t^4)*(x^4)*(w^7),
    (t^3)*(x^5)*(w^7),
    (t^2)*(x^6)*(w^7),
    (t^5)*x*(y^2)*(w^7),
    (t^4)*x*(y^3)*(w^7),
    t*(z^7)*(w^7),
    (t^7)*(x^3)*(w^4),
    (t^6)*(x^3)*(w^5),
    (t^7)*x*y*(w^5),
    (t^6)*x*(y^2)*(w^5),
    (t^5)*x*(y^3)*(w^5),
    t*y*(z^6)*(w^6),
    (t^2)*y*(z^7)*(w^3),
    (t^5)*(x^4)*(w^4),
    (t^2)*y*(z^6)*(w^4),
    (t^2)*(z^7)*(w^4),
    (t^7)*(y^3)*z*w,
    (t^5)*(y^3)*z*(w^3),
    (t^7)*(x^2)*y*w,
    (t^7)*x*(y^2)*w,
    (t^3)*(x^2)*y*(w^5),
    (t^3)*y*(z^6),
    (t^2)*(y^3)*z*(w^4),
    t*(x^2)*z*(w^4),
    (y^2)*(z^2)*(w^3),
    (t^2)*x*y*z*w,
    t*(y^2)*(z^2),
    t*x*(z^2)
  };
  
  cout << "For { ";
  for (auto & t : B) cout << t << ", ";
  cout << " } the numerator is\n";
  auto h = hilbert_numerator_bigatti(B);
  cout << *h;

  delete h;

}