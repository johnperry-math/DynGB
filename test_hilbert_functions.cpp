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

int main() {
  Monomial x0 { 1, 0, 0, 0, 0, 0 };
  WT_TYPE wt6 [6] { 3, 2, 1, 1, 1, 1 };
  list <Monomial> T0 = { x0 };
  cout << "Hilbert numerator of { " << x0 << "}\n";
  Dense_Univariate_Integer_Polynomial * hn = hilbert_numerator_bigatti(T0);
  cout << *hn << endl;
  cout << "Hilbert reduced numerator of { " << x0 << "}\n";
  Dense_Univariate_Integer_Polynomial * hnr = hilbert_second_numerator(6, hn);
  cout << *hnr << endl;
  cout << "Graded Hilbert numerator of { " << x0 << "}\n";
  cout << "\t(with grading (";
  for (unsigned i = 0; i < 6; ++i) cout << wt6[i] << ',';
  cout << ") )\n";
  Dense_Univariate_Integer_Polynomial * hng = hilbert_numerator_bigatti(T0, wt6);
  cout << *hng << endl;
  cout << "Graded Hilbert reduced numerator of { " << x0 << "}\n";
  cout << "\t(with grading (";
  for (unsigned i = 0; i < 6; ++i) cout << wt6[i] << ',';
  cout << ") )\n";
  Dense_Univariate_Integer_Polynomial * hngr = hilbert_second_numerator(6, hng, wt6);
  cout << *hngr << endl;
  cout << "Hilbert polynomial of { " << x0 << "}\n";
  Dense_Univariate_Rational_Polynomial * hp = hilbert_polynomial(6, 0, T0, hn);
  cout << *hp << endl;
  delete hn; delete hnr; delete hng; delete hngr;
  delete hp;
  cout << endl; cout << endl;
  //
  Monomial t { 5, 0 };
  Monomial u { 3, 1 };
  Monomial v { 1, 4 };
  WT_TYPE wt2 [2] { 5, 3 };
  list<Monomial> T { t, u, v };
  cout << "Hilbert numerator of { " << t << ", " << u << ", " << v << "}\n";
  hn = hilbert_numerator_bigatti(T);
  cout << *hn << endl;
  cout << "Hilbert reduced numerator of { " << t << ", " << u << ", " << v << "}\n";
  hnr = hilbert_second_numerator(2, hn);
  cout << *hnr << endl;
  cout << "Graded Hilbert numerator of { " << t << ", " << u << ", " << v << "}\n";
  cout << "\t(with grading (" << wt2[0] << ',' << wt2[1] << ") )\n";
  hng = hilbert_numerator_bigatti(T, wt2);
  cout << *hng << endl;
  cout << "Graded Hilbert reduced numerator of { " << t << ", " << u << ", " << v << "}\n";
  cout << "\t(with grading (" << wt2[0] << ',' << wt2[1] << ") )\n";
  hngr = hilbert_second_numerator(2, hng, wt2);
  cout << *hngr << endl;
  cout << "Hilbert polynomial of { " << t << ", " << u << ", " << v << "}\n";
  hp = hilbert_polynomial(2, 0, T, hn);
  cout << *hp << endl;
  delete hn; delete hnr; delete hng; delete hngr;
  delete hp;
  cout << endl; cout << endl;
  //
  Monomial t1 { 0, 2, 0, 0, 0 };
  Monomial t2 { 0, 1, 1, 0, 0 };
  Monomial t3 { 0, 0, 2, 0, 0 };
  Monomial t4 { 0, 1, 0, 1, 0 };
  Monomial t5 { 0, 0, 1, 1, 0 };
  Monomial t6 { 0, 0, 0, 2, 0 };
  WT_TYPE wt5 [5] { 1, 2, 3, 4, 5 };
  list<Monomial> T1 { t1, t2, t3, t4, t5, t6 };
  cout << "Hilbert numerator of { " << t1 << ',' << t2 << ',' << t3 << ','
       << t4 << ',' << t5 << ',' << t6 << "}\n";
  hn = hilbert_numerator_bigatti(T1);
  cout << *hn << endl;
  cout << "Hilbert reduced numerator of { " << t1 << ',' << t2 << ',' << t3 << ','
       << t4 << ',' << t5 << ',' << t6 << "}\n";
  hnr = hilbert_second_numerator(5, hn);
  cout << *hnr << endl;
  cout << "Graded Hilbert numerator of { " << t1 << ',' << t2 << ',' << t3 << ','
       << t4 << ',' << t5 << ',' << t6 << "}\n";
  cout << "\t(with grading (";
  for (unsigned i = 0; i < 5; ++i) cout << wt5[i] << ',';
  cout << ") )\n";
  hng = hilbert_numerator_bigatti(T1, wt5);
  cout << *hng << endl;
  cout << "Graded Hilbert reduced numerator of { " << t1 << ',' << t2 << ',' << t3 << ','
       << t4 << ',' << t5 << ',' << t6 << "}\n";
  cout << "\t(with grading (";
  for (unsigned i = 0; i < 5; ++i) cout << wt5[i] << ',';
  cout << ") )\n";
  hngr = hilbert_second_numerator(5, hng, wt5);
  cout << *hngr << endl;
  cout << "Hilbert polynomial:\n";
  hp = hilbert_polynomial(5, 0, T1, hn);
  cout << *hp << endl;
  delete hn; delete hnr; delete hng; delete hngr;
  delete hp;
  cout << endl;
  //
  Monomial u1 = { 1, 2, 3, 0, 0 };
  Monomial u2 = { 0, 1, 2, 3, 0 };
  Monomial u3 = { 0, 0, 1, 2, 3 };
  Monomial u4 = { 3, 0, 0, 1, 2 };
  Monomial u5 = { 2, 3, 0, 0, 1 };
  Monomial u6 = { 3, 2, 1, 0, 0 };
  Monomial u7 = { 0, 3, 2, 1, 0 };
  Monomial u8 = { 0, 0, 3, 2, 1 };
  Monomial u9 = { 1, 0, 0, 3, 2 };
  Monomial u10 = { 2, 1, 0, 0, 3 };
  Monomial u11 = { 2, 2, 2, 0, 0 };
  Monomial u12 = { 0, 2, 2, 2, 0 };
  Monomial u13 = { 0, 0, 2, 2, 2 };
  list<Monomial> U1 = { u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13 };
  cout << "Hilbert numerator of { " << u1 << ',' << u2 << ',' << u3 << ','
       << u4 << ',' << u5 << ',' << u6 << ',' << u7 << ',' << u8 << ','
       << u9 << ',' << u10 << ',' << u11 << ',' << u12 << ',' << u13 << "}\n";
  hn = hilbert_numerator_bigatti(U1);
  cout << *hn << endl;
  cout << "Hilbert reduced numerator of { " << u1 << ',' << u2 << ',' << u3 << ','
       << u4 << ',' << u5 << ',' << u6 << ',' << u7 << ',' << u8 << ','
       << u9 << ',' << u10 << ',' << u11 << ',' << u12 << ',' << u13 << "}\n";
  hnr = hilbert_second_numerator(5, hn);
  cout << *hnr << endl;
  cout << "Graded Hilbert numerator of { " << u1 << ',' << u2 << ',' << u3 << ','
       << u4 << ',' << u5 << ',' << u6 << ',' << u7 << ',' << u8 << ','
       << u9 << ',' << u10 << ',' << u11 << ',' << u12 << ',' << u13 << "}\n";
  cout << "\t(with grading (";
  for (unsigned i = 0; i < 5; ++i) cout << wt5[i] << ',';
  cout << ") )\n";
  hng = hilbert_numerator_bigatti(U1, wt5);
  cout << *hng << endl;
  cout << "Graded Hilbert reduced numerator of { " << u1 << ',' << u2 << ',' << u3 << ','
       << u4 << ',' << u5 << ',' << u6 << ',' << u7 << ',' << u8 << ','
       << u9 << ',' << u10 << ',' << u11 << ',' << u12 << ',' << u13 << "}\n";
  cout << "\t(with grading (";
  for (unsigned i = 0; i < 5; ++i) cout << wt5[i] << ',';
  cout << ") )\n";
  hngr = hilbert_second_numerator(5, hng, wt5);
  cout << *hngr << endl;
  cout << "Hilbert polynomial:\n";
  hp = hilbert_polynomial(5, 0, U1, hn);
  cout << *hp << endl;
  delete hn; delete hnr; delete hng; delete hngr;
  delete hp;
  cout << endl;  
}