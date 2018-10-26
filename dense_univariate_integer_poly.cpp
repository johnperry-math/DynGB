#ifndef __DENSE_UNIVARIATE_INTEGER_POLY_CPP_
#define __DENSE_UNIVARIATE_INTEGER_POLY_CPP_

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

#include "dense_univariate_integer_poly.hpp"

using std::cout; using std::endl;

Dense_Univariate_Integer_Polynomial::Dense_Univariate_Integer_Polynomial(
  DEG_TYPE n
) {
  coeffs = new MPZCOEF_TYPE [n];
  for (DEG_TYPE i = 0; i < n; ++i)
    coeffs[i] = 0;
  size = n;
  deg = 0;
}

Dense_Univariate_Integer_Polynomial::Dense_Univariate_Integer_Polynomial(
  const Dense_Univariate_Integer_Polynomial & p
) {
  size = p.size;
  deg = p.deg;
  coeffs = new MPZCOEF_TYPE[size];
  for (DEG_TYPE i = 0; i <= deg; ++i)
    coeffs[i] = p.coeffs[i];
}

Dense_Univariate_Integer_Polynomial::Dense_Univariate_Integer_Polynomial(
  DEG_TYPE n, int64_t * C
) {
  size = n + 1;
  deg = n;
  coeffs = new MPZCOEF_TYPE[size] { 0 };
  for (DEG_TYPE i = 0; i <= deg; ++i) {
    coeffs[i] = (long )(C[i] >> 32);
    coeffs[i] << 32;
    coeffs[i] = coeffs[i] + (long )(C[i] % (((long long)1) << 32));
  }
}

void Dense_Univariate_Integer_Polynomial::expand_poly(DEG_TYPE n) {
  if (n + 1 > size) {
    MPZCOEF_TYPE * new_nums = new MPZCOEF_TYPE [n + 1];
    for (DEG_TYPE i = 0; i < deg + 1; ++i)
      new_nums[i] = coeffs[i];
    delete [] coeffs;
    coeffs = new_nums;
    for (DEG_TYPE i = deg + 1; i < n + 1; ++i)
      coeffs[i] = 0;
    size = n + 1;
  }
}

void Dense_Univariate_Integer_Polynomial::set_coefficient(
  DEG_TYPE k, MPZCOEF_TYPE a
) {
  coeffs[k] = a;
  if (k > deg and a != 0) { deg = k; }
  else if (k == deg and a == 0) {
    while (deg > 0 and coeffs[deg] == 0) { --deg; }
  }
}

void Dense_Univariate_Integer_Polynomial::scale_by(MPZCOEF_TYPE a) {
  for (DEG_TYPE i = 0; i <= deg; ++i)
    if (coeffs[i] != 0)
      set_coefficient(i, coeffs[i] * a);
}

void Dense_Univariate_Integer_Polynomial::multiply_by_monomial_of_degree(
  DEG_TYPE k
) {
  expand_poly(deg + k);
  for (DEG_TYPE i = deg; i > 0; --i) {
    if (coeffs[i] != 0) {
      set_coefficient(i + k, coeffs[i]);
      set_coefficient(i, 0);
    }
  }
  set_coefficient(k, coeffs[0]);
  set_coefficient(0, 0);
}

void Dense_Univariate_Integer_Polynomial::multiply_by(
  const Dense_Univariate_Integer_Polynomial & q
) {
  DEG_TYPE n = deg + q.deg + 1; // add 1 in case deg == q.deg == 0
  DEG_TYPE nq = q.deg;
  n = (n > size) ? n : size;
  MPZCOEF_TYPE * new_nums = new MPZCOEF_TYPE [n] { 0 };
  MPZCOEF_TYPE * b = q.coeffs;
  /*for (DEG_TYPE i = 0; i < n; ++i)
    new_nums[i] = 0;*/
  for (DEG_TYPE i = 0; i < deg + 1; ++i)
    for (DEG_TYPE j = 0; j < nq + 1; ++j)
      if (coeffs[i] != 0 and b[j] != 0)
        new_nums[i + j] += coeffs[i] * b[j];
  delete [] coeffs;
  coeffs = new_nums;
  size = n;
  deg = deg + q.deg;
}

void Dense_Univariate_Integer_Polynomial::negate() {
  for (DEG_TYPE i = 0; i <= deg; ++i)
    if (coeffs[i] != 0)
      coeffs[i] = -coeffs[i];
}

void Dense_Univariate_Integer_Polynomial::add(
  const Dense_Univariate_Integer_Polynomial & q
) {
  DEG_TYPE new_deg = (deg > q.deg) ? deg : q.deg;
  expand_poly(new_deg);
  deg = new_deg;
  for (DEG_TYPE i = 0; i <= q.deg; ++i)
    if (q.coeffs[i] != 0)
      set_coefficient(i, coeffs[i] + q.coeffs[i]);
}

Dense_Univariate_Integer_Polynomial
Dense_Univariate_Integer_Polynomial::operator-(
  const Dense_Univariate_Integer_Polynomial & q
) const {
  DEG_TYPE n = (deg > q.degree()) ? deg : q.degree();
  Dense_Univariate_Integer_Polynomial r(n+1);
  DEG_TYPE i = 0;
  for ( /* already initialized */ ; i <= n; ++i)
    r.set_coefficient(i, coeffs[i] - q.coeff(i));
  // only one of the next two while loops should be executed
  while (i < deg) {
    r.set_coefficient(i, coeffs[i]);
    ++i;
  }
  while (i < q.degree()) {
    r.set_coefficient(i, -q.coeff(i));
    ++i;
  }
  return r;
}

#endif