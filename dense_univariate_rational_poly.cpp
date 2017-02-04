#ifndef __DENSE_UNIVARIATE_RATIONAL_POLY_CPP_
#define __DENSE_UNIVARIATE_RATIONAL_POLY_CPP_

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

#include "dense_univariate_rational_poly.hpp"

template <typename T, typename U>
void divide_by_common_term(T & a, U & b) {
  T c = (a < 0) ? -a : a;
  //U d = (b < 0) ? -b : b;
  U d = b;
  while (d != 0) {
    T r = c % d;
    c = d;
    d = r;
  }
  a /= c;
  b /= c;
}

Dense_Univariate_Rational_Polynomial::Dense_Univariate_Rational_Polynomial(
  DEG_TYPE n
) {
  coeffs = new MPQCOEF_TYPE [n];
  for (DEG_TYPE i = 0; i < n; ++i)
    coeffs[i] = 0;
  size = n;
  deg = 0;
}

Dense_Univariate_Rational_Polynomial::Dense_Univariate_Rational_Polynomial(
  const Dense_Univariate_Rational_Polynomial & other
) {
  size = other.size;
  deg = other.deg;
  coeffs = new MPQCOEF_TYPE [size] { 0 };
  for (DEG_TYPE i = 0; i < size; ++i)
    coeffs[i] = other.coeffs[i];
}

Dense_Univariate_Rational_Polynomial::Dense_Univariate_Rational_Polynomial(
    DEG_TYPE n, int64_t * nums, uint64_t * denoms
) {
  size = n + 1;
  deg = n;
  coeffs = new MPQCOEF_TYPE [size] { 0 };
  for (DEG_TYPE i = 0; i <= deg; ++i)
    coeffs[i] = long(nums[i]) / (unsigned long)(denoms[i]);
}

void Dense_Univariate_Rational_Polynomial::expand_poly(DEG_TYPE n) {
  if (n + 1 > size) {
    MPQCOEF_TYPE * new_coeffs = new MPQCOEF_TYPE [n + 1];
    for (DEG_TYPE i = 0; i < deg + 1; ++i)
      new_coeffs[i] = coeffs[i];
    delete [] coeffs;
    coeffs = new_coeffs;
    for (DEG_TYPE i = deg + 1; i < n + 1; ++i)
      coeffs[i] = 0;
    size = n + 1;
  }
}

void Dense_Univariate_Rational_Polynomial::scale_by(COEF_TYPE a) {
  for (DEG_TYPE i = 0; i <= deg; ++i)
    if (coeffs[i] != 0)
      coeffs[i] *= (long)a;
}

void Dense_Univariate_Rational_Polynomial::scale_by(MPQCOEF_TYPE a) {
  for (DEG_TYPE i = 0; i <= deg; ++i)
    if (coeffs[i] != 0)
      coeffs[i] *= a;
}

void Dense_Univariate_Rational_Polynomial::scale_by(
  COEF_TYPE a, UCOEF_TYPE b
) {
  for (DEG_TYPE i = 0; i <= deg; ++i)
    if (coeffs[i] != 0) {
      coeffs[i] *= (long)a;
      coeffs[i] /= (unsigned long)b;
    }
}

void Dense_Univariate_Rational_Polynomial::multiply_by_monomial_of_degree(
  DEG_TYPE k
) {
  expand_poly(deg + k);
  for (DEG_TYPE i = deg; i > 0; --i) {
    if (coeffs[i] != 0) {
      coeffs[i + k] = coeffs[i];
      coeffs[i] = 0;
    }
  }
  coeffs[k] = coeffs[0];
  coeffs[0] = 0;
}

void Dense_Univariate_Rational_Polynomial::multiply_by(
  const Dense_Univariate_Rational_Polynomial & q
) {
  DEG_TYPE n = deg + q.deg + 1; // add 1 in case deg == q.deg == 0
  n = (n > size) ? n : size;
  MPQCOEF_TYPE * new_coeffs = new MPQCOEF_TYPE [n];
  for (DEG_TYPE i = 0; i < n; ++i)
    new_coeffs[i] = 0;
  for (DEG_TYPE i = 0; i < deg + 1; ++i)
    for (DEG_TYPE j = 0; j < q.deg + 1; ++j) {
      if (coeffs[i] != 0 and q.coeffs[j] != 0)
        new_coeffs[i + j] += coeffs[i] * q.coeffs[j];
    }
  delete [] coeffs;
  coeffs = new_coeffs;
  size = n;
  deg = deg + q.deg;
}

void Dense_Univariate_Rational_Polynomial::negate() {
  for (DEG_TYPE i = 0; i <= deg; ++i)
    if (coeffs[i] != 0)
      coeffs[i] = -coeffs[i];
}

void Dense_Univariate_Rational_Polynomial::add(
  const Dense_Univariate_Rational_Polynomial & q
) {
  DEG_TYPE new_deg = (deg > q.deg) ? deg : q.deg;
  expand_poly(new_deg);
  deg = new_deg;
  for (DEG_TYPE i = 0; i <= q.deg; ++i)
    if (q.coeffs[i] != 0)
      coeffs[i] += q.coeffs[i];
  if (coeffs[deg] == 0) {
    DEG_TYPE i = deg;
    while (i > 0 and coeffs[i] == 0)
      --i;
    deg = i;
  }
}

void Dense_Univariate_Rational_Polynomial::subtract(
  const Dense_Univariate_Rational_Polynomial & q
) {
  DEG_TYPE new_deg = (deg > q.deg) ? deg : q.deg;
  expand_poly(new_deg);
  deg = new_deg;
  for (DEG_TYPE i = 0; i <= q.deg; ++i)
    if (q.coeffs[i] != 0)
      coeffs[i] -= q.coeffs[i];
  if (coeffs[deg] == 0) {
    DEG_TYPE i = deg;
    while (i > 0 and coeffs[i] == 0)
      --i;
    deg = i;
  }
}

Dense_Univariate_Rational_Polynomial
Dense_Univariate_Rational_Polynomial::operator-(
    const Dense_Univariate_Rational_Polynomial & other) const {
  DEG_TYPE m = (deg < other.degree()) ? deg : other.degree();
  DEG_TYPE n = (deg > other.degree()) ? deg : other.degree();
  Dense_Univariate_Rational_Polynomial result(n + 1);
  DEG_TYPE i = 0;
  for ( /* already initialized */; i <= m; ++i)
    result.coeffs[i] = coeffs[i] - other.coeffs[i];
  // only one of the next two loops will be performed
  while (i < deg) {
    result.coeffs[i] = coeffs[i];
    ++i;
  }
  while (i < other.degree()) {
    result.coeffs[i] = -other.coeffs[i];
    ++i;
  }
  return result;
}

#endif