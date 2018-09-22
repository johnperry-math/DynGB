#ifndef __DENSE_UNIVARIATE_RATIONAL_POLY_CPP_
#define __DENSE_UNIVARIATE_RATIONAL_POLY_CPP_

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
  numerators = new COEF_TYPE [n];
  denominators = new UCOEF_TYPE [n];
  for (DEG_TYPE i = 0; i < n; ++i) {
    numerators[i] = 0;
    denominators[i] = 1;
  }
  size = n;
  deg = 0;
}

Dense_Univariate_Rational_Polynomial::Dense_Univariate_Rational_Polynomial(
  const Dense_Univariate_Rational_Polynomial & other
) {
  size = other.size;
  deg = other.deg;
  numerators = new COEF_TYPE [size] { 0 };
  denominators = new UCOEF_TYPE [size] { 1 };
  DEG_TYPE i = 0;
  for (/* */; i <= deg; ++i) {
    numerators[i] = other.numerator(i);
    denominators[i] = other.denominator(i);
  }
  for (/* */; i < size; ++i) {
    numerators[i] = 0;
    denominators[i] = 1;
  }
}

Dense_Univariate_Rational_Polynomial::Dense_Univariate_Rational_Polynomial(
    DEG_TYPE n, int64_t * nums, uint64_t * denoms
) {
  size = n + 1;
  deg = n;
  numerators = new COEF_TYPE [size] { 0 };
  denominators = new UCOEF_TYPE [size] { 1 };
  for (DEG_TYPE i = 0; i <= deg; ++i) {
    numerators[i] = nums[i];
    denominators[i] = denoms[i];
  }
}

void Dense_Univariate_Rational_Polynomial::expand_poly(DEG_TYPE n) {
  if (n + 1 > size) {
    COEF_TYPE * new_nums = new COEF_TYPE [n + 1];
    UCOEF_TYPE * new_dens = new UCOEF_TYPE [n + 1];
    for (DEG_TYPE i = 0; i < deg + 1; ++i) {
      new_nums[i] = numerators[i];
      new_dens[i] = denominators[i];
    }
    delete [] numerators;
    delete [] denominators;
    numerators = new_nums;
    denominators = new_dens;
    for (DEG_TYPE i = deg + 1; i < n + 1; ++i) {
      numerators[i] = 0;
      denominators[i] = 1;
    }
    size = n + 1;
  }
}

void Dense_Univariate_Rational_Polynomial::scale_by(COEF_TYPE a) {
  for (DEG_TYPE i = 0; i <= deg; ++i)
    if (numerators[i] != 0)
      set_coefficient(i, numerators[i] * a, denominators[i]);
}

void Dense_Univariate_Rational_Polynomial::scale_by(
  COEF_TYPE a, UCOEF_TYPE b
) {
  for (DEG_TYPE i = 0; i <= deg; ++i)
    if (numerators[i] != 0) {
      set_coefficient(i, numerators[i] * a, numerators[i] * b);
    }
}

void Dense_Univariate_Rational_Polynomial::multiply_by_monomial_of_degree(
  DEG_TYPE k
) {
  expand_poly(deg + k);
  for (DEG_TYPE i = deg; i > 0; --i) {
    if (numerators[i] != 0) {
      set_coefficient(i + k, numerators[i], denominators[i]);
      set_coefficient(i, 0, 1);
    }
  }
  set_coefficient(k, numerators[0], denominators[0]);
  set_coefficient(0, 0, 1);
}

void Dense_Univariate_Rational_Polynomial::multiply_by(
  const Dense_Univariate_Rational_Polynomial & q
) {
  DEG_TYPE n = deg + q.deg + 1; // add 1 in case deg == q.deg == 0
  n = (n > size) ? n : size;
  COEF_TYPE * new_nums = new COEF_TYPE [n];
  UCOEF_TYPE * new_dens = new UCOEF_TYPE [n];
  for (DEG_TYPE i = 0; i < n; ++i) {
    new_nums[i] = 0;
    new_dens[i] = 1;
  }
  for (DEG_TYPE i = 0; i < deg + 1; ++i)
    for (DEG_TYPE j = 0; j < q.deg + 1; ++j) {
      if (numerators[i] != 0 and q.numerators[j] != 0) {
        COEF_TYPE new_a = numerators[i] * q.numerators[j];
        UCOEF_TYPE new_b = denominators[i] * q.denominators[j];
        new_nums[i + j] = new_nums[i + j] * new_b + new_dens[i + j] * new_a;
        new_dens[i + j] *= new_b;
        divide_by_common_term<COEF_TYPE, UCOEF_TYPE>(new_nums[i + j], new_dens[i + j]);
      }
    }
  delete [] numerators;
  delete [] denominators;
  numerators = new_nums;
  denominators = new_dens;
  size = n;
  deg = deg + q.deg;
}

void Dense_Univariate_Rational_Polynomial::negate() {
  for (DEG_TYPE i = 0; i <= deg; ++i)
    if (numerators[i] != 0)
      numerators[i] = -numerators[i];
}

void Dense_Univariate_Rational_Polynomial::add(
  const Dense_Univariate_Rational_Polynomial & q
) {
  DEG_TYPE new_deg = (deg > q.deg) ? deg : q.deg;
  expand_poly(new_deg);
  deg = new_deg;
  for (DEG_TYPE i = 0; i <= q.deg; ++i) {
    if (q.numerators[i] != 0) {
      numerators[i] = numerators[i] * q.denominators[i]
          + q.numerators[i] * denominators[i];
      denominators[i] *= q.denominators[i];
      if (numerators[i] != 0)
        divide_by_common_term(numerators[i], denominators[i]);
    }
  }
  if (numerators[deg] == 0) {
    DEG_TYPE i = deg;
    while (i > 0 and numerators[i] == 0)
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
  for (DEG_TYPE i = 0; i <= q.deg; ++i) {
    if (q.numerators[i] != 0) {
      numerators[i] = numerators[i] * q.denominators[i]
          - q.numerators[i] * denominators[i];
      denominators[i] *= q.denominators[i];
      if (numerators[i] != 0)
        divide_by_common_term(numerators[i], denominators[i]);
    }
  }
  if (numerators[deg] == 0) {
    DEG_TYPE i = deg;
    while (i > 0 and numerators[i] == 0)
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
    result.set_coefficient(i,
        numerators[i] * other.denominator(i)
            - other.numerator(i) * denominators[i],
        denominators[i] * other.denominator(i)
    );
  // only one of the next two loops will be performed
  while (i < deg) {
    result.set_coefficient(i, numerators[i], denominators[i]);
    ++i;
  }
  while (i < other.degree()) {
    result.set_coefficient(i, -other.numerator(i), other.denominator(i));
    ++i;
  }
  return result;
}

#endif