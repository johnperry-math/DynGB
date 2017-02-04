#ifndef __DENSE_UNIVARIATE_RATIONAL_POLY_HPP_
#define __DENSE_UNIVARIATE_RATIONAL_POLY_HPP_

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

#include <iostream>

#include "system_constants.hpp"

using std::ostream; using std::cout; using std::endl;

/**
  @ingroup polygroup
  @class Dense_Univariate_Rational_Polynomial
  @author John Perry
  @date 2016
  @brief quick-&rsquo;n-dirty Dense_Univariate rational polynomial class
  @warning You must have \c deg smaller than \c size,
    as \c deg indexes the monomial of largest degree,
    which can be at most \c size&nbsp;-&nbsp;1.
*/
class Dense_Univariate_Rational_Polynomial {
public:
  /** @name Construction */
  ///@{
  /** @brief construct with the number of expected terms */
  Dense_Univariate_Rational_Polynomial(DEG_TYPE);
  /** @brief copy constructor */
  Dense_Univariate_Rational_Polynomial(
      const Dense_Univariate_Rational_Polynomial &
  );
  /** @brief construct from given data, with the expected number of terms */
  Dense_Univariate_Rational_Polynomial(DEG_TYPE, int64_t *, uint64_t *);
  ///@}
  /** @name Destruction */
  ///@{
  ~Dense_Univariate_Rational_Polynomial() { delete [] coeffs; }
  ///@}
  /** @name Modification */
  ///@{
  /**
    @brief expand to allow for higher-degree monomials
    @details You should do this at the beginning of any computation that
      could expand the size of the polynomial.
  */
  void expand_poly(DEG_TYPE);
  /** @brief set the coefficient of \f$x^k\f$ to \f$\frac{a}{b}\f$ */
  void set_coefficient(DEG_TYPE k, long a, unsigned long b) {
    expand_poly(k);
    MPZCOEF_TYPE az = a;
    MPZCOEF_TYPE bz = b;
    coeffs[k] = az;
    coeffs[k] /= bz;
    //coeffs[k] = a;
    //coeffs[k] /= b;
    if (k > deg and a != 0) { deg = k; }
    else if (k == deg and a == 0) {
      while (deg > 0 and coeffs[deg] == 0) { --deg; }
    }
  }
  /**
    @brief multiplies every monomial by a constant integer
  */
  void scale_by(COEF_TYPE a);
  /**
    @brief multiplies every monomial by a constant rational
    @param a rational&rsquo;s numerator
    @param b rational&rsquo;s denominator
  */
  void scale_by(COEF_TYPE a, UCOEF_TYPE b);
  /**
    @brief multiplies every monomial by a constant rational
  */
  void scale_by(MPQCOEF_TYPE);
  /**
    @brief a hopefully efficient multiplication algorithm
    @details Moves exponents from higher degrees to even higher degrees,
      freeing up space for smaller exponents.
  */
  void multiply_by_monomial_of_degree(DEG_TYPE);
  /**
    @brief highly inefficient polynomial multiplication (\f$O(mn)\f$)
  */
  void multiply_by(const Dense_Univariate_Rational_Polynomial &);
  /** @brief negates the numerators */
  void negate();
  /** @brief adds \c other to \c this */
  void add(const Dense_Univariate_Rational_Polynomial &);
  /** @brief alias for add() */
  void operator +=(const Dense_Univariate_Rational_Polynomial & other) { add(other); }
  /**
    @brief subtracts \c other from \c this
  */
  void subtract(const Dense_Univariate_Rational_Polynomial &);
  /** @brief alias for subtract() */
  void operator -=(const Dense_Univariate_Rational_Polynomial & other) { subtract(other); }
  ///@}
  /** @name Computation */
  ///@{
  /** @brief returns the difference between \c this and the other */
  Dense_Univariate_Rational_Polynomial operator-(
    const Dense_Univariate_Rational_Polynomial &) const;
  ///@}
  /** @name Basic properties */
  ///@{
  /** @brief indicates whether the polynomial is zero */
  bool is_zero() const {
    bool nonzero = true;
    for (DEG_TYPE i = 0; nonzero and i <= deg; ++i)
      nonzero = (coeffs[i] == 0);
    return nonzero;
  }
  /** @brief returns the \f$k\f$th numerator */
  MPZCOEF_TYPE numerator(DEG_TYPE k) const { return coeffs[k].get_num(); }
  /** @brief returns the \f$k\f$th denominator */
  MPZCOEF_TYPE denominator(DEG_TYPE k) const { return coeffs[k].get_den(); }
  /** @brief returns the polynomial&rsquo;s degree */
  DEG_TYPE degree() const { return deg; }
  ///@}
  /** @name I/O */
  ///@{
  friend ostream & operator<<(
      ostream & os, const Dense_Univariate_Rational_Polynomial & p
  ) {
    if (p.coeffs[p.deg] < 0) { os << '-'; }
    for (DEG_TYPE i = p.deg; i > 0; --i) {
      if (p.coeffs[i].get_num() != 0) {
        if (p.coeffs[i].get_num() == 1) {
          if (p.coeffs[i].get_den() != 1) {
            os << "1 / " << p.coeffs[i].get_den();
          }
        } else if (p.coeffs[i].get_num() == -1) {
          if (p.coeffs[i].get_den() != 1) {
            os << "1 / " << p.coeffs[i].get_den();
          }
        } else {
          if (p.coeffs[i].get_num() > 0) { os << p.coeffs[i].get_num(); }
          else { os << -p.coeffs[i].get_num(); }
          if (p.coeffs[i].get_den() != 1) { os << " / " << p.coeffs[i].get_den(); }
        }
        if (i != 1) { os << "t^" << i; }
        else { os << "t"; }
      }
      if (p.coeffs[i - 1].get_num() < 0) { os << " - "; }
      else if (p.coeffs[i - 1].get_num() > 0) { os << " + "; }
    }
    if (p.coeffs[0].get_num() != 0) {
      if (p.coeffs[0].get_num() > 0) { os << p.coeffs[0].get_num(); }
      else if (p.coeffs[0].get_num() < 0) { os << -p.coeffs[0].get_num(); }
      if (p.coeffs[0].get_den() != 1) { os << " / " << p.coeffs[0].get_den(); }
    }
    return os;
  }
  ///@}
protected:
  /** @brief list of coefficients; index 0 is the constant&rsquo;s */
  MPQCOEF_TYPE * coeffs;
  /** @brief degree of the polynomial (largest nonzero exponent) */
  DEG_TYPE deg;
  /** @brief number of slots for coefficients */
  DEG_TYPE size;
};

#endif