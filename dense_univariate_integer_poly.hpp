#ifndef __DENSE_UNIVARIATE_INTEGER_POLY_HPP_
#define __DENSE_UNIVARIATE_INTEGER_POLY_HPP_

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

#include "system_constants.hpp"

using std::ostream;

/**
  @ingroup polygroup
  @class Dense_Univariate_Integer_Polynomial
  @author John Perry
  @date 2016
  @brief quick-&rsquo;n-dirty Dense_Univariate integer polynomial class
  @warning You must have @c deg smaller than @c size,
    as @c deg indexes the monomial of largest degree,
    which can be at most @c size&nbsp;-&nbsp;1.
*/
class Dense_Univariate_Integer_Polynomial {
public:
  /** @name Construction */
  ///@{
  /** @brief construct with the number of expected terms */
  explicit Dense_Univariate_Integer_Polynomial(DEG_TYPE);
  /** @brief copy constructor */
  Dense_Univariate_Integer_Polynomial(
      const Dense_Univariate_Integer_Polynomial &
  );
  /** @brief from serialized data */
  Dense_Univariate_Integer_Polynomial(DEG_TYPE, int64_t *);
  ///@}
  /** @name Destruction */
  ///@{
  ~Dense_Univariate_Integer_Polynomial() { delete [] coeffs; }
  ///@}
  /** @name Modification */
  ///@{
  /**
    @brief expand to allow for higher-degree monomials
    @details Do this at the beginning of any computation that
      could expand the size of the polynomial.
  */
  void expand_poly(DEG_TYPE);
  /** @brief set the coefficient of @f$x^k@f$ to @f$\frac{a}{b}@f$ */
  void set_coefficient(DEG_TYPE k, MPZCOEF_TYPE a);
  /**
    @brief multiplies every monomial by a constant
  */
  void scale_by(MPZCOEF_TYPE a);
  /**
    @brief a hopefully efficient multiplication algorithm
    @details Moves exponents from higher degrees to even higher degrees,
      freeing up space for smaller exponents.
  */
  void multiply_by_monomial_of_degree(DEG_TYPE);
  /**
    @brief highly inefficient polynomial multiplication (@f$O(mn)@f$)
  */
  void multiply_by(const Dense_Univariate_Integer_Polynomial &);
  /** @brief negates the coefficients */
  void negate();
  /**
    @brief reasonably efficient, given our dense representation
    @param q polynomial to add to @c this
  */
  void add(const Dense_Univariate_Integer_Polynomial & q);
  /**
    @brief @see add()
    @param q polynomial to add to @c this
    @return @c this, after adding @p q
  */
  Dense_Univariate_Integer_Polynomial & operator +=(
      const Dense_Univariate_Integer_Polynomial & q
  ) {
    add(q);
    return *this;
  }
  ///@}
  /** @name Computation */
  ///@{
  /** @brief returns the result of subtracting the other from @c this */
  Dense_Univariate_Integer_Polynomial operator-(
      const Dense_Univariate_Integer_Polynomial &
  ) const;
  ///@}
  /** @name Basic properties */
  ///@{
  /** @brief the @f$k@f$th coefficient */
  MPZCOEF_TYPE coeff(DEG_TYPE k) const { return coeffs[k]; }
  /** @brief synonym for @c coeff(k) */
  MPZCOEF_TYPE operator[](DEG_TYPE k) const { return coeffs[k]; }
  /** @brief the polynomial&rsquo;s degree (exponent of largest nonzero term) */
  DEG_TYPE degree() const { return deg; }
  /** @brief returns @c True if and only if every valid coefficient is zero */
  bool is_zero() const {
    bool result = true;
    for (DEG_TYPE i = 0; result and i <= deg; ++i)
      result = (coeffs[i] == 0);
    return result;
  }
  /** @brief all the coefficients */
  const MPZCOEF_TYPE * coefficient_array() { return coeffs; }
  ///@}
  /** @name I/O */
  ///@{
  friend ostream & operator<<(
      ostream & os, const Dense_Univariate_Integer_Polynomial & p
  ) {
    if (p.coeffs[p.deg] < 0) { os << '-'; }
    for (DEG_TYPE i = p.deg; i > 0; --i) {
      if (p.coeffs[i] != 0) {
        if (p.coeffs[i] != 1 and p.coeffs[i] != -1) {
          if (p.coeffs[i] > 0) { os << p.coeffs[i]; }
          else { os << -p.coeffs[i]; }
          os << '*';
        }
        if (i != 1) { os << "t^" << i; }
        else { os << "t"; }
      }
      if (p.coeffs[i - 1] < 0) { os << " - "; }
      else if (p.coeffs[i - 1] > 0) { os << " + "; }
    }
    if (p.coeffs[0] != 0) {
      if (p.coeffs[0] > 0) { os << p.coeffs[0]; }
      else if (p.coeffs[0] < 0) { os << -p.coeffs[0]; }
    }
    return os;
  }
  ///@}
protected:
  /** @brief list of numerators; index 0 is the constant&rsquo;s */
  MPZCOEF_TYPE * coeffs;
  /** @brief degree of the polynomial (largest nonzero exponent) */
  DEG_TYPE deg;
  /** @brief number of slots for coefficients */
  DEG_TYPE size;
};

#endif