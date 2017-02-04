#ifndef __DENSE_UNIVARIATE_RATIONAL_POLY_HPP_
#define __DENSE_UNIVARIATE_RATIONAL_POLY_HPP_

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
  ~Dense_Univariate_Rational_Polynomial() {
    delete [] numerators;
    delete [] denominators;
  }
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
  void set_coefficient(DEG_TYPE k, COEF_TYPE a, UCOEF_TYPE b) {
    expand_poly(k);
    numerators[k] = a;
    denominators[k] = b;
    if (k > deg and a != 0) { deg = k; }
    else if (k == deg and a == 0) {
      while (deg > 0 and numerators[deg] == 0) { --deg; }
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
      nonzero = (numerators[i] == 0);
    return nonzero;
  }
  /** @brief returns the \f$k\f$th numerator */
  COEF_TYPE numerator(DEG_TYPE k) const { return numerators[k]; }
  /** @brief returns the \f$k\f$th denominator */
  UCOEF_TYPE denominator(DEG_TYPE k) const { return denominators[k]; }
  /** @brief returns the polynomial&rsquo;s degree */
  DEG_TYPE degree() const { return deg; }
  /** @brief all the numerators */
  const COEF_TYPE * numerator_array() { return numerators; }
  /** @brief all the denominators */
  const UCOEF_TYPE * denominator_array() { return denominators; }
  ///@}
  /** @name I/O */
  ///@{
  friend ostream & operator<<(
      ostream & os, Dense_Univariate_Rational_Polynomial & p
  ) {
    if (p.numerators[p.deg] < 0) { os << '-'; }
    for (DEG_TYPE i = p.deg; i > 0; --i) {
      if (p.numerators[i] != 0) {
        if (p.numerators[i] == 1) {
          if (p.denominators[i] != 1) {
            os << "1 / " << p.denominators[i];
          }
        } else if (p.numerators[i] == -1) {
          if (p.denominators[i] != 1) {
            os << "1 / " << p.denominators[i];
          }
        } else {
          if (p.numerators[i] > 0) { os << p.numerators[i]; }
          else { os << -p.numerators[i]; }
          if (p.denominators[i] != 1) { os << " / " << p.denominators[i]; }
        }
        if (i != 1) { os << "t^" << i; }
        else { os << "t"; }
      }
      if (p.numerators[i - 1] < 0) { os << " - "; }
      else if (p.numerators[i - 1] > 0) { os << " + "; }
    }
    if (p.numerators[0] != 0) {
      if (p.numerators[0] > 0) { os << p.numerators[0]; }
      else if (p.numerators[0] < 0) { os << -p.numerators[0]; }
      if (p.denominators[0] != 1) { os << " / " << p.denominators[0]; }
    }
    return os;
  }
  ///@}
protected:
  /** @brief list of numerators; index 0 is the constant&rsquo;s */
  COEF_TYPE * numerators;
  /** @brief list of denominators; index 0 is the constant&rsquo;s */
  UCOEF_TYPE * denominators;
  /** @brief degree of the polynomial (largest nonzero exponent) */
  DEG_TYPE deg;
  /** @brief number of slots for coefficients */
  DEG_TYPE size;
};

#endif