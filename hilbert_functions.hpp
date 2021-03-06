#ifndef __HILBERT_FUNCTIONS_
#define __HILBERT_FUNCTIONS_

#include <set>
#include <list>
#include <cstring>
#include <iostream>

using std::set; using std::list;

#include "system_constants.hpp"

#include "fields.hpp"
#include "monomial.hpp"
#include "polynomial_linked_list.hpp"
#include "dense_univariate_integer_poly.hpp"
#include "dense_univariate_rational_poly.hpp"

/**
  @defgroup commalg Commutative Algebra
  @brief Classes useful for computational Commutative Algebra
*/

/**
  @defgroup utils Utilities
  @brief Useful functions and classes that don't fit in easily elsewhere
*/

/**
  @ingroup utils
  @author John Perry
  @date 2016
  @brief divides out the common term of the two given numbers
  @details Useful for simplifying terms of Dense_Univariate_Rational_Polynomial.
    Uses the (simple) Euclidean Algorithm to determine the gcd.
*/
void divide_by_common_term(COEF_TYPE &, UCOEF_TYPE &);

/**
  @ingroup commalg
  @brief test for the 0-base case @cite Bigatti97
  @return @c true if and only if the zero base case applies
*/
bool is_zero_base_case(const list<Monomial> &);

/**
  @ingroup commalg
  @brief computes Hilbert numerator when the 0-base case applies @cite Bigatti97
  @param T list of generators of monomial ideal
  @param grading use for a graded polynomial ring
  @return the Hilbert polynomial for the 0-base case
  @details If you do not want to use a grading, set it to @c nullptr.
*/
Dense_Univariate_Integer_Polynomial * solve_zero_base_case(
    const list<Monomial> & T, const WT_TYPE * grading
);

/**
  @ingroup commalg
  @brief test for the 1-base case @cite Bigatti97
  @details If result is <c>T.end()</c>, this is not a 1-base case.
  @return index to an element of @p T that makes the 0-base case.
*/
list<Monomial>::const_iterator is_one_base_case(const list<Monomial> &);

/**
  @ingroup commalg
  @brief applies Bigatti&rsquo;s algorithm for the 1-base case
  @param T list of generators of monomial ideal
  @param ti index of the one term in @c T which is \e not a simple power
  @param grading use for a graded polynomial ring
  @return Hilbert numerator for a 1-base case
  @details If you do not want to use a grading, set it to @c nullptr.
*/
Dense_Univariate_Integer_Polynomial * solve_one_base_case(
    const list<Monomial> & T, list<Monomial>::const_iterator ti,
    const WT_TYPE * grading
);

/**
  @ingroup commalg
  @brief test for the &ldquo;splitting case&rdquo; @cite Bigatti97
  @param T list of generators of monomial ideal
  @param U empty list of iterators on T; will be modified
  @param V empty list of iterators on T; will be modified
  @return @c true if and only if @p T can be split into two lists
    of relatively prime monomials
  @details If the result is @c true, the monomials indexed by iterators in @p U
    should be relatively prime to the monomials indexed by iterators in @c V.
  @warning If the result is @c false, disregard the entries of @c U and @c V.
*/
bool is_splitting_case(
    const list<Monomial> & T,
    list< list<Monomial>::const_iterator > & U,
    list< list<Monomial>::const_iterator > & V
);

/**
  @ingroup commalg
  @brief applies Bigatti&rsquo;s algorithm for the 1-base case
  @param T list of generators of monomial ideal
  @param U list of iterators indexing elements of @p T
  @param V list of iterators indexing elements of @p U
  @param grading use for a graded polynomial ring
  @return the Hilbert numerator for the splitting case @cite Bigatti97
  @details If you do not want to use a grading, set it to @c nullptr.
  @warning The entries of @p U @p V need to be relatively prime
    to the entries of @p V; otherwise the result is wrong.
*/
Dense_Univariate_Integer_Polynomial * solve_splitting_case(
    const list<Monomial> & T,
    const list< list<Monomial>::const_iterator > & U,
    const list< list<Monomial>::const_iterator > & V,
    const WT_TYPE * grading
);

/**
  @ingroup commalg
  @brief chooses a pivot for the Bigatti algorithm
  @param T list of generators of monomial ideal
  @return a monomial used to pivot in the Bigatti algorithm @cite Bigatti97
  @details This uses the strategy described by @cite RouneHilbert2010
    and due to @cite Bigatti97 : select @f$x_i^e@f$ such that @f$i@f$ maximizes the number
    of terms divisible by that variable, and @f$e@f$ is the median of the powers
    of @f$x_i@f$ appearing in @f$T@f$.
*/
Monomial choose_hilbert_pivot(const list<Monomial> & T);

/**
  @ingroup commalg
  @brief the Bigatti algorithm to compute the Hilbert numerator @cite Bigatti97
  @param T list of generators of a monomial ideal
  @param grading use for a graded polynomial ring
  @return the basic, unreduced Hilbert numerator of the ideal generated by
      @f$ T @f$ .
  @details This computes the first Hilbert numerator.
      If you do not want to use a grading, set it to @c nullptr.
  @warning @c T needs to be non-redundant; i.e., @f$\forall t,u\in T@f$
    @f$t\not\mid u@f$.
*/
Dense_Univariate_Integer_Polynomial * hilbert_numerator_bigatti(
    const list<Monomial> & T, const WT_TYPE * grading = nullptr
);

/**
  @ingroup commalg
  @brief computes the second Hilbert numerator (after reduction by
    @f$(1-t)^{\deg I}@f$)
  @param n number of variables (used to bound computation)
  @param hn the first Hilbert numerator (before reduction)
  @param grading the grading of the ideal (if @c NULLPTR , the standard grading)
  @return the reduced Hilbert numerator of the specified Hilbert numerator
*/
Dense_Univariate_Integer_Polynomial * hilbert_second_numerator(
  NVAR_TYPE n,
  Dense_Univariate_Integer_Polynomial * hn,
  const WT_TYPE * grading = nullptr
);

/**
  @ingroup commalg
  @brief computes the dimension of the ideal
    by subtracting the Hilbert numerators
  @param nvars number of variables in ring
  @param first_numerator the unreduced Hilbert numerator
  @param second_numerator the reduced Hilbert numerator
  @return dimension of the ideal whose Hilbert numerators are given
*/
unsigned ideal_dimension(
  NVAR_TYPE nvars,
  const Dense_Univariate_Integer_Polynomial *first_numerator,
  const Dense_Univariate_Integer_Polynomial *second_numerator
);

/**
  @ingroup commalg
  @brief computes the number of combinations @f$C(t+a,b)@f$
  @param a an integer
  @param b an integer
  @return the polynomial @f$C(t+a,b)@f$
*/
Dense_Univariate_Rational_Polynomial * polynomial_binomial(
  long long a, long long b
);

/**
  @ingroup commalg
  @brief computes the Hilbert polynomial for an ideal
  @param n number of variables in ideal
  @param dim dimension of ideal (if known)
  @param T ideal&rsquo;s generators
  @param hn Hilbert numerator (if known)
  @param hn2 reduced Hilbert numerator (if known)
  @return the Hilbert polynomial of the ideal generated by @p T
*/
Dense_Univariate_Rational_Polynomial * hilbert_polynomial(
    NVAR_TYPE n,
    unsigned int dim,
    const list<Monomial> T,
    Dense_Univariate_Integer_Polynomial * hn = nullptr,
    Dense_Univariate_Integer_Polynomial * hn2 = nullptr
);

#endif