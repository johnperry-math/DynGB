#ifndef __ARRAY_POLYNOMIALS_H_
#define __ARRAY_POLYNOMIALS_H_

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

#include <cstdlib>
#include <iostream>
#include <list>
using std::list;
#include <vector>
using std::vector;

#include "polynomial.hpp"

class Constant_Polynomial;

/**
  @class Constant_Polynomial_Iterator
  @author John Perry
  @date 2015
  @brief Iterates through a Constant_Polynomial.
  @ingroup IteratorGroup
*/
class Constant_Polynomial_Iterator : public Polynomial_Iterator
{
public:
  /** @name Construction */
  ///@{
  /**
    @brief Creates an iterator for <c>poly</c> and starts at the leading term.
    @param at_end whether to start at the end of the polynomial
    @param q polynomial to iterate over
  */
  Constant_Polynomial_Iterator(const Constant_Polynomial *q, bool at_end=false);
  ///@}
  /** @name Destruction */
  ///@{
  ~Constant_Polynomial_Iterator();
  ///@}
  /** @name Iteration */
  ///@{
  virtual void restart_iteration() override;
  virtual void moveRight() override;
  virtual void moveLeft() override;
  virtual bool canMoveRight() const override;
  virtual bool canMoveLeft() const override;
  virtual bool fellOff() const override;
  ///@}
  /** @name Data access */
  ///@{
  virtual const Monomial & currMonomial() const override;
  virtual const Prime_Field_Element & currCoeff() const override;
  ///@}
protected:
  /** @brief the polynomial we iterate on */
  const Constant_Polynomial * p;
  /** @brief current position in p&rsquo;s array */
  long i;
};

/**
  @class Constant_Polynomial
  @author John Perry
  @date 2015
  @brief A Constant_Polynomial is a polynomial that should not change.
  @ingroup polygroup

  @details We do allow on change at this level: to resort the monomials.
  We do not consider this a real change (it&rsquo;s still the same polynomial),
  and in any case even this is effectively unimplemented right now
  (unless you set up the polynomial with unsorted terms,
  in which case sorting them helps).
*/
class Constant_Polynomial : public Abstract_Polynomial
{
public:
  /** @name Construction */
  ///@{
  /**
    @param length how long @c this should be
    @param R parent ring
    @param mons array of Monomial to populate the terms
    @param coeffs array of coefficients to populate the terms, in same order as
        @p mons
    @param order monomial ordering used to first sort the polynomial;
        @c nullptr gives @c generic_grevlex_ptr
    @details We assume that mons and coeffs both have length n.
    We do not check that the monomials have the same number of variables,
    nor that the coefficients come from the same field,
    nor do we sort the monomials. The client needs to do this!
  */
  Constant_Polynomial(
      unsigned length,
      Polynomial_Ring & R,
      const Monomial *mons,
      const Prime_Field_Element *coeffs,
      const Monomial_Ordering * order = nullptr
  );
  /**
    @param length how long @c this should be
    @param R parent ring
    @param mons vector of pointers to Monomial to populate the terms
    @param coeffs array of coefficients to populate the terms, in same order as
        @p mons
    @param start where in the array to start populating the polynomial's
        coefficients
    @details This initializes a polynomial from @p mons and @p coeffs.
        It is assumed that @p coeffs may have some zero elements;
        we do not add these to the polynomial.
        It is also assumed that the monomials are already in correct order,
        so that no sorting is required.
        The monomial ordering is taken from the first monomial.
  */
  Constant_Polynomial(
      unsigned length,
      Polynomial_Ring & R,
      const vector<Monomial *> mons,
      const COEF_TYPE *coeffs,
      unsigned start
  );
  /**
    @param R parent ring
    @param mons list of Monomial to populate the terms
    @param coeffs list of coefficients to populate the terms, in same order as
        @p mons
    @param order monomial ordering used to first sort the polynomial;
        @c nullptr gives @c generic_grevlex_ptr
    @details We assume that mons and coeffs have the same length.
    We do not check that the monomials have the same number of variables,
    nor that the coefficients come from the same field,
    nor do we sort the monomials. The client needs to do this!
  */
  Constant_Polynomial(
      Polynomial_Ring & R,
      const list<Monomial> & mons,
      const list<Prime_Field_Element> & coeffs,
      const Monomial_Ordering * order = nullptr
  );
  /**
    @param length how long @c this should be
    @param R parent ring
    @param order monomial ordering
    @details This is a pretty useless constructor unless you know what you&rsquo;re doing.
    Unless you&rsquo;re a child of Constant_Polynomial,
    you don&rsquo;t know what you&rsquo;re doing.
    Even then, you <i>probably</i> don&rsquo;t know what you&rsquo;re doing.
  */
  Constant_Polynomial(
      unsigned length,
      Polynomial_Ring & R,
      const Monomial_Ordering * order = generic_grevlex_ptr
  );
  /**
    @brief Creates a constant polynomial copy of p.
  */
  explicit Constant_Polynomial(const Abstract_Polynomial & p);
  /** @brief from serial data */
  Constant_Polynomial(Polynomial_Ring &, const Monomial_Ordering *, uint64_t, uint64_t *);
  ///@}
  /** @name Destruction */
  ///@{
  ~Constant_Polynomial();
  ///@}
  /** @name Basic properties */
  ///@{
  /**
    @brief sets the ordering of monomials in this polynomial
    @param order the (new?) monomial ordering
    @param sort_anew set this to @c true to re-sort the polynomial
  */
  virtual void set_monomial_ordering(
      const Monomial_Ordering * order, bool sort_anew = true
  ) override;
  /**
    @brief sort by order
    @details currently uses insertion sort
  */
  virtual void sort_by_order() override;
  /** @return leading monomial -- call after sort_by_order()! */
  virtual Monomial & leading_monomial() const override;
  /** @return leading coefficient -- call after sort_by_order()! */
  virtual Prime_Field_Element leading_coefficient() const override;
  /** @return number of monomials */
  virtual unsigned length() const override;
  /** @return is this polynomial zero? */
  virtual bool is_zero() const override;
  ///@}
  /** @name Computation */
  ///@{
  /**
    @return zero constant polynomial: one entry, and it&rsquo;s zero!
  */
  virtual Constant_Polynomial * zero_polynomial() const override;
  /**
    @return multiple of @c this and @p u
    @param t a Monomial in the same ring
  */
  virtual Constant_Polynomial * monomial_multiple(const Monomial &t) const override;
  /**
    @return multiple of @c this and @p c
    @param c a scalar in the same prime field
  */
  virtual Constant_Polynomial * scalar_multiple(const Prime_Field_Element &c)
      const override;
  ///@}
  /** @name Iteration */
  ///@{
  /** @brief an iterator that poses no risk of modifying the polynomial */
  virtual Constant_Polynomial_Iterator * new_iterator() const override;
  /** @brief iterator to the first element */
  virtual Polynomial_Iterator * begin() const override;
  /** @brief iterator to the last element */
  virtual Polynomial_Iterator * end() const override;
  ///@}
  /** @name Serialization */
  ///@{
  /**
    @return an array of @c uint64_t containing the polynomial data
      and assigns the length of this array to @c size
    @param size modified to the number of elements in the array
    @details The data is stored in the form
        [ A1 M11 M12 ... M1n A2 M21 M22 ... M2n ... ]
      where Ai is the ith coefficient (lead is first) and Mij is the exponent
      of xj in the ith monomial
  */
  uint64_t * serialized(uint64_t & size);
  ///@}
  /** @brief to iterate without changing @c this */
  friend class Constant_Polynomial_Iterator;
  /** @brief to iterate and possibly change @c this */
  friend class Mutable_Constant_Polynomial_Iterator;
protected:
  /** @brief array of monomials, in one-to-one correspondence with <c>A</c> */
  Monomial *M;
  /**
    @brief array of coefficients, in one-to-one correspondence with <c>M</c>
  */
  Prime_Field_Element *A;
  /** @brief position <b>after</b> last monomial */
  unsigned m;
  /** @brief location of leading term in array (always farther left) */
  unsigned head;
};

/**
  @class Mutable_Constant_Polynomial_Iterator
  @author John Perry
  @date 2015
  @brief An iterator to modify monomials and coefficients of a
    Constant_Polynomial.
  @ingroup IteratorGroup

  @details Sometimes we want to be able to change particular monomials in a
  Constant_Polynomial.
  This a bit naughty, since it violates the spirit of a
  Constant_Polynomial.
  I guess I need a different name for the latter.
  This generally mimics a Constant_Polynomial_Iterator,
  but cannot inherit from it due to issues with the <c>const</c> keyword.
  This may reflect bad design.
*/
class Mutable_Constant_Polynomial_Iterator
    : public Mutable_Polynomial_Iterator
{
public:
  /** @name Construction */
  ///@{
  /**
    @param poly the polynomial @c this will iterate over
    @brief Creates an iterator for <c>poly</c> and starts at its leading term.
  */
  explicit Mutable_Constant_Polynomial_Iterator(Constant_Polynomial * poly);
  ///@}
  /** @name Destruction */
  ///@{
  ~Mutable_Constant_Polynomial_Iterator();
  ///@}
  /** @name Iteration */
  ///@{
  virtual void restart_iteration() override;
  virtual void moveRight() override;
  virtual void moveLeft() override;
  virtual bool fellOff() const override;
  ///@}
  /** @name Data access */
  ///@{
  virtual const Monomial & currMonomial() const override;
  virtual const Prime_Field_Element & currCoeff() const override;
  ///@}
  /** @name Data modification */
  ///@{
  virtual void set_currCoeff(const Prime_Field_Element & a) override;
  virtual void set_currMonomial(const Monomial & t) override;
  ///@}
protected:
  /** @brief the polynomial over which we iterate */
  Constant_Polynomial * p;
  /** @brief the monomial we currently point to */
  long long i;
};

#endif