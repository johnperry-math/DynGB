#ifndef __POLYNOMIAL_H_
#define __POLYNOMIAL_H_

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

/**
  Generic virtual class for polynomials.
  From this we derive others.
*/

#include <cstdlib>
#include <iostream>

using std::cout; using std::endl;

#include "fields.hpp"
#include "monomial.hpp"
#include "strategies.hpp"
#include "polynomial_ring.hpp"

extern Monomial_Ordering * generic_grevlex_ptr;

// forward declarations
class Polynomial;
class Polynomial_Iterator;
class Poly_Strategy_Data;

/**
  @defgroup polygroup Polynomials
  @brief classes related to the structure of polynomials
*/

/**
  @defgroup IteratorGroup Iterators over polynomials
  @brief classes related to iterating over polynomials
*/

/**
  @class Polynomial_Term
  @author John Perry
  @date 2016
  @brief convenience class to help iterate through polynomials
  @ingroup IteratorGroup
*/
class Polynomial_Term {
public:
  /** @name Construction */
  ///@{
  Polynomial_Term(const Monomial & m, const Prime_Field_Element & a) : t(m), c(a)
  { /* already initialized */ }
  ///@}
  /** @name Basic properties */
  ///@{
  /** @brief the monomial of this term */
  const Monomial & monomial() { return t; }
  /** @brief the coefficient for this term */
  const Prime_Field_Element & coefficient() { return c; }
  ///@}
protected:
  /** @brief the monomial part of this term */
  const Monomial & t;
  /** @brief the coefficient of this term */
  const Prime_Field_Element & c;
};

/**
  @class Abstract_Polynomial
  @author John Perry
  @date 2015
  @brief The general class of a polynomial
  @ingroup polygroup

  @details This class encapsulates the minimum feature set necessary for our application.
  Naturally, it must be possible to identify a leading monomial and coefficient.
  Every polynomial must also be able to produce an iterator;
  indicate its its length (number of monomials), and its number of variables;
  describe its ground field; and produce a zero polynomial of the same type
  (this has its uses).

  @warning Monomials should have the same number of variables, and
    coefficients should all come from the same field. Behavior is undefined if
    these assumptions are violated. From an algebraic point of view, it
    doesn&rsquo;t make much sense to violate them, anyway;
    (in/pro)ject into a different ring if you want to screw around like this.
*/
class Abstract_Polynomial
{
public:
  /** @name Construction */
  ///@{
  /**
    @param ring the polynomial ring in which this polynomial will reside
    @param ordering the monomial ordering that first sorts the monomials
  */
  Abstract_Polynomial(Polynomial_Ring & ring, const Monomial_Ordering * ordering)
    : R(ring)
  { }
  ///@}
  /** @name Destruction */
  ///@{
  /**
    @brief deletes the strategy, if there is one
  */
  virtual ~Abstract_Polynomial() { if (strat != nullptr) delete strat; };
  ///@}
  /** @name Basic properties */
  ///@{
  /** @brief ring in which this polynomial resides */
  Polynomial_Ring & base_ring() const;
  /** @brief ground field -- all coefficients should be in this field */
  const Prime_Field & ground_field() const;
  /**
    @brief number of variables -- all monomials should agree with this
      (though it is never tested by the class)
  */
  unsigned number_of_variables() const;
  /** @brief reports leading monomial&rsquo;s monomial ordering */
  const Monomial_Ordering * monomial_ordering() const {
    return leading_monomial().monomial_ordering();
  }
  /** @brief leading monomial -- call after sort_by_order()! */
  virtual Monomial & leading_monomial() const = 0;
  /** @brief leading coefficient -- call after sort_by_order()! */
  virtual Prime_Field_Element leading_coefficient() const = 0;
  /** @brief number of monomials */
  virtual unsigned length() const = 0;
  /** @brief is this polynomial zero? */
  virtual bool is_zero() const = 0;
  /** @brief can <c>this</c> reduce <c>other</c>? */
  virtual bool can_reduce(Abstract_Polynomial &other) const;
  /** @brief strategy related information */
  virtual Poly_Strategy_Data * strategy() const { return strat; }
  /** @brief maximum sum of exponents for any monomial */
  virtual DEG_TYPE standard_degree() const;
  /**
    @return largest weighted sum of exponents for any monomial
    @param w weights to use for weighted degree
    @details Equivalent to standard_degree() if <c>w == nullptr</c>.
  */
  virtual DEG_TYPE weighted_degree(const WT_TYPE * w = nullptr) const;
  ///@}
  /** @name Computation */
  ///@{
  /** @brief new zero polynomial of this same type */
  virtual Abstract_Polynomial * zero_polynomial() const = 0;
  /** @brief multiple of this and @f$u@f$ */
  virtual Abstract_Polynomial * monomial_multiple(const Monomial &) const = 0;
  /** @brief multiple of this and @f$c@f$ */
  virtual Abstract_Polynomial * scalar_multiple(const Prime_Field_Element &)
      const = 0;
  /** @brief sets the polynomial&rsquo;s strategy to @c psd */
  void set_strategy(Poly_Strategy_Data * psd);
  ///@}
  /** @name Ordering monomials */
  ///@{
  /**
    @brief set the monomial ordering and sort the polynomials
      (optionally, but by default)
    @param order new monomial ordering
    @param sort_anew whether to sort the polynomial anew
    @warning In most cases you will want to sort anew \e immediately
      after setting the ordering. Otherwise, the monomials may be in the
      wrong order! That is therefore the default behavior of this function,
      but in case you don&rsquo;t want to sort, that option is provided.
  */
  virtual void set_monomial_ordering(
      const Monomial_Ordering * order, bool sort_anew = true
  ) = 0;
  /**
    @brief sort according to the leading monomial&rsquo;s ordering

    Note that it makes sense to sort even Constant_Polynomial&rsquo;s,
    as it is not the polynomial that changes, only the order of its monomials.
  */
  virtual void sort_by_order() = 0;
  ///@}
  /** @name Iteration */
  ///@{
  /** @brief An iterator that poses no risk of modifying the polynomial */
  virtual Polynomial_Iterator * new_iterator() const = 0;
  /** @brief returns an iterator to the polynomial&rsquo;s leading monomial */
  virtual Polynomial_Iterator * begin() const = 0;
  /** @brief iterator to last monomial */
  virtual Polynomial_Iterator * end() const = 0;
  ///@}
  /** @name I/O */
  ///@{
  /** @brief output */
  friend ostream & operator << (ostream & os,
      const Abstract_Polynomial & p);
  /** @brief prints the polynomial to @p os */
  virtual void print(ostream & os=cout) const;
  /** @brief prints the polynomial to @p os, followed by a carriage return */
  virtual void println(ostream & os=cout) const;
  /** @brief prints the polynomial to @c cout, followed by a carriage return */
  virtual void printlncout() const { println(); }
  ///@}
protected:
  /** @brief data about polynomial ring */
  Polynomial_Ring & R;
  /** @brief data for computational strategies */
  Poly_Strategy_Data * strat = nullptr;
};

/**
  @class Polynomial_Iterator
  @author John Perry
  @date 2015
  @brief Used to iterate through a polynomial.
  @ingroup IteratorGroup

  @details This class encapsulates the basic functionality needed to iterate through
  the monomials of a polynomial. A Polynomial_Iterator cannot modify
  a polynomial, however; to do that, you need a
  Mutable_Polynomial_Iterator.
*/
class Polynomial_Iterator
{
public:
  /** @name Destruction */
  ///@{
  /** @brief needed to avoid undefined behavior when disposing */
  virtual ~Polynomial_Iterator() = 0;
  ///@}
  /** @name Iteration */
  ///@{
  /** @brief This should move the iterator to the leading term. */
  virtual void restart_iteration() = 0;
  /** @brief Can this iterator start at the tail? */
  virtual bool canStartAtTail() { return false; }
  /**
    @brief This should move the iterator to the last term (if this is possible).
  */
  virtual void tail() { restart_iteration(); }
  /** @brief Moves right in the polynomial, to the next smaller monomial. */
  virtual void moveRight() = 0;
  /** @brief Moves right in the polynomial, to the next smaller monomial. */
  const Polynomial_Iterator & operator++() { moveRight(); return *this; }
  /** @brief Can this iterator move right, or would it fall off? */
  virtual bool canMoveRight() const = 0;
  /** @brief Moves left in the polynomial, to the next larger monomial. */
  virtual void moveLeft() = 0;
  /** @return Can this iterator move left, or would it fall off? */
  virtual bool canMoveLeft() const = 0;
  /**
    @return true iff the iterator no longer points to a valid monomial.

    This is NOT the same as pointing to a monomial with coefficient zero;
    this is true when the iterator would probably report inaccurate data.
  */
  virtual bool fellOff() const = 0;
  ///@}
  /** @name Data access */
  ///@{
  /** @brief Reports the polynomial on which @c this is iterating. */
  virtual const Abstract_Polynomial * my_poly() const { return p_base; }
  /** @brief Reports the monomial at the current position. */
  virtual const Monomial & currMonomial() const = 0;
  /** @brief Reports the coefficient at the current position. */
  virtual const Prime_Field_Element & currCoeff() const = 0;
  /** @brief C++11 iteration */
  const Polynomial_Term operator *() const {
    return Polynomial_Term(currMonomial(), currCoeff());
  }
  virtual bool operator !=(const Polynomial_Iterator & other) const {
    return my_poly() != other.my_poly() or currMonomial() != other.currMonomial();
  }
  ///@}
protected:
  /** @brief the polynomial @c this points to */
  const Abstract_Polynomial * p_base;
};

/**
  @class Mutable_Polynomial_Iterator
  @author John Perry
  @date 2015
  @brief A Mutable_Polynomial_Iterator allows one to modify the
    terms of a polynomial.
  @ingroup IteratorGroup
*/
class Mutable_Polynomial_Iterator : public Polynomial_Iterator
{
public:
  /** @name Data modification */
  ///@{
  /** @brief change coefficient in current position */
  virtual void set_currCoeff(const Prime_Field_Element &) = 0;
  /** @brief change monomial in current position */
  virtual void set_currMonomial(const Monomial &) = 0;
  ///@}
};

/**
  @class Mutable_Polynomial
  @author John Perry
  @date 2015
  @brief Polynomials that need arithmetic typically descend from this class.
  @ingroup polygroup
  
  @details This class extends Abstract_Polynomial to allow for basic
  arithmetic of a polynomial.
*/
class Mutable_Polynomial : public Abstract_Polynomial
{
public:
  /** @name Construction */
  ///@{
  /** @brief constructor */
  Mutable_Polynomial(
      Polynomial_Ring & R,
      const Monomial_Ordering * ordering = generic_grevlex_ptr
  ) : Abstract_Polynomial(R, ordering) { }
  ///@}
  /** @name Destruction */
  ///@{
  /** @brief destructor */
  virtual ~Mutable_Polynomial() = 0;
  ///@}
  /** @name Computation */
  ///@{
  /** @brief zero polynomial of this type */
  virtual Mutable_Polynomial * zero_polynomial() const = 0;
  /** @brief add another polynomial */
  virtual Mutable_Polynomial & operator +=(const Abstract_Polynomial &) = 0;
  /** @brief subtract another polynomial */
  virtual Mutable_Polynomial & operator -=(const Abstract_Polynomial &) = 0;
  /** @brief add monomial multiple of other */
  virtual void add_polynomial_multiple(const Prime_Field_Element &,
                                       const Monomial &,
                                       const Abstract_Polynomial &,
                                       bool subtract = false
  ) = 0;
  /** @brief multiply by scalar */
  virtual void multiply_by_scalar(const Prime_Field_Element &a);
  /** @brief multiply by monomial */
  virtual void multiply_by_monomial(const Monomial &t);
  /** @brief reduce by @f$p@f$ until no further reduction possible */
  virtual void reduce_by(const Abstract_Polynomial &p);
  // utility
  /** @brief multiple of this and u */
  virtual Mutable_Polynomial * monomial_multiple(const Monomial &) const = 0;
  /** @brief multiple of this and c */
  virtual Mutable_Polynomial * scalar_multiple(const Prime_Field_Element &)
      const = 0;
  ///@}
  /** @name Iteration */
  ///@{
  /** @brief An iterator that may modify the current position */
  virtual Mutable_Polynomial_Iterator * new_mutable_iterator() = 0;
  ///@}
  /** @name Modification */
  ///@{
  /** @brief Remove and return the head */
  virtual Mutable_Polynomial * detach_head() = 0;
  /** @brief Attach a new monomial to the tail -- check that it belongs at tail! */
  virtual void add_last(const Prime_Field_Element &, const Monomial &) = 0;
  ///@}
};

#endif