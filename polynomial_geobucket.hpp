#ifndef __POLYNOMIAL_GEOBUCKET_H_
#define __POLYNOMIAL_GEOBUCKET_H_

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
#include <stdexcept>

using std::cout; using std::endl;
using std::runtime_error;

#include "fields.hpp"
#include "monomial.hpp"
#include "polynomial.hpp"
#include "polynomial_array.hpp"
#include "polynomial_linked_list.hpp"
#include "polynomial_ring.hpp"

extern Monomial_Ordering * generic_grevlex_ptr;

#define NUM_BUCKETS  32
#define BUCKET_BASE   4
#define BUCKET_SHIFT  2

// forward declaration
class Polynomial_Geobucket;

/**
  @class Geobucket_Iterator
  @author John Perry
  @date 2015
  @brief Iterates through polynomials using a geobucket representation.
  See Mutable_Polynomial_Iterator for details on methods.
  @ingroup IteratorGroup

  @details Implementation based on @cite YanGeobuckets.
  @warning A polynomial in geobucket representation may not be fully simplified.
    Clients must not expect the monomials to be in any sort of order,
    and multiple instances of a monomial are possible.
    Geobuckets collapse only at canonicalization.
*/
class Geobucket_Iterator : public Mutable_Polynomial_Iterator {
public:
  /** @name Construction */
  ///@{
  Geobucket_Iterator(const Polynomial_Geobucket *, bool at_end = false);
  Geobucket_Iterator(Polynomial_Geobucket *, bool at_end = false);
  ///@}
  /** @name Destruction */
  ///@{
  ~Geobucket_Iterator();
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
  /** @name Data modification */
  ///@{
  virtual void set_currCoeff(const Prime_Field_Element & a) override;
  virtual void set_currMonomial(const Monomial & t) override;
  ///@}
protected:
  /**
    @brief to save myself the hassle of
    <c>static_cast<const Polynomial_Geobucket *>(p)</c> every time I turn around
  */
  Polynomial_Geobucket * g;
  /** @brief the bucket number at which we&rsquo;re stopped */
  unsigned bucket_number;
  /** @brief an iterator for the current bucket */
  LLPolynomial_Iterator * pi;
};

/**
  @class Polynomial_Geobucket
  @author John Perry
  @date 2015
  @brief Implementation of geobuckets
  @ingroup polygroup

  @warning After any operation that might modify the leading term
    (such as adding, subtracting polynomials -- <i>not,</i> however,
    scalar or monomial multiplication) you will need to call
    recompute_leading_monomial().
  @warning Exceeding the number of buckets will raise an exception.
    This will occur when a partial sum&rsquo;s size exceeds <c>NUM_BUCKETS</c>.
*/
class Polynomial_Geobucket : public Mutable_Polynomial {
public:
  /** @name Construction */
  ///@{
  /**
    @brief Initializes a polynomial with the given number of variables,
      over the given field.
  */
  Polynomial_Geobucket(Polynomial_Ring & R,
                       Monomial_Ordering * order = generic_grevlex_ptr
  );
  /**
    @brief Initializes a geobucket that is a copy of @f$p@f$.
  */
  Polynomial_Geobucket(Abstract_Polynomial & p);
  ///@}
  /** @name Destruction */
  ///@{
  /**
    @brief Deallocates buckets, if they exist.
  */
  ~Polynomial_Geobucket();
  ///@}
  /** @name Basic properties */
  ///@{
  /**
    @brief sorts each geobucket
    @warning This performs no simplification between buckets.
  */
  virtual void sort_by_order() override;
  /**
    @brief sets the monomial ordering for each bucket
    @warning This performs no simplification bewteen buckets.
    @param order the new monomial ordering
    @param sort_anew whether to re-sort the polynomials
  */
  virtual void set_monomial_ordering(
      const Monomial_Ordering * order, bool sort_anew = true
  ) override;
  /**
    @brief returns the leading monomial
    @return the leading monomial
    @warning Does not recompute the leading monomial first.
      Make sure the leading monomial exists and is valid.
    @warning Descendants that perform any arithemtic
      that might change the leading monomial
      should invoke recompute_leading_monomial() before this.
      The class&rsquo;s own functions for arithmetic do this automatically
      (or should).
  */
  virtual Monomial & leading_monomial() const override;
  /**
    @brief returns the leading coefficient
    @return the leading coefficient
    @warning Does not recompute the leading monomial first.
      Make sure the leading monomial exists and is valid.
    @warning Descendants that perform any arithemtic
      that might change the leading monomial
      should invoke recompute_leading_monomial() before this.
      The class&rsquo;s own functions for arithmetic do this automatically
      (or should).
  */
  virtual Prime_Field_Element leading_coefficient() const override;
  /**
    @brief how long is this polynomial?
    @return number of monomials in this polynomial
    @warning This will not include potential simplification of monomials.
      A geobucket does not simplify non-leading terms until canonicalization.
  */
  virtual unsigned length() const override;
  /**
    @brief creates and returns a geobucket initialized to zero
    @return a geobucket initialized to zero
  */
  virtual Polynomial_Linked_List * zero_polynomial() const override;
  /**
    @brief is this polynomial zero?
    @return <c>true</c> iff all the buckets are <c>nullptr</c> or themselves zero
    @see Fixed_Length_Polynomial
  */
  virtual bool is_zero() const override;
  /**
    @brief Whether <c>this</c> can reduce <c>other</c>.
      A geobucket should not generally be used to reduce other polynomials,
      so avoid this like the plague.
  */
  virtual bool can_reduce(Abstract_Polynomial & other) const override;
  virtual Geobucket_Iterator * new_iterator() const override;
  virtual Polynomial_Iterator * begin() const override;
  virtual Polynomial_Iterator * end() const override;
  virtual Geobucket_Iterator * new_mutable_iterator() override;
  ///@}
  /** @brief to iterate over @c this and possibly change it */
  friend class Geobucket_Iterator;
  /** @name Computation */
  ///@{
  /**
    @brief You will need to call this after every operation
    that might modify the leading term.

    @details The assumption here is that the first bucket no longer contains
    a valid leading term, so we skip it and try to extract from a later one.
    If the other buckets are all zero, we indicate this.
  */
  void recompute_leading_monomial();
  /**
    @brief Returns @f$\texttt{this}\times t@f$.
  */
  virtual Polynomial_Geobucket * monomial_multiple(const Monomial &t) const override;
  /**
    @brief Returns @f$\texttt{this}\times a@f$.
  */
  virtual Polynomial_Geobucket * scalar_multiple(const Prime_Field_Element & a)
      const override;
  /**
    @brief Adds @f$g@f$ to <c>this</c>. Recomputes leading monomial.
  */
  virtual Polynomial_Geobucket & operator +=(const Abstract_Polynomial & g) override;
  /**
    @brief Subtracts @f$g@f$ from <c>this</c>. Recomputes leading monomial.
  */
  virtual Polynomial_Geobucket & operator -=(const Abstract_Polynomial & g) override;
  /**
    @brief Adds @f$bug@f$ to <c>this</c>. Recomputes leading monomial.
  */
  virtual void add_polynomial_multiple(
       const Prime_Field_Element & b, const Monomial & u,
       const Abstract_Polynomial & g, bool subtract=false
      ) override;
  /**
    @brief Detaches the head and recomputes leading monomial.
  */
  virtual Polynomial_Linked_List * detach_head() override;
  /**
    @brief Adds @f$at@f$ as a monomial of <c>this</c>.
      (Not necessarily the last!)
    @param a coefficient of the term to add
    @param t Monomial of the term to add
    @warning In the case of a geobucket, add_last()
    doesn&rsquo;t really make sense
    (the last monomial might not be in order, after all) so this should not
    be used. As in the usual case, we do assume it belongs at the tail;
    we simply add it to a tail bucket (i.e., not <c>bucket[0]</c>).
  */
  virtual void add_last(const Prime_Field_Element & a, const Monomial & t) override;
  /**
    @brief returns a copy of <c>this</c> in a simplified linear form
    @param constant_result whether you want a Constant_Polynomial
    @details If @p constant_result is @c true, the result is a
      @c Constant_Polynomial.
      Otherwise, the result is a @c Polynomial_Linked_List.
      Choose the latter if you anticipate further reduction,
      such as tail reduction.
      (Geobuckets by default do not reduce lower-order terms.)
    @return a copy of <c>this</c> in a simplified linear form
  */
  Abstract_Polynomial * canonicalize(bool constant_result=false);
  ///@}
  /** @name I/O */
  ///@{
  /**
    @brief prints the @f$i@f$th bucket
    @param i which bucket to print
    @param os where to print the bucket
    @details In the context of geobuckets, this makes more sense than printing
    the polynomial. Practically, printing the polynomial does the same thing
    as printing the contents of the buckets, but the more general function
    can be misleading on account of the un-simplified nature of the polynomial.
  */
  void print(unsigned i, ostream & os=cout) const;
  /**
    @brief prints the polynomial with an explicitly bucket form
    @param os where to print the polynomial
    @details Prints the polynomial in the form
    @f[(b_1) + (b_2) + \cdots + (b_\textrm{last}).@f]
    Parentheses separate buckets, whose sums have not (yet) been simplified.
    Uninitiated and nonzero buckets are not printed.
  */
  virtual void print(ostream & os=cout) const override;
  inline friend ostream & operator << (ostream &, const Polynomial_Geobucket &);
  ///@}
protected:
  /**
    @brief Log-length of @f$i@f$, used to select a bucket for a polynomial
      of length @f$i@f$.
    @param i number of terms we need space for
    @return the smallest power of 2 larger than i
  */
  inline unsigned lglen(unsigned i) {
    unsigned result = 1;
    while ((i >>= BUCKET_SHIFT)) { ++result; }
    return result;
  }
  /**
    @brief Array of ptrs to linked list polys; most initialized to <c>nullptr</c>.
  */
  Polynomial_Linked_List ** buckets;
};

#endif
