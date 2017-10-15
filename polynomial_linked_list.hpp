#ifndef __POLYNOMIAL_LINKED_LIST_H_
#define __POLYNOMIAL_LINKED_LIST_H_

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
  Doubly linked list representation of a polynomial.
*/

#include "fields.hpp"
#include "monomial.hpp"
#include "polynomial.hpp"
#include "polynomial_ring.hpp"

extern Monomial_Ordering * generic_grevlex_ptr;

/**
  @class Monomial_Node
  @author John Perry
  @date 2015
  @brief Tool for Polynomial_Linked_List.
  @ingroup polygroup

  @details Each node in a linked list polynomial contains both the coefficient and the
  monomial. 
*/
class Monomial_Node
{
public:
  /** @name Construction */
  ///@{
  /**
    @brief Monomial u with coefficient a. Both are copied, and can be deleted.
  */
  Monomial_Node(const Prime_Field_Element & a, const Monomial & u);
  /**
    @brief Monomial u (copied) with coefficient 1.
  */
  Monomial_Node(Prime_Field & F, const Monomial & u);
  ///@}
  /** @name Destruction */
  ///@{
  ///@}
  /** @name Basic properties */
  ///@{
  /**
    @brief This term&rsquo;s monomial, or power product.
      The coefficient is not included.
  */
  Monomial & monomial();
  /** @brief This term&rsquo;s coefficient. */
  Prime_Field_Element & coefficient();
  ///@}
  /** @name Iteration */
  ///@{
  friend class LLPolynomial_Iterator;
  friend class Polynomial_Linked_List;
  ///@}
  /** @name Memory management */
  ///@{
  /** @brief requests memory form Monomial_Node's Grading_Order_Data_Allocator */
  void * operator new(size_t);
  /** @brief returns data to Monomial_Node's Grading_Order_Data_Allocator */
  void operator delete(void *);
  ///@}
protected:
  /** @brief the monomial in this node */
  Monomial t;
  /** @brief the monomial&rsquo;s coefficient */
  Prime_Field_Element c;
  /** @brief for linking */
  Monomial_Node *next;
  /** @brief for linking */
  Monomial_Node *prev;
};

// forward declaration
class Polynomial_Linked_List;

/**
  @class LLPolynomial_Iterator
  @author John Perry
  @date 2015
  @brief Iterator over linked list polynomials.
  @ingroup IteratorGroup
*/
class LLPolynomial_Iterator : public Mutable_Polynomial_Iterator
{
public:
  /** @name Construction */
  ///@{
  /**
    @brief Initializes at the leading monomial.
  */
  LLPolynomial_Iterator(Polynomial_Linked_List * poly, bool at_end = false);
  /**
    @brief Initializes at the leading monomial.
  */
  LLPolynomial_Iterator(const Polynomial_Linked_List * poly, bool at_end = false);
  ///@}
  /** @name Iteration */
  ///@{
  /**
    @brief Initializes at the leading monomial.
  */
  virtual void restart_iteration() override;
  /**
    @brief Returns the monomial the iterator currently points to. 
  */
  /**
    @brief Moves the iterator right: to the next smaller monomial.
  */
  virtual void moveRight() override { iter_curr = iter_curr->next; }
  /**
    @brief Moves the iterator left: to the next larger monomial.
  */
  virtual void moveLeft() override { iter_curr = iter_curr->prev; }
  /** @brief Can this iterator move left, or would it fall off? */
  virtual bool canMoveLeft() const override;
  /** @brief Can this iterator move right, or would it fall off? */
  virtual bool canMoveRight() const override;
  /**
    @brief true iff the iterator no longer points to a valid monomial.
  */
  virtual bool fellOff() const override;
  ///@}
  /** @name Data access */
  ///@{
  virtual const Monomial & currMonomial() const override;
  /**
    @brief Returns the coefficient of the monomial the iterator currently points to.
  */
  virtual const Prime_Field_Element & currCoeff() const override;
  ///@}
  /** @name Data modification */
  ///@{
  /** @brief change coefficient in current position */
  virtual void set_currCoeff(const Prime_Field_Element & a) override;
  /** @brief change monomial in current position */
  virtual void set_currMonomial(const Monomial & t) override;
  ///@
protected:
  /** @brief the polynomial over which we iterate */
  Polynomial_Linked_List *p;
  /** @brief the node at which we have stopped */
  Monomial_Node *iter_curr;
};

/**
  @class  Polynomial_Linked_List
  @author John Perry
  @date 2015
  @brief Polynomials represented as a doubly linked list.
  @ingroup polygroup
*/
class Polynomial_Linked_List : public Mutable_Polynomial
{
public:
  /** @name Construction */
  ///@{
  /** @brief initialize to zero */
  Polynomial_Linked_List(
      Polynomial_Ring & R,
      const Monomial_Ordering * order = generic_grevlex_ptr
  );
  /** @brief initialize to monomial; monomial is copied */
  Polynomial_Linked_List(
      Polynomial_Ring & R,
      const Monomial & t,
      const Monomial_Ordering * order = nullptr
  );
  /** @brief initialize to monomial and coefficient; monomial is copied */
  Polynomial_Linked_List(
      Polynomial_Ring & R,
      const Prime_Field_Element & c, const Monomial & t,
      const Monomial_Ordering * order = nullptr
  );
  /** @brief initialize to given monomial node: nothing is copied */
  Polynomial_Linked_List(
      Polynomial_Ring & R,
      Monomial_Node * node,
      const Monomial_Ordering * order = nullptr
  );
  /** @brief copy constructor: deep copy of monomials */
  Polynomial_Linked_List(const Polynomial_Linked_List & other);
  /** @brief constructor from abstract polynomial: deep copy of monomials */
  explicit Polynomial_Linked_List(const Abstract_Polynomial & p);
  ///@}
  /** @name Destruction */
  ///@{
  /** @brief deletes all monomial nodes */
  virtual ~Polynomial_Linked_List();
  ///@}
  /** @name Basic properties */
  ///@{
  /**
    @brief true iff the first node in the list is <c>nullptr</c>
      or has zero coeff
  */
  virtual bool is_zero() const override;
  /**
    @brief Returns the leading monomial &mdash;
      call sort_by_order() first!
  */
  virtual Monomial & leading_monomial() const override;
  /**
    @brief Returns the leading coefficient &mdash; call sort_by_order first!
  */
  virtual Prime_Field_Element leading_coefficient() const override;
  /** @brief Returns the number of polynomials in the list. */
  virtual unsigned length() const override;
  virtual void set_monomial_ordering(
      const Monomial_Ordering * ord, bool sort_anew = true
  ) override;
  ///@}
  /** @name Computation */
  ///@{
  /**
    @brief Returns as simple a zero polynomial as I can muster,
    short of this being <c>nullptr</c>.
  */
  virtual Polynomial_Linked_List * zero_polynomial() const override;
  /**
    @brief Returns a new polynomial whose value is @f$\textit{this}\times u@f$.
  */
  virtual Polynomial_Linked_List * monomial_multiple(const Monomial &u) const override;
  /**
    @brief Returns a new polynomial whose value is @f$\textit{this}\times c@f$.
  */
  virtual Polynomial_Linked_List * scalar_multiple(const Prime_Field_Element &c)
      const override;
  /**
    @brief Adds <c>other</c> to <c>this</c>, and returns the result.
  */
  virtual Polynomial_Linked_List & operator +=(const Abstract_Polynomial &other) override;
  /**
    @brief Subtracts <c>other</c> from <c>this</c>, and returns the result.
  */
  virtual Polynomial_Linked_List & operator -=(const Abstract_Polynomial &other) override;
  /**
    @brief "Fast" addition of @f$atq@f$ to <c>this</c>.
    @details If <c>subtract==true</c>, subtract instead.
    @param a coefficient of the term to multiply to @p q
    @param t Monomial of the term to multiply to @p q
    @param q polynomial to add to @c this, after multiplying @p q by @f$ at @f$
    @param subtract whether to add instead of subtract
  */
  virtual void add_polynomial_multiple(
      const Prime_Field_Element &a,
      const Monomial &t,
      const Abstract_Polynomial &q,
      bool subtract
  ) override;
  /**
    @brief Sort by specified weight order.
  */
  virtual void sort_by_order() override;
  ///@}
  /** @name Iteration */
  ///@{
  /** @brief an iterator that poses no risk of modifying the polynomial */
  virtual LLPolynomial_Iterator * new_iterator() const override;
  /** @brief an iterator that poses no risk of modifying the polynomial */
  virtual Polynomial_Iterator * begin() const override;
  /** @brief an iterator that poses no risk of modifying the polynomial */
  virtual Polynomial_Iterator * end() const override;
  /** @brief An iterator that may modify the current position */
  virtual LLPolynomial_Iterator * new_mutable_iterator() override;
  ///@}
  /** @name Modification */
  ///@{
  /** @brief Detach and return leading term */
  virtual Polynomial_Linked_List * detach_head() override;
  /** @brief Add this monomial as the last leading term */
  virtual void add_last(const Prime_Field_Element & c, const Monomial & t) override;
  ///@}
  /** @brief to iterate over @c this and possibly change it */
  friend class LLPolynomial_Iterator;
protected:
  /** @brief first monomial in the list, returned by leading_monomial() */
  Monomial_Node *head;
};

#endif