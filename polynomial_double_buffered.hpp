#ifndef __DBL_BUF_POLYNOMIAL_H_
#define __DBL_BUF_POLYNOMIAL_H_

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
using std::cout; using std::endl;

#include "polynomial.hpp"
#include "polynomial_linked_list.hpp"

extern Monomial_Ordering * generic_grevlex_ptr;

class Double_Buffered_Polynomial;

/**
  @class DB_Polynomial_Iterator
  @author John Perry
  @date 2015
  @brief Iterator over double-buffered polynomials.
  \ingroup IteratorGroup
*/
class DB_Polynomial_Iterator : public Mutable_Polynomial_Iterator {

public:

  /** @name Construction */
  ///@{
  DB_Polynomial_Iterator(const Double_Buffered_Polynomial * f, bool at_end=false);
  ///@}

  /** @name Destruction */
  ///@{
  ~DB_Polynomial_Iterator();
  ///@}

  /** @name Iteration */
  ///@{
  virtual void restart_iteration() override;

  virtual void moveRight() override;

  virtual void moveLeft() override;

  virtual bool fellOff() const override;
  virtual bool canMoveLeft() const override;
  virtual bool canMoveRight() const override;
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

  /** @brief pointer to polynomial */
  const Double_Buffered_Polynomial * p;
  /** @brief pointer to coefficients */
  Prime_Field_Element * A;
  /** @brief pointer to monomials */
  Monomial * T;
  /** @brief position of first monomial */
  unsigned head;
  /** @brief position of last monomial */
  unsigned tail;
  /** @brief current position in the list */
  long long current_position;

};

/**
  @class Double_Buffered_Polynomial
  @author John Perry
  @date 2015
  @brief Polynomials implemented using double buffers.
  \ingroup polygroup

  A double-buffered polynomial maintains at all times two arrays to store
  its terms. (Technically, it retains four arrays: two for the monomials,
  and two for the coefficients.) Any operation that might change the
  \e length of the polynomials reads the data in one buffer and writes
  the result to the other buffer. The goal of this approach is to avoid the
  penalties associated with allocating, deallocating, and traversing the
  nodes of a linked list.
*/

class Double_Buffered_Polynomial : public Mutable_Polynomial {

public:

  /** @name Construction */
  ///@{

  Double_Buffered_Polynomial(
    Polynomial_Ring & R,
    Monomial_Ordering * order = generic_grevlex_ptr
  );

  explicit Double_Buffered_Polynomial(Abstract_Polynomial const & p);

  ///@}
  /** @name Destruction */
  ///@{

  ~Double_Buffered_Polynomial();

  ///@}

  /** @name Basic properties */
  ///@{

  virtual Monomial & leading_monomial() const override;

  virtual Prime_Field_Element leading_coefficient() const override;

  virtual unsigned length() const override;

  virtual bool is_zero() const override;

  virtual bool can_reduce(Abstract_Polynomial &other) const override;

  virtual Double_Buffered_Polynomial * zero_polynomial() const override;

  virtual void set_monomial_ordering(
      const Monomial_Ordering * order, bool sort_anew=true
  ) override;

  ///@}
  /** @name Computation */
  ///@{

  virtual Double_Buffered_Polynomial * monomial_multiple(
      const Monomial & t
  ) const override;

  virtual Double_Buffered_Polynomial * scalar_multiple(
      const Prime_Field_Element & c)
  const override;

  virtual Mutable_Polynomial & operator += (const Abstract_Polynomial & p) override;

  virtual Mutable_Polynomial & operator -= (const Abstract_Polynomial & p) override;

  virtual void add_polynomial_multiple(
      const Prime_Field_Element & a,
      const Monomial & u,
      const Abstract_Polynomial & p,
      bool subtract
    ) override;

  virtual void sort_by_order() override;

  ///@}
  /** @name Iteration */
  ///@{

  virtual DB_Polynomial_Iterator * new_iterator() const override;

  virtual DB_Polynomial_Iterator * new_mutable_iterator() override;

  virtual Polynomial_Iterator * begin() const override;
  virtual Polynomial_Iterator * end() const override;

  ///@}
  /** @name Modification */
  ///@{

  virtual void add_last(const Prime_Field_Element & a, const Monomial & t) override;

  virtual Polynomial_Linked_List * detach_head() override;

  ///@}

protected:

  /** @name Buffering */
  ///@{

  /**
    @brief @c true iff buffer @p b has space to hold @p n elements;
    expand other buffer if not.
  */
  inline bool test_buffer(unsigned b, unsigned n) { return sizes[b] >= n; }

  /**
    @brief Expand buffer @p b to hold @f$ 2n @f$ elements.
    @param b the buffer to expand
    @param n number of monomials you&rsquo;d like it to hold now
    @details This does \e not test to see if the space is already available.
  */
  void expand_buffer(unsigned b, unsigned n) {
    if (mons[b] != nullptr) {
      for (unsigned k = 0; k < sizes[b]; ++k)
        mons[b][k].deinitialize();
      free(mons[b]);
      free(coeffs[b]);
    }
    unsigned new_size = 2*n;
    sizes[b] = new_size;
    mons[b] = static_cast<Monomial *>(malloc(new_size*sizeof(Monomial)));
    for (unsigned k = 0; k < new_size; ++k) {
      mons[b][k].initialize_exponents(number_of_variables());
      mons[b][k].clear_ordering_data();
    }
    coeffs[b] = static_cast<Prime_Field_Element *>(
      malloc(new_size*sizeof(Prime_Field_Element))
    );
  }

  ///@}

  /** @brief to iterate over @c this and possibly change it */
  friend class DB_Polynomial_Iterator;

private:

  Monomial * mons[2];
  Prime_Field_Element * coeffs[2];
  unsigned sizes[2];
  unsigned active_buffer;
  unsigned head;
  unsigned tail;

};

#endif
