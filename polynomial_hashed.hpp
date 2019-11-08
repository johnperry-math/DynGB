#ifndef __POLYNOMIAL_HASHED_H_
#define __POLYNOMIAL_HASHED_H_

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

#include "f4_hash.hpp"
#include "monomial.hpp"
#include "polynomial.hpp"

class F4_Hash;

class Polynomial_Hashed;

class Hashed_Polynomial_Iterator : public Polynomial_Iterator {
public:
  /** @name Construction */
  ///@{
  Hashed_Polynomial_Iterator(const Polynomial_Hashed *p, bool at_end=false);
  ///@}
  /** @name Destruction */
  ///@{
  ~Hashed_Polynomial_Iterator();
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
  const Polynomial_Hashed * p;
  long i;
  bool started_at_beginning;
};

class Polynomial_Hashed : public Abstract_Polynomial {
public:
  /** @name Construction */
  ///@{
  Polynomial_Hashed(
      Polynomial_Ring & R,
      vector< Monomial * > & monomials,
      F4_Hash & monomial_hash,
      const vector< pair< unsigned, COEF_TYPE > > & from_row,
      const vector< Monomial * > & row_monomials,
      Monomial_Ordering * mord
  );
  Polynomial_Hashed(
      Abstract_Polynomial &, vector< Monomial * > &, F4_Hash &,
      Monomial_Ordering *
  );
  ///@}
  /** @name Destruction */
  ///@{
  ~Polynomial_Hashed();
  ///@}
  /** @name Basic properties */
  ///@{
  virtual void set_monomial_ordering(
      const Monomial_Ordering *, bool = true
  ) override;
  virtual void sort_by_order() override;
  virtual Monomial & leading_monomial() const override;
  virtual Prime_Field_Element leading_coefficient() const override;
  virtual unsigned length() const override;
  virtual bool is_zero() const override;
  ///@}
  /** @name Computation */
  ///@{
  virtual Polynomial_Hashed * zero_polynomial() const override;
  virtual Polynomial_Hashed * monomial_multiple(const Monomial &t) const override;
  virtual Polynomial_Hashed * scalar_multiple(const Prime_Field_Element &) const override;
  ///@}
  /** @name Iteration */
  ///@{
  virtual Hashed_Polynomial_Iterator * new_iterator() const override;
  virtual Hashed_Polynomial_Iterator * begin() const override;
  virtual Hashed_Polynomial_Iterator * end() const override;
  ///@}
  friend class Hashed_Polynomial_Iterator;
protected:
  bool sort_indirect(
      pair< unsigned, Prime_Field_Element > & first,
      pair< unsigned, Prime_Field_Element > & second
  );
  Poly_Strategy_Data * strat = nullptr;
  vector< pair< unsigned, Prime_Field_Element > > terms;
  F4_Hash & hash;
  vector< Monomial * > & M;
  Monomial_Ordering * mord;
};

#endif