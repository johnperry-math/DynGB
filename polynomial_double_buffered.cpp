#ifndef __POLYNOMIAL_DOUBLE_BUFFERED_CPP_
#define __POLYNOMIAL_DOUBLE_BUFFERED_CPP_

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

#include "polynomial_double_buffered.hpp"

DB_Polynomial_Iterator::~DB_Polynomial_Iterator() { /* no action needed */ }

void DB_Polynomial_Iterator::restart_iteration() { current_position = head; }

void DB_Polynomial_Iterator::moveRight() { ++current_position; }

void DB_Polynomial_Iterator::moveLeft() { --current_position; }

bool DB_Polynomial_Iterator::fellOff() const {
  return current_position < head or current_position >= tail;
}

bool DB_Polynomial_Iterator::canMoveLeft() const {
  return current_position - 1 < head;
}

bool DB_Polynomial_Iterator::canMoveRight() const {
   return current_position + 1 < tail;
}

const Monomial & DB_Polynomial_Iterator::currMonomial() const {
  return T[current_position];
}

const Prime_Field_Element & DB_Polynomial_Iterator::currCoeff() const {
  return A[current_position];
}

void DB_Polynomial_Iterator::set_currCoeff(const Prime_Field_Element & a) {
  A[current_position] = a;
}

void DB_Polynomial_Iterator::set_currMonomial(const Monomial & t) {
  T[current_position] = t;
}

DB_Polynomial_Iterator::DB_Polynomial_Iterator(
    const Double_Buffered_Polynomial * f, bool at_end
) {
  p_base = p = f;
  A = f->coeffs[f->active_buffer];
  T = f->mons[f->active_buffer];
  if (at_end) current_position = f->tail;
  else current_position = 0;
  head = f->head;
  tail = f->tail;
}

Double_Buffered_Polynomial::Double_Buffered_Polynomial(Polynomial_Ring & R,
                           Monomial_Ordering * order
) : Mutable_Polynomial(R, order) {
  mons[0] = mons[1] = nullptr;
  coeffs[0] = coeffs[1] = nullptr;
  sizes[0] = sizes[1] = 0;
  active_buffer = 0;
  head = 0;
  tail = 0;
}

Double_Buffered_Polynomial::Double_Buffered_Polynomial(
    Abstract_Polynomial const & p
) : Mutable_Polynomial(p.base_ring(), p.monomial_ordering()) {
  mons[0] = mons[1] = nullptr;
  coeffs[0] = coeffs[1] = nullptr;
  sizes[0] = sizes[1] = 0;
  active_buffer = 0;
  head = 0;
  tail = p.length();
  expand_buffer(0, tail);
  Polynomial_Iterator * pi = p.new_iterator();
  unsigned i = 0;
  while (not pi->fellOff()) {
    mons[0][i] = pi->currMonomial();
    coeffs[0][i] = pi->currCoeff();
    pi->moveRight();
    ++i;
  }
  delete pi;
}

Double_Buffered_Polynomial::~Double_Buffered_Polynomial() {
  if (mons[0] != nullptr) {
    for (unsigned i = 0; i < sizes[0]; ++i)
      mons[0][i].deinitialize();
    free(mons[0]);
    free(coeffs[0]);
  }
  if (mons[1] != nullptr) {
    for (unsigned i = 0; i < sizes[1]; ++i)
      mons[1][i].deinitialize();
    free(mons[1]);
    free(coeffs[1]);
  }
}

Monomial & Double_Buffered_Polynomial::leading_monomial() const {
  return mons[active_buffer][head];
}

Prime_Field_Element Double_Buffered_Polynomial::leading_coefficient()
const {
  return coeffs[active_buffer][head];
}

unsigned Double_Buffered_Polynomial::length() const {
  return tail - head;
}

bool Double_Buffered_Polynomial::is_zero() const { return head == tail; }

bool Double_Buffered_Polynomial::can_reduce(Abstract_Polynomial &other)
const {
  return mons[active_buffer][head] | other.leading_monomial();
}

Double_Buffered_Polynomial * Double_Buffered_Polynomial::zero_polynomial(
) const {
  Double_Buffered_Polynomial * result
    = new Double_Buffered_Polynomial(base_ring());
  return result;
}

void Double_Buffered_Polynomial::set_monomial_ordering(
    const Monomial_Ordering * order, bool sort_anew
) {
  for (unsigned i = head; i < tail; ++i)
    mons[active_buffer][i].set_monomial_ordering(order);
  if (sort_anew)
    sort_by_order();
}

Double_Buffered_Polynomial * Double_Buffered_Polynomial::monomial_multiple(
    const Monomial & t
) const {
  Double_Buffered_Polynomial * p = new Double_Buffered_Polynomial(*this);
  unsigned i = p->active_buffer;
  Monomial * M = p->mons[i];
  for (unsigned j = head; j < tail; ++j)
    M[j] *= t;
  return p;
}

Double_Buffered_Polynomial * Double_Buffered_Polynomial::scalar_multiple(
    const Prime_Field_Element & c)
const {
  Double_Buffered_Polynomial * p = new Double_Buffered_Polynomial(*this);
  unsigned i = p->active_buffer;
  Prime_Field_Element * C = p->coeffs[i];
  for (unsigned j = head; j < tail; ++j)
    C[j] *= c;
  return p;
}

Mutable_Polynomial & Double_Buffered_Polynomial::operator += (
  const Abstract_Polynomial & p)
{
  Polynomial_Iterator * pi = p.new_iterator();
  unsigned i = active_buffer;
  unsigned j = 1 - active_buffer;
  // make sure buffer can contain sum -- this is an overestimate
  unsigned m = length() + p.length();
  if (not test_buffer(j, m))
    expand_buffer(j, m);
  // insertion will start from point 0 in new buffer
  unsigned k = 0;
  unsigned l = head;
  while (not pi->fellOff() and l < tail) {
    if (mons[i][l].is_like(pi->currMonomial())) {
      mons[j][k] = mons[i][l];
      coeffs[j][k] = coeffs[i][l] + pi->currCoeff();
      ++l;
      pi->moveRight();
    } else if (pi->currMonomial() > mons[i][l]) {
      mons[j][k] = pi->currMonomial();
      coeffs[j][k] = pi->currCoeff();
      pi->moveRight();
    } else {
      mons[j][k] = mons[i][l];
      coeffs[j][k] = coeffs[i][l];
      ++l;
    }
    if (not coeffs[j][k].is_zero())
      ++k;
  }
  while (not pi->fellOff()) {
    mons[j][k] = pi->currMonomial();
    coeffs[j][k] = pi->currCoeff();
    ++k;
    pi->moveRight();
  }
  while (l < tail) {
    mons[j][k] = mons[i][l];
    coeffs[j][k] = coeffs[i][l];
    ++k;
    ++l;
  }
  delete pi;
  // adjust to new buffer
  active_buffer = j;
  head = 0;
  tail = k;
  return *this;
}

Mutable_Polynomial & Double_Buffered_Polynomial::operator -= (
    const Abstract_Polynomial & p
) {
  Polynomial_Iterator * pi = p.new_iterator();
  unsigned i = active_buffer;
  unsigned j = 1 - active_buffer;
  // make sure buffer can contain sum -- this is an overestimate
  unsigned m = length() + p.length();
  if (not test_buffer(j, m))
    expand_buffer(j, m);
  // insertion will start from point 0 in new buffer
  unsigned k = 0;
  unsigned l = head;
  while (not pi->fellOff() and l < tail) {
    if (mons[i][l].is_like(pi->currMonomial())) {
      mons[j][k] = mons[i][l];
      coeffs[j][k] = coeffs[i][l] - pi->currCoeff();
      ++l;
      pi->moveRight();
    } else if (pi->currMonomial() > mons[i][l]) {
      mons[j][k] = pi->currMonomial();
      coeffs[j][k] = -(pi->currCoeff());
      pi->moveRight();
    } else {
      mons[j][k] = mons[i][l];
      coeffs[j][k] = coeffs[i][l];
      ++l;
    }
    if (not coeffs[j][k].is_zero())
      ++k;
  }
  while (not pi->fellOff()) {
    mons[j][k] = pi->currMonomial();
    coeffs[j][k] = -(pi->currCoeff());
    ++k;
    pi->moveRight();
  }
  while (l < tail) {
    mons[j][k] = mons[i][l];
    coeffs[j][k] = coeffs[i][l];
    ++k;
    ++l;
  }
  delete pi;
  // adjust to new buffer
  active_buffer = j;
  head = 0;
  tail = k;
  return *this;
}

void Double_Buffered_Polynomial::add_polynomial_multiple(
    const Prime_Field_Element & a,
    const Monomial & u,
    const Abstract_Polynomial & p,
    bool subtract
  )
{
  Polynomial_Iterator * pi = p.new_iterator();
  bool true_multiple = not u.is_one();
  unsigned i = active_buffer;
  unsigned j = 1 - active_buffer;
  // make sure buffer can contain sum -- this is an overestimate
  unsigned m = length() + p.length();
  if (not test_buffer(j, m))
    expand_buffer(j, m);
  // insertion will start from point 0 in new buffer
  unsigned k = 0;
  unsigned l = head;
  while (not pi->fellOff() and l < tail) {
    if (mons[i][l].like_multiple(pi->currMonomial(), u)) {
      mons[j][k] = mons[i][l];
      coeffs[j][k] = pi->currCoeff();
      coeffs[j][k] *= a;
      if (subtract)
        coeffs[j][k].negate();
      coeffs[j][k] += coeffs[i][l];
      ++l;
      pi->moveRight();
    } else if (mons[i][l].larger_than_multiple(pi->currMonomial(), u)) {
      mons[j][k] = mons[i][l];
      coeffs[j][k] = coeffs[i][l];
      ++l;
    } else {
      mons[j][k] = pi->currMonomial();
      if (true_multiple)
        mons[j][k] *= u;
      coeffs[j][k] = pi->currCoeff();
      coeffs[j][k] *= a;
      if (subtract)
        coeffs[j][k].negate();
      pi->moveRight();
    }
    if (not coeffs[j][k].is_zero())
      ++k;
  }
  while (not pi->fellOff()) {
    mons[j][k] = pi->currMonomial();
    mons[j][k] *= u;
    coeffs[j][k] = pi->currCoeff();
    coeffs[j][k] *= a;
    if (subtract)
      coeffs[j][k].negate();
    ++k;
    pi->moveRight();
  }
  while (l < tail) {
    mons[j][k] = mons[i][l];
    coeffs[j][k] = coeffs[i][l];
    ++k;
    ++l;
  }
  delete pi;
  // adjust to new buffer
  active_buffer = j;
  head = 0;
  tail = k;
}

void Double_Buffered_Polynomial::sort_by_order()
{
  Monomial * M = mons[active_buffer];
  Prime_Field_Element * A = coeffs[active_buffer];
  // if (weights == nullptr)
  if (true) {
    // insertion sort, as we don't expect worst case
    for (int i = 0; i < tail; ++i) {
      Monomial t(M[i]);
      Prime_Field_Element c(A[i]);
      int j = i;
      for (/* */; j > 0 and M[i] > M[j-1]; --j) {
        A[j] = A[j-1];
        M[j] = M[j-1];
      }
      if (i != j) {
        A[j] = c;
        M[j] = t;
      }
    }
  }
}

DB_Polynomial_Iterator * Double_Buffered_Polynomial::new_iterator() const {
  DB_Polynomial_Iterator * result = new DB_Polynomial_Iterator(this);
  return result;
}

Polynomial_Iterator * Double_Buffered_Polynomial::begin() const {
  return new DB_Polynomial_Iterator(this);
}

Polynomial_Iterator * Double_Buffered_Polynomial::end() const {
  return new DB_Polynomial_Iterator(this, true);
}

DB_Polynomial_Iterator *
    Double_Buffered_Polynomial::new_mutable_iterator()
{
  DB_Polynomial_Iterator * result = new DB_Polynomial_Iterator(this);
  return result;
}

void Double_Buffered_Polynomial::add_last(
    const Prime_Field_Element & a, const Monomial & t
) {
  unsigned i = active_buffer;
  // can we store in current buffer?
  if (tail + 1 >= sizes[i]) {
    // no: store in other buffer
    unsigned j = 1 - i;
    if (not test_buffer(j, length() + 1))
      expand_buffer(j, length() + 1);
    // copy data from current buffer to other buffer
    Monomial * T = mons[j];
    Monomial * U = mons[i];
    Prime_Field_Element * A = coeffs[j];
    Prime_Field_Element * B = coeffs[i];
    unsigned k = head = 0;
    for (unsigned l = head; l < tail; ++l, ++k) {
      T[k] = U[l];
      A[k] = B[l];
    }
    tail = k;
    active_buffer = i = j;
  }
  // add new term
  coeffs[i][tail] = a;
  mons[i][tail] = t;
  ++tail;
}

Polynomial_Linked_List * Double_Buffered_Polynomial::detach_head() {
  Polynomial_Linked_List * result = new Polynomial_Linked_List(
      base_ring(), coeffs[active_buffer][head], mons[active_buffer][head]);
  ++head;
  return result;
}

#endif