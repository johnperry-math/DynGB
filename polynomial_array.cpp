#ifndef __POLYNOMIAL_ARRAY_CPP_
#define __POLYNOMIAL_ARRAY_CPP_

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

#include "polynomial_array.hpp"

Constant_Polynomial_Iterator::Constant_Polynomial_Iterator(
      const Constant_Polynomial * poly, bool at_end
) {
  p = poly;
  if (at_end) i = p->m - 1;
  else i = p->head;
}

Constant_Polynomial_Iterator::~Constant_Polynomial_Iterator() { }

const Monomial & Constant_Polynomial_Iterator::currMonomial()
    const
{ return p->M[i];}

const Prime_Field_Element &
    Constant_Polynomial_Iterator::currCoeff() const
{ return p->A[i]; }

void Constant_Polynomial_Iterator::restart_iteration() { i = p->head; }

void Constant_Polynomial_Iterator::moveRight() { ++i; }

void Constant_Polynomial_Iterator::moveLeft() { --i; }

bool Constant_Polynomial_Iterator::canMoveRight() const { return i + 1 < p->m; }

bool Constant_Polynomial_Iterator::canMoveLeft() const { return i - 1 < p->head; }

bool Constant_Polynomial_Iterator::fellOff() const
{ return i < p->head or i == p->m; }

Constant_Polynomial::Constant_Polynomial(
    unsigned length,
    Polynomial_Ring & R,
    const Monomial *mons,
    const Prime_Field_Element *coeffs,
    const Monomial_Ordering * order
) : Abstract_Polynomial(R, order)
{
  if (order == nullptr) {
    if (mons[0].monomial_ordering() == nullptr)
      order = generic_grevlex_ptr;
    else
      order = mons[0].monomial_ordering();
  }
  m = length;
  M = (Monomial *)calloc(m, sizeof(Monomial));
  A = (Prime_Field_Element *)calloc(m, sizeof(Prime_Field_Element));
  for (unsigned i = 0; i < m; ++i) {
    M[i].common_initialization();
    M[i].initialize_exponents(mons[i].num_vars());
    M[i] = mons[i];
    if (order != M[i].monomial_ordering())
      M[i].set_monomial_ordering(order);
    A[i] = coeffs[i];
  }
  head = 0;
}

Constant_Polynomial::Constant_Polynomial(
    unsigned length,
    Polynomial_Ring & R,
    const vector<Monomial *> mons,
    const COEF_TYPE *coeffs,
    unsigned start
) : Abstract_Polynomial(R, mons[0]->monomial_ordering())
{
  m = length;
  M = (Monomial *)calloc(m, sizeof(Monomial));
  A = (Prime_Field_Element *)calloc(m, sizeof(Prime_Field_Element));
  auto F = R.ground_field();
  auto n = R.number_of_variables();
  Prime_Field_Element scale(F.inverse(coeffs[start]), &F);
  unsigned j = 0;
  for (unsigned i = start; j < m and i < mons.size(); ++i) {
    if (coeffs[i] != 0) {
      M[j].common_initialization();
      M[j].initialize_exponents(n);
      M[j] = *(mons[i]);
      A[j] = scale * coeffs[i];
      ++j;
    }
  }
  head = 0;
}

Constant_Polynomial::Constant_Polynomial(
    Polynomial_Ring & R,
    const list<Monomial> & mons,
    const list<Prime_Field_Element> & coeffs,
    const Monomial_Ordering * order
) : Abstract_Polynomial(R, order)
{
  if (order == nullptr) {
    if (mons.front().monomial_ordering() == nullptr)
      order = generic_grevlex_ptr;
    else
      order = mons.front().monomial_ordering();
  }
  m = mons.size();
  M = (Monomial *)calloc(m, sizeof(Monomial));
  A = (Prime_Field_Element *)calloc(m, sizeof(Prime_Field_Element));
  auto Mi = mons.begin();
  auto Ai = coeffs.begin();
  unsigned i = 0;
  while (Mi != mons.end()) {
    M[i].common_initialization();
    M[i].initialize_exponents(Mi->num_vars());
    M[i] = *Mi;
    if (order != Mi->monomial_ordering())
      M[i].set_monomial_ordering(order);
    A[i] = *Ai;
    ++Mi; ++Ai; ++i;
  }
  head = 0;
}

Constant_Polynomial::Constant_Polynomial(
    unsigned length,
    Polynomial_Ring & R,
    const Monomial_Ordering * order
) : Abstract_Polynomial(R, order)
{
  m = length;
  M = (Monomial *)malloc(m*sizeof(Monomial));
  A = (Prime_Field_Element *)malloc(m*sizeof(Prime_Field_Element));
  for (unsigned i = 0; i < m; ++i) {
    M[i].common_initialization();
    M[i].initialize_exponents(number_of_variables());
    M[i].set_monomial_ordering(order);
  }
  head = 0;
}

Constant_Polynomial::Constant_Polynomial(
  Polynomial_Ring & R,
  const Monomial_Ordering * order,
  uint64_t size,
  uint64_t * AM
) : Abstract_Polynomial(R, order) {
  m = size;
  const Prime_Field & F = R.ground_field();
  M = (Monomial *)malloc(m*sizeof(Monomial));
  A = (Prime_Field_Element *)malloc(m*sizeof(Prime_Field_Element));
  unsigned j = 0;
  for (unsigned i = 0; i < m; ++i) {
    A[i].assign(AM[j++], &F);
    M[i].common_initialization();
    M[i].initialize_exponents(number_of_variables());
    for (NVAR_TYPE k = 0; k < number_of_variables(); ++k)
      M[i].set_exponent(k, AM[j++]);
    M[i].set_monomial_ordering(order);
  }
  head = 0;
}

uint64_t * Constant_Polynomial::serialized(uint64_t & size) {
  NVAR_TYPE n = number_of_variables();
  /*uint64_t * result = (uint64_t *)malloc(
      sizeof(uint64_t)*m*(n+1)
  );*/
  uint64_t * result = new uint64_t[m*(n+1)];
  unsigned j = 0;
  for (unsigned i = 0; i < m; ++i) {
    result[j++] = A[i].value();
    for (NVAR_TYPE k = 0; k < n; ++k)
      result[j++] = M[i][k];
  }
  size = m*(n+1);
  return result;
}

Constant_Polynomial::~Constant_Polynomial() {
  for (unsigned i = head; i < m; ++i)
    M[i].deinitialize();
  free(M); free(A); M = nullptr; A = nullptr;
}

void Constant_Polynomial::set_monomial_ordering(
    const Monomial_Ordering * order, bool sort_anew
) {
  for (int i = head; i < m; ++i)
    M[i].set_monomial_ordering(order);
  if (sort_anew)
    sort_by_order();
}

void Constant_Polynomial::sort_by_order()
{
  // insertion sort, as we don't expect worst case
  for (int i = head + 1; i < m; ++i) {
    int j = i;
    Monomial t(M[i]);
    Prime_Field_Element c(A[i]);
    for (/* */; j > head and t > M[j-1]; --j) {
      A[j] = A[j-1];
      M[j] = M[j-1];
    }
    if (i != j) {
      A[j] = c;
      M[j] = t;
    }
  }
}

Monomial & Constant_Polynomial::leading_monomial() const { return M[head]; }

Prime_Field_Element Constant_Polynomial::leading_coefficient() const {
  return A[head];
}

unsigned Constant_Polynomial::length() const { return m - head; }

bool Constant_Polynomial::is_zero() const {
  return m == head or A[head].is_zero();
}

Constant_Polynomial * Constant_Polynomial::zero_polynomial() const {
  Monomial empty(M[0].num_vars());
  Prime_Field_Element zero(0, A[0].field());
  Monomial newM [] { empty };
  Prime_Field_Element newA [] { zero };
  return new Constant_Polynomial(1, R, newM, newA, monomial_ordering());
}

Constant_Polynomial * Constant_Polynomial::monomial_multiple(const Monomial &t)
const {
  Constant_Polynomial *result = new Constant_Polynomial(*this);
  for (unsigned i = head; i < m + 1; ++i)
    result->M[i] *= t;
  return result;
}

Constant_Polynomial * Constant_Polynomial::scalar_multiple(
    const Prime_Field_Element &c
) const {
  Constant_Polynomial *result = new Constant_Polynomial(*this);
  for (unsigned i = head; i < m; ++i)
    result->A[i] *= c;
  return result;
}

Constant_Polynomial_Iterator *
Constant_Polynomial::new_iterator() const
{ return new Constant_Polynomial_Iterator(this); }

Polynomial_Iterator * Constant_Polynomial::begin() const
{ return new Constant_Polynomial_Iterator(this); }

Polynomial_Iterator * Constant_Polynomial::end() const
{ return new Constant_Polynomial_Iterator(this, true); }

Constant_Polynomial::Constant_Polynomial(const Abstract_Polynomial & p)
: Abstract_Polynomial(p.base_ring(), p.monomial_ordering())
{
  Polynomial_Iterator * pi = p.new_iterator();
  m = p.length();
  M = (Monomial *)malloc(m*sizeof(Monomial));
  A = (Prime_Field_Element *)malloc(m*sizeof(Prime_Field_Element));
  for (unsigned i = 0; i < m and !pi->fellOff(); pi->moveRight())
  {
    M[i].common_initialization();
    M[i].initialize_exponents(pi->currMonomial().num_vars());
    M[i] = pi->currMonomial();
    A[i] = pi->currCoeff();
    ++i;
  }
  head = 0;
  delete pi;
}

Mutable_Constant_Polynomial_Iterator::Mutable_Constant_Polynomial_Iterator(
    Constant_Polynomial * poly)
{
  p = poly;
  i = p->head;
}

Mutable_Constant_Polynomial_Iterator::~Mutable_Constant_Polynomial_Iterator() { }

void Mutable_Constant_Polynomial_Iterator::restart_iteration() { i = p->head; }

void Mutable_Constant_Polynomial_Iterator::moveRight() { ++i; }

void Mutable_Constant_Polynomial_Iterator::moveLeft() { --i; }

bool Mutable_Constant_Polynomial_Iterator::fellOff() const {
  return i < p->head or i >= p->m;
}

const Monomial & Mutable_Constant_Polynomial_Iterator::currMonomial() const {
  return p->M[i];
}

const Prime_Field_Element & Mutable_Constant_Polynomial_Iterator::currCoeff() const {
  return p->A[i];
}

void Mutable_Constant_Polynomial_Iterator::set_currCoeff(const Prime_Field_Element & a) {
  p->A[i] = a;
}

void Mutable_Constant_Polynomial_Iterator::set_currMonomial(const Monomial & t) {
  p->M[i] = t;
}

#endif