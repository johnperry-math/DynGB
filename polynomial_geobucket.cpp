#ifndef __POLYNOMIAL_GEOBUCKET_CPP_
#define __POLYNOMIAL_GEOBUCKET_CPP_

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

#include "polynomial_geobucket.hpp"

bool Geobucket_Iterator::canMoveLeft() const {
  return pi != nullptr and pi->canMoveLeft();
}

bool Geobucket_Iterator::canMoveRight() const {
  return pi != nullptr and pi->canMoveRight();
}

bool Geobucket_Iterator::fellOff() const {
  return pi == nullptr or pi->fellOff();
}

const Monomial & Geobucket_Iterator::currMonomial() const {
  return pi->currMonomial();
}

const Prime_Field_Element & Geobucket_Iterator::currCoeff() const {
  return pi->currCoeff();
}

void Geobucket_Iterator::set_currCoeff(const Prime_Field_Element & a)
{
  pi->set_currCoeff(a);
}

void Geobucket_Iterator::set_currMonomial(const Monomial & t) {
  pi->set_currMonomial(t);
}

Geobucket_Iterator * Polynomial_Geobucket::new_iterator() const {
  return new Geobucket_Iterator(this);
}

Polynomial_Iterator * Polynomial_Geobucket::begin() const {
  return new Geobucket_Iterator(this);
}

Polynomial_Iterator * Polynomial_Geobucket::end() const {
  return new Geobucket_Iterator(this, true);
}

Geobucket_Iterator * Polynomial_Geobucket::new_mutable_iterator() {
  return new Geobucket_Iterator(this);
}

Geobucket_Iterator::Geobucket_Iterator(
  const Polynomial_Geobucket * poly, bool at_end
) {
  p = g = const_cast<Polynomial_Geobucket *>(poly);
  bucket_number = 0;
  if (at_end)
    for (unsigned i = 0; i < NUM_BUCKETS; ++i)
      if (poly->buckets[i] != nullptr)
        bucket_number = i;
  pi = new LLPolynomial_Iterator(poly->buckets[bucket_number]);
}

Geobucket_Iterator::Geobucket_Iterator(Polynomial_Geobucket * poly, bool at_end)
{
  p = g = poly;
  bucket_number = 0;
  if (at_end)
    for (unsigned i = 0; i < NUM_BUCKETS; ++i)
      if (poly->buckets[i] != nullptr)
        bucket_number = i;
  pi = new LLPolynomial_Iterator(poly->buckets[bucket_number]);
}

Geobucket_Iterator::~Geobucket_Iterator() {
  delete pi;
}

void Geobucket_Iterator::restart_iteration() {
  if (bucket_number != 0) {
    delete pi;
    bucket_number = 0;
    pi = new LLPolynomial_Iterator(g->buckets[0]);
  }
}

void Geobucket_Iterator::moveRight() {
  if (pi != nullptr) {
    pi->moveRight();
    if (pi->fellOff()) {
      delete pi;
      pi = nullptr;
      do
        ++bucket_number;
      while (bucket_number < NUM_BUCKETS
             and g->buckets[bucket_number] == nullptr);
      if (bucket_number < NUM_BUCKETS and g->buckets[bucket_number] != nullptr)
        pi = new LLPolynomial_Iterator(g->buckets[bucket_number]);
    }
  }
}

void Geobucket_Iterator::moveLeft() {
  if (pi != nullptr) {
    pi->moveLeft();
    if (pi->fellOff()) {
      delete pi;
      pi = nullptr;
      // unlike moving right, we can always fall back on buckets[0]
      do
        --bucket_number;
      while (g->buckets[bucket_number] == nullptr);
      pi = new LLPolynomial_Iterator(g->buckets[bucket_number]);
    }
  }
}

Polynomial_Geobucket::Polynomial_Geobucket(
    Polynomial_Ring & R, Monomial_Ordering * order
) : Mutable_Polynomial(R, order) {
  // I don't expect to need NUM_BUCKETS, but the mere pointers to potential
  // buckets can't be entirely bad (I hope...)
  buckets = (Polynomial_Linked_List **)malloc(
      NUM_BUCKETS*sizeof(Polynomial_Linked_List *));
  buckets[0] = new Polynomial_Linked_List(R, order);
  for (unsigned i = 1; i < NUM_BUCKETS; ++i)
    buckets[i] = nullptr;
}

Polynomial_Geobucket::Polynomial_Geobucket(Abstract_Polynomial & p)
: Mutable_Polynomial(p.base_ring(), p.monomial_ordering())
{
  unsigned n = p.number_of_variables();
  buckets = (Polynomial_Linked_List **)malloc(
      NUM_BUCKETS*sizeof(Polynomial_Linked_List *));
  for (unsigned i = 1; i < NUM_BUCKETS; ++i)
    buckets[i] = nullptr;
  unsigned i = lglen(p.length());
  buckets[i] = new Polynomial_Linked_List(p);
  if (i > 0) {
    buckets[0] = buckets[i]->detach_head();
  }
}

Polynomial_Geobucket::~Polynomial_Geobucket() {
  for (unsigned i = 0; i < NUM_BUCKETS; ++i)
    if (buckets[i] != nullptr)
      delete buckets[i];
  free(buckets);
}

void Polynomial_Geobucket::sort_by_order()
{
  for (unsigned i = 0; i < NUM_BUCKETS; ++i)
    buckets[i]->sort_by_order();
}

void Polynomial_Geobucket::set_monomial_ordering(
    const Monomial_Ordering * order, bool sort_anew
) {
  for (unsigned i = 0; i < NUM_BUCKETS; ++i)
    buckets[i]->set_monomial_ordering(order, sort_anew);
}

Monomial & Polynomial_Geobucket::leading_monomial() const {
  return buckets[0]->leading_monomial();
}

Prime_Field_Element Polynomial_Geobucket::leading_coefficient() const {
  return buckets[0]->leading_coefficient();
}

unsigned Polynomial_Geobucket::length() const {
  unsigned result = 0;
  for (unsigned i = 0; i < NUM_BUCKETS; ++i)
    if (buckets[i] != 0)
      result += buckets[i]->length();
  return result;
}

Polynomial_Linked_List * Polynomial_Geobucket::zero_polynomial() const {
  return new Polynomial_Linked_List(R);
}

bool Polynomial_Geobucket::is_zero() const {
  unsigned result = true;
  for (unsigned i = 0; result and i < NUM_BUCKETS; ++i)
    result = buckets[i] == nullptr or buckets[i]->is_zero();
  return result;
}

bool Polynomial_Geobucket::can_reduce(Abstract_Polynomial & other) const {
  // I'm not sure why you'd want to do this,
  // but we need it for class consistency.
  return leading_monomial() | other.leading_monomial();
}

void Polynomial_Geobucket::recompute_leading_monomial() {
  unsigned i = 0;
  do {
    i = 0;
    for (unsigned j = 1; j < NUM_BUCKETS; ++j)
      if (i != j and buckets[j] != nullptr and (not buckets[j]->is_zero())) {
        if (buckets[i]->is_zero() or
            buckets[j]->leading_monomial() > buckets[i]->leading_monomial())
          i = j;
        else {
          if (buckets[i]->leading_monomial().is_like(
                buckets[j]->leading_monomial())) {
            Polynomial_Linked_List * q = buckets[j]->detach_head();
            *buckets[i] += *q;
            delete q;
            j = 0; // need to restart from bucket 1
            if (buckets[i]->is_zero()) {
              delete buckets[i];
              buckets[i] = nullptr;
              i = 0;
            }
          }
        }
      }
  } while (i != 0 and buckets[i]->is_zero());
  if (i == 0)
    delete buckets[0]->detach_head();
  else {
    if (buckets[0] != nullptr)
      delete buckets[0];
    buckets[0] = buckets[i]->detach_head();
  }
  // condensation does not seem to help
  /*++aops_performed;
  if (aops_performed > 50) {
    aops_performed = 0;
    condense_buckets();
  }*/
}

Polynomial_Geobucket * Polynomial_Geobucket::monomial_multiple(const Monomial &t) const {
  Polynomial_Geobucket * result = new Polynomial_Geobucket(base_ring());
  for (unsigned i = 0; i < NUM_BUCKETS; ++i)
  if (buckets[i] != nullptr)
    result->buckets[i] = buckets[i]->monomial_multiple(t);
  return result;
}

Polynomial_Geobucket * Polynomial_Geobucket::scalar_multiple(const Prime_Field_Element & a)
const {
  Polynomial_Geobucket * result = new Polynomial_Geobucket(base_ring());
  for (unsigned i = 0; i < NUM_BUCKETS; ++i)
  if (buckets[i] != nullptr)
    result->buckets[i] = buckets[i]->scalar_multiple(a);
  return result;
}

Polynomial_Geobucket & Polynomial_Geobucket::operator +=(
    const Abstract_Polynomial & g
) {
  Polynomial_Linked_List * f = buckets[0];
  buckets[0] = new Polynomial_Linked_List(base_ring());
  (*f) += g;
  unsigned i = lglen(f->length());
  unsigned m = i;
  bool placed = false;
  while (not placed and i < NUM_BUCKETS) {
    if (buckets[i] != nullptr and not buckets[i]->is_zero()) {
      (*f) += *(buckets[i]);
      m = lglen(f->length());
    }
    if (m < (BUCKET_BASE << BUCKET_SHIFT * i)) {
      placed = true;
      if (not f->is_zero()) {
        if (buckets[i] != nullptr)
          delete buckets[i];
        buckets[i] = f;
      }
    }
  }
  recompute_leading_monomial();
  return *this;
}

Polynomial_Geobucket & Polynomial_Geobucket::operator -=(const Abstract_Polynomial & g) {
  Polynomial_Linked_List * f = buckets[0];
  buckets[0] = new Polynomial_Linked_List(base_ring());
  (*f) += g;
  unsigned i = lglen(f->length());
  bool placed = false;
  unsigned m;
  while (not placed and i < NUM_BUCKETS) {
    if (buckets[i] != nullptr and not buckets[i]->is_zero()) {
      (*f) -= *(buckets[i]);
      m = lglen(f->length());
    }
    if (m < (BUCKET_BASE << BUCKET_SHIFT * i)) {
      placed = true;
      if (not f->is_zero()) {
        if (buckets[i] != nullptr)
          delete buckets[i];
        buckets[i] = f;
      }
    }
  }
  recompute_leading_monomial();
  return *this;
}

void Polynomial_Geobucket::add_polynomial_multiple(
     const Prime_Field_Element & b, const Monomial & u,
     const Abstract_Polynomial & g, bool subtract
    )
{
  unsigned i = lglen(g.length() + 1);
  if (i >= NUM_BUCKETS) throw(new runtime_error("too large for buckets"));
  Polynomial_Linked_List * f = buckets[0];
  f->add_polynomial_multiple(b, u, g, subtract);
  buckets[0] = new Polynomial_Linked_List(base_ring());
  if (buckets[i] == nullptr)
    buckets[i] = f;
  else if (buckets[i]->is_zero()) {
    delete buckets[i];
    buckets[i] = f;
  } else {
    (*f) += *(buckets[i]);
    delete buckets[i];
    buckets[i] = nullptr;
    while (i < NUM_BUCKETS and lglen(f->length()) > i)
    {
      ++i;
      if (buckets[i] != nullptr) {
        (*f) += *(buckets[i]);
        delete buckets[i];
        buckets[i] = new Polynomial_Linked_List(base_ring());
      }
    }
    if (i < NUM_BUCKETS) {
      if (buckets[i] != nullptr)
        delete buckets[i];
      buckets[i] = f;
    }
    else
      throw(new runtime_error("too large for buckets"));
  }
  recompute_leading_monomial();
}

Polynomial_Linked_List * Polynomial_Geobucket::detach_head() {
  Polynomial_Linked_List * result = buckets[0];
  buckets[0] = new Polynomial_Linked_List(base_ring());
  recompute_leading_monomial();
  return result;
}

void Polynomial_Geobucket::add_last(
    const Prime_Field_Element & a, const Monomial & t
) {
  unsigned i = 1;
  while (i < NUM_BUCKETS and buckets[i] != nullptr
         and buckets[i]->length() + 1 == BUCKET_BASE << BUCKET_SHIFT * i)
    ++i;
  if (i >= NUM_BUCKETS)
    throw(new runtime_error("Buckets are full"));
  if (buckets[i] == nullptr)
    buckets[i] = new Polynomial_Linked_List(R);
  buckets[i]->add_last(a, t);
  if (buckets[0]->is_zero()) recompute_leading_monomial();
}

Abstract_Polynomial * Polynomial_Geobucket::canonicalize(bool constant_result) {
  Polynomial_Linked_List * f;
  if (is_zero())
    f = new Polynomial_Linked_List(R);
  else
    f = new Polynomial_Linked_List(*(buckets[0]));
  for (unsigned i = 1; i < NUM_BUCKETS; ++i)
    if (buckets[i] != nullptr and (not buckets[i]->is_zero())) {
      *f += *(buckets[i]);
      //delete buckets[i]; // removed to avoid deallocation errors
    }
  Abstract_Polynomial * result;
  if (not constant_result)
    result = f;
  else {
    result = new Constant_Polynomial(*f);
    delete f;
  }
  return result;
}

void Polynomial_Geobucket::condense_buckets() {
  if (not is_zero()) {
    Polynomial_Linked_List * f;
    unsigned i = 1;
    while (i < NUM_BUCKETS and (buckets[i] == nullptr or buckets[i]->is_zero()))
      ++i;
    if (i < NUM_BUCKETS) {
      f = new Polynomial_Linked_List(*(buckets[i]));
      buckets[i] = new Polynomial_Linked_List(base_ring());
      ++i;
      while (i < NUM_BUCKETS) {
        if (buckets[i] != nullptr and (not buckets[i]->is_zero())) {
          *f += *(buckets[i]);
          delete buckets[i];
          buckets[i] = new Polynomial_Linked_List(base_ring());
        }
        ++i;
      }
      unsigned i = lglen(f->length() + 1);
      if (i >= NUM_BUCKETS) throw(new runtime_error("too large for buckets"));
      delete buckets[i];
      buckets[i] = f;
    }
  }
}

void Polynomial_Geobucket::print(unsigned i, ostream & os) const {
  os << *(buckets[i]) << endl;
}

void Polynomial_Geobucket::print(ostream & os) const {
  for (unsigned i = 0; i < NUM_BUCKETS; ++i)
    if (buckets[i] != nullptr and (not buckets[i]->is_zero())) {
      if (i > 0)
        os << " + ";
      os << '(' << *(buckets[i]) << ')';
    }
  os << endl;
}

/**
  \brief prints the geobucket in an explicitly geobucket form
  \see print(ostream &) const
*/
ostream & operator << (ostream & os, const Polynomial_Geobucket & p) {
  for (unsigned i = 0; i < NUM_BUCKETS; ++i)
    if (p.buckets[i] != nullptr and (not p.buckets[i]->is_zero())) {
      if (i > 0)
        os << " + ";
      os << '(' << *(p.buckets[i]) << ')';
    }
  return os;
}

#endif