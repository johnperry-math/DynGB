#ifndef __POLYNOMIAL_LINKED_LIST_
#define __POLYNOMIAL_LINKED_LIST_

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

#include "polynomial_linked_list.hpp"

/**
  @brief memory manager for Monomial_Node
  @ingroup memorygroup
  @details Automatically initialized, but clients need to call the destructor
    when finished.
*/
Grading_Order_Data_Allocator<Monomial_Node> * monododa = nullptr;

void * Monomial_Node::operator new(size_t size) {
  if (monododa == nullptr) monododa = new Grading_Order_Data_Allocator<Monomial_Node>(size);
  Monomial_Node * result = monododa->get_new_block();
  return result;
}

void Monomial_Node::operator delete(void *t) {
  monododa->return_used_block(static_cast<Monomial_Node *>(t));
}

Monomial_Node::Monomial_Node(const Prime_Field_Element & a, const Monomial & u)
    : t(u), c(a)
{ }

Monomial_Node::Monomial_Node(Prime_Field & F, const Monomial & u)
: t(u), c(F.unity())
{}

Monomial & Monomial_Node::monomial() { return t; }

Prime_Field_Element & Monomial_Node::coefficient() { return c; }

bool LLPolynomial_Iterator::canMoveRight() const {
  return iter_curr != nullptr and iter_curr->next != nullptr;
}

bool LLPolynomial_Iterator::canMoveLeft() const {
  return iter_curr != nullptr and iter_curr->prev != nullptr;
}

bool LLPolynomial_Iterator::fellOff() const { return iter_curr == nullptr; }

const Monomial & LLPolynomial_Iterator::currMonomial() const {
  return iter_curr->t;
}

const Prime_Field_Element & LLPolynomial_Iterator::currCoeff() const {
  return iter_curr->c;
}

void LLPolynomial_Iterator::set_currCoeff(const Prime_Field_Element & a) {
  if (!fellOff()) iter_curr->c = a;
}

void LLPolynomial_Iterator::set_currMonomial(const Monomial & t) {
  if (!fellOff()) iter_curr->t = t;
}


LLPolynomial_Iterator::LLPolynomial_Iterator(
    Polynomial_Linked_List * poly, bool at_end
) {
  p = poly;
  iter_curr = p->head;
  if (at_end and iter_curr != nullptr)
    while (iter_curr->next != nullptr)
      iter_curr = iter_curr->next;
}

LLPolynomial_Iterator::LLPolynomial_Iterator(
    const Polynomial_Linked_List * poly, bool at_end)
{
  p = const_cast<Polynomial_Linked_List *>(poly);
  iter_curr = p->head;
  if (at_end and iter_curr != nullptr)
    while (iter_curr->next != nullptr)
      iter_curr = iter_curr->next;
}

void LLPolynomial_Iterator::restart_iteration() {
  iter_curr = p->head;
}

Polynomial_Linked_List::Polynomial_Linked_List(
    Polynomial_Ring & R,
    const Monomial_Ordering * order
) : Mutable_Polynomial(R, order)
{
  head = nullptr;
}

Polynomial_Linked_List::Polynomial_Linked_List(
    Polynomial_Ring & R,
    const Monomial & t,
    const Monomial_Ordering * order
) : Mutable_Polynomial(R, order)
{
  if (order == nullptr) {
    if (t.monomial_ordering() == nullptr)
      order = generic_grevlex_ptr;
    else
      order = t.monomial_ordering();
  }
  head = new Monomial_Node(R.ground_field(), t);
  if (order != t.monomial_ordering())
    head->t.set_monomial_ordering(order);
  head->next = nullptr;
}

Polynomial_Linked_List::Polynomial_Linked_List(
    Polynomial_Ring & R,
    const Prime_Field_Element & c, const Monomial & t,
    const Monomial_Ordering * order
) : Mutable_Polynomial(R, order), head(new Monomial_Node(c, t))
{
  if (order == nullptr) {
    if (t.monomial_ordering() == nullptr)
      order = generic_grevlex_ptr;
    else
      order = t.monomial_ordering();
  }
  if (order != t.monomial_ordering())
    head->t.set_monomial_ordering(order);
  head->next = nullptr;
}

Polynomial_Linked_List::Polynomial_Linked_List(
    Polynomial_Ring & R,
    Monomial_Node * node,
    const Monomial_Ordering * order
) : Mutable_Polynomial(R, order)
{
  if (order == nullptr) {
    if (node == nullptr or node->t.monomial_ordering() == nullptr)
      order = generic_grevlex_ptr;
    else
      order = node->t.monomial_ordering();
  }
  head = node;
  if (node != nullptr) {
    node->next = node->prev = nullptr; // to be safe
    head->t.set_monomial_ordering(order);
  }
}

Polynomial_Linked_List::Polynomial_Linked_List(
    const Polynomial_Linked_List & other
) : Mutable_Polynomial(other.base_ring(), other.monomial_ordering())
{
  head = nullptr;
  Monomial_Node *curr = nullptr;
  Monomial_Node *other_curr = other.head;
  while (other_curr != nullptr) {
    Monomial_Node * next = new Monomial_Node(*other_curr);
    if (head == nullptr) head = curr = next;
    else
    {
      curr->next = next;
      curr->next->prev = curr;
      curr = next;
    }
    other_curr = other_curr->next;
  }
  curr->next = nullptr;
}

Polynomial_Linked_List::Polynomial_Linked_List(const Abstract_Polynomial & p)
: Mutable_Polynomial(p.base_ring(), p.monomial_ordering())
{
  Polynomial_Iterator *pi = p.new_iterator();
  Monomial_Node *curr = nullptr;
  while (!pi->fellOff()) {
    if (curr == nullptr) {
      head = curr = new Monomial_Node(pi->currCoeff(), pi->currMonomial());
      curr->prev = nullptr;
    }
    else {
      curr->next = new Monomial_Node(pi->currCoeff(), pi->currMonomial());
      curr->next->prev = curr;
      curr = curr->next;
    }
    pi->moveRight();
  }
  delete pi;
  curr->next = nullptr;
}

Polynomial_Linked_List::~Polynomial_Linked_List() {
  while (head != nullptr) {
    Monomial_Node * tmp = head->next;
    delete head;
    head = tmp;
  }
}

bool Polynomial_Linked_List::is_zero() const {
  return head == nullptr or head->c.is_zero();
}

Monomial & Polynomial_Linked_List::leading_monomial() const { return head->t; }

Prime_Field_Element Polynomial_Linked_List::leading_coefficient() const {
  return head->c;
}

unsigned Polynomial_Linked_List::length() const {
  unsigned i = 0;
  for (Monomial_Node *curr = head; curr != nullptr; ++i) { curr = curr->next; }
  return i;
}

void Polynomial_Linked_List::set_monomial_ordering(
    const Monomial_Ordering * ord, bool sort_anew
) {
  for (Monomial_Node *curr = head; curr != nullptr; curr = curr->next)
    curr->t.set_monomial_ordering(ord);
  if (sort_anew)
    sort_by_order();
}

Polynomial_Linked_List * Polynomial_Linked_List::zero_polynomial() const {
  return new Polynomial_Linked_List(R);
}

Polynomial_Linked_List * Polynomial_Linked_List::monomial_multiple(
    const Monomial &u
) const {
  Polynomial_Linked_List * p = new Polynomial_Linked_List(*this);
  for (Monomial_Node *curr = p->head; curr != nullptr and !(curr->c.is_zero());
         curr = curr->next)
    curr->t *= u;
  return p;
}

Polynomial_Linked_List * Polynomial_Linked_List::scalar_multiple(
    const Prime_Field_Element &c
) const {
  Polynomial_Linked_List * p = new Polynomial_Linked_List(*this);
  for (Monomial_Node *curr = p->head; curr != nullptr and !(curr->c.is_zero());
         curr = curr->next)
    curr->c *= c;
  return p;
}

Polynomial_Linked_List & Polynomial_Linked_List::operator +=(
      const Abstract_Polynomial &other
) {
  // need an iterator
  Polynomial_Iterator * oi = other.new_iterator();
  // add until one of the polynomials is exhausted
  for (Monomial_Node *p = head;
       (!oi->fellOff() and !(oi->currCoeff().is_zero()))
          and (p != nullptr and !(p->c.is_zero()));
      )
  {
    DEG_TYPE a = p->t.ordering_degree();
    while (not oi->fellOff() and a < oi->currMonomial().ordering_degree()){
      Monomial_Node *r = new Monomial_Node(oi->currCoeff(), oi->currMonomial());
      r->prev = p->prev;
      if (p->prev != nullptr) p->prev->next = r;
      if (head == p) head = r;
      p->prev = r;
      r->next = p;
      oi->moveRight();
    }
    if (not oi->fellOff()) {
      if (p->t.is_like(oi->currMonomial())) {
        p->c += oi->currCoeff();
        oi->moveRight();
        if (not p->c.is_zero())
          p = p->next;
        else {
          Monomial_Node *r = p->next;
          if (p->prev != nullptr) p->prev->next = r;
          if (r != nullptr) r->prev = p->prev;
          if (head == p) head = r;
          delete p;
          p = r;
        }
      } else if (p->t > oi->currMonomial())
        p = p->next;
      else {
        Monomial_Node *r = new Monomial_Node(oi->currCoeff(), oi->currMonomial());
        r->prev = p->prev;
        if (p->prev != nullptr) p->prev->next = r;
        if (head == p) head = r;
        p->prev = r;
        r->next = p;
        oi->moveRight();
      }
    }
  }
  // if this polynomial is exhausted, we need to add what remains of other
  if (!oi->fellOff() and !(oi->currCoeff().is_zero())) {
    if (head == nullptr) {
      head = new Monomial_Node(oi->currCoeff(), oi->currMonomial());
      head->next = head->prev = nullptr;
      oi->moveRight();
    }
    Monomial_Node *p = head;
    while (p->next != nullptr) p = p->next;
    while (!oi->fellOff() and !(oi->currCoeff().is_zero())) {
      p->next = new Monomial_Node(oi->currCoeff(), oi->currMonomial());
      p->next->prev = p;
      p = p->next;
      oi->moveRight();
    }
    p->next = nullptr;
  }
  delete oi;
  return *this;
}

Polynomial_Linked_List & Polynomial_Linked_List::operator -=(
    const Abstract_Polynomial &other
) {
  // need an iterator
  Polynomial_Iterator * oi = other.new_iterator();
  // add until one of the polynomials is exhausted
  for (Monomial_Node *p = head;
       (!oi->fellOff() and !(oi->currCoeff().is_zero()))
          and (p != nullptr and !(p->c.is_zero()));
      )
  {
    if (p->t.is_like(oi->currMonomial())) {
      p->c += oi->currCoeff();
      if (p->c.is_zero()) {
        Monomial_Node *r = p->next;
        if (p->prev != nullptr) p->prev->next = r;
        if (r != nullptr) r->prev = p->prev;
        if (head == p) head = r;
        delete p;
        p = r;
        oi->moveRight();
      }
    } else if (p->t > oi->currMonomial())
      p = p->next;
    else {
      Monomial_Node *r = new Monomial_Node(oi->currCoeff(), oi->currMonomial());
      r->prev = p->prev;
      if (p->prev != nullptr) p->prev->next = r;
      if (head == p) head = r;
      p->prev = r;
      r->next = p;
      oi->moveRight();
    }
  }
  // if this polynomial is exhausted, we need to add what remains of other
  if (!oi->fellOff() and !(oi->currCoeff().is_zero())) {
    if (head == nullptr) {
      head = new Monomial_Node(-(oi->currCoeff()), oi->currMonomial());
      head->next = head->prev = nullptr;
      oi->moveRight();
    }
    Monomial_Node *p = head;
    while (p->next != nullptr) p = p->next;
    while (!oi->fellOff() and !(oi->currCoeff().is_zero())) {
      p->next = new Monomial_Node(oi->currCoeff(), oi->currMonomial());
      p->next->prev = p;
      p = p->next;
      oi->moveRight();
    }
    p->next = nullptr;
  }
  delete oi;
  return *this;
}

void Polynomial_Linked_List::add_polynomial_multiple(
    const Prime_Field_Element &a,
    const Monomial &t,
    const Abstract_Polynomial &q,
    bool subtract
) {
  // cout << "\tadding " << a << " * " << t << " * (" << q << ')' << endl;
  // need an iterator
  Polynomial_Iterator * oi = q.new_iterator();
  // first add monomials until this or q is exhausted
  for (Monomial_Node *curr = head;
       curr != nullptr and !(curr->c.is_zero()) and
          !oi->fellOff() and !(oi->currCoeff().is_zero());
      )
  {
    if (curr->t.like_multiple(oi->currMonomial(), t)) {
      if (subtract) curr->c -= a * oi->currCoeff();
      else curr->c += a * oi->currCoeff();
      if ((curr->c).is_zero()) {
        if (curr->prev != nullptr) curr->prev->next = curr->next;
        if (curr->next != nullptr) curr->next->prev = curr->prev;
        Monomial_Node *tmp = curr->next;
        if (head == curr) head = tmp;
        delete curr;
        curr = tmp;
      }
      else
        curr = curr->next;
      oi->moveRight();
    } else if (curr->t.larger_than_multiple(oi->currMonomial(), t)) {
      curr = curr->next;
    } else { // smaller
      Prime_Field_Element c(oi->currCoeff());
      if (subtract) c = -c;
      Monomial_Node *new_node = new Monomial_Node(c, oi->currMonomial());
      new_node->t *= t;
      new_node->c *= a;
      if (new_node->c.is_zero())
        delete new_node;
      else {
        new_node->prev = curr->prev;
        if (curr->prev != nullptr) curr->prev->next = new_node;
        if (head == curr) head = new_node;
        curr->prev = new_node;
        new_node->next = curr;
        oi->moveRight();
      }
    }
  }
  // if this is exhausted but q is not, add the remaining monomials
  if (!oi->fellOff() and !(oi->currCoeff().is_zero())) {
    if (head == nullptr) {
      head = new Monomial_Node(-(oi->currCoeff()), oi->currMonomial());
      head->t *= t;
      head->c *= a;
      head->next = head->prev = nullptr;
      oi->moveRight();
    }
    Monomial_Node * curr = head;
    while (curr->next != nullptr) curr = curr->next;
    while (!oi->fellOff() and !(oi->currCoeff().is_zero())) {
      if (subtract)
        curr->next = new Monomial_Node(-(oi->currCoeff()), oi->currMonomial());
      else
        curr->next = new Monomial_Node(oi->currCoeff(), oi->currMonomial());
      curr->next->prev = curr;
      curr = curr->next;
      curr->t *= t;
      curr->c *= a;
      oi->moveRight();
    }
    curr->next = nullptr;
  }
  delete oi;
  // cout << "\tresulted in " << *this << endl;
}

void Polynomial_Linked_List::sort_by_order()
{
  // insertion sort, as we don't expect worst case
  for (Monomial_Node *curr = head->next; curr != nullptr; /* */) {
    while (curr->t > curr->prev->t) {
      Monomial_Node *tmp = curr->next;
      curr->prev->next = tmp;
      if (tmp != nullptr) tmp->prev = curr->prev;
      curr->next = curr->prev;
      curr->prev = curr->prev->prev;
      curr->prev->prev = curr;
      curr = tmp;
    }
  }
}

LLPolynomial_Iterator * Polynomial_Linked_List::new_iterator() const {
  return new LLPolynomial_Iterator(this);
}

Polynomial_Iterator * Polynomial_Linked_List::begin() const {
  return new LLPolynomial_Iterator(this);
}

Polynomial_Iterator * Polynomial_Linked_List::end() const {
  return new LLPolynomial_Iterator(this, true);
}

LLPolynomial_Iterator * Polynomial_Linked_List::new_mutable_iterator() {
  return new LLPolynomial_Iterator(this);
}

Polynomial_Linked_List * Polynomial_Linked_List::detach_head() {
  Monomial_Node *tmp = head;
  if (head != nullptr) {
    head = head->next;
    if (head != nullptr) head->prev = nullptr;
  }
  Polynomial_Linked_List * result = new Polynomial_Linked_List(R, tmp);
  return result;
}

void Polynomial_Linked_List::add_last(
    const Prime_Field_Element & c,
    const Monomial & t
) {
  if (head == nullptr) {
    head = new Monomial_Node(c, t);
    head->prev = head->next = nullptr;
  } else {
    Monomial_Node * tail = head;
    while (tail != nullptr and tail->next != nullptr) { tail = tail->next; }
    tail->next = new Monomial_Node(c, t);
    tail->next->prev = tail;
    tail->next->next = nullptr;
  }
}

#endif