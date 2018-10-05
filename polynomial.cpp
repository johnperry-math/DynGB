#ifndef __POLYNOMIAL_CPP_
#define __POLYNOMIAL_CPP_

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

#include "polynomial.hpp"

Polynomial_Iterator::~Polynomial_Iterator() {}

Polynomial_Ring & Abstract_Polynomial::base_ring() const { return R; }

const Prime_Field & Abstract_Polynomial::ground_field() const {
  return R.ground_field();
}

unsigned Abstract_Polynomial::number_of_variables() const {
  return R.number_of_variables();
}

bool Abstract_Polynomial::can_reduce(Abstract_Polynomial &other) const
{ return leading_monomial() | other.leading_monomial(); }

void Abstract_Polynomial::set_strategy(Poly_Strategy_Data * psd) { strat = psd; }

void Abstract_Polynomial::print(ostream & os) const { os << *this; }

void Abstract_Polynomial::println(ostream & os) const {
  print(); os << endl;
}

ostream & operator << (ostream & os, const Abstract_Polynomial & p) {
  if (p.is_zero())
    os << "0 ";
  else {
    Polynomial_Iterator * pi = p.new_iterator();
    while (!pi->fellOff()) {
      if (pi->currMonomial().is_one() or not pi->currCoeff().is_one())
        os << pi->currCoeff() << " * ";
      //os << pi->currMonomial();
      if (not pi->currMonomial().is_one())
        pi->currMonomial().print(false, os, p.base_ring().name_list());
      pi->moveRight();
      if (!pi->fellOff()) os << " + ";
    }
    delete pi;
  }
  return os;
}

Mutable_Polynomial::~Mutable_Polynomial() { }

void Mutable_Polynomial::multiply_by_scalar(const Prime_Field_Element &a) {
  Mutable_Polynomial_Iterator * pi = new_mutable_iterator();
  while (!pi->fellOff()) {
    Prime_Field_Element b = pi->currCoeff() * a;
    pi->set_currCoeff(b);
    pi->moveRight();
  }
  delete pi;
}

void Mutable_Polynomial::multiply_by_monomial(const Monomial &t) {
  Mutable_Polynomial_Iterator * pi = new_mutable_iterator();
  while (!pi->fellOff()) {
    Monomial b = pi->currMonomial();
    b *= t;
    pi->set_currMonomial(b);
    pi->moveRight();
  }
  delete pi;
}

void Mutable_Polynomial::reduce_by(const Abstract_Polynomial &p) {
  while (p.can_reduce(*this)) {
    Monomial u = leading_monomial();
    u /= p.leading_monomial();
    Prime_Field_Element a = leading_coefficient();
    a *= p.leading_coefficient().inverse();
    add_polynomial_multiple(a, u, p, true);
  }
}

DEG_TYPE Abstract_Polynomial::standard_degree() const {
  DEG_TYPE d = 0;
  Polynomial_Iterator * pi = new_iterator();
  while (not pi->fellOff()) {
    if (pi->currMonomial().total_degree() > d)
      d = pi->currMonomial().total_degree();
    pi->moveRight();
  }
  delete pi;
  return d;
}

DEG_TYPE Abstract_Polynomial::weighted_degree(const WT_TYPE * w) const {
  if (w == nullptr) return standard_degree();
  else {
    DEG_TYPE d = 0;
    Polynomial_Iterator * pi = new_iterator();
    while (not pi->fellOff()) {
      if (pi->currMonomial().weighted_degree(w) > d)
        d = pi->currMonomial().weighted_degree(w);
      pi->moveRight();
    }
    delete pi;
    return d;
  }
}

#endif
