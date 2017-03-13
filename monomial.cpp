#ifndef __MONOMIAL_CPP_
#define __MONOMIAL_CPP_

/*****************************************************************************\
* This file is part of DynGB.                                                 *
*                                                                             *
* DynGB is free software: you can redistribute it and/or modify               *
* it under the terms of the GNU General Public License as published by        *
* the Free Software Foundation, either version 2 of the License, or           *
* (at your option) any later version.                                         *
*                                                                             *
* Foobar is distributed in the hope that it will be useful,                   *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with DynGB. If not, see <http://www.gnu.org/licenses/>.               *
\*****************************************************************************/

#include <cstring>
#include <bitset>

#include "monomial.hpp"

#if EPACK
  #define MASKALL 0b1111111111111111111111111111111111111111111111111111111111111111
  const uint64_t MASK [] = {
    0b0000000000000000000000000000000000000000000000000000000011111111,
    0b0000000000000000000000000000000000000000000000001111111100000000,
    0b0000000000000000000000000000000000000000111111110000000000000000,
    0b0000000000000000000000000000000011111111000000000000000000000000,
    0b0000000000000000000000001111111100000000000000000000000000000000,
    0b0000000000000000111111110000000000000000000000000000000000000000,
    0b0000000011111111000000000000000000000000000000000000000000000000,
    0b1111111100000000000000000000000000000000000000000000000000000000
  };
#endif

/**
  @brief memory manager for monomial exponents
  @ingroup memorygroup
  @details Automatically initialized, but clients need to call the destructor
    when finished.
*/
Grading_Order_Data_Allocator<EXP_TYPE> * moda = nullptr;
/**
  @brief memory manager for monomials (not their exponents; see moda for that).
  @ingroup memorygroup
  @details Automatically initialized, but clients need to call the destructor
    when finished.
*/
Grading_Order_Data_Allocator<Monomial> * monoda = nullptr;

void * Monomial::operator new(size_t size) {
  if (monoda == nullptr) monoda = new Grading_Order_Data_Allocator<Monomial>(size);
  Monomial * result = monoda->get_new_block();
  return result;
}

void Monomial::operator delete(void *t) {
  monoda->return_used_block(static_cast<Monomial *>(t));
}

void Monomial::initialize_exponents(NVAR_TYPE number_of_vars) {
  n = number_of_vars;
  if (moda == nullptr) moda = new Grading_Order_Data_Allocator<EXP_TYPE>(n);
  exponents = moda->get_new_block();
  for (NVAR_TYPE i = 0; i < n; ++i) exponents[i] = 0;
  mask = 0;
}

void Monomial::deinitialize() {
  moda->return_used_block(exponents);
  if (ordering_data != nullptr) delete ordering_data;
}

void Monomial::set_exponent(NVAR_TYPE i, DEG_TYPE e) {
  exponents[i] = e;
  if (e == 0) {
    mask &= ~ (1 << i);
    #if EPACK
    if (i < NPACK) emask &= ~MASK[i];
    #endif
  }
  else {
    mask |= (1 << i);
    #if EPACK
    if (i < NPACK) emask |= (MASK[i] & (((uint64_t )((PACKSIZE )e)) << (i*NPACK)));
    #endif
  }
}

Monomial::Monomial(
    NVAR_TYPE number_of_vars, const Monomial_Ordering * order
) {
  common_initialization(order);
  initialize_exponents(number_of_vars);
  ordering->set_data(*this);
  mask = 0;
  #if EPACK
  emask = 0;
  #endif
}

Monomial::Monomial(const Monomial &other) {
  common_initialization(other.ordering);
  n = other.n;
  mask = other.mask;
  #if EPACK
  emask = other.emask;
  #endif
  exponents = moda->get_new_block();
  memcpy(exponents, other.exponents, n*sizeof(EXP_TYPE));
  ordering_data = (other.ordering_data == nullptr) ? nullptr
      : ordering_data = other.monomial_ordering_data()->clone();
  ord_degree = other.ord_degree;
}

Monomial::Monomial(
    initializer_list<EXP_TYPE> powers, const Monomial_Ordering * order
) {
  common_initialization(order);
  mask = 0;
  #if EPACK
  emask = 0;
  #endif
  initialize_exponents(powers.size());
  NVAR_TYPE i = 0;
  for (
       auto pi = powers.begin();
       pi != powers.end();
       ++pi
  ) {
    exponents[i] = *pi;
    if (*pi != 0) {
      mask |= (1 << i);
      #if EPACK
      if (i < NPACK) emask |= (MASK[i] & (((uint64_t )((PACKSIZE )exponents[i])) << (i*NPACK)));
      #endif
    }
    ++i;
  }
  ordering->set_data(*this);
}

Monomial::Monomial(
    NVAR_TYPE size, const EXP_TYPE *powers, const Monomial_Ordering * order
) {
  common_initialization(order);
  n = size;
  mask = 0;
  #if EPACK
  emask = 0;
  #endif
  if (moda == nullptr) moda = new Grading_Order_Data_Allocator<EXP_TYPE>(n);
  exponents = moda->get_new_block();
  for (NVAR_TYPE i = 0; i < n; ++i) {
    exponents[i] = powers[i];
    if (exponents[i] != 0) {
      mask |= (1 << i);
      #if EPACK
      if (i < NPACK) emask |= (MASK[i] & (((uint64_t )((PACKSIZE )exponents[i])) << (i*NPACK)));
      #endif
    }
  }
  ordering->set_data(*this);
}

Monomial::~Monomial() {
  moda->return_used_block(exponents);
  if (ordering_data != nullptr)
    delete ordering_data; 
}

bool Monomial::is_one() const {
  bool result = true;
  for (NVAR_TYPE i = 0; result and i < n; ++i)
    result = result and exponents[i] == 0;
  return result;
}

DEG_TYPE Monomial::degree(NVAR_TYPE i) const {
  //if (i < n) return exponents[i]; else return 0;
  return exponents[i];
}

DEG_TYPE Monomial::operator[](NVAR_TYPE i) const { return degree(i); }

DEG_TYPE Monomial::total_degree(NVAR_TYPE m) const {
  DEG_TYPE result = 0;
  m = ( m == 0 or m > n ) ? n : m;
  if (exponents != nullptr)
    for (NVAR_TYPE i = 0; i < m; ++i)
      result += exponents[i];
  return result;
};

DEG_TYPE Monomial::weighted_degree(const WT_TYPE *weights, NVAR_TYPE m)
const {
  DEG_TYPE result = 0;
  m = ( m == 0 or m > n ) ? n : m;
  if (weights == nullptr)
    result = total_degree();
  else if (exponents != nullptr)
    for (NVAR_TYPE i = 0; i < m; ++i)
      result += weights[i]*exponents[i];
  return result;
}

void Monomial::set_monomial_ordering(const Monomial_Ordering * mord) {
  if (ordering != mord) {
    ordering = mord;
    ordering->set_data(*this);
  }
}

void Monomial::set_ordering_data(Monomial_Order_Data * mordat) {
  if (ordering_data != nullptr)
    delete ordering_data;
  ordering_data = mordat;
}

bool Monomial::operator ==(const Monomial &other) const {
  bool result = (n == other.n) and (mask == other.mask)
      #if EPACK
      and (emask == other.emask)
      #endif
  ;
  #if !EPACK
  for (NVAR_TYPE i = 0; result and i < n; ++i)
  #else
  for (NVAR_TYPE i = NPACK; result and i < n; ++i)
  #endif
    result = result and exponents[i] == other.exponents[i];
  return result;
}

bool Monomial::operator !=(const Monomial &other) const {
  return not (*this == other);
}

bool Monomial::is_coprime(const Monomial &other) const {
  return (mask & other.mask) == 0;
}

bool Monomial::is_like(const Monomial &other) const { return *this == other; }

bool Monomial::like_multiple(EXP_TYPE * e, const Monomial & v) const {
  bool result = (n == v.n);
  for (NVAR_TYPE i = 0; result and i < n; ++i)
    result = result and exponents[i] == e[i] + v.exponents[i];
  return result;
}

bool Monomial::like_multiple(const Monomial &u, const Monomial &v) const {
  bool result = (n == u.n and n == v.n) and (mask == (u.mask | v.mask))
      #if EPACK
      and (emask == (u.emask + v.emask))
      #endif
      ;
  #if !EPACK
  for (NVAR_TYPE i = 0; result and i < n; ++i)
  #else
  for (NVAR_TYPE i = NPACK; result and i < n; ++i)
  #endif
    result = result and exponents[i] == u.exponents[i] + v.exponents[i];
  return result;
}

bool Monomial::operator >(const Monomial & u) const { return larger_than(u); }

bool Monomial::operator >=(const Monomial & u) const
{ return ordering->first_larger_or_equal(*this, u); }

bool Monomial::operator <(const Monomial & u) const { return !(*this >= u); }

bool Monomial::operator <=(const Monomial & u) const { return !(*this > u); }

bool Monomial::larger_than_multiple(
    const Monomial & u, const Monomial &v
) const {
  return ordering->first_larger_than_multiple(*this, u, v);
}

bool Monomial::divisible_by(const Monomial &other) const {
  bool result = (n == other.n) and (((~ mask) & other.mask) == 0);
  #if !EPACK
  for (NVAR_TYPE i = 0; result and i < n; ++i)
    result = result and exponents[i] >= other.exponents[i];
  #else
  for (NVAR_TYPE i = 0; result and i < NPACK; ++i)
    result = result
  and ((emask & MASK[i]) >= (other.emask & MASK[i]));
  for (NVAR_TYPE i = NPACK; result and i < n; ++i)
    result = result and exponents[i] >= other.exponents[i];
  #endif
  return result;
}

bool Monomial::operator |(const Monomial &other) const {
  return other.divisible_by(*this);
}

Monomial & Monomial::operator =(const Monomial &other) {
  if (this != &other)
  {
    mask = other.mask;
    #if EPACK
    emask = other.emask;
    #endif
    if (other.n > n) {
      delete [] exponents;
      exponents = new EXP_TYPE [other.n];
    }
    n = other.n;
    for (NVAR_TYPE i = 0; i < n; ++i)
      exponents[i] = other.exponents[i];
    ordering = other.ordering;
    if (ordering_data != nullptr) {
      delete ordering_data;
      ordering_data = nullptr;
    }
    ordering->set_data(*this);
  }
  return *this;
}

Monomial & Monomial::operator *=(const Monomial &other) {
  mask |= other.mask;
  #if EPACK
  emask += other.emask;
  #endif
  NVAR_TYPE i;
  for (i = 0; i + 1 < n; i += 2) {
    exponents[i] += other.exponents[i];
    exponents[i + 1] += other.exponents[i + 1];
  }
  if (n % 2) exponents[i] += other.exponents[i];
  ordering->set_data(*this);
  return *this;
}

Monomial Monomial::operator *(const Monomial & other) const {
  Monomial u(n);
  for (NVAR_TYPE i = 0; i < n; ++i)
    u.set_exponent(i, exponents[i] + other.exponents[i]);
  return u;
}

bool Monomial::operator /=(const Monomial &other) {
  bool result = (n == other.n);
  #if EPACK
  emask -= other.emask;
  #endif
  if (result)
    for (NVAR_TYPE i = 0; i < n; ++i) {
      exponents[i] -= other.exponents[i];
      if (exponents[i] == 0) mask &= ~(1 << i);
    }
  ordering->set_data(*this);
  return result;
}

Monomial Monomial::lcm(const Monomial & u) const {
   Monomial result(n, monomial_ordering());
   result.mask = mask | u.mask;
   #if EPACK
   result.emask = emask;
   #endif
   result.set_monomial_ordering(monomial_ordering());
   for (NVAR_TYPE i = 0; i < n; ++i)
     if (exponents[i] >= u.exponents[i]) result.exponents[i] = exponents[i];
     else {
       result.exponents[i] = u.exponents[i];
       #if EPACK
       if (i < NPACK) {
         result.emask &= ~MASK[i];
         result.emask |= ((((uint64_t )((PACKSIZE )u.exponents[i])) << (i*NPACK)) & MASK[i]);
       }
       #endif
     }
   ordering->set_data(result);
   return result;
}

Monomial Monomial::gcd(const Monomial & u) const {
  Monomial result(n, monomial_ordering());
  result.mask = mask & u.mask;
  #if EPACK
  result.emask = emask;
  #endif
  for (NVAR_TYPE i = 0; i < n; ++i)
    if (exponents[i] <= u.exponents[i]) { result.exponents[i] = exponents[i]; }
    else {
      result.exponents[i] = u.exponents[i];
      #if EPACK
      if (i < NPACK) {
        result.emask &= ~MASK[i];
        result.emask |= ((((uint64_t )((PACKSIZE )u.exponents[i])) << (i*NPACK)) & MASK[i]);
      }
      #endif
    }
  ordering->set_data(result);
  return result;
}

Monomial Monomial::colon(const Monomial & u) const {
  Monomial result(n, monomial_ordering());
  result.mask = 0;
  #if EPACK
  result.emask = 0;
  #endif
  for (NVAR_TYPE i = 0; i < n; ++i)
    if (exponents[i] <= u.exponents[i])
      result.exponents[i] = 0;
    else {
      result.exponents[i] = exponents[i] - u.exponents[i];
      result.mask += (1 << i);
      #if EPACK
      if (i < NPACK) result.emask |= (((uint64_t )((PACKSIZE )exponents[i])) << (i*NPACK)) & MASK[i];
      #endif
    }
  ordering->set_data(result);
  return result;
}

ostream & operator << (ostream & os, const Monomial &t) {
  t.print(false, os);
  return os;
}

bool Monomial::larger_than(const Monomial & u) const {
  return ordering->first_larger(*this, u);
}

#define SHOW_MASK false
#define SHOW_EPACK false
void Monomial::print(bool add_newline, ostream & os, const string * names) const
{
  if (is_one()) cout << "1 ";
  else
    for (NVAR_TYPE i = 0; i < n; ++i)
      if (exponents[i] != 0)
      {
        if (names == nullptr)
          os << "x" << i;
        else
          os << names[i];
        if (exponents[i] != 1)
          os << "^" << exponents[i];
        os << ' ';
      }
  #if SHOW_MASK
  os << '(' << mask << ')';
  #endif
  #if SHOW_EPACK
  os << '(';
  for (NVAR_TYPE i = 0; i < NPACK; ++i)
    os << std::bitset<NPACK>(((emask & MASK[i]) >> (i*NPACK))) << ' ';
  os << ')';
  #endif
  if (add_newline)
    os << endl;
}

void Monomial_Ordering::set_data(Monomial & t) const {
  t.set_ordering_data(nullptr); 
}

bool Monomial_Ordering::first_larger_or_equal(
    const Monomial & t, const Monomial & u
) const {
  return t.is_like(u) or first_larger(t, u);
}

bool Monomial_Ordering::first_larger_or_equal_than_multiple(
    const Monomial & t, const Monomial & u, const Monomial & v
) const {
  return t.like_multiple(u, v) or first_larger_than_multiple(t, u, v);
}

#endif