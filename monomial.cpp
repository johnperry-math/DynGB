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
* DynGB is distributed in the hope that it will be useful,                    *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with DynGB. If not, see <http://www.gnu.org/licenses/>.               *
\*****************************************************************************/

#include <cstring>
#include <bitset>
#include <mutex>
using std::mutex;

#include "monomial.hpp"

/**
  @brief memory manager for monomial exponents
  @ingroup memorygroup
  @details Automatically initialized, but clients need to call the destructor
    when finished.
*/
Grading_Order_Data_Allocator<EXP_TYPE> * moda = nullptr;
mutex moda_mutex;
/**
  @brief memory manager for monomials (not their exponents; see moda for that).
  @ingroup memorygroup
  @details Automatically initialized, but clients need to call the destructor
    when finished.
*/
Grading_Order_Data_Allocator<Monomial> * monoda = nullptr;
mutex monoda_mutex;

void * Monomial::operator new(size_t size) {
  //monoda_mutex.lock();
  if (monoda == nullptr) monoda = new Grading_Order_Data_Allocator<Monomial>(size, "monoda");
  Monomial * result = monoda->get_new_block();
  //monoda_mutex.unlock();
  return result;
}

void Monomial::operator delete(void *t) {
  //monoda_mutex.lock();
  monoda->return_used_block(static_cast<Monomial *>(t));
  //monoda_mutex.unlock();
}

void Monomial::initialize_exponents(NVAR_TYPE number_of_vars) {
  n = number_of_vars;
  make_last_invalid();
  //moda_mutex.lock();
  if (moda == nullptr) moda = new Grading_Order_Data_Allocator<EXP_TYPE>(2*n, "moda");
  exponents = moda->get_new_block();
  //moda_mutex.unlock();
}

void Monomial::deinitialize() {
  //moda_mutex.lock();
  moda->return_used_block(exponents);
  //moda_mutex.unlock();
  if (ordering_data != nullptr) delete ordering_data;
}

Exponent_Location Monomial::find_exponent(NVAR_TYPE i) const {
  Exponent_Location result;
  if (valid_exponents()) {
    result.loc = last;
    result.already_set = false;
    NVAR_TYPE k = 0;
    while (k < last) {
      if (exponents[k] < i) k += 2;
      else {
        result.loc = k;
        result.already_set = (exponents[k] == i);
        break;
      }
    }
  }
  return result;
}

void Monomial::shift_exponents_right(NVAR_TYPE i) {
  for (NVAR_TYPE k = last; k > i; k -= 2) {
    exponents[k] = exponents[k - 2];
    exponents[k + 1] = exponents[k - 1];
  }
  last += 2;
}

void Monomial::compress() {
  NVAR_TYPE d = 0; // distance to move each entry
  NVAR_TYPE i = 0;
  for (/* already initialized */; i < last; /* */) {
    if (d != 0) {
      exponents[i] = exponents[i + d];
      exponents[i + 1] = exponents[i + d + 1];
    }
    if (exponents[i + 1] != 0)
      i += 2;
    else {
      d += 2;
      last -= 2;
    }
  }
  if (i == 0) make_last_invalid();
  else last = i;
}

void Monomial::set_exponent(NVAR_TYPE i, DEG_TYPE e) {
  //if (e != 0) exponent_mask |= 1 << i; else exponent_mask &= !(1 << i);
  if (not valid_exponents()) {
    last = 2;
    exponents[0] = i;
    exponents[1] = e;
  } else {
    auto pos = find_exponent(i);
    if (pos.already_set) { exponents[pos.loc] = e; }
    else {
      auto k = pos.loc;
      shift_exponents_right(k);
      exponents[k] = i;
      exponents[k + 1] = e;
    }
  }
}

void Monomial::fix_mask() const {
  exponent_mask = 0;
  for (NVAR_TYPE i = 0; i < last; i += 2)
    if (exponents[i + 1] != 0) exponent_mask |= ( 1 << exponents[i] );
    else exponent_mask &= !( 1 << exponents[i] );
}

bool Monomial::same_mask(const Monomial & t) const {
  return exponent_mask == t.exponent_mask;
}

bool Monomial::same_mask(const Monomial & t, const Monomial & u) const {
  return exponent_mask == ( t.exponent_mask | u.exponent_mask );
}

Monomial::Monomial(
    NVAR_TYPE number_of_vars, const Monomial_Ordering * order
) {
  common_initialization(order);
  initialize_exponents(number_of_vars);
  //fix_mask();
  ordering->set_data(*this);
}

void Monomial::make_product_or_quotient(
    const Monomial & t, const Monomial & u, bool product
) {
  n = t.n;
  last = 0;
  auto te = t.exponents, ue = u.exponents;
  NVAR_TYPE i = 0, j = 0;
  if (product) {
    for (/* already initialized */; i < t.last and j < u.last; last += 2) {
      if (t.exponents[i] == u.exponents[j]) {
        exponents[last] = t.exponents[i];
        exponents[last + 1] = t.exponents[i + 1] + u.exponents[j + 1];
        i += 2;
        j += 2;
      } else if (t.exponents[i] < u.exponents[j]) {
        exponents[last] = t.exponents[i];
        exponents[last + 1] = t.exponents[i + 1];
        i += 2;
      }
      else {
        exponents[last] = u.exponents[j];
        exponents[last + 1] = u.exponents[j + 1];
        j += 2;
      }
    }
    for (/* */; i < t.last; last += 2) {
      exponents[last] = t.exponents[i];
      exponents[last + 1] = t.exponents[i + 1];
      i += 2;
    }
    for (/* */; j < u.last; last += 2) {
      exponents[last] = u.exponents[j];
      exponents[last + 1] = u.exponents[j + 1];
      j += 2;
    }
  } else {
    for (/* already initialized */; i < t.last and j < u.last; last += 2) {
      if (t.exponents[i] == u.exponents[j]) {
        exponents[last] = t.exponents[i];
        exponents[last + 1] = t.exponents[i + 1] - u.exponents[j + 1];
        i += 2;
        j += 2;
      } else if (t.exponents[i] < u.exponents[j]) {
        exponents[last] = t.exponents[i];
        exponents[last + 1] = t.exponents[i + 1];
        i += 2;
      }
      else {
        /* this should never happen! */
        exponents[last] = u.exponents[j];
        exponents[last + 1] = -u.exponents[j + 1];
        j += 2;
      }
    }
    for (/* */; i < t.last; last += 2) {
      exponents[last] = t.exponents[i];
      exponents[last + 1] = t.exponents[i + 1];
      i += 2;
    }
    for (/* */; j < u.last; last += 2) {
      /* this loop should never happen! */
      exponents[last] = u.exponents[j];
      exponents[last + 1] = -u.exponents[j + 1];
      j += 2;
    }
    compress();
  }
  //fix_mask();
  ordering->set_data(*this); 
}

Monomial::Monomial(const Monomial & t, const Monomial & u, bool product) {
  common_initialization(t.ordering);
  exponents = moda->get_new_block();
  make_product_or_quotient(t, u, product);
}

Monomial::Monomial(const Monomial &other) {
  common_initialization(other.ordering);
  n = other.n;
  last = other.last;
  //exponent_mask = other.exponent_mask;
  //moda_mutex.lock();
  exponents = moda->get_new_block();
  //moda_mutex.unlock();
  memcpy(exponents, other.exponents, last*sizeof(EXP_TYPE));
  ordering_data = (other.ordering_data == nullptr) ? nullptr
      : ordering_data = other.monomial_ordering_data()->clone();
  ord_degree = other.ord_degree;
}

Monomial::Monomial(
    initializer_list<EXP_TYPE> powers, const Monomial_Ordering * order
) {
  common_initialization(order);
  initialize_exponents(powers.size());
  NVAR_TYPE i = 0;
  last = 0;
  for (
       auto pi = powers.begin();
       pi != powers.end();
       ++pi
  ) {
    if (*pi != 0) {
      exponents[last] = i;
      exponents[last + 1] = *pi;
      last += 2;
    }
    ++i;
  }
  //fix_mask();
  ordering->set_data(*this);
}

Monomial::Monomial(
    NVAR_TYPE size, const EXP_TYPE *powers, const Monomial_Ordering * order
) {
  common_initialization(order);
  n = size;
  last = 0;
  //moda_mutex.lock();
  if (moda == nullptr) moda = new Grading_Order_Data_Allocator<EXP_TYPE>(2*n, "moda");
  exponents = moda->get_new_block();
  //moda_mutex.unlock();
  for (NVAR_TYPE i = 0; i < n; ++i) {
    if (powers[i] != 0) {
      exponents[last] = i;
      exponents[last + 1] = powers[i];
      last += 2;
    }
  }
  //fix_mask();
  ordering->set_data(*this);
}

Monomial::~Monomial() {
  //moda_mutex.lock();
  moda->return_used_block(exponents);
  //moda_mutex.unlock();
  if (ordering_data != nullptr)
    delete ordering_data; 
}

bool Monomial::is_one() const {
  return not valid_exponents();
}

DEG_TYPE Monomial::degree(NVAR_TYPE i) const {
  auto pos = find_exponent(i);
  return (pos.already_set) ? exponents[pos.loc + 1] : 0;
}

DEG_TYPE Monomial::operator[](NVAR_TYPE i) const { return degree(i); }

DEG_TYPE Monomial::total_degree(NVAR_TYPE m) const {
  DEG_TYPE result = 0;
  if (m == 0) m = n;
  for (NVAR_TYPE i = 0; i < last and exponents[i] < m; i += 2)
    result += exponents[i + 1];
  return result;
};

DEG_TYPE Monomial::weighted_degree(const WT_TYPE *weights, NVAR_TYPE m) const {
  DEG_TYPE result = 0;
  if (m == 0) m = n;
  if (weights == nullptr)
    result = total_degree();
  else
    for (NVAR_TYPE i = 0; i < last and exponents[i] < m; i += 2)
      result += weights[exponents[i]]*exponents[i + 1];
  return result;
}

unsigned monomial_cache_hits = 0;
unsigned monomial_cache_misses = 0;

#include <chrono>
using std::chrono::duration; using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;

double caching_time = 0;

DEG_TYPE Monomial::cached_weighted_degree(const WT_TYPE *weights) const {
  //high_resolution_clock::time_point start = high_resolution_clock::now();
  DEG_TYPE result = 0;
  if ( cached_signature ) {
    result = cached_degree;
    ++monomial_cache_hits;
  } else {
    ++monomial_cache_misses;
    if (weights == nullptr)
      result = total_degree();
    else {
      for (NVAR_TYPE i = 0; i < last; i += 2)
        result += weights[exponents[i]]*exponents[i + 1];
      //result += total_degree();
    }
    cached_degree = result;
    cached_signature = true;
  }
  //high_resolution_clock::time_point stop = high_resolution_clock::now();
  //caching_time += duration_cast<duration<double> >(stop - start).count();
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
  //if (
  //    ( cached_signature != 0 ) &&
  //    ( cached_signature == other.cached_signature ) &&
  //    ( cached_degree != other.cached_degree )
  //)
  //  return false;
  //if (exponent_mask != other.exponent_mask) return false;
  bool result = (n == other.n) and (last == other.last);
  for (NVAR_TYPE i = 0; result and i < last; i += 2)
    result = exponents[i] == other.exponents[i]
        and exponents[i+1] == other.exponents[i+1];
  return result;
}

bool Monomial::operator !=(const Monomial &other) const {
  return not (*this == other);
}

bool Monomial::is_coprime(const Monomial &other) const {
  //if ( exponent_mask & other.exponent_mask ) return false;
  bool result = true;
  NVAR_TYPE i = 0, j = 0;
  auto oe = other.exponents;
  for (/* already initialized */; result and i < last and j < other.last; /* */) {
    auto a = exponents[i], b = oe[j];
    result = a != b;
    if (result) {
      if (a < b) i += 2;
      else if (a > b) j += 2;
    }
  }
  return result;
}

bool Monomial::is_like(const EXP_TYPE * e) const {
  bool result = true;
  for (unsigned i = 0; result and i < n; ++i) result = exponents[i] == e[i];
  return result;
}

bool Monomial::is_like(const Monomial &other) const {
  bool result = (n == other.n) and (last == other.last);
  //bool result = exponent_mask == other.exponent_mask;
  auto & oe = other.exponents;
  for (NVAR_TYPE i = 0; result and i < last; i += 2)
    result = exponents[i] == oe[i] and exponents[i+1] == oe[i+1];
  return result;
}

bool Monomial::like_multiple(const EXP_TYPE * e, const Monomial & v) const {
  bool result = (n == v.n) and (last == v.last);
  auto & ve = v.exponents;
  for (NVAR_TYPE i = 0; result and i < last; i += 2)
    result = exponents[i] == ve[i]
        and exponents[i+1] == ve[i+1] + e[exponents[i]];
  return result;
}

bool Monomial::like_multiple(const Monomial &u, const Monomial &v) const {
  //if (
  //    ( cached_signature != 0 ) && ( cached_signature == u.cached_signature ) &&
  //    ( cached_signature == v.cached_signature ) &&
  //    ( cached_degree != u.cached_degree + v.cached_degree )
  //   )
  //  return false;
  //if ( exponent_mask != ( u.exponent_mask | v.exponent_mask ) ) return false;
  bool result = (n == u.n and n == v.n);
  auto ue = u.exponents, ve = v.exponents;
  NVAR_TYPE i = 0, j = 0, k = 0;
  for (
        /* already initialized */;
        result and i < last and j < u.last and k < v.last;
        /* */
  ) {
    if (exponents[i] == ue[j]) {
      if (exponents[i] != ve[k]) {
        result = exponents[i+1] == ue[j+1];
        i += 2; j += 2;
      } else {
        result = exponents[i+1] == ue[j+1] + ve[k+1];
        i += 2; j += 2; k += 2;
      }
    } else if (exponents[i] == ve[k]) {
      result = exponents[i+1] == ve[k+1];
      i += 2; k += 2;
    } else
      result = false;
  }
  for (/* */; result and i < last and j < u.last; /* */) {
    result = (exponents[i] == ue[j]) and (exponents[i+1] == ue[j+1]);
    i += 2; j += 2;
  }
  for (/* */; result and i < last and k < v.last; /* */) {
    result = (exponents[i] == ve[k]) and (exponents[i+1] == ve[k+1]);
    i += 2; k += 2;
  }
  result = result and (i == last) and (j == u.last) and (k == v.last);
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
  //if ( (!exponent_mask) & other.exponent_mask ) return false; 
  bool result = (n == other.n);
  NVAR_TYPE i = 0, j = 0;
  auto oe = other.exponents;
  for (/* already initialized */; result and i < last and j < other.last; /* */)
  {
    if (exponents[i] == oe[j]) {
      result = exponents[i+1] >= oe[j+1];
      i += 2;
      j += 2;
    } else if (exponents[i] < oe[j])
      i += 2;
    else
      result = false;
  }
  result = result and j == other.last;
  return result;
}

bool Monomial::operator |(const Monomial &other) const {
  return other.divisible_by(*this);
}

bool Monomial::divisible_by_power(const Monomial &other, unsigned power) const {
  //if ( (!exponent_mask) & other.exponent_mask ) return false;
  bool result = (n == other.n);
  NVAR_TYPE i = 0, j = 0;
  auto oe = other.exponents;
  for (/* already initialized */; result and i < last and j < other.last; /* */)
  {
    if (exponents[i] == oe[j]) {
      result = exponents[i+1] >= oe[j+1]*power;
      i += 2;
      j += 2;
    } else if (exponents[i] < oe[j])
      i += 2;
    else
      result = false;
  }
  result = result and j == other.last;
  return result;
}

bool Monomial::divides_lcm(const Monomial &t, const Monomial &u) const {
  //if ( ( !exponent_mask ) & ( t.exponent_mask | u.exponent_mask ) ) return false;
  bool result = (n == t.n) and (n == u.n);
  NVAR_TYPE i = 0, j = 0, k = 0;
  auto te = t.exponents, ue = u.exponents;
  for (/* already initialized */; result and i < last and j < t.last and k < u.last; /* */) {
    NVAR_TYPE index;
    EXP_TYPE value;
    int which;
    // first figure out next power
    if (te[j] == ue[k]) {
      index = te[j];
      value = (te[j+1] < ue[k+1]) ? ue[k+1] : te[j+1];
      which = 0;
    } else if (te[j] < ue[k]) {
      index = te[j];
      value = te[j+1];
      which = 1;
    } else {
      index = ue[k];
      value = ue[k+1];
      which = 2;
    }
    if (exponents[i] == index) {
      result = exponents[i+1] <= value;
      i += 2;
      switch (which) {
      case 0: j += 2; k += 2; break;
      case 1: j += 2; break;
      case 2: k += 2; break;
      }
    } else if (exponents[i] < index)
      result = false;
    else {
      switch (which) {
      case 0: j += 2; k += 2; break;
      case 1: j += 2; break;
      case 2: k += 2; break;
      }
    }
  }
  if (i < last and exponents[i] < t.last) {
    for (/* */; result and i < last and j < t.last; /* */) {
      if (exponents[i] < te[j])
        result = false;
      else if (exponents[i] > te[j])
        j += 2;
      else {
        result = exponents[i+1] <= te[j+1];
        i += 2; j += 2;
      }
    }
  }
  if (i < last and exponents[i] < u.last) {
    for (/* */; result and i < last and k < u.last; /* */) {
      if (exponents[i] < ue[k])
        result = false;
      else if (exponents[i] > ue[k])
        k += 2;
      else {
        result = exponents[i+1] <= ue[k+1];
        i += 2; k += 2;
      }
    }
  }
  result = result and i == last;
  return result;
}

bool Monomial::like_lcm(const Monomial &t, const Monomial &u) const {
  //if ( exponent_mask != ( t.exponent_mask | u.exponent_mask ) ) return false;
  bool result = (n == t.n) and (n == u.n);
  NVAR_TYPE i = 0, j = 0, k = 0;
  auto te = t.exponents, ue = u.exponents;
  for (/* already initialized */; result and i < last and j < t.last and k < u.last; /* */) {
    NVAR_TYPE index;
    EXP_TYPE value;
    int which;
    // first figure out next power
    if (te[j] == ue[k]) {
      index = te[j];
      value = (te[j+1] < ue[k+1]) ? ue[k+1] : te[j+1];
      which = 0;
    } else if (te[j] < ue[k]) {
      index = te[j];
      value = te[j+1];
      which = 1;
    } else {
      index = ue[k];
      value = ue[k+1];
      which = 2;
    }
    if (exponents[i] == index) {
      result = exponents[i+1] == value;
      i += 2;
      switch (which) {
      case 0: j += 2; k += 2; break;
      case 1: j += 2; break;
      case 2: k += 2; break;
      }
    } else
      result = false;
  }
  if (i < last and exponents[i] < t.last) {
    for (/* */; result and i < last and j < t.last; /* */) {
      if (exponents[i] != te[j])
        result = false;
      else {
        result = exponents[i+1] == te[j+1];
        i += 2; j += 2;
      }
    }
  }
  if (i < last and exponents[i] < u.last) {
    for (/* */; result and i < last and k < u.last; /* */) {
      if (exponents[i] != ue[k])
        result = false;
      else {
        result = exponents[i+1] == ue[k+1];
        i += 2; k += 2;
      }
    }
  }
  result = result and i == last;
  return result;
}

Monomial & Monomial::operator =(const Monomial &other) {
  if (this != &other)
  {
    //exponent_mask = other.exponent_mask;
    last = other.last;
    if (other.n > n) {
      moda->return_used_block(exponents);
      //exponents = new EXP_TYPE [other.n];
      exponents = moda->get_new_block();
      n = other.n;
    }
    auto oe = other.exponents;
    for (NVAR_TYPE i = 0; i < last; i += 2) {
      exponents[i] = oe[i];
      exponents[i+1] = oe[i+1];
    }
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
  //exponent_mask |= other.exponent_mask;
  static EXP_TYPE * buffer = nullptr;
  if (buffer == nullptr) buffer = moda->get_new_block();
  NVAR_TYPE i = 0, j = 0, k = 0;
  auto oe = other.exponents;
  for (/* already initialized */; i < last and j < other.last; /* */) {
    if (exponents[i] == oe[j]) {
      buffer[k] = exponents[i];
      buffer[k+1] = exponents[i+1] + oe[j+1];
      i += 2;
      j += 2;
      k += 2;
    } else if (exponents[i] < oe[j]) {
      buffer[k] = exponents[i];
      buffer[k+1] = exponents[i+1];
      i += 2;
      k += 2;
    } else {
      buffer[k] = oe[j];
      buffer[k+1] = oe[j+1];
      j += 2;
      k += 2;
    }
  }
  for (/* */; i < last; /* */) {
    buffer[k] = exponents[i];
    buffer[k+1] = exponents[i+1];
    i += 2;
    k += 2;
  }
  for (/* */; j < other.last; /* */) {
    buffer[k] = oe[j];
    buffer[k+1] = oe[j+1];
    j += 2;
    k += 2;
  }
  auto tmp = exponents;
  exponents = buffer;
  buffer = tmp;
  last = k;
  ordering->set_data(*this);
  return *this;
}

Monomial Monomial::operator *(const Monomial & other) const {
  Monomial u(n);
  //u.exponent_mask = exponent_mask | other.exponent_mask;
  NVAR_TYPE i = 0, j = 0, k = 0;
  auto oe = other.exponents, ue = u.exponents;
  for (/* already initialized */; i < last and j < other.last; /* */) {
    if (exponents[i] == oe[j]) {
      ue[k] = exponents[i];
      ue[k+1] = exponents[i+1] + oe[j+1];
      i += 2;
      j += 2;
      k += 2;
    } else if (exponents[i] < oe[j]) {
      ue[k] = exponents[i];
      ue[k+1] = exponents[i+1];
      i += 2;
      k += 2;
    } else {
      ue[k] = oe[j];
      ue[k+1] = oe[j+1];
      j += 2;
      k += 2;
    }
  }
  for (/* */; i < last; i += 2) {
    ue[k] = exponents[i];
    ue[k+1] = exponents[i+1];
    i += 2;
    k += 2;
  }
  for (/* */; j < other.last; j += 2) {
    ue[k] = oe[j];
    ue[k+1] = oe[j+1];
    j += 2;
    k += 2;
  }
  u.last = k;
  ordering->set_data(u);
  return u;
}

Monomial Monomial::operator *(const Indeterminate & other) const {
  Monomial u(n);
  auto i = other.index_in_ring();
  //u.exponent_mask = exponent_mask | ( 1 << i );
  auto ue = u.exponents;
  NVAR_TYPE k = 0;
  for (/* already initialized */; k < last and exponents[k] < i; k += 2) {
    ue[k] = exponents[k];
    ue[k+1] = exponents[k+1];
  }
  if (k < last and exponents[k] == i)
    ++ue[k+1];
  else {
    ue[k] = i;
    ue[k+1] = 1;
  }
  k += 2;
  for (/* */; k - 2 < last; k += 2) {
    ue[k] = exponents[k - 2];
    ue[k + 1] = exponents[k - 1];
  }
  u.last = k;
  u.ordering->set_data(u);
  return u;
}

bool Monomial::operator /=(const Monomial &other) {
  bool result = (n == other.n);
  if (result) {
    NVAR_TYPE i = 0, j = 0;
    auto oe = other.exponents;
    for (/* already initialized */; j < other.last; /* */) {
      if (exponents[i] == oe[j]) {
        exponents[i+1] -= oe[j+1];
        i += 2;
        j += 2;
      } else // should only have exponents[i] < oe[j] here
        i += 2;
    }
  }
  compress();
  ordering->set_data(*this);
  //fix_mask();
  return result;
}

Monomial Monomial::lcm(const Monomial & u) const {
   Monomial result(n, monomial_ordering());
   //result.exponent_mask = exponent_mask | u.exponent_mask;
   result.set_monomial_ordering(monomial_ordering());
   NVAR_TYPE i = 0, j = 0, k = 0;
   auto ue = u.exponents, re = result.exponents;
   for (/* already initialized */; i < last and j < u.last; /* */)
     if (exponents[i] == ue[j]) {
       re[k] = exponents[i];
       re[k+1] = (exponents[i+1] > ue[j+1]) ? exponents[i+1] : ue[j+1];
       i += 2;
       j += 2;
       k += 2;
     } else if (exponents[i] < ue[j]) {
       re[k] = exponents[i];
       re[k+1] = exponents[i+1];
       i += 2;
       k += 2;
     } else {
       re[k] = ue[j];
       re[k+1] = ue[j+1];
       j += 2;
       k += 2;
     }
   for (/* */; i < last; i += 2) {
     re[k] = exponents[i];
     re[k+1] = exponents[i+1];
     k += 2;
   }
   for (/* */; j < u.last; j += 2) {
     re[k] = ue[j];
     re[k+1] = ue[j+1];
     k += 2;
   }
   result.last = k;
   ordering->set_data(result);
   return result;
}

Monomial Monomial::gcd(const Monomial & u) const {
   Monomial result(n, monomial_ordering());
   //result.exponent_mask = exponent_mask & u.exponent_mask;
   result.set_monomial_ordering(monomial_ordering());
   NVAR_TYPE i = 0, j = 0, k = 0;
   auto ue = u.exponents, re = result.exponents;
   for (/* already initialized */; i < last and j < u.last; /* */)
   if (exponents[i] == ue[j]) {
     re[k] = exponents[i];
     re[k+1] = (exponents[i+1] < ue[j+1]) ? exponents[i+1] : ue[j+1];
     i += 2;
     j += 2;
     k += 2;
   } else if (exponents[i] < ue[j])
     i += 2;
   else
     j += 2;
   result.last = k;
   ordering->set_data(result);
   return result;
}

DEG_TYPE Monomial::gcd_degree(const Monomial & u) const {
   DEG_TYPE result;
   NVAR_TYPE i = 0, j = 0, k = 0;
   auto ue = u.exponents;
   for (/* already initialized */; i < last and j < u.last; /* */)
   if (exponents[i] == ue[j]) {
     result += (exponents[i+1] < ue[j+1]) ? exponents[i+1] : ue[j+1];
     i += 2;
     j += 2;
     k += 2;
   } else if (exponents[i] < ue[j])
     i += 2;
   else
     j += 2;
   return result;
}

Monomial Monomial::colon(const Monomial & u) const {
  Monomial result(n, monomial_ordering());
  NVAR_TYPE i = 0, j = 0, k = 0;
  auto ue = u.exponents, re = result.exponents;
  for (/* already initialized */; i < last and j < u.last; /* */)
    if (exponents[i] == ue[j]) {
      if (exponents[i+1] > ue[j+1]) {
        re[k] = exponents[i];
        re[k+1] = exponents[i+1] - ue[j+1];
        k += 2;
      }
      i += 2;
      j += 2;
    } else if (exponents[i] > ue[j]) {
      j += 2;
    } else {
      re[k] = exponents[i];
      re[k+1] = exponents[i+1];
      i += 2;
      k += 2;
    }
  for (/* */; i < last; i += 2) {
    re[k] = exponents[i];
    re[k+1] = exponents[i+1];
    k += 2;
  }
  result.last = k;
  ordering->set_data(result);
  //result.fix_mask();
  return result;
}

ostream & operator << (ostream & os, const Monomial &t) {
  t.print(false, os);
  return os;
}

bool Monomial::larger_than(const Monomial & u) const {
  return ordering->first_larger(*this, u);
}

void Monomial::print(bool add_newline, ostream & os, const string * names) const
{
  if (is_one()) os << "1 ";
  else {
    bool first_printed = false;
    for (NVAR_TYPE i = 0; i < last; i += 2)
      if (exponents[i+1] != 0)
      {
        if (first_printed)
          os << "* ";
        else
          first_printed = true;
        if (names == nullptr)
          os << 'x' << exponents[i];
        else
          os << names[exponents[i]];
        if (exponents[i+1] != 1)
          os << '^' << exponents[i+1];
        os << ' ';
      }
  }
  //os << " << " << exponent_mask << " >> ";
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
