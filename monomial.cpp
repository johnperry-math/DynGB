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
  if (not valid_exponents()) {
    // no exponents have been initialized, so it is not already set
    result.loc = 0;
    result.already_set = false;
  } else {
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
  /*  NVAR_TYPE start = 0, stop = (valid_exponents()) ? last : start;
    NVAR_TYPE k = (start + stop) / 2;
    if (k % 2 == 1) k -= 1;
    bool searching = start != stop;
    while (searching) {
      if (exponents[k] == i) {
        start = stop = k;
        searching = false;
      } else if (exponents[k] < i) {
        start = k;
        k += (stop - k) / 2;
        if (k % 2 == 1) k -= 1;
        searching = k != start;
      } else {
        stop = k;
        k /= 2;
        if (k % 2 == 1) k -= 1;
        searching = k != start;
      }
    }
    if (k == last or exponents[k] == i) {
      result.loc = k;
      result.already_set = k != last;
    } else {
      result.loc = k + 2;
      result.already_set = false;
    }
  }*/
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

Monomial::Monomial(
    NVAR_TYPE number_of_vars, const Monomial_Ordering * order
) {
  common_initialization(order);
  initialize_exponents(number_of_vars);
  ordering->set_data(*this);
}

void Monomial::make_product_or_quotient(
    const Monomial & t, const Monomial & u, bool product
) {
  n = t.n;
  last = 0;
  exponents = moda->get_new_block();
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
  ordering->set_data(*this); 
}

Monomial::Monomial(const Monomial & t, const Monomial & u, bool product) {
  common_initialization(t.ordering);
  make_product_or_quotient(t, u, product);
}

Monomial::Monomial(const Monomial &other) {
  common_initialization(other.ordering);
  n = other.n;
  last = other.last;
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

DEG_TYPE Monomial::weighted_degree(const WT_TYPE *weights, NVAR_TYPE m)
const {
  DEG_TYPE result = 0;
  if (m == 0) m = n;
  if (weights == nullptr)
    result = total_degree();
  else
    for (NVAR_TYPE i = 0; i < last and exponents[i] < m; i += 2)
      result += weights[exponents[i]]*exponents[i + 1];
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

bool Monomial::is_like(const Monomial &other) const { return *this == other; }

bool Monomial::like_multiple(const EXP_TYPE * e, const Monomial & v) const {
  bool result = (n == v.n) and (last == v.last);
  auto ve = v.exponents;
  for (NVAR_TYPE i = 0; result and i < last; i += 2)
    result = exponents[i] == ve[i]
        and exponents[i+1] == ve[i+1] + e[exponents[i]];
  return result;
}

bool Monomial::like_multiple(const Monomial &u, const Monomial &v) const {
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
        result = exponents[i+1] + ue[j+1] + ve[k+1];
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

Monomial & Monomial::operator =(const Monomial &other) {
  if (this != &other)
  {
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
  return result;
}

Monomial Monomial::lcm(const Monomial & u) const {
   Monomial result(n, monomial_ordering());
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

Monomial Monomial::colon(const Monomial & u) const {
  Monomial result(n, monomial_ordering());
  result.set_monomial_ordering(monomial_ordering());
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
  if (is_one()) cout << "1 ";
  else {
    bool first_printed = false;
    for (NVAR_TYPE i = 0; i < last; i += 2)
      if (exponents[i+1] != 0)
      {
        if (first_printed)
          cout << "* ";
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