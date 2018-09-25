#ifndef __PARTICULAR_ORDERINGS_CPP_
#define __PARTICULAR_ORDERINGS_CPP_

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

#include "system_constants.hpp"

#include "goda.hpp"
#include "particular_orderings.hpp"

/**
  @ingroup memorygroup
  @brief Memory manager for graded ordering data.
  @details Automatically initialized, but clients need to call the destructor
    when finished.
*/
Grading_Order_Data_Allocator<WT_TYPE> * goda = nullptr;

/**
  @ingroup memorygroup
  @brief Memory manager for graded orderings (not their data; see goda for that).
  @details Automatically initialized, but clients need to call the destructor
    when finished.
*/
Grading_Order_Data_Allocator<WGrevlex_Order_Data> * woda = nullptr;

void * WGrevlex_Order_Data::operator new(size_t size) {
  if (woda == nullptr) woda = new Grading_Order_Data_Allocator<WGrevlex_Order_Data>(size, "woda");
  WGrevlex_Order_Data * result = woda->get_new_block();
  return result;
}

void WGrevlex_Order_Data::operator delete(void *t) {
  woda->return_used_block(static_cast<WGrevlex_Order_Data *>(t));
}

void Generic_Grevlex::set_data(Monomial & t) const {
  t.set_ordering_data(nullptr);
  t.set_ordering_degree(t.total_degree());
}

bool Generic_Grevlex::first_larger(
    const Monomial & t, const Monomial & u
) const {
  DEG_TYPE dtk = t.ordering_degree();
  DEG_TYPE duk = u.ordering_degree();
  bool searching = dtk == duk;
  if (searching) {
    auto a = t.packed_log();
    auto b = u.packed_log();
    auto i = t.packed_size();
    auto j = u.packed_size();
    for (/* */; searching and i > 0 and j > 0; /* */) {
      i -= 2; j -= 2;
      if (a[i] >= b[j])
        dtk -= a[i+1];
      if (b[j] >= a[i])
        duk -= b[j+1];
      searching = dtk == duk;
    }
  }
  return dtk > duk;
}

bool Generic_Grevlex::first_smaller(
    const Monomial & t, const Monomial & u
) const {
  DEG_TYPE dtk = t.ordering_degree();
  DEG_TYPE duk = u.ordering_degree();
  bool searching = dtk == duk;
  auto a = t.packed_log();
  auto b = u.packed_log();
  auto i = t.packed_size();
  auto j = u.packed_size();
  if (searching) {
    for (/* */; searching and i > 0 and j > 0; /* */) {
      i -= 2; j -= 2;
      if (a[i] >= b[j])
        dtk -= a[i+1];
      if (b[j] >= a[i])
        duk -= b[j+1];
      searching = dtk == duk;
    }
  }
  return dtk < duk;
}

bool Generic_Grevlex::first_larger_than_multiple(
    const Monomial & t, const Monomial & u, const Monomial & v
) const {
  DEG_TYPE dtk = t.ordering_degree();
  DEG_TYPE duk = u.ordering_degree();
  DEG_TYPE dvk = v.ordering_degree();
  bool searching = dtk == duk + dvk;
  if (searching) {
    if (t.is_one()) {
      if (u.is_one()) {
        if (v.is_one()) {
          auto a = t.packed_log();
          auto b = u.packed_log();
          auto c = v.packed_log();
          auto i = t.packed_size();
          auto j = u.packed_size();
          auto k = v.packed_size();
          i -= 2; j -= 2; k -= 2;
          for (/* */; searching and i > 0 and j > 0 and k > 0; /* */) {
            auto max_index = (a[i] > b[j]) ? a[i] : b[j];
            max_index = (max_index > c[k]) ? max_index : c[k];
            if (max_index == a[i]) {
              dtk -= a[i+1];
              i -= 2;
            }
            if (max_index == b[j]) {
              duk -= b[j+1];
              j -= 2;
            }
            if (max_index == c[k]) {
              dvk -= c[k+1];
              k -= 2;
            }
            searching = dtk == duk + dvk;
          }
        } else
          return first_larger(t, u);
      } else if (v.is_one())
        return first_larger(t, v);
      else
        return not t.is_one();
    } else
      return not (u.is_one() or v.is_one());
  }
  return dtk > duk + dvk;
}

Generic_Grevlex generic_grevlex;
Monomial_Ordering * generic_grevlex_ptr = &generic_grevlex;

WGrevlex::WGrevlex(NVAR_TYPE num_vars, WT_TYPE * wts, bool thorough) :
    n(num_vars), weights(wts), thorough_weighting(thorough)
{}

WGrevlex::WGrevlex(Ray r, bool thorough) :
    n(r.get_dimension()), thorough_weighting(thorough)
{
  weights = new WT_TYPE[n];
  for (NVAR_TYPE i = 0; i < n; ++i)
    weights[i] = r[i];
}

void WGrevlex::set_data(Monomial & t) const {
  t.set_ordering_data(nullptr);
  t.set_ordering_degree(t.weighted_degree(weights));
}

bool WGrevlex::first_larger(
    const Monomial & t, const Monomial & u
) const {
  DEG_TYPE dtk = t.ordering_degree();
  DEG_TYPE duk = u.ordering_degree();
  bool searching = dtk == duk;
  if (searching) {
    auto a = t.packed_log();
    auto b = u.packed_log();
    auto i = t.packed_size();
    auto j = u.packed_size();
    for (/* */; searching and i > 0 and j > 0; /* */) {
      i -= 2; j -= 2;
      if (a[i] >= b[j])
        dtk -= a[i+1] * weights[a[i]];
      if (b[j] >= a[i])
        duk -= b[j+1] * weights[b[j]];
      searching = dtk == duk;
    }
  }
  return dtk > duk;
}

bool WGrevlex::first_smaller(
    const Monomial & t, const Monomial & u
) const {
  DEG_TYPE dtk = t.ordering_degree();
  DEG_TYPE duk = u.ordering_degree();
  bool searching = dtk == duk;
  if (searching) {
    auto a = t.packed_log();
    auto b = u.packed_log();
    auto i = t.packed_size();
    auto j = u.packed_size();
    for (/* */; searching and i > 0 and j > 0; /* */) {
      i -= 2; j -= 2;
      if (a[i] >= b[j])
        dtk -= a[i+1] * weights[a[i]];
      if (b[j] >= a[i])
        duk -= b[j+1] * weights[b[j]];
      searching = dtk == duk;
    }
  }
  return dtk < duk;
}

bool WGrevlex::first_larger_than_multiple(
    const Monomial & t, const Monomial & u, const Monomial & v
) const {
  DEG_TYPE dtk = t.ordering_degree();
  DEG_TYPE duk = u.ordering_degree();
  DEG_TYPE dvk = v.ordering_degree();
  bool searching = dtk == duk + dvk;
  if (searching) {
    if (t.is_one()) {
      if (u.is_one()) {
        if (v.is_one()) {
          auto a = t.packed_log();
          auto b = u.packed_log();
          auto c = v.packed_log();
          auto i = t.packed_size();
          auto j = u.packed_size();
          auto k = v.packed_size();
          i -= 2; j -= 2; k -= 2;
          for (/* */; searching and i > 0 and j > 0 and k > 0; /* */) {
            auto max_index = (a[i] > b[j]) ? a[i] : b[j];
            max_index = (max_index > c[k]) ? max_index : c[k];
            if (max_index == a[i]) {
              dtk -= a[i+1] * weights[a[i]];
              i -= 2;
            }
            if (max_index == b[j]) {
              duk -= b[j+1] * weights[b[j]];
              j -= 2;
            }
            if (max_index == c[k]) {
              dvk -= c[k+1] * weights[c[k]];
              k -= 2;
            }
            searching = dtk == duk + dvk;
          }
        } else
          return first_larger(t, u);
      } else if (v.is_one())
        return first_larger(t, v);
      else
        return not t.is_one();
    } else
      return not (u.is_one() or v.is_one());
  }
  return dtk > duk + dvk;
}

void Grevlex_Order_Data::assign_gradings(const Monomial & t) {
  DEG_TYPE value = 0;
  for (NVAR_TYPE l = 0; l < number_of_gradings; ++l) {
    value += t.degree(l);
    gradings[number_of_gradings-l-1] = value;
  }
}

Grevlex_Order_Data::Grevlex_Order_Data(const Monomial & t)
  : number_of_gradings(t.num_vars())
{
  const NVAR_TYPE n = number_of_gradings;
  if (goda == nullptr) goda = new Grading_Order_Data_Allocator<WT_TYPE>(n, "goda");
  gradings = goda->get_new_block();
  assign_gradings(t);
}

Grevlex_Order_Data::Grevlex_Order_Data(const Grevlex_Order_Data & other)
  : number_of_gradings(other.number_of_gradings)
{
  gradings = goda->get_new_block();
  memcpy(gradings, other.gradings, number_of_gradings*sizeof(DEG_TYPE));
}

Grevlex_Order_Data * Grevlex_Order_Data::clone() {
  return new Grevlex_Order_Data(*this);
}

Grevlex_Order_Data::~Grevlex_Order_Data() { goda->return_used_block(gradings); }

DEG_TYPE Grevlex_Order_Data::operator [] (NVAR_TYPE i) const {
  return gradings[i];
}

Grevlex_Ordering::Grevlex_Ordering(NVAR_TYPE number_of_variables)
  : n(number_of_variables)
{}

bool Grevlex_Ordering::first_larger(const Monomial & t, const Monomial & u) const
{
  bool still_tied = true;
  bool result = false;
  for (NVAR_TYPE k = 0; still_tied and k < n; ++k) {
    DEG_TYPE dtk = partial_degree(t, k);
    DEG_TYPE duk = partial_degree(u, k);
    if (dtk < duk)
      still_tied = false;
    else if (dtk > duk) {
      result = true;
      still_tied = false;
    }
  }
  return result;
}

bool Grevlex_Ordering::first_smaller(const Monomial & t, const Monomial & u) const
{
  bool still_tied = true;
  bool result = false;
  for (NVAR_TYPE k = 0; still_tied and k < n; ++k) {
    DEG_TYPE dtk = partial_degree(t, k);
    DEG_TYPE duk = partial_degree(u, k);
    if (dtk > duk)
      still_tied = false;
    else if (dtk < duk) {
      result = true;
      still_tied = false;
    }
  }
  return result;
}

bool Grevlex_Ordering::first_larger_than_multiple(
      const Monomial & t, const Monomial & u, const Monomial & v
) const {
  bool still_tied = true;
  bool result = false;
  for (NVAR_TYPE k = 0; still_tied and k < n; ++k) {
    DEG_TYPE dtk = partial_degree(t, k);
    DEG_TYPE duk = partial_degree(u, k);
    DEG_TYPE dvk = partial_degree(v, k);
    if (dtk < duk + dvk)
      still_tied = false;
    else if (dtk > duk + dvk) {
      result = true;
      still_tied = false;
    }
  }
  return result;
}

DEG_TYPE Grevlex_Ordering::partial_degree(
    const Monomial & t, NVAR_TYPE i
) const {
  return (*(static_cast<Grevlex_Order_Data *>(t.monomial_ordering_data())))[i];
}

DEG_TYPE Grevlex_Ordering::compute_ith_weight(
    const Monomial & t, NVAR_TYPE i
) const {
  DEG_TYPE result = 0;
  for (NVAR_TYPE k = 0; k < n - i; ++k)
    result += t.degree(k);
  return result;
}

void Grevlex_Ordering::set_data(Monomial & t) const {
  if (t.monomial_ordering_data() == nullptr)
    t.set_ordering_data(new Grevlex_Order_Data(t));
  else {
    (static_cast<Grevlex_Order_Data *>(t.monomial_ordering_data()))->assign_gradings(t);
  }
}

Lex_Ordering::Lex_Ordering(NVAR_TYPE number_of_variables) : n(number_of_variables)
{}

bool Lex_Ordering::first_larger(const Monomial & t, const Monomial & u) const {
  bool still_tied = true;
  bool result = false;
  for (NVAR_TYPE k = 0; still_tied and k < n; ++k) {
    if (t[k] < u[k])
      still_tied = false;
    else if (t[k] > u[k]) {
      result = true;
      still_tied = false;
    }
  }
  return result;
}

bool Lex_Ordering::first_smaller(const Monomial & t, const Monomial & u) const {
  bool still_tied = true;
  bool result = false;
  for (NVAR_TYPE k = 0; still_tied and k < n; ++k) {
    if (t[k] > u[k])
      still_tied = false;
    else if (t[k] < u[k]) {
      result = true;
      still_tied = false;
    }
  }
  return result;
}

bool Lex_Ordering::first_larger_than_multiple(
    const Monomial & t, const Monomial & u, const Monomial & v
) const {
  bool still_tied = true;
  bool result = false;
  for (NVAR_TYPE k = 0; still_tied and k < n; ++k) {
    if (t[k] < u[k] + v[k])
      still_tied = false;
    else if (t[k] > u[k] + v[k]) {
      result = true;
      still_tied = false;
    }
  }
  return result;
}

void WGrevlex_Order_Data::assign_gradings(Monomial &t) {
  DEG_TYPE value = 0;
  const WT_TYPE * w
    = (
        static_cast<const CachedWGrevlex_Ordering *>(t.monomial_ordering())
      )->order_weights();
  for (NVAR_TYPE l = 0; l < number_of_gradings; ++l) {
    gradings[number_of_gradings-l-1] = t.weighted_degree(w, l + 1);
  }
  t.set_ordering_degree(gradings[0]);
}

WGrevlex_Order_Data::WGrevlex_Order_Data(Monomial & t)
  : number_of_gradings(t.num_vars())
{
  const NVAR_TYPE n = number_of_gradings;
  if (goda == nullptr) goda = new Grading_Order_Data_Allocator<WT_TYPE>(n, "goda");
  gradings = goda->get_new_block();
  assign_gradings(t);
}

WGrevlex_Order_Data::WGrevlex_Order_Data(const WGrevlex_Order_Data & other)
  : number_of_gradings(other.number_of_gradings)
{
  gradings = goda->get_new_block();
  memcpy(gradings, other.gradings, number_of_gradings*sizeof(DEG_TYPE));
}

WGrevlex_Order_Data * WGrevlex_Order_Data::clone() {
  return new WGrevlex_Order_Data(*this);
}

WGrevlex_Order_Data::~WGrevlex_Order_Data() { goda->return_used_block(gradings); }

DEG_TYPE WGrevlex_Order_Data::operator [] (NVAR_TYPE i) const {
  return gradings[i];
}

CachedWGrevlex_Ordering::CachedWGrevlex_Ordering(
    NVAR_TYPE number_of_variables, WT_TYPE * w, bool thorough
) : n(number_of_variables), weights(w), fully_apply(thorough)
{}

int CachedWGrevlex_Ordering::cmp(const Monomial & t, const Monomial & u) const {
  int result;
  long long cmp = t.ordering_degree() - u.ordering_degree();
  if (cmp < 0) result = -1;
  else if (cmp > 0) result = 1;
  for (NVAR_TYPE k = 1; result == 0 and k < n; ++k) {
    DEG_TYPE dtk = partial_degree(t, k);
    DEG_TYPE duk = partial_degree(u, k);
    result = dtk - duk;
  }
  return result;
}

bool CachedWGrevlex_Ordering::first_larger(const Monomial & t, const Monomial & u)
const {
  bool still_tied = true;
  bool result = false;
  if (t.ordering_degree() > u.ordering_degree()) {
    result = true;
    still_tied = false;
  } else if (t.ordering_degree() < u.ordering_degree()) {
    result = false;
    still_tied = false;
  }
  for (NVAR_TYPE k = 1; still_tied and k < n; ++k) {
    DEG_TYPE dtk = partial_degree(t, k);
    DEG_TYPE duk = partial_degree(u, k);
    if (dtk < duk)
      still_tied = false;
    else if (dtk > duk) {
      result = true;
      still_tied = false;
    }
  }
  return result;
}

bool CachedWGrevlex_Ordering::first_smaller(const Monomial & t, const Monomial & u)
const {
  bool still_tied = true;
  bool result = false;
  if (t.ordering_degree() < u.ordering_degree()) {
    result = true;
    still_tied = false;
  } else if (t.ordering_degree() > u.ordering_degree()) {
    result = false;
    still_tied = false;
  }
  for (NVAR_TYPE k = 1; still_tied and k < n; ++k) {
    DEG_TYPE dtk = partial_degree(t, k);
    DEG_TYPE duk = partial_degree(u, k);
    if (dtk > duk)
      still_tied = false;
    else if (dtk < duk) {
      result = true;
      still_tied = false;
    }
  }
  return result;
}

const WT_TYPE * CachedWGrevlex_Ordering::order_weights() const {
  return weights;
}

bool CachedWGrevlex_Ordering::first_larger_than_multiple(
    const Monomial & t, const Monomial & u, const Monomial & v
) const {
  bool still_tied = true;
  bool result = false;
  if (t.ordering_degree() > u.ordering_degree() + v.ordering_degree()) {
    result = true;
    still_tied = false;
  } else if (t.ordering_degree() < u.ordering_degree() + v.ordering_degree()) {
    result = false;
    still_tied = false;
  }
  for (NVAR_TYPE k = 1; still_tied and k < n; ++k) {
    DEG_TYPE dtk = partial_degree(t, k);
    DEG_TYPE duk = partial_degree(u, k);
    DEG_TYPE dvk = partial_degree(v, k);
    if (dtk < duk + dvk)
      still_tied = false;
    else if (dtk > duk + dvk) {
      result = true;
      still_tied = false;
    }
  }
  return result;
}

DEG_TYPE CachedWGrevlex_Ordering::compute_ith_weight(
    const Monomial & t, NVAR_TYPE i
) const {
  DEG_TYPE result = 0;
  for (NVAR_TYPE k = 0; k < n - i; ++k)
    if (fully_apply)
      result += weights[k]*t[k];
    else
      result += t[k];
  return result;
}

void CachedWGrevlex_Ordering::set_data(Monomial & t) const {
  if (t.monomial_ordering_data() == nullptr) {
    WGrevlex_Order_Data * new_data = new WGrevlex_Order_Data(t);
    t.set_ordering_data(new_data);
  }
  else
    (static_cast<WGrevlex_Order_Data *>(t.monomial_ordering_data()))
        ->assign_gradings(t);
}

const char * Nonsingular_Matrix_Ordering_Exception::what() const throw() {
  return "Nonsingular matrix supplied for matrix ordering";
}

/**
  @ingroup orderinggroup
  @author John Perry
  @date 2016
  @brief verifies that a matrix supplied for an ordering is not nonsingular
  @details For the time being, this reduces the matrix to upper triangular form
    (if possible) and then checks that the diagonal is nonzero.
  @param m number of rows in the matrix referenced by @p A
  @param n number of columns in the matrix referenced by @p A
  @param A array of at least \f$ m\times n \f$ elements
  @return @c true if and only if the matrix referenced by @p A is nonsingular
*/
bool nonsingular(NVAR_TYPE m, NVAR_TYPE n, const WT_TYPE **A) {
  bool possibly = true;
  // first copy A
  long long ** M = new long long * [m];
  for (NVAR_TYPE i = 0; i < m; ++i) {
    M[i] = new long long [n];
    for (NVAR_TYPE j = 0; j < n; ++j)
      M[i][j] = A[i][j];
  }
  std::cout << m << ',' << n << std::endl;
  for (NVAR_TYPE i = 0; i < m; ++i) {
    for (NVAR_TYPE j = 0; j < n; ++j)
      std::cout << M[i][j] << ' ';
    std::cout << std::endl;
  }
  std::cout << std::endl;
  for (NVAR_TYPE i = 0; possibly and i < m; ++i) {
    // first find a nonzero element in this column
    bool searching = true;
    NVAR_TYPE j = i;
    for (/* */; searching and j < m; ++j)
      if (M[j][i] != 0) {
        searching = false;
        --j;
      }
    possibly = not searching;
    if (possibly and j < m) {
      // check if we need to swap rows
      if (j != i) {
        for (NVAR_TYPE k = i; k < n; ++k) {
          long long temp = M[j][k];
          M[j][k] = M[i][k];
          M[i][k] = temp;
        }
      }
      // clear out the rest of the column
      for (j = i + 1; j < m; ++j) {
        WT_TYPE header = M[j][i];
        if (header != 0) {
          // delete header & adjust row
          for (NVAR_TYPE k = i; k < n; ++k)
            M[j][k] = M[j][k] * M[i][i] - M[i][k] * header;
        }
      }
    }
  }
  std::cout << m << ',' << n << std::endl;
  for (NVAR_TYPE i = 0; i < m; ++i) {
    for (NVAR_TYPE j = 0; j < n; ++j)
      std::cout << M[i][j] << ' ';
    std::cout << std::endl;
  }
  std::cout << std::endl;
  // free memory
  for (NVAR_TYPE i = 0; i < m; ++i)
    delete [] M[i];
  delete [] M;
  // should be no need to check main diagonal explicitly;
  // obtaining zero in a main diagonal should mean possibly is now false
  return possibly;
}

Matrix_Ordering::Matrix_Ordering(
    NVAR_TYPE rows, NVAR_TYPE cols, const WT_TYPE **data
) : m(rows), n(cols), W(data) {
  if (not nonsingular(m, n, W))
    throw new Nonsingular_Matrix_Ordering_Exception();
}

bool Matrix_Ordering::first_larger(const Monomial & t, const Monomial & u) const
{
  bool result = false;
  bool searching = true;
  for (NVAR_TYPE i = 0; searching and i < m; ++i) {
    DEG_TYPE wt = 0;
    DEG_TYPE wu = 0;
    for (NVAR_TYPE j = 0; j < n; ++j) {
      wt += W[i][j] * t[j];
      wu += W[i][j] * u[j];
    }
    if (wt < wu)
      searching = false;
    else if (wt > wu) {
      searching = false;
      result = true;
    }
  }
  return result;
}

bool Matrix_Ordering::first_smaller(const Monomial & t, const Monomial & u) const
{
  bool result = false;
  bool searching = true;
  for (NVAR_TYPE i = 0; searching and i < m; ++i) {
    DEG_TYPE wt = 0;
    DEG_TYPE wu = 0;
    for (NVAR_TYPE j = 0; j < n; ++j) {
      wt += W[i][j] * t[j];
      wu += W[i][j] * u[j];
    }
    if (wt > wu)
      searching = false;
    else if (wt < wu) {
      searching = false;
      result = true;
    }
  }
  return result;
}

bool Matrix_Ordering::first_larger_than_multiple(
    const Monomial & t, const Monomial & u, const Monomial & v
) const {
  bool result = false;
  bool searching = true;
  for (NVAR_TYPE i = 0; searching and i < m; ++i) {
    DEG_TYPE wt = 0;
    DEG_TYPE wu = 0;
    DEG_TYPE wv = 0;
    for (NVAR_TYPE j = 0; j < n; ++j) {
      wt += W[i][j] * t[j];
      wu += W[i][j] * u[j];
      wv += W[i][j] * v[j];
    }
    if (wt < wu + wv)
      searching = false;
    else if (wt > wu + wv) {
      searching = false;
      result = true;
    }
  }
  return result;
}

#endif