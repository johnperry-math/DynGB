#ifndef __HILBERT_FUNCTIONS_CPP_
#define __HILBERT_FUNCTIONS_CPP_

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

#include <set>
#include <list>
#include <vector>
#include <cstring>
#include <iostream>
#include <algorithm>

using std::set; using std::list; using std::vector;

#include "system_constants.hpp"

#include "fields.hpp"
#include "monomial.hpp"
#include "polynomial_linked_list.hpp"

#include "hilbert_functions.hpp"

extern list<Monomial> colon_ideal_without_ideals(
    const list<Monomial> &,
    const Monomial &
);

bool is_zero_base_case(const list<Monomial> & T) {
  bool result = true; // innocent until proven guilty
  const Monomial & t = T.front();
  NVAR_TYPE n = t.num_vars();
  bool * d = new bool [n] {false};
  for (const Monomial & t : T) {
    bool found = false;
    for (NVAR_TYPE k = 0; result and k < n; ++k) {
      if (t.degree(k) != 0) {
        if (found) { result = false; }
        else {
          found = true;
          if (d[k]) { result = false; }
          else { d[k] = true; }
        }
      }
    }
  }
  delete [] d;
  return result;
}

Dense_Univariate_Integer_Polynomial * solve_one_monomial_case(
    const list<Monomial> & T,
    const WT_TYPE * grading
) {
  DEG_TYPE n = T.front().weighted_degree(grading) + 1;
  DEG_TYPE d = n - 1;
  Dense_Univariate_Integer_Polynomial * result
      = new Dense_Univariate_Integer_Polynomial(n);
  result->set_coefficient(0, 1);
  result->set_coefficient(d, -1);
  return result;
}

Dense_Univariate_Integer_Polynomial * solve_zero_base_case(
    const list<Monomial> & T,
    const WT_TYPE * grading
) {
  DEG_TYPE n = T.front().weighted_degree(grading) + 1;
  DEG_TYPE d = n;
  list<Monomial>::const_iterator ti = T.begin();
  for (++ti; ti != T.end(); ++ti)
  {
    n += ti->weighted_degree(grading);
    if (ti->weighted_degree(grading) + 1 > d)
      d = ti->weighted_degree(grading) + 1;
  }
  Dense_Univariate_Integer_Polynomial * result
      = new Dense_Univariate_Integer_Polynomial(n);
  result->set_coefficient(0, 1);
  Dense_Univariate_Integer_Polynomial * intermediate
      = new Dense_Univariate_Integer_Polynomial(d);
  intermediate->set_coefficient(0, 1);
  for (ti = T.begin(); ti != T.end(); ++ti) {
    intermediate->set_coefficient(ti->weighted_degree(grading), -1);
    result->multiply_by(*intermediate);
    intermediate->set_coefficient(ti->weighted_degree(grading), 0);
  }
  delete intermediate;
  return result;
}

list<Monomial>::const_iterator is_one_base_case(const list<Monomial> & T) {
  list<Monomial>::const_iterator result = T.end();
  for (
    list<Monomial>::const_iterator ti = T.begin();
    (result != T.end()) and ti != T.end();
    ++ti
   ) {
    list<Monomial> U;
    for (const Monomial & u : T)
      if (not u.is_like(*ti))
        U.push_back(u);
    if (is_zero_base_case(U))
      result = ti;
  }
  return result;
}

Dense_Univariate_Integer_Polynomial * solve_one_base_case(
    const list<Monomial> & T, list<Monomial>::const_iterator ui,
    const WT_TYPE * grading
) {
  // copy T and remove the monomial that isn't a simple power
  list<Monomial> U(T);
  Monomial p = *ui;
  U.erase(ui);
  Dense_Univariate_Integer_Polynomial * result = solve_zero_base_case(U, grading);
  DEG_TYPE d = 0;
  DEG_TYPE max_d = 0;
  for (const Monomial & u : U)
  {
    bool found_index = false;
    unsigned j = 0;
    for (/* already initialized */ ; not found_index and j < u.num_vars(); ++j)
      if (u.degree(j) != 0)
        found_index = true;
    --j; // to compensate for a forced increment at the end of the loop
    DEG_TYPE e = (u.degree(j) > p.degree(j)) ?
        u.degree(j) - p.degree(j) : 0;
    d += (grading == nullptr) ? e : grading[j] * e;
    if (max_d < e)
      max_d = e;
  }
  Dense_Univariate_Integer_Polynomial * prod
      = new Dense_Univariate_Integer_Polynomial(d + p.total_degree() + 1);
  Dense_Univariate_Integer_Polynomial * fac
      = new Dense_Univariate_Integer_Polynomial(max_d + 1);
  fac->set_coefficient(0, 1);
  prod->set_coefficient(0, 1);
  for (const Monomial & u : U)
  {
    NVAR_TYPE k = 0;
    for (/* */; k < u.num_vars() and u.degree(k) == 0; ++k) {
      /* do nothing */
    }
    //cout << "one base case: " << *ui << " and " << p << endl;
    if (grading == nullptr) {
      fac->set_coefficient(u.degree(k) - p.degree(k), -1);
      prod->multiply_by(*fac);
      fac->set_coefficient(u.degree(k) - p.degree(k), 0);
    } else {
      WT_TYPE d = grading[k];
      fac->set_coefficient(d*(u.degree(k) - p.degree(k)), -1);
      prod->multiply_by(*fac);
      fac->set_coefficient(d*(u.degree(k) - p.degree(k)), 0);
    }
  }
  delete fac;
  prod->multiply_by_monomial_of_degree(p.weighted_degree(grading));
  prod->negate();
  result->add(*prod);
  delete(prod);
  return result;
}

list<Monomial>::const_iterator is_splitting_case(const list<Monomial> & T) {
  list<Monomial>::const_iterator result = T.end(); // guilty until proven innocent
  for (
    list<Monomial>::const_iterator ti = T.begin();
    result == T.end() and ti != T.end();
    ++ti
  ) {
    // check whether other monomials are relatively prime to t
    bool relprime = true;
    for (
      list<Monomial>::const_iterator ui = T.begin();
      relprime and ui != T.end();
      ++ui
    ) {
      if (ti != ui)
        relprime = ti->is_coprime(*ui);
    }
    if (relprime)
      result = ti;
  }
  return result;
}

Dense_Univariate_Integer_Polynomial * solve_splitting_case(
  const list<Monomial> & T, list<Monomial>::const_iterator ui,
  const WT_TYPE * grading
) {
  list<Monomial> U, V;
  for (list<Monomial>::const_iterator ti = T.begin(); ti != T.end(); ++ti)
    if (ti != ui)
      U.push_back(*ti);
    else
      V.push_back(*ti);
  Dense_Univariate_Integer_Polynomial * result = hilbert_numerator_bigatti(U, grading);
  Dense_Univariate_Integer_Polynomial * other = solve_one_monomial_case(V, grading);
  result->multiply_by(*other);
  delete other;
  return result;
}

Monomial choose_hilbert_pivot(const list<Monomial> & T) {
  NVAR_TYPE n = T.front().num_vars();
  unsigned * xcount = new unsigned [n] {0};
  for (const Monomial & t : T)
    for (NVAR_TYPE k = 0; k < n; ++k)
      if (t.degree(k) != 0)
        ++xcount[k];
  int i = 0;
  for (NVAR_TYPE k = 1; k < n; ++k)
    if (xcount[k] > xcount[i])
      i = k;
  delete [] xcount;
  unsigned * td = new unsigned [T.size()];
  unsigned j = 0;
  for (const Monomial & t : T)
    if (t.degree(i) != 0)
      td[j++] = t.degree(i);
  unsigned m = j;
  std::sort(td, td + m);
  j = (m == 2) ? 0 : m / 2;
  Monomial result(n);
  result.set_exponent(i, td[j]);
  delete [] td;
  return result;
}

Dense_Univariate_Integer_Polynomial * hilbert_numerator_bigatti(
    const list<Monomial> & T, const WT_TYPE * grading
) {
  bool verbose = false;
  if (verbose) {
    cout << "T = { ";
    for (const Monomial & t : T) cout << t << " ,";
    cout << "}\n";
    if (grading != nullptr) {
      cout << "w = { ";
      for (NVAR_TYPE i = 0; i < T.front().num_vars(); ++i)
        cout << grading[i] << ", ";
      cout << "}\n";
    }
  }
  Dense_Univariate_Integer_Polynomial * result;
  int i;
  list<Monomial>::const_iterator ti;
  if (T.size() == 1)
    result = solve_one_monomial_case(T, grading);
  else if (is_zero_base_case(T))
    result = solve_zero_base_case(T, grading);
  else if ((ti = is_one_base_case(T)) != T.end())
    result = solve_one_base_case(T, ti, grading);
  else if ((ti = is_splitting_case(T)) != T.end())
    result = solve_splitting_case(T, ti, grading);
  else {
    Monomial p = choose_hilbert_pivot(T);
    if (verbose) cout << "pivot = " << p << endl;
    // find pivot precisely
    unsigned pi = 0;
    while (p.degree(pi) == 0) { ++pi; }
    list<Monomial> U, V;
    for (const Monomial & t : T) {
      if (not t.divisible_by(p))
        U.push_back(t);
    }
    V = colon_ideal_without_ideals(T, p);
    U.push_back(p);
    if (verbose) {
      cout << "U = { ";
      for (const Monomial & u : U) cout << u << " ,";
      cout << "}\n";
    }
    if (verbose) {
      cout << "V = { ";
      for (const Monomial & v : V) cout << v << " ,";
      cout << "}\n";
    }
    result = hilbert_numerator_bigatti(V, grading);
    if (verbose) cout << "result from V: " << *result << endl;
    result->multiply_by_monomial_of_degree(p.weighted_degree(grading));
    Dense_Univariate_Integer_Polynomial * other
        = hilbert_numerator_bigatti(U, grading);
    if (verbose) cout << "result from U: " << *other << endl;
    *result += *other;
    delete other;
  }
  return result;
}

Dense_Univariate_Integer_Polynomial * hilbert_second_numerator(
  NVAR_TYPE n,
  Dense_Univariate_Integer_Polynomial * first,
  const WT_TYPE * grading
) {
  Dense_Univariate_Integer_Polynomial * hn
    = new Dense_Univariate_Integer_Polynomial(*first);
  unsigned r = 0;
  if (grading == nullptr) {
    NVAR_TYPE i = 0;
    // perform synthetic division on the Hilbert numerator
    for (/* */; r == 0 and i < n; ++i) {
      COEF_TYPE a = hn->coeff(hn->degree());
      Dense_Univariate_Integer_Polynomial * tmp
        = new Dense_Univariate_Integer_Polynomial(*hn);
      for (DEG_TYPE j = hn->degree(); j != 0; --j) {
        COEF_TYPE b = hn->coeff(j - 1);
        a += b;
        hn->set_coefficient(j - 1, a);
      }
      if ((r = hn->coeff(0)) != 0) {
        delete hn;
        hn = tmp;
      } else {
        delete tmp;
        for (unsigned j = 0; j < hn->degree(); ++j)
          hn->set_coefficient(j, hn->coeff(j+1));
        hn->set_coefficient(hn->degree(), 0);
      }
    }
  } else {
    // perform long division (somewhat optimized) with max unused grading
    // first find max grading
    WT_TYPE curr_grad = grading[0];
    NVAR_TYPE curr_grad_index = 0;
    for (unsigned k = 1; k < n; ++k)
      if (grading[k] > curr_grad) {
        curr_grad = grading[k];
        curr_grad_index = k;
      }
    // divide using max grading
    for (NVAR_TYPE i = 1; hn->degree() != 0 and i < n; ++i) {
      Dense_Univariate_Integer_Polynomial * tmp
        = new Dense_Univariate_Integer_Polynomial(*hn);
      Dense_Univariate_Integer_Polynomial * q
        = new Dense_Univariate_Integer_Polynomial(hn->degree());
      for (DEG_TYPE j = hn->degree(); j > curr_grad - 1; --j) {
        COEF_TYPE a = hn->coeff(j);
        q->set_coefficient(j - curr_grad, a);
        hn->set_coefficient(j, 0);
        hn->set_coefficient(j - curr_grad, hn->coeff(j - curr_grad) + a);
      }
      if (not hn->is_zero()) {
        delete q; delete hn;
        hn = tmp;
      } else {
        delete tmp; delete hn;
        hn = q;
      }
      // find next larger gradient (or repeat if same occurs)
      WT_TYPE tmp_grad = grading[0];
      NVAR_TYPE tmp_grad_index = 0;
      bool searching = true;
      for (unsigned k = 0; searching and k < n; ++k) {
        if (grading[k] == curr_grad and k > curr_grad_index) {
          curr_grad_index = k;
          searching = false;
        }
        else if (grading[k] > tmp_grad and grading[k] < curr_grad)
        {
          tmp_grad = grading[k];
          tmp_grad_index = k;
        }
      }
      curr_grad = tmp_grad;
      curr_grad_index = tmp_grad_index;
    }
  }
  return hn;
}

unsigned ideal_dimension(
  NVAR_TYPE n,
  const Dense_Univariate_Integer_Polynomial * h1,
  const Dense_Univariate_Integer_Polynomial * h2
) {
  return n - (h1->degree() - h2->degree());
}

Dense_Univariate_Rational_Polynomial * polynomial_binomial(
  COEF_TYPE a, COEF_TYPE b
) {
  Dense_Univariate_Rational_Polynomial * p
    = new Dense_Univariate_Rational_Polynomial((b > 0) ? b : 1);
  if (b == 0)
    // p = 1
    p->set_coefficient(0, 1, 1);
  else {
    // p = (t + a) / b
    p->set_coefficient(0, a, b);
    p->set_coefficient(1, 1, b);
    // p = p * (t + a - i) / i for i = 1, ..., b
    Dense_Univariate_Rational_Polynomial * q
      = new Dense_Univariate_Rational_Polynomial(2);
    for (unsigned i = 1; i < b; ++i) {
      q->set_coefficient(0, a - i, i);
      q->set_coefficient(1, 1, i);
      //cout << "a - i: " << a - i << " i: " << i << endl;
      //cout << "p: " << *p << endl;
      //cout << "q: " << *q << endl;
      p->multiply_by(*q);
      //cout << "p*q: " << *p << endl;
    }
    delete q;
  }
  return p;
}

Dense_Univariate_Rational_Polynomial * hilbert_polynomial(
    NVAR_TYPE n,
    unsigned int pole_order,
    const list<Monomial> T,
    Dense_Univariate_Integer_Polynomial * hn,
    Dense_Univariate_Integer_Polynomial * hn2
) {
  bool own_hn1 = false;
  bool own_hn2 = false;
  if (hn == nullptr) {
    if (T.size() == 0)
      return nullptr;
    else {
      own_hn1 = true;
      hn = hilbert_numerator_bigatti(T);
    }
  }
  if (hn2 == nullptr) {
    own_hn2 = true;
    hn2 = hilbert_second_numerator(n, hn);
  }
  if (pole_order == 0)
    pole_order = n - (hn->degree() - hn2->degree());
  Dense_Univariate_Rational_Polynomial * hp
    = new Dense_Univariate_Rational_Polynomial(pole_order + 1);
  if (pole_order != 0) {
    DEG_TYPE d = hn2->degree();
    COEF_TYPE d1 = pole_order - 1;
    for (DEG_TYPE i = 0; i <= d; ++i) {
      //q = polynomial_binomial(d1 - d + i, d1);
      Dense_Univariate_Rational_Polynomial * q = polynomial_binomial(d1 - i, d1);
      //q->scale_by(hn2->coeff(d - i));
      q->scale_by(hn2->coeff(i));
      hp->add(*q);
      delete q;
    }
  }
  if (pole_order % 2 != n % 2)
    hp->negate();
  if (own_hn1) delete hn;
  if (own_hn2) delete hn2;
  return hp;
}


#endif