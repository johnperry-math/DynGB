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
#include <bitset>

using std::set; using std::list; using std::vector;
using std::bitset;

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
  for (auto ti = T.cbegin(); result and ti != T.cend(); ++ti) {
    const auto & t = *ti;
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
  if (T.size() == 1 and T.begin()->total_degree() == 0) {
    return new Dense_Univariate_Integer_Polynomial(1);
  }
  DEG_TYPE n = T.front().weighted_degree(grading) + 1;
  DEG_TYPE d = n;
  // determine the size of the result
  list<Monomial>::const_iterator ti = T.begin();
  for (++ti; ti != T.end(); ++ti)
  {
    n += ti->weighted_degree(grading);
    if (ti->weighted_degree(grading) + 1 > d)
      d = ti->weighted_degree(grading) + 1;
  }
  Dense_Univariate_Integer_Polynomial * result =
      new Dense_Univariate_Integer_Polynomial(n);
  result->set_coefficient(0, 1);
  Dense_Univariate_Integer_Polynomial * intermediate =
      new Dense_Univariate_Integer_Polynomial(d);
  intermediate->set_coefficient(0, 1);
  unsigned e = 0;
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
  EXP_TYPE n = T.begin()->num_vars();
  unsigned non_simple_powers = 0;
  for (auto ti = T.begin(); non_simple_powers < 2 and ti != T.end(); ++ti) {
    auto & t = *ti;
    unsigned variables_used = 0;
    for (EXP_TYPE i = 0; variables_used < 2 and i < n; ++i)
      if (t.degree(i) != 0) {
        ++variables_used;
        if (variables_used == 2)
          ++non_simple_powers;
      }
    if (variables_used > 1)
      // this happens only when there are more than two non-simple powers
      result = ti;
  }
  if (non_simple_powers > 1) result = T.end();
  return result;
}

Dense_Univariate_Integer_Polynomial * solve_one_base_case(
    const list<Monomial> & T, list<Monomial>::const_iterator ui,
    const WT_TYPE * grading
) {
  // copy T and remove the monomial that isn't a simple power
  list<Monomial> U(T);
  Monomial p = *ui;
  U.remove(*ui);
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
    if (grading == nullptr) {
      fac->set_coefficient(u.degree(k) - p.degree(k), -1);
      prod->multiply_by(*fac);
      fac->set_coefficient(u.degree(k) - p.degree(k), 0);
    } else {
      WT_TYPE d = grading[k];
      WT_TYPE e = d*(u.degree(k) - p.degree(k));
      fac->expand_poly(e);
      fac->set_coefficient(e, -1);
      prod->multiply_by(*fac);
      fac->set_coefficient(e, 0);
    }
  }
  delete fac;
  prod->multiply_by_monomial_of_degree(p.weighted_degree(grading));
  prod->negate();
  result->add(*prod);
  delete(prod);
  return result;
}

bool is_splitting_case(
    const list<Monomial> & T,
    list<std::list<Monomial>::const_iterator> & U,
    list<std::list<Monomial>::const_iterator> & V
) {
  bitset<MASK_SIZE> U_mask, V_mask;
  auto ti = T.begin();
  U.push_back(ti);
  U_mask |= ti->mask();
  for (++ti; ti != T.end(); ++ti) {
    auto & t = *ti;
    auto t_mask = t.mask();
    if ((t_mask & U_mask) == 0) {
      V.push_back(ti);
      V_mask |= t_mask;
    } else {
      U.push_back(ti);
      U_mask |= t_mask;
      if ((V_mask & t_mask) != 0) {
        bitset<MASK_SIZE> new_V_mask;
        for (auto vii = V.begin(); vii != V.end(); /* */) {
          auto v_mask = (*vii)->mask();
          if ((U_mask & v_mask) == 0) {
            new_V_mask |= v_mask;
            ++vii;
          } else {
            U.push_back(*vii);
            V.erase(vii);
            U_mask |= v_mask;
            vii = V.begin();
          }
        }
        V_mask = new_V_mask;
      }
    }
  }
  return V.size() != 0;
}

Dense_Univariate_Integer_Polynomial * solve_splitting_case(
    const list<Monomial> & T,
    const list< list<Monomial>::const_iterator > & U,
    const list< list<Monomial>::const_iterator > & V,
    const WT_TYPE * grading
) {
  list<Monomial> U_mons, V_mons;
  for (auto & ui : U) U_mons.push_back(*ui);
  for (auto & vi : V) V_mons.push_back(*vi);
  auto first  = hilbert_numerator_bigatti(U_mons, grading);
  auto second = hilbert_numerator_bigatti(V_mons, grading);
  first->multiply_by(*second);
  delete second;
  return first;
}

Monomial choose_hilbert_pivot(const list<Monomial> & T) {
  NVAR_TYPE n = T.front().num_vars();
  unsigned * xcount = new unsigned [n] {0};
  vector< vector<EXP_TYPE> > tds(n, vector<EXP_TYPE> (T.size(), 0));
  vector< unsigned > lasts(n, 0);
  NVAR_TYPE i = T.size();
  EXP_TYPE max = 0;
  for (const Monomial & t : T) {
    for (NVAR_TYPE k = 0; k < n; ++k) {
      auto td = t.degree(k);
      if (td != 0) {
        tds[k][lasts[k]] = td;
        ++xcount[k];
        ++lasts[k];
        if (xcount[k] > max) {
          i = k;
          max = xcount[k];
        }
      }
    }
  }
  delete [] xcount;
  auto & td = tds[i];
  auto m = lasts[i];
  td.resize(m);
  std::sort(td.begin(), td.end());
  auto j = (m == 2) ? 0 : m / 2;
  Monomial result(n, T.begin()->monomial_ordering());
  result.set_exponent(i, td[j]);
  return result;
}

Dense_Univariate_Integer_Polynomial * hilbert_numerator_bigatti(
    const list<Monomial> & T, const WT_TYPE * grading
) {
  static unsigned iterations = 0, level = 0;
  ++iterations;
  ++level;
  auto this_level = level;
  bool verbose = false;
  if (verbose) {
    cout << "level " << this_level << " T = { ";
    for (const Monomial & t : T) cout << '(' << t << " ," << t.monomial_ordering() << ") ";
    cout << "}\n";
    if (grading != nullptr) {
      cout << "w = { ";
      for (NVAR_TYPE i = 0; i < T.front().num_vars(); ++i)
        cout << grading[i] << ", ";
      cout << "}\n";
    }
  }
  Dense_Univariate_Integer_Polynomial * result;
  list<Monomial>::const_iterator ti;
  list<std::list<Monomial>::const_iterator> U_split, V_split;
  if (T.size() == 1) {
    if (verbose) cout << "level " << this_level << " one monomial\n";
    result = solve_one_monomial_case(T, grading);
  } else if (is_zero_base_case(T)) {
    if (verbose) cout << "level " << this_level << " 0-base\n";
    result = solve_zero_base_case(T, grading);
  } else if ((ti = is_one_base_case(T)) != T.end()) {
    if (verbose) cout << "level " << this_level << " 1-base " << *ti << endl;
    result = solve_one_base_case(T, ti, grading);
  } else if ((is_splitting_case(T, U_split, V_split))) {
    if (verbose) {
      cout << "level " << this_level << " splitting\n";
      for (auto t : T) { cout << t << ", "; } cout << endl << "into:\n";
      for (auto ui : U_split) { cout << *ui << ", "; } cout << endl;
      for (auto vi : V_split) { cout << *vi << ", "; } cout << endl;
    }
    result = solve_splitting_case(T, U_split, V_split, grading);
  } else {
    Monomial p = choose_hilbert_pivot(T);
    if (verbose) cout << "level " << this_level << " pivot = " << p << endl;
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
      cout << "U(" << this_level << ") = { ";
      for (const Monomial & u : U) cout << u << " ,";
      cout << "}\n";
    }
    if (verbose) {
      cout << "V(" << this_level << ") = { ";
      for (const Monomial & v : V) cout << v << " ,";
      cout << "}\n";
    }
    result = hilbert_numerator_bigatti(V, grading);
    if (verbose) cout << "result from V(" << this_level << "): " << *result << endl;
    result->multiply_by_monomial_of_degree(p.weighted_degree(grading));
    Dense_Univariate_Integer_Polynomial * other
        = hilbert_numerator_bigatti(U, grading);
    if (verbose) cout << "result from U(" << this_level << "): " << *other << endl;
    *result += *other;
    delete other;
  }
  if (verbose) cout << "result from level " << this_level << *result << endl;
  verbose = false;
  return result;
}

Dense_Univariate_Integer_Polynomial * hilbert_second_numerator(
  NVAR_TYPE n,
  Dense_Univariate_Integer_Polynomial * first,
  const WT_TYPE * grading
) {
  Dense_Univariate_Integer_Polynomial * hn
    = new Dense_Univariate_Integer_Polynomial(*first);
  if (grading == nullptr) {
    MPZCOEF_TYPE r = 0;
    NVAR_TYPE i = 0;
    // perform synthetic division on the Hilbert numerator
    for (/* */; r == 0 and i < n; ++i) {
      MPZCOEF_TYPE a = hn->coeff(hn->degree());
      Dense_Univariate_Integer_Polynomial * tmp
        = new Dense_Univariate_Integer_Polynomial(*hn);
      for (DEG_TYPE j = hn->degree(); j != 0; --j) {
        MPZCOEF_TYPE b = hn->coeff(j - 1);
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
        MPZCOEF_TYPE a = hn->coeff(j);
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