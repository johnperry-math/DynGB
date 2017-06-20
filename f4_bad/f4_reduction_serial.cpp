#ifndef __F4_REDUCTION_SERIAL_CPP__
#define __F4_REDUCTION_SERIAL_CPP__

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

#include "f4_reduction_serial.hpp"
#include "algorithm_buchberger_basic.hpp"

void find_position(
  list<Monomial *>::iterator & ti,
  list<Monomial *>::iterator stop,
  list<Abstract_Polynomial *>::iterator & pi,
  const Monomial & u
) {
  while (ti != stop and **ti > u) {
    ++ti;
    ++pi;
  }
}

F4_Reduction_Data::F4_Reduction_Data(
    const list<Critical_Pair_Basic *> & P,
    const list<Abstract_Polynomial *> & B
) : G(B), Rx(P.front()->first()->base_ring()) {
  cols = 0;
  A = nullptr;
  R_build.clear();
  M_build.clear();
  strategies.clear();
  auto mi = M_build.begin();
  auto ri = R_build.begin();
  mord = P.front()->first()->monomial_ordering();
  for (auto p : P) {
    strategies.push_back(new Poly_Sugar_Data(p->first()));
    if (strategies.back() != nullptr) {
      strategies.back()->at_generation_tasks(p->first_multiplier());
      if (p->second() != nullptr)
        strategies.back()->pre_reduction_tasks(p->second_multiplier().log(), *(p->second()));
    }
    add_monomials(mi, ri, p->first(), p->first_multiplier(), true);
    if (p->second() != nullptr) {
      *ri = const_cast<Abstract_Polynomial *>(p->second());
      ++ri; ++mi;
      add_monomials(mi, ri, p->second(), p->second_multiplier());
    }
  }
  // for each monomial, find an appropriate reducer
  mi = M_build.begin();
  ri = R_build.begin();
  while (mi != M_build.end()) {
    auto g = G.begin();
    bool found = ((*ri) != nullptr);
    while (not found and g != G.end()) {
      if ((**mi).divisible_by((*g)->leading_monomial())) {
        found = true;
        *ri = *g;
        Monomial u(**mi);
        u /= (*g)->leading_monomial();
        add_monomials(mi, ri, *g, u);
        g = G.end();
      }
      ++g;
    }
    ++mi;
    ++ri;
  }
  initialize_many(P);
}

void F4_Reduction_Data::initialize_many(const list<Critical_Pair_Basic *> & P) {
  cols = M_build.size();
  Abstract_Polynomial * p = const_cast<Abstract_Polynomial *>(P.front()->first());
  Prime_Field & F = const_cast<Prime_Field &>(p->ground_field());
  rows = P.size();
  heads.resize(rows);
  nonzero_entries.resize(rows);
  for (unsigned i = 0; i < rows; ++i)
    heads[i] = 0;
  A = (Prime_Field_Element *)malloc(
      sizeof(Prime_Field_Element)*cols*rows*P.size()
  );
  for (unsigned i = 0; i < cols*rows; ++i)
    A[i].assign(0, &F);
  M.clear();
  for (auto mi = M_build.begin(); mi != M_build.end(); ++mi)
    M.push_back(*mi);
  strategies.resize(P.size());
  for (unsigned i = 0; i < strategies.size(); ++i)
    strategies[i] = nullptr;
  unsigned row = 0;
  for (auto cp : P) {
    auto p = cp->first();
    auto t = cp->first_multiplier();
    auto pi = p->new_iterator();
    nonzero_entries[row] = 0;
    auto mi = M_build.begin();
    bool head_found = false;
    for (unsigned i = 0; i < cols; ++i) {
      if ((not pi->fellOff()) and M[i]->like_multiple(pi->currMonomial(), t)) {
        *(A + row*cols + i) = pi->currCoeff();
        pi->moveRight();
        ++nonzero_entries[row];
        if (not head_found) {
          heads[row] = i;
          head_found = true;
        }
      }
      ++mi;
    }
    delete pi;
    strategies[row] = new Poly_Sugar_Data(cp->first());
    delete cp;
    ++row;
  }
  R.resize(R_build.size());
  unsigned i = 0;
  for (Abstract_Polynomial * r : R_build)
    R[i++] = r;
}

void F4_Reduction_Data::add_monomials(
    list<Monomial *>::iterator & t1,
    list<Abstract_Polynomial *>::iterator & r1,
    const Abstract_Polynomial *g,
    const Monomial & u,
    bool new_row
) {
  NVAR_TYPE n = g->leading_monomial().num_vars();
  if (new_row) {
    t1 = M_build.begin();
    r1 = R_build.begin();
    Monomial * t = new Monomial(g->leading_monomial());
    (*t) *= u;
    find_position(t1, M_build.end(), r1, *t);
    if (t1 == M_build.end()) {
      M_build.push_back(t);
      R_build.push_back(nullptr);
      t1 = M_build.begin();
      r1 = R_build.begin();
      auto t2(t1);
      ++t2;
      while (t2 != M_build.end()) {
        ++t1; ++r1;
        ++t2;
      }
    } else if (**t1 != *t) {
      M_build.insert(t1, t);
      R_build.insert(r1, nullptr);
      --t1; --r1;
    }
  }
  auto t2(t1); auto r2(r1);
  Polynomial_Iterator * pi = g->new_iterator();
  pi->moveRight();
  while (not (pi->fellOff())) {
    Monomial * t = new Monomial(pi->currMonomial());
    (*t) *= u;
    find_position(t2, M_build.end(), r2, *t);
    if (t2 == M_build.end()) {
      M_build.push_back(t);
      R_build.push_back(nullptr);
    } else if (**t2 != *t) {
      M_build.insert(t2, t);
      R_build.insert(r2, nullptr);
    }
    pi->moveRight();
  }
  delete pi;
}

F4_Reduction_Data::~F4_Reduction_Data() {
  free(A);
  for (auto strat : strategies) {
    if (strat != nullptr)
      delete strat;
  }
}

void F4_Reduction_Data::print_matrix(bool show_data) {
  if (show_data) {
    for (auto m : M)
      cout << *m << ", ";
    cout << endl;
  }
  for (unsigned i = 0; i < rows; ++i) {
    cout << "( ";
    for (unsigned j = 0; j < cols; ++j) {
      cout << A[i*cols + j];
      if (j < cols - 1)
        cout << ", ";
    }
    cout << ")\n";
  }
}

void F4_Reduction_Data::list_reducers() {
  for (unsigned i = 0; i < cols; ++i) {
    cout << *(M[i]) << " to be reduced by ";
    if (R[i] == nullptr)
      cout << "none\n";
    else
      cout << R[i]->leading_monomial() << endl;
  }
}

bool F4_Reduction_Data::is_zero() {
  bool is_zero_so_far = true;
  for (unsigned i = 0; is_zero_so_far and i < rows; ++i)
    is_zero_so_far = is_zero_so_far and (nonzero_entries[i] == 0);
  return is_zero_so_far;
}

void F4_Reduction_Data::reduce_by_old() {
  NVAR_TYPE n = M[0]->num_vars();
  EXP_TYPE * u = new EXP_TYPE[n]; // exponents of multiplier
  // loop through rows
  for (unsigned k = 0; k < rows; ++k) {
    // loop through our terms
    for (unsigned i = heads[k]; i < cols; ++i) {
      // do we need to reduce this term, and can we?
      if ((not A[k*cols + i].is_zero()) and R[i] != nullptr) {
        // get reducer for this monomial
        const Abstract_Polynomial * g = R[i];
        Polynomial_Iterator * gi = g->new_iterator();
        // determine multiplier
        for (NVAR_TYPE k = 0; k < n; ++k)
          u[k] = (*M[i])[k] - (gi->currMonomial())[k];
        // determine reduction coefficient
        Prime_Field_Element a(A[k*cols + i]*gi->currCoeff().inverse());
        // loop through g's terms
        // by construction, monomials of u*g[i] should already appear in M,
        // so the line marked *** SHOULD NOT need a guard
        unsigned j = i;
        while (not gi->fellOff()) {
          const Monomial & t = gi->currMonomial();
          while (not M[j]->like_multiple(u, t)) ++j;
          bool was_zero = A[k*cols + j].is_zero();
          A[k*cols + j] -= a*gi->currCoeff();
          if (was_zero and not A[k*cols + j].is_zero())
            ++nonzero_entries[k];
          else if (not was_zero and A[k*cols + j].is_zero())
            --nonzero_entries[k];
          ++j;
          gi->moveRight();
        }
        if (strategies[k] != nullptr)
          strategies[k]->pre_reduction_tasks(u, *g);
        advance_head(k);
        delete gi;
      }
    }
  }
  delete [] u;
}

void F4_Reduction_Data::reduce_by_new() {
  for (unsigned i = 0; i < rows; ++i) {
    if (nonzero_entries[i] > 0) {
      for (unsigned j = 0; j < rows; ++j) {
        if (j != i) {
          auto c = heads[i];
          if (not A[j*cols + c].is_zero()) {
            auto a = A[j*cols + c];
            a *= -A[i*cols + c].inverse();
            while (c < cols) {
              if (not A[i*cols + c].is_zero()) {
                bool was_zero = A[j*cols + c].is_zero();
                A[j*cols + c] += a*A[i*cols + c];
                if (was_zero and not A[j*cols + c].is_zero())
                  ++nonzero_entries[j];
                else if (not was_zero and A[j*cols + c].is_zero())
                  --nonzero_entries[j];
              }
              ++c;
            }
            if (heads[i] == heads[j])
              advance_head(j);
          }
        }
      }
    }
  }
}

vector<Constant_Polynomial *> F4_Reduction_Data::finalize() {
  vector<Constant_Polynomial *> result;
  NVAR_TYPE n = M[0]->num_vars();
  for (unsigned i = 0; i < rows; ++i) {
    if (nonzero_entries[i] == 0) {
      delete strategies[i];
      strategies[i] = nullptr;
    } else {
      Monomial * M_final = (Monomial *)malloc(sizeof(Monomial)*nonzero_entries[i]);
      Prime_Field_Element * A_final
          = (Prime_Field_Element *)malloc(sizeof(Prime_Field_Element)*nonzero_entries[i]);
      unsigned k = 0;
      Prime_Field_Element scale(A[i*cols + heads[i]].inverse(), A[i*cols + heads[i]].field());
      for (unsigned j = heads[i]; k < nonzero_entries[i] and j < cols; ++j) {
        if (not A[i*cols + j].is_zero()) {
          A_final[k] = A[i*cols + j]*scale;
          M_final[k].common_initialization();
          M_final[k].initialize_exponents(n);
          M_final[k].set_monomial_ordering(mord);
          M_final[k] = *M[j];
          ++k;
        }
      }
      result.push_back(new Constant_Polynomial(
          nonzero_entries[i],
          Rx,
          M_final, A_final,
          mord
      ));
      result.back()->set_strategy(strategies[i]);
      strategies[i] = nullptr;
      free(M_final);
      free(A_final);
    }
  }
  return result;
}

list<Constant_Polynomial *> f4_control(const list<Abstract_Polynomial *> &F) {
  unsigned number_of_spolys = 0;
  list<Abstract_Polynomial *> G;
  list<Critical_Pair_Basic *> P;
  // set up basis with generators
  for (Abstract_Polynomial * fo : F)
  {
    Constant_Polynomial * f = new Constant_Polynomial(*fo);
    f->set_strategy(new Poly_Sugar_Data(f));
    f->strategy()->at_generation_tasks();
    P.push_back(new Critical_Pair_Basic(f, StrategyFlags::SUGAR_STRATEGY));
  }
  // main loop
  bool verbose = false;
  bool very_verbose = false;
  list<Critical_Pair_Basic *> Pnew;
  while (not P.empty()) {
    sort_pairs_by_strategy(P);
    report_critical_pairs(P);
    Critical_Pair_Basic * p = P.front();
    report_front_pair(p, StrategyFlags::SUGAR_STRATEGY);
    cout << "\tdegree: " << p->lcm().total_degree(0) << endl;
    P.pop_front();
    Pnew.push_back(p);
    DEG_TYPE mindeg = p->lcm().total_degree();
    for (auto pi = P.begin(); pi != P.end(); ++pi) { 
      p = *pi;
      if (p->lcm().total_degree() < mindeg) {
        for (auto qi = Pnew.begin(); qi != Pnew.end(); ++qi)
          P.push_front(*qi);
        Pnew.clear();
        mindeg = p->lcm().total_degree();
        report_front_pair(p, StrategyFlags::SUGAR_STRATEGY);
        cout << "\tdegree: " << p->lcm().total_degree(0) << endl;
        Pnew.push_back(p);
        auto qi = pi;
        ++qi;
        P.erase(pi);
        pi = qi;
      }
      else if (p->lcm().total_degree() == mindeg) {
        report_front_pair(p, StrategyFlags::SUGAR_STRATEGY);
        Pnew.push_back(p);
        P.erase(pi);
      }
    }
    // make s-poly
    F4_Reduction_Data s(Pnew, G);
    number_of_spolys += Pnew.size();
    Pnew.clear();
    if (not s.is_zero()) {
      s.reduce_by_old();
      s.reduce_by_new();
    }
    if (s.is_zero()) {
      cout << "\tmatrix reduced to zero\n";
      // delete s;
    } else {
      vector<Constant_Polynomial *> R = s.finalize();
      for (auto r : R) {
        cout << "\tadded " << r->leading_monomial() << endl;
        very_verbose = false;
        if (very_verbose) { cout << "\tadded "; r->println(); }
        gm_update(P, G, r, StrategyFlags::SUGAR_STRATEGY);
      }
    }
  }
  cout << number_of_spolys << " s-polynomials computed and reduced\n";
  // cleanup
  cout << G.size() << " polynomials before interreduction\n";
  //check_correctness(G, strategy);
  G = reduce_basis(G);
  cout << G.size() << " polynomials after interreduction\n";
  list<Constant_Polynomial *> B;
  for (Abstract_Polynomial * g : G)
    B.push_back(new Constant_Polynomial(*g));
  return B;
}

#endif