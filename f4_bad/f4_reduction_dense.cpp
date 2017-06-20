#ifndef __F4_REDUCTION_CPP__
#define __F4_REDUCTION_CPP__

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

#include "f4_reduction.hpp"
#include "algorithm_buchberger_basic.hpp"

#include <thread>
using std::thread;

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
  num_cols = 0;
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

std::mutex cout_mutex;

void F4_Reduction_Data::initialize_my_rows(
  unsigned row,
  const set<Critical_Pair_Basic *> & my_pairs
) {
  for (auto cp : my_pairs) {
    const Monomial & t = cp->first_multiplier();
    Constant_Polynomial_Iterator pi(static_cast<const Constant_Polynomial *>(cp->first()));
    nonzero_entries[row] = 0;
    auto mi = M_build.begin();
    bool head_found = false;
    for (unsigned i = 0; i < num_cols; ++i) {
      if ((not pi.fellOff()) and M[i]->like_multiple(pi.currMonomial(), t)) {
        A[row*num_cols + i] = pi.currCoeff().value();
        pi.moveRight();
        ++nonzero_entries[row];
        if (not head_found) {
          heads[row] = i;
          head_found = true;
        }
      }
      ++mi;
    }
    strategies[row] = new Poly_Sugar_Data(cp->first());
    delete cp;
    ++row;
  }
}

void F4_Reduction_Data::initialize_many(const list<Critical_Pair_Basic *> & P) {
  num_cols = M_build.size();
  Abstract_Polynomial * p = const_cast<Abstract_Polynomial *>(P.front()->first());
  Prime_Field & F = const_cast<Prime_Field &>(p->ground_field());
  num_rows = P.size();
  heads.resize(num_rows);
  nonzero_entries.resize(num_rows);
  for (unsigned i = 0; i < num_rows; ++i)
    heads[i] = 0;
  A = (COEF_TYPE *)calloc(
      num_cols*num_rows, sizeof(COEF_TYPE)
  );
  M.clear();
  for (auto mi = M_build.begin(); mi != M_build.end(); ++mi)
    M.push_back(*mi);
  strategies.resize(P.size());
  for (unsigned i = 0; i < strategies.size(); ++i)
    strategies[i] = nullptr;
  R.resize(R_build.size());
  unsigned i = 0;
  for (Abstract_Polynomial * r : R_build)
    R[i++] = r;
  //
  unsigned cores = std::thread::hardware_concurrency() * 2;
  unsigned num_threads = (cores < num_rows) ? cores : num_rows;
  set<Critical_Pair_Basic *> * thread_pairs = new set<Critical_Pair_Basic *>[num_threads];
  vector<unsigned> thread_rows(num_threads);
  // loop through num_rows
  unsigned k = 0;
  for (auto cp : P) {
    thread_pairs[k % num_threads].insert(cp);
    ++k;
  }
  unsigned start_row = 0;
  for (unsigned j = 0; j < num_threads; ++j) {
    thread_rows[j] = start_row;
    start_row += thread_pairs[j].size();
  }
  thread * workers = new thread[num_threads];
  cout << "initializing with " << num_threads << " threads\n";
  for (unsigned c = 0; c < num_threads; ++c)
    workers[c] = thread(
        &F4_Reduction_Data::initialize_my_rows, this, thread_rows[c], std::cref(thread_pairs[c])
    );
  for (unsigned c = 0; c < num_threads; ++c)
    workers[c].join();
  delete [] workers;
  delete [] thread_pairs;
  //
  unsigned nonzero_els = 0;
  for (unsigned i = 0; i < num_rows; ++i)
    nonzero_els += nonzero_entries[i];
  cout << "sparsity: " << (double)nonzero_els * 100. / double(num_rows * num_cols) << endl;
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
  for (unsigned i = 0; i < num_rows; ++i) {
    cout << "( ";
    for (unsigned j = 0; j < num_cols; ++j) {
      cout << A[i*num_cols + j];
      if (j < num_cols - 1)
        cout << ", ";
    }
    cout << ")\n";
  }
}

void F4_Reduction_Data::list_reducers() {
  for (unsigned i = 0; i < num_cols; ++i) {
    cout << *(M[i]) << " to be reduced by ";
    if (R[i] == nullptr)
      cout << "none\n";
    else
      cout << R[i]->leading_monomial() << endl;
  }
}

bool F4_Reduction_Data::is_zero() {
  bool is_zero_so_far = true;
  for (unsigned i = 0; is_zero_so_far and i < num_rows; ++i)
    is_zero_so_far = is_zero_so_far and (nonzero_entries[i] == 0);
  return is_zero_so_far;
}

void F4_Reduction_Data::reduce_my_rows(const set<unsigned> & my_rows) {
  const NVAR_TYPE & n = Rx.number_of_variables();
  const Prime_Field & F = Rx.ground_field();
  const UCOEF_TYPE & mod = F.modulus();
  EXP_TYPE * u = new EXP_TYPE[n]; // exponents of multiplier
  EXP_TYPE * v = new EXP_TYPE[n];
  for (unsigned k : my_rows) {
    // loop through our terms
    const unsigned k0 = k*num_cols;
    for (unsigned i = heads[k]; i < num_cols; ++i) {
      // do we need to reduce this term, and can we?
      if ((not (A[k0 + i] == 0)) and R[i] != nullptr) {
        // get reducer for this monomial
        const Abstract_Polynomial * g = R[i];
        Polynomial_Iterator * gi = g->new_iterator();
        // determine multiplier
        for (NVAR_TYPE k = 0; k < n; ++k)
          u[k] = (*M[i])[k] - (gi->currMonomial())[k];
        // determine reduction coefficient
        COEF_TYPE a = (A[k0 + i] * F.inverse(gi->currCoeff().value())) % mod;
        if (a < 0) a += mod;
        //Prime_Field_Element a(A[k*num_cols + i]*gi->currCoeff().inverse());
        // loop through g's terms
        // by construction, monomials of u*g should already appear in M,
        unsigned j = i;
        while (not gi->fellOff()) {
          const EXP_TYPE * t = gi->currMonomial().log();
          for (unsigned l = 0; l < n; ++l) v[l] = t[l] + u[l];
          while (not M[j]->is_like(v)) ++j;
          unsigned pos = k0 + j;
          bool was_zero = (A[pos] == 0);
          A[pos] -= a*(gi->currCoeff().value());
          while (A[pos] < 0)
            A[pos] += mod;
          if (was_zero and (A[pos] != 0))
            ++nonzero_entries[k];
          else if (not was_zero and (A[pos] == 0))
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
  delete [] u; delete [] v;
}

void F4_Reduction_Data::reduce_by_old() {
  unsigned cores = std::thread::hardware_concurrency() * 2;
  unsigned num_threads = (cores < num_rows) ? cores : num_rows;
  set<unsigned> * thread_rows = new set<unsigned>[num_threads];
  // loop through num_rows
  for (unsigned k = 0; k < num_rows; ++k)
    thread_rows[k % num_threads].insert(k);
  thread * workers = new thread[num_threads];
  for (unsigned c = 0; c < num_threads; ++c)
    workers[c] = thread(
        &F4_Reduction_Data::reduce_my_rows, this, std::cref(thread_rows[c])
    );
  for (unsigned c = 0; c < num_threads; ++c)
    workers[c].join();
  delete [] workers;
  delete [] thread_rows;
}

void F4_Reduction_Data::reduce_some_by_one_new(
  unsigned i,
  const Prime_Field & F,
  const set<unsigned> & to_reduce
) {
  auto mod = F.modulus();
  auto i0 = i*num_cols;
  auto ihead = heads[i];
  for (auto j : to_reduce) {
    auto c = ihead;
    auto j0 = j*num_cols;
    if (A[j0 + c] != 0) {
      auto a = A[j0 + c];
      a *= -F.inverse(A[i0 + c]);
      while (c < num_cols) {
        if (not (A[i0 + c] == 0)) {
          bool was_zero = (A[j0 + c] == 0);
          A[j0 + c] += a*A[i0 + c];
          A[j0 + c] %= mod;
          while (A[j0 + c] < 0)
            A[j0 + c] += mod;
          if (was_zero and not (A[j0 + c] == 0))
            ++nonzero_entries[j];
          else if (not was_zero and A[j0 + c] == 0)
            --nonzero_entries[j];
        }
        ++c;
      }
      if (heads[i] == heads[j])
        advance_head(j);
    }
  }
}

void F4_Reduction_Data::reduce_by_new() {
  auto F = Rx.ground_field();
  unsigned cores = std::thread::hardware_concurrency() * 2;
  unsigned num_threads = (cores < num_rows) ? cores : num_rows;
  set<unsigned> * thread_rows = new set<unsigned>[num_threads];
  thread * workers = new thread[num_threads];
  for (unsigned i = 0; i < num_rows; ++i) {
    if (nonzero_entries[i] > 0) {
      for (unsigned j = 0; j < num_threads; ++j)
        thread_rows[j].clear();
      for (unsigned j = 0; j < num_rows; ++j) {
        if (j != i) {
          thread_rows[j % num_threads].insert(j);
        }
      }
      for (unsigned c = 0; c < num_threads; ++c)
        workers[c] = thread(
            &F4_Reduction_Data::reduce_some_by_one_new, this,
            i, std::cref(F), std::cref(thread_rows[c])
        );
      for (unsigned c = 0; c < num_threads; ++c)
        workers[c].join();
    }
  }
  delete [] workers;
  delete [] thread_rows;
}

vector<Constant_Polynomial *> F4_Reduction_Data::finalize() {
  unsigned nonzero_els = 0;
  for (unsigned i = 0; i < num_rows; ++i) {
    nonzero_els += nonzero_entries[i];
  }
  cout << "sparsity: " << (double)nonzero_els * 100. / double(num_rows * num_cols) << endl;
  auto F = Rx.ground_field();
  vector<Constant_Polynomial *> result;
  NVAR_TYPE n = M[0]->num_vars();
  for (unsigned i = 0; i < num_rows; ++i) {
    if (nonzero_entries[i] == 0) {
      delete strategies[i];
      strategies[i] = nullptr;
    } else {
      Monomial * M_final = (Monomial *)malloc(sizeof(Monomial)*nonzero_entries[i]);
      Prime_Field_Element * A_final
          = (Prime_Field_Element *)malloc(sizeof(Prime_Field_Element)*nonzero_entries[i]);
      unsigned k = 0;
      Prime_Field_Element scale(F.inverse(A[i*num_cols + heads[i]]), &F);
      for (unsigned j = heads[i]; k < nonzero_entries[i] and j < num_cols; ++j) {
        if (not (A[i*num_cols + j] == 0)) {
          A_final[k] = scale * A[i*num_cols + j];
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
  time_t start_time, end_time;
  double total_time = 0;
  time_t start_f4 = time(nullptr);
  cout << "computation started at " << asctime(localtime(&start_f4)) << endl;
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
      start_time = time(nullptr);
      s.reduce_by_old();
      end_time = time(nullptr);
      total_time += difftime(end_time, start_time);
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
  time_t end_f4 = time(nullptr);
  double duration = difftime(end_f4, start_f4);
  cout << "computation ended at " << asctime(localtime(&end_f4)) << endl;
  cout << "parallel f4 took " << duration << " seconds\n";
  cout << "parallel f4 spent " << total_time << " seconds in timed section\n";
  return B;
}

#endif