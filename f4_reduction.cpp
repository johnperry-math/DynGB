#ifndef __F4_REDUCTION_CPP__
#define __F4_REDUCTION_CPP__

double copy_time = 0;

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

#include "f4_reduction.hpp"
#include "algorithm_buchberger_basic.hpp"

#include <thread>
using std::thread;
#include <mutex>
using std::mutex;

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
) :
    G(B), Rx(P.front()->first()->base_ring()),
    M_table(P.front()->first()->base_ring().number_of_variables())
{
  num_cols = 0;
  A.clear();
  head.clear();
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
  unsigned els = 0;
  for (auto Ak : A) els += Ak.size();
  cout << "saved " << (num_rows*num_cols - els)*100 / (num_rows*num_cols) << "%\n";
}

void F4_Reduction_Data::initialize_some_rows(
    const list<Critical_Pair_Basic *> & P, unsigned row
) {
  const unsigned num_cols = M_build.size();
  const COEF_TYPE F0 = 0;
  for (auto cp : P) {
    auto p = cp->first();
    const Monomial & t = cp->first_multiplier();
    auto pi = p->new_iterator();
    vector<COEF_TYPE> & Arow = A[row];
    nonzero_entries[row] = 0;
    unsigned i = M_table.lookup_product(pi->currMonomial(), t);
    Arow.resize(num_cols - i, F0);
    Arow[0] = pi->currCoeff().value();
    offset[row] = i;
    head[row] = 0;
    nonzero_entries[row] = 1;
    pi->moveRight();
    for (/* */; not pi->fellOff(); ++i) {
      i = M_table.lookup_product(pi->currMonomial(), t);
      Arow[i - offset[row]] = pi->currCoeff().value();
      pi->moveRight();
      ++nonzero_entries[row];
    }
    delete pi;
    strategies[row] = new Poly_Sugar_Data(cp->first());
    delete cp;
    ++row;
  }
}

void F4_Reduction_Data::initialize_many(const list<Critical_Pair_Basic *> & P) {
  num_cols = M_build.size();
  num_rows = P.size();
  strategies.resize(num_rows);
  nonzero_entries.resize(num_rows);
  R_built.resize(num_cols);
  num_readers.assign(num_cols, 0);
  M.clear();
  size_t m = 0;
  for (auto mi = M_build.begin(); mi != M_build.end(); ++mi) {
    M.push_back(*mi);
    M_table.add_monomial(*mi, m);
    ++m;
  }
  for (unsigned i = 0; i < strategies.size(); ++i)
    strategies[i] = nullptr;
  A.resize(P.size());
  head.resize(P.size());
  offset.resize(P.size());
  R.resize(R_build.size());
  unsigned i = 0;
  for (Abstract_Polynomial * r : R_build)
    R[i++] = r;
  //
  unsigned cores = std::thread::hardware_concurrency();
  unsigned num_threads = (cores < num_rows) ? cores : num_rows;
  list<Critical_Pair_Basic *> * thread_rows
    = new list<Critical_Pair_Basic *>[num_threads];
  // loop through num_rows
  unsigned k = 0;
  for (auto pi : P) {
    thread_rows[k % num_threads].push_back(pi);
    ++k;
  }
  unsigned start_work[num_threads];
  unsigned start_row = 0;
  for (unsigned i = 0; i < num_threads; ++i) {
    start_work[i] = start_row;
    cout << "thread " << i << " starts at " << start_row << endl;
    start_row += thread_rows[i].size();
  }
  thread * workers = new thread[num_threads];
  for (unsigned c = 0; c < num_threads; ++c) {
    workers[c] = thread(
        &F4_Reduction_Data::initialize_some_rows, this, std::cref(thread_rows[c]),
        start_work[c]
    );
  }
  for (unsigned c = 0; c < num_threads; ++c)
    workers[c].join();
  delete [] workers;
  delete [] thread_rows;
  //
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
    cout << "A[" << i << "]: ( ";
    unsigned j;
    for (j = 0; j < offset[i] + head[i]; ++j)
      cout << "0, ";
    for (/* */; j < offset[i] + A[i].size(); ++j)
      cout << A[i][j - offset[i]] << ", ";
    for (/* */; j < num_cols; ++j)
      cout << "0, ";
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

double reduction_time = 0;
mutex red_mutex;

void F4_Reduction_Data::reduce_my_rows(
  const set<unsigned> & my_rows
) {
  NVAR_TYPE n = Rx.number_of_variables();
  const Prime_Field & F = Rx.ground_field();
  red_mutex.lock();
  Monomial u(n);
  red_mutex.unlock();
  UCOEF_TYPE mod = F.modulus();
  unsigned new_nonzero_entries;
  for (unsigned mi = 0; mi < num_cols; ++mi) {
    // is this monomial reducible?
    if (R[mi] != nullptr) {
      for (unsigned k : my_rows) {
        vector<COEF_TYPE> & Ak = A[k];
        // do we need to reduce this poly?
        if (offset[k] <= mi and Ak[mi - offset[k]] != 0) {
          unsigned i = mi - offset[k];
          // get reducer for this monomial
          const Abstract_Polynomial * g = R[mi];
          Polynomial_Iterator * gi = g->new_iterator();
          // determine multiplier
          Monomial & t = *M[mi];
          for (NVAR_TYPE l = 0; l < n; ++l)
            u.set_exponent(l, t[l] - (gi->currMonomial())[l]);
          red_mutex.lock();
          vector<pair<unsigned, COEF_TYPE> > & r = R_built[mi];
          if (r.size() == 0) { // need to create reducer
            size_t j = mi;
            // loop through g's terms
            while (not gi->fellOff()) {
              const Monomial & t = gi->currMonomial();
              j = M_table.lookup_product(u, t);
              r.emplace_back(j, gi->currCoeff().value());
              gi->moveRight();
            }
          }
          delete gi;
          red_mutex.unlock();
          // prepare for reduction
          new_nonzero_entries = nonzero_entries[k];
          auto sk = strategies[k];
          sk->pre_reduction_tasks(u, *g);
          // determine reduction coefficient
          COEF_TYPE a(mod - ((Ak[i]*F.inverse(r[0].second)) % mod));
          // loop through g's terms
          for (auto ri : r) {
            unsigned j = ri.first;
            unsigned l = j - offset[k];
            bool was_zero = Ak[l] == 0;
            Ak[l] += a*ri.second; Ak[l] %= mod;
            if (was_zero)
              ++new_nonzero_entries;
            else if (Ak[l] == 0)
              --new_nonzero_entries;
            // advance
          }
          unsigned & hk = head[k];
          while (hk < Ak.size() and Ak[hk] == 0) ++hk;
          nonzero_entries[k] = new_nonzero_entries;
        }
      }
    }
  }
}

void F4_Reduction_Data::reduce_by_old() {
  unsigned cores = std::thread::hardware_concurrency() * 2;
  unsigned num_threads = (cores < num_rows) ? cores : num_rows;
  set<unsigned> * thread_rows = new set<unsigned>[num_threads];
  // loop through num_rows
  for (unsigned k = 0; k < num_rows; ++k)
    thread_rows[k % num_threads].insert(k);
  thread * workers = new thread[num_threads];
  for (unsigned c = 0; c < num_threads; ++c) {
    workers[c] = thread(
        &F4_Reduction_Data::reduce_my_rows, this, std::cref(thread_rows[c])
    );
  }
  for (unsigned c = 0; c < num_threads; ++c)
    workers[c].join();
  //reduce_my_rows(thread_rows[0], buffer[0], indices[0]);
  delete [] workers;
  delete [] thread_rows;
}

void F4_Reduction_Data::reduce_my_new_rows(
  unsigned i,
  const Prime_Field & F,
  const set<unsigned> & to_reduce
) {
  auto mod = F.modulus();
  const auto & Ai = A[i];
  for (auto j : to_reduce) {
    auto ci = head[i];
    auto & Aj = A[j];
    auto c = ci + offset[i] - offset[j]; // pos in A[j] of A[i]'s head
    auto a = Aj[c];
    unsigned ops = 0;
    a *= mod - F.inverse(Ai[head[i]]);
    while (c < Aj.size() and ops < nonzero_entries[i]) {
      if (Ai[ci] != 0) {
        bool was_zero = (Aj[c] == 0);
        Aj[c] += a*Ai[ci]; Aj[c] %= mod;
        if (was_zero)
          ++nonzero_entries[j];
        else if (Aj[c] == 0)
          --nonzero_entries[j];
        ++ops;
      }
      ++c; ++ci;
    }
    unsigned & hj = head[j];
    while (hj < Aj.size() and Aj[hj] == 0) ++hj;
  }
}

void F4_Reduction_Data::reduce_by_new() {
  Prime_Field F(Rx.ground_field());
  unsigned cores = std::thread::hardware_concurrency() * 2;
  unsigned num_threads = (cores < num_rows) ? cores : num_rows;
  set<unsigned> * thread_rows = new set<unsigned>[num_threads];
  thread * workers = new thread[num_threads];
  for (unsigned i = 0; i < num_rows; ++i) {
    if (nonzero_entries[i] > 0) {
      COEF_TYPE Ai0 = head[i] + offset[i]; // abs pos of A[i]'s head
      for (unsigned j = 0; j < num_threads; ++j)
        thread_rows[j].clear();
      unsigned k = 0;
      for (unsigned j = 0; j < num_rows; ++j) {
        if (
            j != i and Ai0 > offset[j]
            and (A[j][Ai0 - offset[j]] != 0)
        ) {
          thread_rows[k].insert(j);
          ++k; k %= num_threads;
        }
      }
      for (unsigned c = 0; c < num_threads; ++c)
        workers[c] = thread(
            &F4_Reduction_Data::reduce_my_new_rows, this,
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
  cout << "spent " << copy_time << " seconds copying\n";
  vector<Constant_Polynomial *> result;
  const Prime_Field & F = Rx.ground_field();
  UCOEF_TYPE mod = F.modulus();
  NVAR_TYPE n = M[0]->num_vars();
  for (unsigned i = 0; i < num_rows; ++i) {
    if (nonzero_entries[i] == 0) {
      delete strategies[i];
      strategies[i] = nullptr;
    } else {
      vector<COEF_TYPE> & Ai = A[i];
      Monomial * M_final = static_cast<Monomial *>(
          malloc(sizeof(Monomial)*nonzero_entries[i])
      );
      Prime_Field_Element * A_final
          = static_cast<Prime_Field_Element *>(
                malloc(sizeof(Prime_Field_Element)*nonzero_entries[i])
            );
      COEF_TYPE scale = F.inverse(Ai[head[i]]);
      unsigned k = 0;
      for (unsigned j = head[i]; k < nonzero_entries[i]; ++j) {
        if (Ai[j] != 0) {
          A_final[k].assign(Ai[j]*scale % mod, &F);
          M_final[k].common_initialization();
          M_final[k].initialize_exponents(n);
          M_final[k].set_monomial_ordering(mord);
          M_final[k] = *M[offset[i] + j];
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
  time_t start_f4 = time(nullptr);
  cout << "computation started at " << asctime(localtime(&start_f4)) << endl;
  unsigned number_of_spolys = 0;
  double total_time = 0;
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
    for (auto pi = P.begin(); pi != P.end(); /* */) { 
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
      } else if (p->lcm().total_degree() == mindeg) {
        report_front_pair(p, StrategyFlags::SUGAR_STRATEGY);
        Pnew.push_back(p);
        auto qi = pi;
        ++qi;
        P.erase(pi);
        pi = qi;
      } else
        ++pi;
    }
    // make s-poly
    //time_t start_time = time(nullptr);
    F4_Reduction_Data s(Pnew, G); // cyc8h ~20sec
    //time_t end_time = time(nullptr);
    //total_time += difftime(end_time, start_time);
    number_of_spolys += Pnew.size();
    Pnew.clear();
    if (not s.is_zero()) {
      s.reduce_by_old(); // cyc8h ~38sec
      s.reduce_by_new(); // cyc8h ~4sec
    }
    if (s.is_zero()) {
      cout << "\tmatrix reduced to zero\n";
      // delete s;
    } else {
      time_t start_time = time(nullptr);
      vector<Constant_Polynomial *> R = s.finalize();
      time_t end_time = time(nullptr);
      total_time += difftime(end_time, start_time);
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