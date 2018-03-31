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
#include <algorithm>
using std::fill;

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
  for (auto cp : P) {
    auto p = cp->first();
    const Monomial & t = cp->first_multiplier();
    auto pi = p->new_iterator();
    vector<pair<unsigned, COEF_TYPE> > & Arow = A[row];
    unsigned i = M_table.lookup_product(pi->currMonomial(), t);
    Arow.emplace_back(i, pi->currCoeff().value());
    pi->moveRight();
    for (/* */; not pi->fellOff(); ++i) {
      i = M_table.lookup_product(pi->currMonomial(), t);
      Arow.emplace_back(i, pi->currCoeff().value());
      pi->moveRight();
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
    } else delete t;
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
    } else delete t;
    pi->moveRight();
  }
  delete pi;
}

F4_Reduction_Data::~F4_Reduction_Data() {
  for (auto strat : strategies) {
    if (strat != nullptr)
      delete strat;
  }
  for (auto t : M_build) delete t;
}

void F4_Reduction_Data::print_row(unsigned i)  {
  for (auto term : A[i])
    cout << term.second << ' ' << *M[term.first] << " + ";
  cout << endl;
}

void F4_Reduction_Data::print_matrix(bool show_data) {
  if (show_data) {
    for (auto m : M)
      cout << *m << ", ";
    cout << endl;
  }
  for (unsigned i = 0; i < num_rows; ++i) {
    cout << "A[" << i << "]: ( ";
    auto & Ai = A[i];
    unsigned j = 0, l = 0;
    while (j < Ai[l].first) {
      cout << "0, ";
      ++j;
    }
    while (l < Ai.size()) {
      cout << Ai[l].second << ", ";
      ++l; ++j;
      if (l < Ai.size())
        while (j < Ai[l].first) {
          cout << "0, ";
          ++j;
        }
    }
    while (j < num_cols) {
      cout << "0, ";
      ++j;
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
    is_zero_so_far = is_zero_so_far and (A[i].size() == 0);
  return is_zero_so_far;
}

double reduction_time = 0;
mutex red_mutex;

void F4_Reduction_Data::densify_row(
    vector<COEF_TYPE> & buffer,
    const vector<pair<unsigned, COEF_TYPE> > & source
) {
  fill(buffer.begin(), buffer.end(), 0);
  for (auto & term : source) {
    buffer[term.first] = term.second;
  }
}

void F4_Reduction_Data::reduce_buffer(
  vector<COEF_TYPE> & buffer,
  const vector<pair<unsigned, COEF_TYPE> > & source,
  COEF_TYPE mod
) {
  auto a = buffer[source[0].first];
  for (auto & term : source) {
    buffer[term.first] += mod - (a * term.second % mod);
    buffer[term.first] %= mod;
  }
}

void F4_Reduction_Data::sparsify_row(
    const vector<COEF_TYPE> & buffer, vector<pair<unsigned, COEF_TYPE> > & source
) {
  source.clear();
  for (unsigned k = 0; k < num_cols; ++k)
    if (buffer[k] != 0)
      source.emplace_back(k, buffer[k]);
}

void F4_Reduction_Data::add_reductor(
    vector<COEF_TYPE> & buffer, COEF_TYPE a,
    const vector<pair<unsigned, COEF_TYPE> > & reductor, COEF_TYPE mod
) {
  for (unsigned i = 0; i < reductor.size(); ++i) {
    auto & ri = reductor[i];
    buffer[ri.first] += mod - (a * ri.second % mod);
    buffer[ri.first] %= mod;
  }
}

void F4_Reduction_Data::build_reducer(vector<COEF_TYPE> & buffer, unsigned mi) {
  NVAR_TYPE n = Rx.number_of_variables();
  const Prime_Field & F = Rx.ground_field();
  UCOEF_TYPE mod = F.modulus();
  const auto g = R[mi];
  Polynomial_Iterator * gi = g->new_iterator();
  EXP_TYPE v[n];
  const Monomial & t = *M[mi];
  const Monomial & u = g->leading_monomial();
  for (NVAR_TYPE k = 0; k < n; ++k)
    v[k] = t[k] - u[k];
  // build needed subreducers first
  // we do this recursively, building only from right-to-left, to save space
  while (gi->canMoveRight()) gi->moveRight(); // don't try to reduce by self
  while (gi->currMonomial() != g->leading_monomial()) {
    auto j = M_table.lookup_product(gi->currMonomial(), v);
    red_lock[j].lock();
    if (R[j] != nullptr and R_built[j].size() == 0)
      build_reducer(buffer, j);
    red_lock[j].unlock();
    gi->moveLeft();
  }
  // having build all necessary reducers, we can now populate & reduce
  buffer.resize(M.size());
  fill(buffer.begin(), buffer.end(), 0);
  gi->restart_iteration();
  while (not gi->fellOff()) {
    auto j = M_table.lookup_product(gi->currMonomial(), v);
    if (R[j] == nullptr)
      buffer[j] += gi->currCoeff().value();
    else
      add_reductor(buffer, gi->currCoeff().value(), R_built[j], mod);
    gi->moveRight();
  }
  delete gi;
  // finalize reducer
  sparsify_row(buffer, R_built[mi]);
}

void F4_Reduction_Data::reduce_my_rows(
  vector<COEF_TYPE> & A_buffer, vector<COEF_TYPE> & R_buffer,
  const set<unsigned> & my_rows
) {
  NVAR_TYPE n = Rx.number_of_variables();
  const Prime_Field & F = Rx.ground_field();
  EXP_TYPE v[n];
  UCOEF_TYPE mod = F.modulus();
  A_buffer.resize(M.size());
  for (unsigned i : my_rows) {
    auto & Ai = A[i];
    fill(A_buffer.begin(), A_buffer.end(), 0);
    // loop through terms of Ai
    for (auto term : Ai) {
      auto mi = term.first; // monomial index
      auto a = term.second; // coefficient
      // is this monomial reducible?
      if (R[mi] == nullptr) {
        A_buffer[mi] += a;
        A_buffer[mi] %= mod;
      } else {
        // can we reduce this monomial?
        red_lock[mi].lock();
        auto & r = R_built[mi];
        if (r.size() == 0)
          build_reducer(R_buffer, mi);
        red_lock[mi].unlock();
        add_reductor(A_buffer, a, r, mod);
        auto & t = *M[mi];
        auto & u = R[mi]->leading_monomial();
        for (unsigned k = 0; k < n; ++k)
          v[k] = t[k] - u[k];
        strategies[i]->pre_reduction_tasks(v, *R[mi]);
      }
    }
    sparsify_row(A_buffer, Ai);
  }
}

void F4_Reduction_Data::reduce_by_old() {
  if (red_lock.size() < num_cols) {
    vector<mutex> new_list(3*num_cols / 2);
    red_lock.swap(new_list);
  }
  unsigned cores = std::thread::hardware_concurrency() * 2;
  unsigned num_threads = (cores < num_rows) ? cores : num_rows;
  set<unsigned> * thread_rows = new set<unsigned>[num_threads];
  A_buf.resize(num_threads);
  R_buf.resize(num_threads);
  // loop through num_rows
  for (unsigned k = 0; k < num_rows; ++k)
    thread_rows[k % num_threads].insert(k);
  thread * workers = new thread[num_threads];
  for (unsigned c = 0; c < num_threads; ++c) {
    workers[c] = thread(
        &F4_Reduction_Data::reduce_my_rows, this, std::ref(A_buf[c]),
            std::ref(R_buf[c]), std::cref(thread_rows[c])
    );
  }
  for (unsigned c = 0; c < num_threads; ++c)
    workers[c].join();
  delete [] workers;
  delete [] thread_rows;
}

bool row_has_term_of_index(
    const vector<pair<unsigned, COEF_TYPE> > & reducer, unsigned index
) {
  for (const auto & term : reducer)
    if (term.first == index)
      return true;
  return false;
}

void F4_Reduction_Data::reduce_my_new_rows(
  vector<COEF_TYPE> & buffer,
  unsigned i,
  const Prime_Field & F,
  const set<unsigned> & to_reduce
) {
  auto mod = F.modulus();
  auto & Ai = A[i];
  unsigned k = Ai[0].first;
  for (auto j : to_reduce) {
    if (i != j) {
      auto & Aj = A[j];
      if (row_has_term_of_index(Aj, k)) {
        densify_row(buffer, Aj);
        reduce_buffer(buffer, Ai, mod);
        sparsify_row(buffer, Aj);
      }
    }
  }
}

void F4_Reduction_Data::reduce_by_new() {
  Prime_Field F(Rx.ground_field());
  auto mod = F.modulus();
  unsigned cores = std::thread::hardware_concurrency() * 2;
  unsigned num_threads = (cores < num_rows) ? cores : num_rows;
  set<unsigned> * thread_rows = new set<unsigned>[num_threads];
  thread * workers = new thread[num_threads];
  for (unsigned i = 0; i < num_rows; ++i) {
    auto & Ai = A[i];
    if (Ai.size() > 0) {
      COEF_TYPE a = Ai[0].second;
      COEF_TYPE b = F.inverse(a);
      for (auto & term : Ai) {
        term.second *= b;
        term.second %= mod;
      }
      for (unsigned j = 0; j < num_threads; ++j)
        thread_rows[j].clear();
      unsigned k = 0;
      for (unsigned j = 0; j < num_rows; ++j) {
        if (j != i and row_has_term_of_index(A[j], Ai[0].first)) {
          thread_rows[k].insert(j);
          ++k; k %= num_threads;
        }
      }
      for (unsigned c = 0; c < num_threads; ++c)
        workers[c] = thread(
            &F4_Reduction_Data::reduce_my_new_rows, this,
            std::ref(A_buf[c]), i, std::cref(F), std::cref(thread_rows[c])
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
  NVAR_TYPE n = M[0]->num_vars();
  for (unsigned i = 0; i < num_rows; ++i) {
    if (A[i].size() == 0) {
      delete strategies[i];
      strategies[i] = nullptr;
    } else {
      const vector<pair<unsigned, COEF_TYPE> > & Ai = A[i];
      Monomial * M_final = static_cast<Monomial *>(
          malloc(sizeof(Monomial)*Ai.size())
      );
      Prime_Field_Element * A_final
          = static_cast<Prime_Field_Element *>(
                malloc(sizeof(Prime_Field_Element)*Ai.size())
            );
      unsigned k = 0;
      for (auto term : Ai) {
          A_final[k].assign(term.second, &F);
          M_final[k].common_initialization();
          M_final[k].initialize_exponents(n);
          M_final[k].set_monomial_ordering(mord);
          M_final[k] = *M[term.first];
          ++k;
      }
      result.push_back(new Constant_Polynomial(
          A[i].size(),
          Rx,
          M_final, A_final,
          mord
      ));
      result.back()->set_strategy(strategies[i]);
      strategies[i] = nullptr;
      for (k = 0; k < A[i].size(); ++k)
        M_final[k].deinitialize();
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
      time_t start_time = time(nullptr);
      s.reduce_by_old(); // cyc8h ~38sec
      time_t end_time = time(nullptr);
      total_time += difftime(end_time, start_time);
      s.reduce_by_new(); // cyc8h ~4sec
    }
    if (s.is_zero()) {
      cout << "\tmatrix reduced to zero\n";
      // delete s;
    } else {
      //time_t start_time = time(nullptr);
      vector<Constant_Polynomial *> R = s.finalize();
      //time_t end_time = time(nullptr);
      //total_time += difftime(end_time, start_time);
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
  unsigned int num_mons = 0;
  for (Abstract_Polynomial * g : G) {
    B.push_back(new Constant_Polynomial(*g));
    num_mons += g->length();
    delete g;
  }
  cout << num_mons << " monomials in basis (possibly counting multiple times)\n";
  time_t end_f4 = time(nullptr);
  double duration = difftime(end_f4, start_f4);
  cout << "computation ended at " << asctime(localtime(&end_f4)) << endl;
  cout << "parallel f4 took " << duration << " seconds\n";
  cout << "parallel f4 spent " << total_time << " seconds in timed section\n";
  return B;
}

#endif