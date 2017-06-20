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
#include <atomic>

Grading_Order_Data_Allocator<Entry> * eoda = nullptr;

class spin_lock {
private: std::atomic_flag l = ATOMIC_FLAG_INIT;
public:
  void lock() { while (l.test_and_set(std::memory_order_acquire)) { } }
  void unlock() { l.clear(std::memory_order_release); }
};

spin_lock eoda_mutex;

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
  if (eoda == nullptr) eoda = new Grading_Order_Data_Allocator<Entry_struct>(2);
  num_cols = 0;
  R_build.clear();
  M_build.clear();
  strategies.clear();
  auto mi = M_build.begin();
  auto ri = R_build.begin();
  mord = P.front()->first()->monomial_ordering();
  for (auto p : P) {
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

void F4_Reduction_Data::initialize_my_rows(
  unsigned row,
  const set<Critical_Pair_Basic *> & my_pairs
) {
  for (auto cp : my_pairs) {
    auto & ai = A[row];
    const Monomial & t = cp->first_multiplier();
    Polynomial_Iterator * pi = cp->first()->new_iterator();
    nonzero_entries[row] = 0;
    auto mi = M_build.begin();
    bool head_found = false;
    for (unsigned i = 0; i < num_cols; ++i) {
      if ((not pi->fellOff()) and M[i]->like_multiple(pi->currMonomial(), t)) {
        eoda_mutex.lock();
        Entry * e = eoda->get_new_block();
        eoda_mutex.unlock();
        e->idx = i; e->val = pi->currCoeff().value();
        ai.push_back(e);
        pi->moveRight();
        ++nonzero_entries[row];
      }
      ++mi;
    }
    delete pi;
    strategies[row] = new Poly_Sugar_Data(cp->first());
    if (cp->second() == nullptr) delete cp->first();
    delete cp;
    ++row;
  }
}

void F4_Reduction_Data::initialize_many(const list<Critical_Pair_Basic *> & P) {
  num_cols = M_build.size();
  Abstract_Polynomial * p = const_cast<Abstract_Polynomial *>(P.front()->first());
  Prime_Field & F = const_cast<Prime_Field &>(p->ground_field());
  num_rows = P.size();
  nonzero_entries.resize(num_rows);
  A.resize(P.size());
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
  for (auto Ak : A) {
    for (auto a : Ak) {
      eoda_mutex.lock();
      eoda->return_used_block(a);
      eoda_mutex.unlock();
    }
  }
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
    unsigned j = 0;
    for (auto a : A[i]) {
      while (j < a->idx) {
        cout << "0, ";
        ++j;
      }
      cout << a->val << ", ";
      ++j;
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
    is_zero_so_far = is_zero_so_far and (nonzero_entries[i] == 0);
  return is_zero_so_far;
}

void F4_Reduction_Data::reduce_my_rows(
  const set<unsigned> & my_rows,
  const list<Entry *> & r,
  vector<list<Entry *>::iterator> & rj
) {
  const Prime_Field & F = Rx.ground_field();
  const UCOEF_TYPE & mod = F.modulus();
  for (unsigned k : my_rows) {
    // loop through our terms
    auto & Ak = A[k];
    //auto ak = Ak.begin();
    auto ak(rj[k]);
    unsigned k0 = (*ak)->idx;
    auto ri = r.begin();
    // determine reduction coefficient
    while ((*ak)->idx != (*ri)->idx) ++ak;
    COEF_TYPE a = ((mod - (*ak)->val) * F.inverse((*ri)->val)) % mod;
    if (a < 0) a += mod;
    // loop through remaining terms
    while (ak != Ak.end() and ri != r.end()) {
      // loop through g's terms
      // by construction, monomials of u*g should already appear in M,
      while (ak != Ak.end() and (*ak)->idx < (*ri)->idx)
        ++ak;
      auto aki = *ak;
      auto rki = *ri;
      while (ri != r.end() and (ak == Ak.end() or (*ri)->idx < aki->idx)) {
        eoda_mutex.lock();
        auto e = eoda->get_new_block();
        eoda_mutex.unlock();
        e->idx = (*ri)->idx; e->val = (a*(*ri)->val) % mod;
        Ak.insert(ak, e);
        ++nonzero_entries[k];
        ++ri;
        rki = *ri;
      }
      rki = *ri;
      if (ri != r.end() and aki->idx == rki->idx) {
        (*ak)->val += a*(*ri)->val; (*ak)->val %= mod;
        if ((*ak)->val != 0)
          ++ak;
        else {
          auto al = ak; ++al;
          eoda_mutex.lock();
          eoda->return_used_block(*ak);
          eoda_mutex.unlock();
          Ak.erase(ak);
          ak = al;
          --nonzero_entries[k];
        }
        ++ri;
      }
    }
    while (ri != r.end()) {
      eoda_mutex.lock();
      auto e = eoda->get_new_block();
      eoda_mutex.unlock();
      e->idx = (*ri)->idx; e->val = (a*(*ri)->val) % mod;
      Ak.push_back(e);
      ++nonzero_entries[k];
      ++ri;
    }
    // recompute rj[k]
    for (ak = Ak.begin(); ak != Ak.end() and (*ak)->idx < k0; ++ak) { }
    rj[k] = ak;
  }
}

void F4_Reduction_Data::reduce_by_old() {
  static double this_time = 0.0;
  NVAR_TYPE n = Rx.number_of_variables();
  list<Entry *> r;
  unsigned cores = std::thread::hardware_concurrency() * 2;
  thread * workers = new thread[cores];
  set<unsigned> * thread_rows = new set<unsigned>[cores];
  unsigned * v = new unsigned[n];
  unsigned * u = new unsigned[n];
  vector<list<Entry *>::iterator> rj(num_rows);
  for (unsigned i = 0; i < num_rows; ++i) rj[i] = A[i].begin();
  for (unsigned j = 0; j < num_cols; ++j) {
    // how many rows need reduction by R[i]?
    if (R[j] != nullptr) {
      list<unsigned> unreduced_rows;
      for (unsigned i = 0; i < num_rows; ++i) {
        if (nonzero_entries[i] != 0) {
          auto ai(rj[i]);
          while (ai != A[i].end() and (*ai)->idx < j)
            ++ai;
          if (ai != A[i].end() and (*ai)->idx == j) {
            unreduced_rows.push_back(i);
          }
          rj[i] = ai;
        }
      }
      if (unreduced_rows.size() > 0) {
        r.clear();
        for (unsigned c = 0; c < cores; ++c)
          thread_rows[c].clear();
        unsigned num_threads = (cores < unreduced_rows.size()) ?
            cores : unreduced_rows.size();
        auto g = R[j];
        auto t = R[j]->leading_monomial();
        for (unsigned k = 0; k < n; ++k)
          v[k] = (*M[j])[k] - t[k];
        unsigned l = j;
        auto gj = R[j]->new_iterator();
        for (/* */; not gj->fellOff(); gj->moveRight()) {
          const Monomial & t = gj->currMonomial();
          for (unsigned k = 0; k < n; ++k)
            u[k] = v[k] + t[k];
          while (not M[l]->is_like(u)) ++l;
          eoda_mutex.lock();
          auto e = eoda->get_new_block();
          eoda_mutex.unlock();
          e->idx = l; e->val = gj->currCoeff().value();
          r.push_back(e);
        }
        delete gj;
        // loop through num_rows
        unsigned ti = 0;
        for (unsigned k : unreduced_rows) {
          thread_rows[ti].insert(k);
          ti += 1; ti %= num_threads;
          strategies[k]->pre_reduction_tasks(v, *R[j]);
        }
        time_t start_time = time(nullptr);
        for (unsigned c = 0; c < num_threads; ++c)
          workers[c] = thread(
              &F4_Reduction_Data::reduce_my_rows, this,
              std::cref(thread_rows[c]), std::cref(r), std::ref(rj)
          );
        for (unsigned c = 0; c < num_threads; ++c)
          workers[c].join();
        time_t end_time = time(nullptr);
        this_time += difftime(end_time, start_time);
      }
    }
  }
  delete [] v;
  delete [] u;
  delete [] workers;
  delete [] thread_rows;
  cout << "total time in reduction: " << this_time << endl;
}

void F4_Reduction_Data::reduce_some_by_one_new(
  unsigned i,
  const Prime_Field & F,
  const set<unsigned> & to_reduce
) {
  auto mod = F.modulus();
  for (auto j : to_reduce) {
    auto ai = A[i].begin();
    auto aj = A[j].begin();
    while ((*aj)->idx < (*ai)->idx) ++aj;
    // determine reduction coefficient
    COEF_TYPE a = ((mod - (*aj)->val) * F.inverse((*ai)->val)) % mod;
    if (a < 0) a += mod;
    while (ai != A[i].end() and aj != A[j].end()) {
      // loop through g's terms
      // by construction, monomials of u*g should already appear in M,
      while (aj != A[j].end() and (*aj)->idx < (*ai)->idx)
        ++aj;
      if (aj != A[j].end()) {
        auto aik = *ai; auto ajk = *aj;
        while (ai != A[i].end() and (*ai)->idx < (*aj)->idx) {
          eoda_mutex.lock();
          auto e = eoda->get_new_block();
          eoda_mutex.unlock();
          e->idx = (*ai)->idx; e->val = (a*(*ai)->val) % mod;
          A[j].insert(aj, e);
          ++ai;
          ++nonzero_entries[j];
        }
        aik = *ai; ajk = *aj;
        if (ai != A[i].end() and (*ai)->idx == (*aj)->idx) {
          (*aj)->val += a*(*ai)->val; (*aj)->val %= mod;
          if ((*aj)->val != 0)
            ++aj;
          else {
            auto ak = aj; ++ak;
            eoda_mutex.lock();
            eoda->return_used_block(*aj);
            eoda_mutex.unlock();
            A[j].erase(aj);
            aj = ak;
            --nonzero_entries[j];
          }
          ++ai;
        }
      }
    }
    while (ai != A[i].end()) {
      eoda_mutex.lock();
      auto e = eoda->get_new_block();
      eoda_mutex.unlock();
      e->idx = (*ai)->idx; e->val = (a*(*ai)->val) % mod;
      A[j].push_back(e);
      ++nonzero_entries[j];
      ++ai;
    }
  }
}

void F4_Reduction_Data::reduce_by_new() {
  auto F = Rx.ground_field();
  unsigned cores = std::thread::hardware_concurrency() * 2;
  thread * workers = new thread[cores];
  set<unsigned> * thread_rows = new set<unsigned>[cores];
  vector<list<Entry *>::iterator> aj(A.size());
  for (unsigned i = 0; i < num_rows; ++i)
    aj[i] = A[i].begin();
  //
  for (unsigned i = 0; i < num_rows; ++i) {
    if (nonzero_entries[i] > 0) {       
      unsigned i0 = A[i].front()->idx;
      list<unsigned> unreduced_rows;
      for (unsigned j = 0; j < num_rows; ++j) {
        if (j != i and nonzero_entries[j] != 0) {
          auto aj = A[j].begin();
          while (aj != A[j].end() and (*aj)->idx < i0) ++aj;
          if (aj != A[j].end() and (*aj)->idx == i0)
            unreduced_rows.push_back(j);
        }
      }
      if (unreduced_rows.size() > 0) {
        for (unsigned c = 0; c < cores; ++c)
          thread_rows[c].clear();
        unsigned num_threads = (cores < unreduced_rows.size()) ?
            cores : unreduced_rows.size();
        unsigned ti = 0;
        for (unsigned k : unreduced_rows) {
          thread_rows[ti].insert(k);
          ti += 1; ti %= num_threads;
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
  auto & F = Rx.ground_field();
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
      Prime_Field_Element scale(F.inverse(A[i].front()->val), &F);
      for (auto p : A[i]) {
        A_final[k] = scale * p->val;
        M_final[k].common_initialization();
        M_final[k].initialize_exponents(n);
        M_final[k].set_monomial_ordering(mord);
        M_final[k] = *M[p->idx];
        ++k;
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
        P.splice(P.begin(), Pnew, Pnew.begin(), Pnew.end());
        Pnew.clear();
        mindeg = p->lcm().total_degree();
        report_front_pair(p, StrategyFlags::SUGAR_STRATEGY);
        cout << "\tdegree: " << p->lcm().total_degree(0) << endl;
        auto qi = pi;
        ++qi;
        Pnew.splice(Pnew.begin(), P, pi);
        pi = qi;
        --pi;
      }
      else if (p->lcm().total_degree() == mindeg) {
        report_front_pair(p, StrategyFlags::SUGAR_STRATEGY);
        auto qi = pi;
        ++qi;
        Pnew.splice(Pnew.begin(), P, pi);
        pi = qi;
        --pi;
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
      cout << "total timed section " << total_time << endl;
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
    B.push_back(static_cast<Constant_Polynomial *>(g));
  time_t end_f4 = time(nullptr);
  double duration = difftime(end_f4, start_f4);
  cout << "computation ended at " << asctime(localtime(&end_f4)) << endl;
  cout << "parallel f4 took " << duration << " seconds\n";
  cout << "parallel f4 spent " << total_time << " seconds in timed section\n";
  return B;
}

#endif