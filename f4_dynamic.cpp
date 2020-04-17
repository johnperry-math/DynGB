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
* DynGB is distributed in the hope that it will be useful,                    *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with DynGB. If not, see <http://www.gnu.org/licenses/>.               *
\*****************************************************************************/

#include "f4_dynamic.hpp"
#include "algorithm_buchberger_basic.hpp"

using std::cerr;

#include <fstream>
using std::ofstream;

#include <future>
using std::future; using std::async;

#include <thread>
using std::thread;
#include <mutex>
using std::mutex;
#include <bitset>
using std::bitset;
#include <string>
using std::to_string;
#include <algorithm>
using std::fill;
#include <cstdlib>
using std::srand; using std::rand;

#include "lp_solver.hpp"
using LP_Solvers::LP_Solver;
#include "skeleton.hpp"
using LP_Solvers::Skeleton;
#include "glpk_solver.hpp"
using LP_Solvers::GLPK_Solver;
#include "ppl_solver.hpp"
using LP_Solvers::PPL_Solver;

using Dynamic_Engine::less_by_hilbert_then_degree;
using Dynamic_Engine::less_by_hilbert;
using Dynamic_Engine::less_by_degree_then_hilbert;
using Dynamic_Engine::less_by_grad_hilbert_then_degree;
using Dynamic_Engine::less_by_smoothest_degrees;
using Dynamic_Engine::less_by_largest_max_component;
using Dynamic_Engine::less_by_num_crit_pairs;
using Dynamic_Engine::less_by_betti;
using Dynamic_Engine::less_by_grad_betti;
using Dynamic_Engine::less_by_degree_then_grad_hilbert;
using Dynamic_Engine::less_by_random;
using Dynamic_Engine::verify_and_modify_if_necessary;

extern list<Abstract_Polynomial *> reduce_basis(list<Abstract_Polynomial *>G);
extern void report_front_pair(Critical_Pair_Basic *p, StrategyFlags strategy);
extern void gm_update_dynamic(
    list<Critical_Pair_Dynamic *> & P,
    list<Abstract_Polynomial *> & G,
    Abstract_Polynomial * r,
    StrategyFlags strategy,
    ORDERING_TYPE * ordering
);

extern template void report_critical_pairs<Critical_Pair_Dynamic>(
    const list<Critical_Pair_Dynamic *>, bool
);

extern template void sort_pairs_by_strategy<Critical_Pair_Basic>(
    list<Critical_Pair_Basic *> &
);

F4_Reduction_Data::F4_Reduction_Data(
    WGrevlex * curr_ord,
    const list<Critical_Pair_Dynamic *> & P,
    const list<Abstract_Polynomial *> & B,
    Dynamic_Heuristic method
) :
    G(B), Rx(P.front()->first()->base_ring()), heur(method),
    M_table(P.front()->first()->base_ring().number_of_variables()),
    mord(curr_ord)
{
  static double overall_time = 0;
  static double adding_time = 0;
  static double initializing_time = 0;
  time_t ostart = time(nullptr);
  auto mod = Rx.ground_field().modulus();
  num_cols = 0;
  // set up the heuristic
  switch(heur) {
    case Dynamic_Heuristic::ORD_HILBERT_THEN_DEG:
      heuristic_judges_smaller = less_by_hilbert_then_degree; break;
    case Dynamic_Heuristic::ORD_HILBERT_THEN_LEX:
      heuristic_judges_smaller = less_by_hilbert; break;
    case Dynamic_Heuristic::DEG_THEN_ORD_HILBERT:
      heuristic_judges_smaller = less_by_degree_then_hilbert; break;
    case Dynamic_Heuristic::GRAD_HILB_THEN_DEG:
      heuristic_judges_smaller = less_by_grad_hilbert_then_degree; break;
    case Dynamic_Heuristic::SMOOTHEST_DEGREES:
      heuristic_judges_smaller = less_by_smoothest_degrees; break;
    case Dynamic_Heuristic::LARGEST_MAX_COMPONENT:
      heuristic_judges_smaller = less_by_largest_max_component; break;
    case Dynamic_Heuristic::MIN_CRIT_PAIRS:
      heuristic_judges_smaller = less_by_num_crit_pairs; break;
    case Dynamic_Heuristic::BETTI_HILBERT_DEG:
      heuristic_judges_smaller = less_by_betti; break;
    case Dynamic_Heuristic::GRAD_BETTI_HILBERT_DEG:
      heuristic_judges_smaller = less_by_grad_betti; break;
    default:
      heuristic_judges_smaller = less_by_hilbert; break;
  }
  // set up the matrix
  A.clear();
  strategies.clear();
  NVAR_TYPE n = Rx.number_of_variables();
  for (auto p : P) {
    Poly_Sugar_Data * new_sugar = new Poly_Sugar_Data(p->first());
    strategies.push_back(new_sugar);
    if (strategies.back() != nullptr) {
      strategies.back()->at_generation_tasks(p->first_multiplier());
      if (p->second() != nullptr) {
        auto p2log = p->second_multiplier().log();
        strategies.back()->pre_reduction_tasks(p2log, *(p->second()));
        delete [] p2log;
      }
    }
    add_monomials(curr_ord, p->first(), p->first_multiplier(), true);
    if (p->second() != nullptr) {
      add_monomials(curr_ord, p->second(), p->second_multiplier());
      M_builder[const_cast<Monomial *>(&(p->lcm()))] = const_cast<Abstract_Polynomial *>(p->second());
    }
  }
  // for each monomial, find an appropriate reducer
  for (auto mi = M_builder.rbegin(); mi != M_builder.rend(); ++mi) {
    auto g = G.begin();
    bool found = mi->second != nullptr;
    //cout << "for " << *mi->first;
    while (not found and g != G.end()) {
      if (mi->first->divisible_by((*g)->leading_monomial())) {
        found = true;
        //cout << " selected " << (*g)->leading_monomial() << endl;
        mi->second = *g;
        Monomial u(*(mi->first), (*g)->leading_monomial(), false);
        time_t astart = time(nullptr);
        add_monomials(curr_ord, *g, u);
        time_t aend = time(nullptr);
        adding_time += difftime(aend, astart);
        g = G.end();
      }
      ++g;
    }
    //if (not found) cout << " no reducer found\n";
  }
  cout << "adding monomials time " << adding_time << endl;
  time_t istart = time(nullptr);
  initialize_many(P);
  time_t iend = time(nullptr);
  initializing_time += difftime(iend, istart);
  cout << "initializing time " << initializing_time << endl;
  unsigned els = 0;
  for (auto Ak : A) els += Ak.size();
  cout << "saved " << (num_rows*num_cols - els)*100 / (num_rows*num_cols) << "%\n";
  time_t oend = time(nullptr);
  overall_time += difftime(oend, ostart);
  cout << "overall time in setup " << overall_time << endl;
}

mutex print_lock;

bool sort_by_first(pair<unsigned, COEF_TYPE> a, pair<unsigned, COEF_TYPE> b) {
  return a.first < b.first;
}

void F4_Reduction_Data::initialize_some_rows(
    const list<Critical_Pair_Dynamic *> & P, unsigned row
) {
  const unsigned num_cols = M_builder.size();
  const COEF_TYPE F0 = 0;
  for (auto cp : P) {
    auto p = cp->first();
    nonzero_entries[row] = p->length();
    //print_lock.lock();
    //cout << "initializing row " << row << " for " << cp->lcm() << " with poly " << *p << " and poly ";
    //if (cp->second() == nullptr) cout << "0\n"; else cout << *cp->second() << endl;
    //print_lock.unlock();
    const Monomial & t = cp->first_multiplier();
    unsigned j = 0;
    auto & Arow = A[row];
    Arow.resize(p->length());
    auto pi = p->new_iterator();
    for (/* */; not pi->fellOff(); pi->moveRight()) {
      auto i = M_table.lookup_product(pi->currMonomial(), t);
      Arow[j].first = i;
      Arow[j].second = pi->currCoeff().value();
      ++j;
    }
    sort(Arow.begin(), Arow.end(), sort_by_first);
    delete pi;
    delete cp;
    ++row;
  }
}

void F4_Reduction_Data::initialize_many(const list<Critical_Pair_Dynamic *> & P) {
  num_cols = M_builder.size();
  num_rows = P.size();
  cout << "Initializing " << num_rows << " x " << num_cols << " basic matrix\n";
  dirty.resize(num_rows, true);
  last_compatible_ordering.resize(num_rows, nullptr);
  strategies.resize(num_rows);
  nonzero_entries.resize(num_rows);
  R_built.resize(num_cols);
  for (auto m : M) delete m;
  M.clear(); M.resize(M_builder.size());
  R.clear(); R.resize(M_builder.size());
  vector<mutex> new_red_mutex(M_builder.size());
  red_mutex.swap(new_red_mutex);
  size_t m = 0;
  for (auto mi = M_builder.rbegin(); mi != M_builder.rend(); ++mi) {
    M[m] = mi->first;
    R[m] = mi->second;
    M_table.update_location(mi->first, m);
    ++m;
  }
  A.resize(P.size());
  pref_head.resize(P.size());
  pp_weights.resize(M.size());
  compatible_pps.resize(P.size());
  potential_ideals.resize(P.size());
  unsigned cores = std::thread::hardware_concurrency();
  unsigned num_threads = (cores < num_rows) ? cores : num_rows;
  list<Critical_Pair_Dynamic *> * thread_rows
    = new list<Critical_Pair_Dynamic *>[num_threads];
  // loop through num_rows
  unsigned k = 0, rows_added = 0;
  srand(time(NULL));
  for (auto pi : P) {
    thread_rows[k].push_back(pi);
    ++rows_added;
    if (rows_added > P.size() / num_threads + 1) { ++k; rows_added = 0; }
    //thread_rows[rand() % num_threads].push_back(pi);
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
}

void F4_Reduction_Data::print_builder() {
  cout << "[ ";
  for (auto m : M_builder)
    cout << "( " << *(m.first) << ", " << m.second << " ) ";
  cout << "]\n";
}

double emplace_time;

void F4_Reduction_Data::add_monomials(
    const WGrevlex * curr_ord,
    const Abstract_Polynomial *g,
    const Monomial & u,
    bool new_row
) {
  Polynomial_Iterator * pi = g->new_iterator();
  if (not new_row) pi->moveRight();
  while (not (pi->fellOff())) {
    bool already_there = M_table.contains_product(pi->currMonomial(), u);
    if (not already_there) {
      Monomial * t = new Monomial(pi->currMonomial(), u);
      t->set_monomial_ordering(curr_ord);
      M_table.add_monomial(t);
      M_builder.emplace(t, nullptr);
    }
    pi->moveRight();
  }
  //cout << "processed " << monomials_processed << " monomials\n";
  delete pi;
}

F4_Reduction_Data::~F4_Reduction_Data() {
  for (auto strat : strategies) {
    if (strat != nullptr)
      delete strat;
  }
  for (auto t : M_builder) delete t.first;
  cout << "there were at most " << M_table.max_size << " monomials in any list of hash table\n";
  cout << "we spend " << emplace_time << " seconds emplacing\n";
}

void F4_Reduction_Data::print_row(unsigned i, bool as_poly) {
  auto & Ai = A[i];
  for (unsigned j = 0; j < Ai.size(); ++j) {
    if (as_poly) {
      cout << " + " << Ai[j].second << " " << *M[Ai[j].first];
    } else {
      cout << Ai[j].second << " (" << Ai[j].first << "), ";
    }
  }
  cout << endl;
}

void F4_Reduction_Data::print_matrix(bool show_data) {
  if (show_data) { // print monomials
    for (auto m : M)
      cout << *m << ", ";
    cout << endl;
  }
  for (unsigned i = 0; i < num_rows; ++i) { // print entries
    cout << "A[" << i << "]: ( ";
    unsigned j = 0;
    auto & Ai = A[i];
    for (unsigned l = 0; l < Ai.size(); ++l) {
      for (/* */; j < Ai[l].first; ++j) cout << "0, ";
      cout << Ai[l].second << ", ";
      ++j;
    }
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

void F4_Reduction_Data::build_reducer(unsigned mi) {
  const auto g = R[mi];
  Polynomial_Iterator * gi = g->new_iterator();
  auto & r = R_built[mi];
  r.resize(g->length());
  Monomial u(*M[mi], g->leading_monomial(), false);
  unsigned k = 0;
  while (not gi->fellOff()) {
    const Monomial & t = gi->currMonomial();
    r[k].first = M_table.lookup_product(u, t);;
    r[k].second = gi->currCoeff().value();
    ++k;
    gi->moveRight();
  }
  delete gi;
  sort(r.begin(), r.end(), sort_by_first);
}

void expand(
    vector< pair< unsigned, COEF_TYPE > > & row,
    vector< COEF_TYPE > & B,
    vector< unsigned > & prev,
    vector< unsigned > & next
) {
  auto k = row[0].first;
  B[k] = row[0].second;
  prev[k] = B.size();
  next[k] = row[1].first;
  unsigned i = 1, n = row.size() - 1;
  while (i < n) {
    k = row[i].first;
    B[k] = row[i].second;
    next[k] = row[i+1].first;
    prev[k] = row[i-1].first;
    ++i;
  }
  if (i < row.size()) { 
    B[row[i].first] = row[i].second;
    prev[row[i].first] = k;
    next[row[i].first] = B.size();
  } else { // happens only if B.size() == 1
    next[k] = B.size();
  }
}

// always reduces the monomial at start
unsigned reduce_monomial(
    vector<COEF_TYPE> & B,
    const vector< pair< unsigned, COEF_TYPE> > & r,
    COEF_TYPE a,
    COEF_TYPE mod,
    unsigned start, unsigned & head,
    vector<unsigned> & prev, vector<unsigned> & next,
    unsigned & nonzero_entries
) {
  auto num_cols = B.size();
  unsigned i = head, j = 0;
  start = next[start];
  while (i < r[j].first) i = next[i];
  unsigned prev_i = prev[i];
  // add until we run out of monomials in reductee
  while (j < r.size()) {
    auto k = r[j].first;
    if (i < k) {
      prev_i = i; 
      i = next[i];
    } else if (i == k) {
      B[i] -= a*r[j].second; B[i] %= mod;
      if (B[i] < 0) B[i] += mod;
      if (B[i] != 0)
        prev_i = i;
      else {
        if (i == start) start = next[start];
        --nonzero_entries;
        if (i == head) head = next[i];
        if (next[i] < num_cols) prev[next[i]] = prev[i];
        if (prev[i] < num_cols) next[prev[i]] = next[i];
        prev_i = prev[i];
      }
      i = next[i]; ++j;
    } else {
      ++nonzero_entries;
      if (k < start) start = k;
      B[k] = (- r[j].second * a) % mod + mod;
      if (head > k) { // inserting new head
        prev[k] = num_cols;
        next[k] = head;
        prev[head] = k;
        i = head;
        head = k;
      } else { // head < k < i
        if (i < num_cols) {
          prev[k] = prev[i];
          if (prev[i] < num_cols) next[prev[i]] = k;
          prev[i] = k;
        } else {
          prev[k] = prev_i;
          next[prev_i] = k;
        }
        prev_i = k;
        next[k] = i;
      }
      ++j;
    }
  }
  return start;
}

void condense(
    vector< pair< unsigned, COEF_TYPE > > & row,
    unsigned head,
    const vector< COEF_TYPE > & B, const vector< unsigned > & next,
    unsigned nonzero_entries
) {
  if (row.size() != nonzero_entries) row.resize(nonzero_entries);
  unsigned i = 0;
  unsigned n = B.size();
  for (unsigned k = head; k < n; k = next[k]) {
    row[i].first = k;
    row[i].second = B[k];
    ++i;
  }
}

unsigned long long reduced_by_old = 0;

void F4_Reduction_Data::reduce_my_rows(
    const vector<int> & my_rows, vector<COEF_TYPE> & B,
    vector<unsigned> & prev, vector<unsigned> & next
) {
  NVAR_TYPE n = Rx.number_of_variables();
  const Prime_Field & F = Rx.ground_field();
  auto mod = F.modulus();
  Monomial u(n);
  B.resize(num_cols);
  prev.resize(num_cols);
  next.resize(num_cols);
  // expand, reduce, condense each row
  for (unsigned i : my_rows) {
    expand(A[i], B, prev, next);
    unsigned head = A[i][0].first;
    for (unsigned j = head; j < num_cols; /* */) {
      const Abstract_Polynomial * g = R[j];
      if (g == nullptr)
        j = next[j];
      else {
        auto a = B[j];
        red_mutex[j].lock();
        auto & r = R_built[j];
        if (r.size() == 0) build_reducer(j);
        red_mutex[j].unlock();
        auto si = strategies[i];
        si->pre_reduction_tasks(u, *g);
        j = reduce_monomial(
            B, r, a, mod, j, head, prev, next, nonzero_entries[i]
        );
        reduced_by_old += 1;
      }
    }
    condense(A[i], head, B, next, nonzero_entries[i]);
    print_lock.lock(); cout << "row " << i << " completed\n"; print_lock.unlock();
  }
}

void F4_Reduction_Data::reduce_by_old() {
  /*cout << "before reduction\n";
  for (unsigned k = 0; k < num_rows; ++k)
    check_consistency(k);*/
  if (red_lock.size() < num_cols) {
    vector<mutex> new_list(3*num_cols / 2);
    red_lock.swap(new_list);
  }
  unsigned cores = std::thread::hardware_concurrency();
  unsigned num_threads = (cores < num_rows) ? cores : num_rows;
  vector<COEF_TYPE> buffer[num_threads];
  vector<unsigned> prev[num_threads], next[num_threads];
  vector<int> * thread_rows = new vector<int>[num_threads];
  // loop through num_rows
  for (unsigned k = 0; k < num_rows; ++k)
    thread_rows[k % num_threads].push_back(k);
  thread * workers = new thread[num_threads];
  for (unsigned c = 0; c < num_threads; ++c) {
    workers[c] = thread(
        &F4_Reduction_Data::reduce_my_rows, this, std::cref(thread_rows[c]),
        std::ref(buffer[c]), std::ref(prev[c]), std::ref(next[c])
    );
  }
  for (unsigned c = 0; c < num_threads; ++c)
    workers[c].join();
  delete [] workers;
  delete [] thread_rows;
  /*cout << "after reduction\n";
  for (unsigned k = 0; k < num_rows; ++k)
    check_consistency(k);*/
}

unsigned location_of_monomial_index(
    vector< pair< unsigned, COEF_TYPE > > & row, unsigned i
) {
  if (i < row[0].first or i > row[row.size() - 1].first) return row.size(); 
  unsigned j = 0, k = row.size() / 2, l = row.size();
  while (row[k].first != i and j != k and k != l) {
    auto tmp = k;
    if (row[k].first < i) {
      k = (k + l) / 2;
      j = tmp;
    } else {
      k = (j + k) / 2;
      l = tmp;
    }
  }
  if (row[k].first == i) return k;
  if (row[l].first == i) return l;
  return row.size();
}

unsigned long long reduced_by_new = 0;

void F4_Reduction_Data::reduce_my_new_rows(
    unsigned i,
    unsigned lhead_i,
    vector< COEF_TYPE > & B,
    vector< unsigned > & prev, vector< unsigned > & next,
    const set<unsigned> & to_reduce,
    unsigned mod
) {
  auto & Ai = A[i];
  B.resize(num_cols);
  prev.resize(num_cols);
  next.resize(num_cols);
  for (auto j : to_reduce) {
    auto & Aj = A[j];
    expand(Aj, B, prev, next);
    unsigned k = location_of_monomial_index(Aj, lhead_i);
    COEF_TYPE a = Aj[k].second;
    unsigned head = Aj[0].first;
    auto start = Ai[0].first;
    if (B[start] == 0) {
      auto new_start = head;
      while (next[new_start] < start) new_start = next[new_start];
      start = new_start;
    }
    reduce_monomial(B, Ai, a, mod, start, head, prev, next, nonzero_entries[j]);
    reduced_by_new += 1;
    condense(Aj, head, B, next, nonzero_entries[j]);
  }
}

void F4_Reduction_Data::reduce_by_new(
    unsigned i, unsigned lhead_i, const set<unsigned> & unprocessed
) {
  //cout << "pre reduction:\n";
  fill(dirty.begin(), dirty.end(), false);
  //for (auto b: dirty) cout << b << ' '; cout << endl;
  //print_matrix(false);
  auto & F = Rx.ground_field();
  auto mod = F.modulus();
  unsigned cores = std::thread::hardware_concurrency();
  unsigned num_threads = (cores < num_rows) ? cores : num_rows;
  set<unsigned> * thread_rows = new set<unsigned>[num_threads];
  vector<COEF_TYPE> B[num_threads];
  vector<unsigned> prev[num_threads], next[num_threads];
  thread * workers = new thread[num_threads];
  for (unsigned j = 0; j < num_threads; ++j)
    thread_rows[j].clear();
  unsigned k = 0;
  unsigned num_to_reduce = 0;
  for (unsigned j = 0; j < num_rows; ++j) {
    if (j != i and nonzero_entries[j] != 0
        and location_of_monomial_index(A[j], lhead_i) != A[j].size()
    ) {
      dirty[j] = true;
      thread_rows[k].insert(j);
      ++k; k %= num_threads;
      ++num_to_reduce;
    }
  }
  cout << "row " << i << " reduces " << num_to_reduce << " rows\n";
  for (unsigned c = 0; c < num_threads; ++c)
    workers[c] = thread(
        &F4_Reduction_Data::reduce_my_new_rows, this,
        i, lhead_i,
        std::ref(B[c]), std::ref(prev[c]), std::ref(next[c]),
        std::cref(thread_rows[c]), mod
    );
  for (unsigned c = 0; c < num_threads; ++c)
    workers[c].join();
  delete [] workers;
  delete [] thread_rows;
  //cout << "post reduction:\n";
  //for (auto b: dirty) cout << b << ' '; cout << endl;
  //print_matrix(false);
}

Polynomial_Hashed * F4_Reduction_Data::finalize(
    unsigned i, vector< Monomial * > & final_monomials, F4_Hash & final_hash
) {
  Polynomial_Hashed * result;
  const Prime_Field & F = Rx.ground_field();
  UCOEF_TYPE mod = F.modulus();
  NVAR_TYPE n = M[0]->num_vars();
  result = new Polynomial_Hashed(Rx, final_monomials, final_hash, A[i], M, mord);
  result->set_strategy(strategies[i]);
  strategies[i] = nullptr;
  return result;
}

void F4_Reduction_Data::monomials_in_row(unsigned i, list<int> & result) const {
  unsigned processed = 0;
  auto & Ai = A[i];
  for (unsigned k = 0; k < Ai.size(); ++k)
    result.push_back(Ai[k].first);
}

void F4_Reduction_Data::simplify_identical_rows(set<unsigned> & in_use) {
  const Prime_Field & F = Rx.ground_field();
  UCOEF_TYPE mod = F.modulus();
  set<unsigned> removed;
  // identify redundants
  for (auto i : in_use) {
    auto & Ai = A[i];
    if (nonzero_entries[i] != 0) {
      for (unsigned j = i + 1; j < number_of_rows(); ++j) {
        auto & Aj = A[j];
        if (
            nonzero_entries[j] == nonzero_entries[i]
            and Ai[0].first == Aj[0].first
        ) {
          auto a = Aj[0].second * F.inverse(Ai[0].second) % mod;
          unsigned k = 1;
          for (
               /* already initialized */ ;
               k < Ai.size() and Ai[k].first == Aj[k].first
               and ((a*Ai[k].second % mod) == Aj[k].second) ;
               ++k
          ) {
            /* already handled */
          }
          if (k == Ai.size()) {
            nonzero_entries[j] = 0;
            removed.insert(j);
          }
        }
      }
    }
  }
  cout << "identified " << removed.size() << " redundant rows\n";
  // now remove
  for (auto i : removed) in_use.erase(i);
}

unsigned divisible_incompatible = 0;
unsigned old_divisibile_incompatible = 0;

void divisibility_tests(
  list<int> & allPPs,
  F4_Reduction_Data & F4,
  bool & stop
) {
  bool verbose = false;
  if (verbose) {
    cout << "checking divisibility criterion for ";
    for (auto i : allPPs) cout << *F4.M[i] << ", ";
    cout << endl;
  }
  auto it = allPPs.begin();
  auto n = F4.M.front()->num_vars();
  while ((not stop) and it != allPPs.end()) {
    bool incompatible = false, new_incompatible = false;
    auto & u = F4.M[*it];
    set<Monomial *> T;
    Monomial Tt(n, F4.M.front()->monomial_ordering());
    for (auto j : allPPs) {
      if (stop or incompatible) break;
      if (j != *it) {
        auto & t = F4.M[j];
        if (t->divisible_by(*u)) {
          incompatible = true;
          ++old_divisibile_incompatible;
          if (verbose)
            cout << "detected incompatible monomial " << *u
                 << " through divisibility by " << *t << "\n";
        }
        else if (not t->is_coprime(*u)) {
          T.insert(t);
          Tt *= *t;
          if (Tt.divisible_by_power(*u, T.size())) {
            new_incompatible = incompatible = true;
            if (verbose)
              cout << "detected incompatible monomial " << *u
                   << " through divisibility by " << Tt << "(" << T.size() << ")\n";
          }
        }
      }
    }
    while ((not incompatible) and T.size() > 1) {
      int j = 0, d = (*u)[0]*T.size() - Tt[0];
      for (unsigned k = 1; k < n; ++k)
        if ((*u)[k]*T.size() > Tt[k] + d) j = k;
      auto ti = T.begin();
      d = u->gcd_degree(**ti);
      auto remove_me = ti; ++ti;
      for (/* */; ti != T.end(); ++ti) {
        if (((**ti)[j] < (**remove_me)[j]) or u->gcd_degree(**ti) < d)
          remove_me = ti;
      }
      Tt /= **remove_me;
      T.erase(remove_me);
      if (Tt.divisible_by_power(*u, T.size())) {
        if (verbose) {
          cout << "detected incompatible monomial " << *u
               << " through divisibility by " << Tt << " (" << T.size() << ")\n";
          cout << "\t[ "; for (auto t : T) cout << *t << ", "; cout << " ]\n";
        }
        new_incompatible = incompatible = true;
      }
    }
    auto next_it(it); ++next_it;
    if (incompatible) {
      ++divisible_incompatible;
      allPPs.erase(it);
    }
    it = next_it;
  }
}

void divisibility_tests_new(
  list<int> & allPPs,
  F4_Reduction_Data & F4,
  bool & stop
) {
  bool verbose = false;
  bool very_verbose = false;
  if (very_verbose) {
    for (auto i : allPPs) cout << *F4.M[i] << ", ";
    cout << endl;
  }
  // build monomial for multi-divisibility criterion
  auto n = F4.M.front()->num_vars();
  Monomial Ttall(n, F4.M.front()->monomial_ordering());
  for (auto j : allPPs) {
    if (stop) break;
    auto & t = F4.M[j];
    Ttall *= *t;
  }
  auto it = allPPs.begin();
  while ((not stop) and it != allPPs.end()) {
    // check for simple divisibility
    bool incompatible = false, new_incompatible = false;
    auto & u = F4.M[*it];
    for (auto j : allPPs) {
      if (stop or incompatible) break;
      if (j != *it) {
        auto & t = F4.M[j];
        if (t->divisible_by(*u)) {
          incompatible = true;
          if (verbose)
            cout << "detected incompatible monomial " << *u
                 << " through old divisibility by " << *t << "\n";
        }
      }
    }
    // now check for multi-divisibility
    if (not (stop or incompatible)) {
      set<Monomial *> T;
      for (auto set_it = allPPs.begin(); set_it != allPPs.end(); ++set_it) {
        if (it != set_it) {
          T.insert(F4.M[*set_it]);
        }
      }
      Monomial Tt(Ttall); Tt /= *u;
      for (auto new_it = T.begin(); new_it != T.end(); ) {
        if (not u->is_coprime(**new_it)) {
          ++new_it;
        } else {
          Tt /= **new_it;
          auto next_it = new_it; ++next_it; T.erase(new_it); new_it = next_it;
        } 
      }
      if (Tt.divisible_by_power(*u, T.size())) {
        if (verbose) {
          cout << "detected incompatible monomial " << *u
               << " through divisibility by " << Tt << " (" << T.size() << ")\n";
          cout << "\t[ "; for (auto t : T) cout << *t << ", "; cout << " ]\n";
        }
        new_incompatible = incompatible = true;
      } else {
        if (very_verbose)
          cout << "( " << *u << " )^" << T.size() << " does not divide " << Tt << endl;
      }
      while ((not stop) and (not incompatible) and T.size() > 2) {
        int j = 0, d = (*u)[0]*T.size() - Tt[0];
        for (unsigned k = 1; k < n; ++k) {
          if ((*u)[k]*T.size() > Tt[k] + d) j = k;
        }
        auto ti = T.begin();
        auto remove_me = T.end();
        for (/* */; ti != T.end(); ++ti) {
          if (remove_me == T.end()) remove_me = ti;
          else if ((**ti)[j] < (**remove_me)[j] /*or (u->gcd(**ti).total_degree() < u->gcd(**remove_me).total_degree())*/)
            remove_me = ti;
        }
        if (very_verbose) cout << "divide " << Tt << " by " << **remove_me << " to get ";
        Tt /= **remove_me;
        T.erase(remove_me);
        if (very_verbose) cout << Tt << endl;
        if (Tt.divisible_by_power(*u, T.size())) {
          if (verbose) {
            cout << "detected incompatible monomial " << *u
                 << " through divisibility by " << Tt << " (" << T.size() << ")\n";
            cout << "\t[ "; for (auto t : T) cout << *t << ", "; cout << " ]\n";
          }
          new_incompatible = incompatible = true;
        } else {
          if (very_verbose)                    
            cout << "( " << *u << " )^" << T.size() << " does not divide " << Tt << endl;
        }
      }
    }
    auto next_it(it); ++next_it;
    if (incompatible) {
      ++divisible_incompatible;
      Ttall /= *F4.M[*it];
      allPPs.erase(it);
    }
    it = next_it;
  }
}

/**
  @ingroup GBComputation
  @author John Perry
  @date 2019
  @brief Compute the compatible leading monomials of a polynomial.
  @details This differs from the more general case in that monomials are
    indexed by @c M, rather than making copies of monomials.
    In addition, we keep a cache of the monomials for each row,
    so that we don't have to recompute the compatible monomials on each pass.
  @param my_row row to process
  @param F4 matrix structure
  @param skel existing skeleton that defines currently-compatible orderings
  @param stop signal bit whether to stop processing early (indeterminism may result)
  @param completed (out only) whether we complete our given rows
*/
void compatible_pp(
  int my_row,
  F4_Reduction_Data & F4,
  WGrevlex * curr_ord,
  const LP_Solver * skel,
  bool & stop,
  vector<bool> & completed
)
{

  F4.last_compatible_ordering[my_row] = curr_ord;

  int currentLPP_index = F4.A[my_row][0].first;

  // get the exponent vector of the current LPP, insert it
  const Monomial & currentLPP(*F4.M[currentLPP_index]);
  NVAR_TYPE n = currentLPP.num_vars();
  list<int> initial_candidates;
  initial_candidates.push_back(currentLPP_index);

  // compare other monomials with LPP
  list<int> allPPs;

  if (not stop) {

    F4.monomials_in_row(my_row, allPPs);
    //divisibility_tests(allPPs, F4, stop);
    for (const int b : allPPs) {
      if (stop) break;
      const Monomial & u(*F4.M[b]);
      auto & F4b = F4.pp_weights[b];
      auto m = F4b.size();
      //if (currentLPP_index != b and skel->makes_consistent_constraint(u, currentLPP))
      if (currentLPP_index != b) {
        auto & F4c = F4.pp_weights[currentLPP_index];
        for (unsigned i = 0; i < m; ++i)
          if (F4b[i] > F4c[i]) {
            initial_candidates.push_back(b);
            break;
          }
      }
    }

    /*cout << initial_candidates.size() << " initial candidates: ";
    for (auto b : initial_candidates) cout << *F4.M[b] << ", ";
    cout << endl;*/

    list<int> & result = F4.compatible_pps[my_row];

    if (not stop) {
  
      for (int b : initial_candidates)
      {
        if (stop) break;
        auto & F4b = F4.pp_weights[b];
        auto m = F4b.size();
        const Monomial & u(*F4.M[b]);
        bool good_constraints = true;
        for (int c : initial_candidates) {
          if (b != c) {
            auto & F4c = F4.pp_weights[c];
            bool found_one = false;
            //const Monomial & v(*F4.M[c]);
            //if (not skel->makes_consistent_constraint(u, v))
            for (unsigned i = 0; i < m; ++i)
              if (F4b[i] > F4c[i]) {
              found_one = true;
              break;
            }
            if (not found_one) {
              good_constraints = false;
              break;
            }
            //good_constraints = false;
          }
        }
        if (good_constraints)
        {
          result.push_back(b);
          //cout << "\tadded " << u << endl;
        }
    
      }

      if (not stop) completed[my_row] = true;
      if (result.size() == 1) stop = true;

    }

    //if (not stop)
    //  divisibility_tests(result, F4, stop);

  }

}

void F4_Reduction_Data::constraints_for_new_pp(
  const PP_With_Ideal &I,
  const set<int> &monomials_for_comparison,
  vector<Constraint> &result
)
{
  //static time_t all_time = 0;
  //time_t start = time(nullptr);
  // setup
  cout << "adding constraints for " << I.get_pp() << endl;
  NVAR_TYPE n = I.get_ideal().number_of_variables();
  const EXP_TYPE * a, * b;
  CONSTR_TYPE *c = new CONSTR_TYPE[n];  // space for coefficients of constraint
  a = I.get_pp().log();  // exponent vector of candidate
  // loop through exponent vectors of other 
  for (const int ti : monomials_for_comparison)
  {
    const Monomial & t{*M[ti]};
    // insert only different PPs (since I->t should also be in that set)
    if (t != I.get_pp())
    {
      //cout << "\tadding constraint " << I.get_pp() << " > " << t << endl;
      b = t.log();
      for (NVAR_TYPE i = 0; i < n; ++i) c[i] = a[i] - b[i];
      delete [] b;
      result.push_back(Constraint(n,c));
      //cout << "New constraint: " << result.back() << endl;
    }
  }
  delete [] c;
  delete [] a;
  //time_t stop = time(nullptr);
  //all_time += difftime(stop, start);
  //cout << "time spent constructing constraints: " << all_time << endl;
}

pair< bool, LP_Solver * > F4_Reduction_Data::refine(
    unsigned my_row,
    LP_Solver * currSkel,
    const list<Abstract_Polynomial *> & CurrentPolys
) {

  static time_t presolve_time = 0;
  time_t refine_start = time(nullptr);

  static time_t verify_time = 0;

  bool success = true;

  LP_Solver * result = currSkel;

  Skeleton * src_skel    = dynamic_cast<Skeleton *>    (currSkel);
  GLPK_Solver * src_GLPK = dynamic_cast<GLPK_Solver *> (currSkel);
  PPL_Solver * src_PPL   = dynamic_cast<PPL_Solver *>  (currSkel);

  list<PP_With_Ideal> & possibleIdealsBasic = potential_ideals[my_row];

  if (possibleIdealsBasic.size() != 1)
  {
    // test each combination of LPPs for consistency
    // one of them must work (current LPP, if nothing else) 
    set<int> PPunion;
    for (const int t : compatible_pps[my_row]) PPunion.insert(t);
    PP_With_Ideal & I = possibleIdealsBasic.front();
    //cout << currSkel << endl;
    LP_Solver * newSkeleton;
    if (src_skel != nullptr)
      newSkeleton = new Skeleton(*src_skel);
    else if (src_GLPK != nullptr)
      newSkeleton = new GLPK_Solver(*src_GLPK);
    else if (src_PPL != nullptr)
      newSkeleton = new PPL_Solver(*src_PPL);
    vector<Constraint> newvecs;
    //cout << "testing " << I.get_pp() << endl;
    /*cout << '\t' << *I.get_hilbert_polynomial() << endl;
    cout << '\t' << *I.get_hilbert_numerator() << endl;*/
    constraints_for_new_pp(I, PPunion, newvecs);
    if (newSkeleton->solve(newvecs))
    {
      //cout << I.get_pp() << " is consistent\n";
      time_t start_verify = time(nullptr);
      success = verify_and_modify_if_necessary(newSkeleton, CurrentPolys);
      time_t stop_verify = time(nullptr);
      verify_time += difftime(stop_verify, start_verify);
      if (not success) {
        possibleIdealsBasic.pop_front();
        delete newSkeleton;
      } else {
        result = newSkeleton;
      }
    }
    else
    {
      cout << I.get_pp() << "is inconsistent\n";
      // cout << newSkeleton;
      // this monomial is not, in fact, compatible
      //compatible_pps[my_row].erase(M_table.get_index(I.get_pp()));
      success = false;
      possibleIdealsBasic.pop_front();
      delete newSkeleton;
    }
  }
  time_t refine_stop = time(nullptr);
  presolve_time += difftime(refine_stop, refine_start);
  cout << "time spent in presolve: " << presolve_time << endl;

  return pair< bool, LP_Solver * >(success, result);

}

void F4_Reduction_Data::reassign(
    unsigned my_row,
    const Ray & w,
    Dense_Univariate_Integer_Polynomial * & current_hilbert_numerator,
    LP_Solver * currSkel,
    const PP_With_Ideal & winner,
    bool & ordering_changed
) {

  ordering_changed = false;
  pref_head[my_row] = M_table[winner.get_pp()];
  if (current_hilbert_numerator != nullptr) delete current_hilbert_numerator;
  current_hilbert_numerator
      = new Dense_Univariate_Integer_Polynomial(*(winner.get_hilbert_numerator()));
  Ray new_weight = ray_sum(currSkel->get_rays());
  new_weight.simplify_ray();
  for (
        int i = 0;
        not ordering_changed and i < (int )(new_weight.get_dimension());
        ++i
  ) {
    // w is the old ordering
    ordering_changed = ordering_changed or (new_weight[i] != w[i]);
  }
  cout << "ordering changed? " << ordering_changed << endl;

  // signal that weight cache is invalid (avoid recomputation unless necessary)
  if (ordering_changed) pp_weights[0][0] = 0;

}

void subset_analysis(
    F4_Reduction_Data & s,
    const set<unsigned> & rows,
    const list<Monomial> & T,
    const list<Abstract_Polynomial *> & G,
    const Dense_Univariate_Integer_Polynomial * h,
    const list <Critical_Pair_Dynamic *> & P,
    const Ray & w,
    Dynamic_Heuristic method
) {
  for (auto i : rows) {
    s.create_and_sort_ideals(i, T, h, G, P, w, method);
    cout << "sorted " << i << endl;
  }
}

void initial_ideal_analysis(
    F4_Reduction_Data & s,
    const set<unsigned> & rows,
    const list<Monomial> & T,
    const list<Abstract_Polynomial *> & G,
    const Dense_Univariate_Integer_Polynomial * h,
    const list <Critical_Pair_Dynamic *> & P,
    const Ray & w,
    Dynamic_Heuristic method
) {

  if (rows.size() < 5) {
    subset_analysis(s, rows, T, G, h, P, w, method);
  } else {
    unsigned cores = std::thread::hardware_concurrency();
    unsigned threads = cores < rows.size() ? cores : rows.size();
    set<unsigned> subsets[threads];
    auto row = rows.begin();
    auto n = rows.size();
    for (int i = 0; i < threads - 1; ++i) {
      for (int j = 0; j < n / threads; ++j) {
        subsets[i].insert(*row); ++row;
      }
    }
    while (row != rows.end()) {
      subsets[threads - 1].insert(*row); ++row;
    }
    vector< future<void> > tasks;
    for (unsigned i = 0; i < threads; ++i) {
      tasks.emplace_back(
        async(
            subset_analysis,
            std::ref(s), std::cref(subsets[i]), std::cref(T), std::cref(G),
            std::ref(h), std::cref(P), std::cref(w), method
        )
      );
      cout << "launched set " << i << endl;
    }
    for (unsigned i = 0; i < threads; ++i)
      tasks[i].get();
  }

}

void F4_Reduction_Data::create_and_sort_ideals(
    int my_row,
    const list<Monomial> & CurrentLPPs,
    const Dense_Univariate_Integer_Polynomial * current_hilbert_numerator,
    const list<Abstract_Polynomial *> & CurrentPolys,
    const list<Critical_Pair_Dynamic *> & crit_pairs,
    const Ray & w,
    Dynamic_Heuristic method
) {

  //static time_t create_time = 0, create_sort_time = 0;
  //time_t start = time(nullptr);

  //cout << "Have ray " << w << endl;
  vector<WT_TYPE> ord(w.get_dimension());
  for (NVAR_TYPE i = 0; i < w.get_dimension(); ++i) { ord.push_back(w[i]); }
  // list possible future ideals, sort by Hilbert Function
  list<PP_With_Ideal> & possibleIdealsBasic = potential_ideals[my_row];
  list<int> & row_compatibles = compatible_pps[my_row];
  possibleIdealsBasic.clear();
  for (const int ti : row_compatibles)
  {
    PP_With_Ideal newIdeal(*M[ti], CurrentLPPs, w, crit_pairs, current_hilbert_numerator);
    possibleIdealsBasic.push_back(newIdeal);
    /*cout << "pushed back " << *M[ti]
         << " with hilbert poly " << *newIdeal.get_hilbert_polynomial()
         << " and hilbert num " << *newIdeal.get_hilbert_numerator()
         << endl;*/
  }
  /*time_t stop = time(nullptr);
  create_time += difftime(stop, start);
  cout << "time spent in creating ideals: " << create_time << endl;*/
  //cout << "heuristic: " << method << endl;
  //static time_t presolve_time = 0;
  //time_t presolve_start = time(nullptr);
  switch(method)
  {
    case Dynamic_Heuristic::ORD_HILBERT_THEN_LEX:
      possibleIdealsBasic.sort(less_by_hilbert);
      break;
    case Dynamic_Heuristic::ORD_HILBERT_THEN_DEG:
      possibleIdealsBasic.sort(less_by_hilbert_then_degree);
      break;
    case Dynamic_Heuristic::DEG_THEN_ORD_HILBERT:
      possibleIdealsBasic.sort(less_by_degree_then_hilbert);
      break;
    case Dynamic_Heuristic::GRAD_HILB_THEN_DEG:
      possibleIdealsBasic.sort(less_by_grad_hilbert_then_degree);
      break;
    case Dynamic_Heuristic::DEG_THEN_GRAD_HILB:
      possibleIdealsBasic.sort(less_by_degree_then_grad_hilbert);
      break;
    case Dynamic_Heuristic::SMOOTHEST_DEGREES:
      possibleIdealsBasic.sort(less_by_smoothest_degrees);
      break;
    case Dynamic_Heuristic::LARGEST_MAX_COMPONENT:
      possibleIdealsBasic.sort(less_by_largest_max_component);
      break;
    case Dynamic_Heuristic::MIN_CRIT_PAIRS:
      possibleIdealsBasic.sort(less_by_num_crit_pairs);
      break;
    case Dynamic_Heuristic::BETTI_HILBERT_DEG:
      possibleIdealsBasic.sort(less_by_betti);
      break;
    case Dynamic_Heuristic::GRAD_BETTI_HILBERT_DEG:
      possibleIdealsBasic.sort(less_by_grad_betti);
      break;
    case Dynamic_Heuristic::EVIL_RANDOM:
      possibleIdealsBasic.sort(less_by_random);
      break;
    default: possibleIdealsBasic.sort(less_by_hilbert);
  }

  /*stop = time(nullptr);
  create_sort_time += difftime(stop, start);
  cout << "time spent in creating and sorting ideals: " << create_sort_time
       << endl;*/

}

void F4_Reduction_Data::print_cached_weights() {
  for (auto & cache : pp_weights) {
    cout << "( ";
    for (auto w : cache)
      cout << w << ", ";
    cout << ") ";
  }
  cout << endl;
}

void F4_Reduction_Data::recache_weights(LP_Solver * skel) {

  time_t recache_time = 0;
  time_t start_recache = time(nullptr);

  auto ray_set = skel->get_rays();
  vector<Ray> rays;
  for (auto & r : ray_set) rays.emplace_back(r);

  if (rays.size() != pp_weights[0].size()) {
    for (auto & cache : pp_weights) {
      cache.resize(rays.size());
    }
  }

  auto n = Rx.number_of_variables();

  for (unsigned i = 0; i < pp_weights.size(); ++i) {
    auto & cache = pp_weights[i];
    auto & t = *M[i];
    for (unsigned j = 0; j < rays.size(); ++j) {
      cache[j] = 0;
      for (unsigned k = 0; k < n; ++k)
        cache[j] += rays[j][k] * t[k];
    }
  }

  time_t stop_recache = time(nullptr);
  recache_time += difftime(stop_recache, start_recache);
  cout << "recaching time: " << recache_time << endl;
 
}

unsigned F4_Reduction_Data::select_dynamic_single(
    set<unsigned> & unprocessed,
    list<Monomial> & U,
    const list<Abstract_Polynomial *> & G,
    const list<Critical_Pair_Dynamic *> & P,
    WGrevlex * curr_ord,
    LP_Solver * & skel,
    Monomial * & expected_result,
    const Analysis & style
) {
  bool ordering_changed = false;
  PP_With_Ideal * newideal = nullptr;
  //cout << "rays:\n";
  //for (auto & r : skel->get_rays()) cout << '\t' << r << endl;
  Dense_Univariate_Integer_Polynomial *hn = nullptr;
  // select most advantageous unprocessed poly, reduce others
  simplify_identical_rows(unprocessed);

  unsigned winning_row = num_rows;
  LP_Solver * winning_skel = skel;

  Ray w = ray_sum(skel->get_rays());
  Dense_Univariate_Integer_Polynomial * current_hilbert_numerator = nullptr;

  set<unsigned> processed;
  bool found_single = false;

  static time_t compat_time = 0;
  time_t start_compat = time(nullptr);

  list< future<void> > waiters;
  vector<bool> completed(num_rows, false);

  if (pp_weights[0].size() == 0 or pp_weights[0][0] == 0) recache_weights(skel);

  while (nonzero_entries[*unprocessed.begin()] == 0)
    unprocessed.erase(unprocessed.begin());

  if (style == Analysis::row_sequential) {

    cout << "sequential analysis\n";

    unsigned first_row = *unprocessed.begin();
    compatible_pp(first_row, *this, curr_ord, skel, found_single, completed);
  
    create_and_sort_ideals(
        first_row, U, current_hilbert_numerator, G, P, w, heur
    );
    if (compatible_pps[first_row].size() == 1) {
      winning_row = first_row;
    } else {
      // next loop terminates b/c at least one term is truly compatible
      while (winning_row == num_rows) {
        auto refinement_result = refine(first_row, skel, G);
        if (refinement_result.first) {
          winning_row = first_row;
          winning_skel = refinement_result.second;
          break;
        }
      }
    }

  } else {

    cout << "whole analysis\n";

    for (unsigned i: unprocessed) {
      if (nonzero_entries[i] > 0) {
        /*if (
            (compatible_pps[i].size() > 0) and
            (not dirty[i]) and
            (last_compatible_ordering[i] == curr_ord)
        ) {
          completed[i] = true;
        } else {*/
          compatible_pps[i].clear();
          //compatible_pp(i, *this, skel, found_single, completed);
          waiters.push_back( std::move(
              async(
                  compatible_pp, i, std::ref(*this), curr_ord, skel,
                  std::ref(found_single), std::ref(completed)
              )
          ) );
        //}
      }
    }
  
    for (auto & fut : waiters) {
      fut.get();
    }

    time_t stop_compat = time(nullptr);
    compat_time += difftime(stop_compat, start_compat);
    cout << "time spent in compatible_pp: " << compat_time << endl;
  
    for (unsigned i: unprocessed) {
      auto size = compatible_pps[i].size();
      if (size > 0 and completed[i]) {
        //list<int> & row_compatibles = compatible_pps[i];
        //cout << row_compatibles.size() << " compatible at " << i << endl;
        if (size > 1)
          processed.insert(i);
        else {
          processed.clear();
          processed.insert(i);
          break;
        }
      }
    }
  
    cout << "analyzed " << processed.size() << " rows: ";
    for (auto i : processed) cout << i << " (" << compatible_pps[i].size() << "), "; cout << endl;
  
    static time_t sort_time = 0;
    time_t start_sort = time(nullptr);
  
    initial_ideal_analysis(
        *this, processed, U, G, current_hilbert_numerator, P, w, heur
    );
  
    time_t stop_sort = time(nullptr);
    sort_time += difftime(stop_sort, start_sort);
    cout << "time spent creating and sorting ideals: " << sort_time << endl;
  
    // initialize result
    winning_row = number_of_rows();
    winning_skel = skel;
  
    // check each row
    for (unsigned i: processed) {
  
        list<int> & row_compatibles = compatible_pps[i];
  
        if (row_compatibles.size() == 1) {
          winning_row = i; winning_skel = skel;
          break;
        } else {
          while (winning_row == number_of_rows() or
              heuristic_judges_smaller(
                potential_ideals[i].front(), potential_ideals[winning_row].front()
          )) {
            auto refinement_result = refine(i, skel, G);
            if (refinement_result.first) {
              winning_row = i;
              if (winning_skel != skel) delete winning_skel;
              winning_skel = refinement_result.second;
              break;
            }
          }
        }
  
    }

  }

  if (winning_skel != skel) {
    delete skel;
    skel = winning_skel;
    reassign(
        winning_row,
        w, current_hilbert_numerator, skel,
        potential_ideals[winning_row].front(), ordering_changed
    );
  }

  unsigned j = M_table[potential_ideals[winning_row].front().get_pp()];
  pref_head[winning_row] = j;
  expected_result = M[j];
  cout << "selected " << *M[j] << " from row " << winning_row << endl;
  static double new_reduction_time = 0;
  time_t start_time = time(nullptr);
  auto & Ai = A[winning_row];
  auto & F = Rx.ground_field();
  COEF_TYPE a = F.inverse(Ai[location_of_monomial_index(Ai, j)].second);
  auto mod = F.modulus();
  for (auto & term : Ai) {
    term.second *= a; term.second %= mod;
  }
  reduce_by_new(winning_row, j, unprocessed);
  time_t end_time = time(nullptr);
  new_reduction_time += difftime(end_time, start_time);
  cout << "spent " << new_reduction_time << " seconds in reducing by new polys\n";
  U.push_back(*M[j]);
  unprocessed.erase(winning_row);
  for (unsigned i = 0; i < number_of_rows(); ++i)
    if (unprocessed.count(i) > 0)
      if (nonzero_entries[i] == 0)
        unprocessed.erase(i);
  delete newideal;
  delete current_hilbert_numerator;
  return winning_row;
}

extern Grading_Order_Data_Allocator<EXP_TYPE> * moda;
extern Grading_Order_Data_Allocator<Monomial> * monoda;

void hash_integrity_check(
  F4_Hash & hash,
  vector< Monomial * > & mons
) {
  cout << "integrity check\n";
  for (auto mp : mons) {
    if (not hash.contains(mp)) {
      cout << "hash lost " << *mp << " from " << hash.get_index(*mp) << endl;
      hash.vomit(*mp);
    }
  }
}

list<Abstract_Polynomial *> f4_control(
    const list<Abstract_Polynomial *> &F,
    vector< Monomial * > & finalized_monomials,
    F4_Hash & finalized_hash,
    const bool static_algorithm,
    const unsigned max_refinements,
    const Analysis style
) {
  //ofstream log_file("debugging report", std::ofstream::out);
  // butcher85 goes terribly awry if we use the normal strategy
  // so set this to sugar strategy until we add an option to choose strategy
  // with inputs
  StrategyFlags this_strategy = StrategyFlags::SUGAR_STRATEGY;
  list<Monomial> T;
  Dense_Univariate_Integer_Polynomial * hn = nullptr;
  NVAR_TYPE n = F.front()->number_of_variables();
  //LP_Solver * skel = new Skeleton(n);
  LP_Solver * skel = new PPL_Solver(n);
  //LP_Solver * skel = new GLPK_Solver(n);
  time_t start_f4 = time(nullptr);
  Dynamic_Heuristic heur = Dynamic_Heuristic::ORD_HILBERT_THEN_DEG;
  //Dynamic_Heuristic heur = Dynamic_Heuristic::BETTI_HILBERT_DEG;
  cout << "computation started at " << asctime(localtime(&start_f4)) << endl;
  unsigned number_of_spolys = 0;
  double reduce_old_time = 0;
  double dynamic_time = 0;
  double creation_time = 0;
  double gm_time = 0;
  list<Abstract_Polynomial *> G;
  list<Critical_Pair_Dynamic *> P;
  // set up basis with generators
  WT_TYPE * wts = new WT_TYPE[n];
  for (NVAR_TYPE i = 0; i < n; ++i) wts[i] = 1;
  WGrevlex * curr_ord = new ORDERING_TYPE(n, wts);
  list<ORDERING_TYPE *> all_orderings_used; // so we can free them at the end
  list<Abstract_Polynomial *> Ftemp;
  DEG_TYPE operating_degree = 1000000;
  for (Abstract_Polynomial * fo : F)
  {
    Polynomial_Hashed * f = new Polynomial_Hashed(*fo, finalized_monomials, finalized_hash, curr_ord);
    f->set_strategy(new Poly_Sugar_Data(f));
    f->strategy()->at_generation_tasks();
    auto * p = new Critical_Pair_Dynamic(f, this_strategy, curr_ord);
    if (p->lcm().total_degree() < operating_degree)
      operating_degree = p->lcm().total_degree();
    P.push_back(p);
    Ftemp.push_back(f);
  }
  // main loop
  bool verbose = false;
  bool very_verbose = false;
  list<Critical_Pair_Dynamic *> Pnew;
  while (not P.empty()) {
    sort_pairs_by_strategy(P);
    report_critical_pairs(P);
    Critical_Pair_Dynamic * p = P.front();
    operating_degree = p->lcm().total_degree();
    cout << "\tdegree: " << operating_degree << endl;
    while (Pnew.empty() and not P.empty()) {
      for (auto pi = P.begin(); pi != P.end(); /* */) {
        p = *pi;
        if (p->lcm().total_degree() <= operating_degree) {
          report_front_pair(p, this_strategy);
          Pnew.push_back(p);
          auto qi = pi;
          ++qi;
          P.erase(pi);
          pi = qi;
        } else
          ++pi;
      }
      ++operating_degree;
      /*list< Critical_Pair_Dynamic * > Prestore;
      for (auto pi = P.begin(); pi != P.end(); ) { 
        p = *pi;
        if (p->lcm().total_degree() < operating_degree) {
          operating_degree = p->lcm().total_degree();
          cout << "\tno, degree: " << operating_degree << endl;
          report_front_pair(p, StrategyFlags::SUGAR_STRATEGY);
          for (auto q : Pnew) { Prestore.push_back(q); }
          Pnew.clear();
          auto qi = pi; ++qi;
          P.erase(pi);
          Pnew.push_back(p);
          pi = qi;
        } else if (p->lcm().total_degree() == operating_degree) {
          report_front_pair(p, StrategyFlags::SUGAR_STRATEGY);
          Pnew.push_back(p);
          auto qi = pi;
          ++qi;
          P.erase(pi);
          pi = qi;
        } else
          ++pi;
      }
      for (auto q : Prestore) {
        P.push_back(q);
      }*/
    }
    // make s-poly
    time_t start_time = time(nullptr);
    F4_Reduction_Data s(curr_ord, Pnew, G, heur);
    time_t end_time = time(nullptr);
    creation_time += difftime(end_time, start_time);
    number_of_spolys += Pnew.size();
    Pnew.clear();
    if (not s.is_zero()) {
      time_t start_time = time(nullptr);
      s.reduce_by_old();
      time_t end_time = time(nullptr);
      reduce_old_time += difftime(end_time, start_time);
    }
    if (s.is_zero()) {
      cout << "\tmatrix reduced to zero\n";
    } else {
      set<unsigned> unprocessed;
      for (unsigned i = 0; i < s.number_of_rows(); ++i)
        if (s.number_of_nonzero_entries(i) != 0)
          unprocessed.insert(i);
      set<unsigned> all_completed_rows;
      bool ordering_changed = false;
      if (not static_algorithm) {
        const unsigned max_comparisons
            = (max_refinements != 0) ? max_refinements : s.number_of_rows();
        cout << max_comparisons << " comparisons allowed\n";
        int comparisons = 0;
        while (unprocessed.size() != 0 and comparisons < max_comparisons) {
          time_t start_time = time(nullptr);
          //LP_Solver * old_skel = skel;
          Monomial * expected_result;
          unsigned completed_row = s.select_dynamic_single(
              unprocessed, T, G, P, curr_ord, skel, expected_result, style
          );
          time_t end_time = time(nullptr);
          dynamic_time += difftime(end_time, start_time);
          Polynomial_Hashed * r = s.finalize(completed_row, finalized_monomials, finalized_hash);
          if (s.number_of_compatibles(completed_row) > 1) ++comparisons;
          cout << comparisons << " refinements\n";
          all_completed_rows.insert(completed_row);
          Ray w(ray_sum(skel->get_rays()));
          bool ordering_changed_now = false;
          for (unsigned i = 0; (not ordering_changed_now) and (i < w.get_dimension()); ++i)
            ordering_changed_now = w[i] != curr_ord->order_weights()[i];
          if (ordering_changed_now) {
            ordering_changed = true;
            w.simplify_ray();
            WGrevlex * new_ord = new WGrevlex(w);
            cout << "new ordering: " << *new_ord << endl;
            //cout << "skeleton:\n" << *skel << endl;
            all_orderings_used.push_front(curr_ord);
            curr_ord = new_ord;
          }
          if (ordering_changed) {
            r->set_monomial_ordering(curr_ord);
            if (r->leading_coefficient().value() != 1) {
              cout << "ERROR HERE\n";
            }
          }
          cout << "\tadded " << r->leading_coefficient() << " " << r->leading_monomial() << " from row " << completed_row << endl;
          if (r->leading_monomial() != *expected_result) {
            cout << "DID NOT ORDER PROPERLY\n";
          }
          very_verbose = false;
          if (very_verbose) { cout << "\tadded "; r->println(); }
          start_time = time(nullptr);
          gm_update_dynamic(P, G, r, this_strategy, curr_ord);
          end_time = time(nullptr);
          gm_time += difftime(end_time, start_time);
        }
        if (comparisons == max_comparisons) cout << "refinements halted early\n";
      }
      // process remaining pairs statically
      list< unsigned > unfinalized;
      while (unprocessed.size() != 0) {
        unsigned winning_row = *unprocessed.begin();
        if (s.number_of_nonzero_entries(winning_row) != 0) {
          unsigned winning_lm = s.head_monomial_index(winning_row, static_algorithm);
          for (auto i: unprocessed) {
            auto nz = s.number_of_nonzero_entries(i);
            if ( (nz != 0)
                 and ( (s.head_monomial_index(i, static_algorithm) > winning_lm)
                       or (s.head_monomial_index(i, static_algorithm) == winning_lm
                           and nz < s.number_of_nonzero_entries(winning_row)) )
            ) {
              winning_row = i;
              winning_lm = s.head_monomial_index(i, static_algorithm);
            }
          }
cout << "row " << winning_row << " selects " << s.monomial(winning_lm) << endl;
          s.normalize(winning_row, winning_lm);
          s.reduce_by_new(winning_row, winning_lm, unprocessed);
          auto loc = unfinalized.begin();
          while (loc != unfinalized.end() and s.head_monomial_index(*loc, static_algorithm) > winning_lm) ++loc;
          if (loc == unfinalized.end()) unfinalized.push_back(winning_row);
          else unfinalized.insert(loc, winning_row);
        }
        unprocessed.erase(winning_row);
      }
      for (auto row : unfinalized) {
        all_completed_rows.insert(row);
        Polynomial_Hashed * r = s.finalize(row, finalized_monomials, finalized_hash);
        //log_file << " new poly ( " << operating_degree << " , " << r->leading_monomial().total_degree() << " ): " << *r << endl;// r->println(log_file);
        if (ordering_changed) {
          r->set_monomial_ordering(curr_ord);
          if (r->leading_coefficient().value() != 1) {
            cout << "ERROR HERE\n";
          }
        }
        cout << "\tadded " << r->leading_coefficient() << " " << r->leading_monomial() << " ( " << r->leading_monomial().total_degree() << " ) from row " << row << endl;
        very_verbose = false;
        if (very_verbose) { cout << "\tadded "; r->println(); }
        start_time = time(nullptr);
        gm_update_dynamic(P, G, r, this_strategy, curr_ord);
        end_time = time(nullptr);
        gm_time += difftime(end_time, start_time);
      }
      if (ordering_changed) {
        for (auto p : P)
          p->change_ordering(curr_ord);
        for (auto & t : T)
          t.set_monomial_ordering(curr_ord);
      }
    }
    /*list<Polynomial_Hashed *> B;
    cout << "basis of degree " << mindeg << endl;
    for (auto g : G) {
      static_cast<Polynomial_Hashed *>(g)->set_monomial_ordering(curr_ord);
      B.push_back(static_cast<Polynomial_Hashed *>(g));
      //cout << '\t' << *g << ',' << endl;
    }
    check_correctness(B, StrategyFlags::SUGAR_STRATEGY, mindeg);*/
    //cout << "completed degree " << mindeg << " with " << T.size() << " polynomials:\n";
    //for (auto & t : T) { cout << t << ", "; }
    //cout << endl;
  }
  delete skel;
  cout << number_of_spolys << " s-polynomials computed and reduced\n";
  // cleanup
  cout << G.size() << " polynomials before interreduction\n";
  //check_correctness(G);
  G = reduce_basis(G);
  cout << G.size() << " polynomials after interreduction\n";
  for (auto f : Ftemp) delete f;
  list<Abstract_Polynomial *> B;
  unsigned int num_mons = 0;
  for (Abstract_Polynomial * g : G) {
    g->set_monomial_ordering(curr_ord);
    B.push_back(
        new Polynomial_Hashed(*g, finalized_monomials, finalized_hash, curr_ord)
    );
    num_mons += g->length();
    delete g;
  }
  cout << num_mons << " monomials in basis (possibly counting multiple times)\n";
  cout << "final ordering: " << *curr_ord << endl;
  // first one should be curr_ord; we do not want to delete that!
  while (all_orderings_used.size() != 0) {
    ORDERING_TYPE * bye_bye_ordering = all_orderings_used.front();
    all_orderings_used.pop_front();
    delete [] bye_bye_ordering->order_weights();
    delete bye_bye_ordering;
  }
  time_t end_f4 = time(nullptr);
  double duration = difftime(end_f4, start_f4);
  cout << "computation ended at " << asctime(localtime(&end_f4)) << endl;
  cout << "parallel f4 took " << duration << " seconds\n";
  cout << "parallel f4 spent " << reduce_old_time << " seconds in reducing by old polynomials\n";
  cout << reduced_by_new + reduced_by_old << " row reductions performed\n";
  cout << dynamic_time << " seconds were spent in dynamic overhead\n";
  cout << creation_time << " seconds were spent creating the matrices\n";
  cout << gm_time << " seconds were spent analyzing critical pairs\n";
  cout << divisible_incompatible << " monomials detected as incompatible via divisibility\n";
  cout << '(' << old_divisibile_incompatible << " of these were simple divisibility)\n";
  return B;
}

#endif
