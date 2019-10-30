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

#include "f4_dynamic.hpp"
#include "algorithm_buchberger_basic.hpp"

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
    const WGrevlex * curr_ord,
    const list<Critical_Pair_Dynamic *> & P,
    const list<Abstract_Polynomial *> & B,
    Dynamic_Heuristic method
) :
    G(B), Rx(P.front()->first()->base_ring()), heur(method),
    M_table(P.front()->first()->base_ring().number_of_variables())
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
  head.clear();
  strategies.clear();
  NVAR_TYPE n = Rx.number_of_variables();
  mord = curr_ord;
  // sugar data for each row of the matrix
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
      /*cout << "for " << **mi << " selected " << p->second()->leading_monomial() << endl;
      cout << '\t' << p->lcm() << endl;
      cout << '\t' << p->first()->leading_monomial() << ", " << p->first_multiplier() << endl;
      cout << '\t' << p->second()->leading_monomial() << ", " << p->second_multiplier() << endl;*/
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
        Monomial u(*(mi->first));
        u /= (*g)->leading_monomial();
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
  // populate
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

void F4_Reduction_Data::initialize_some_rows(
    const list<Critical_Pair_Dynamic *> & P, unsigned row
) {
  const unsigned num_cols = M_builder.size();
  const COEF_TYPE F0 = 0;
  for (auto cp : P) {
    auto p = cp->first();
    //print_lock.lock();
    //cout << "initializing row " << row << " for " << cp->lcm() << " with poly " << *p << " and poly ";
    //if (cp->second() == nullptr) cout << "0\n"; else cout << *cp->second() << endl;
    //print_lock.unlock();
    const Monomial & t = cp->first_multiplier();
    auto pi = p->new_iterator();
    vector<COEF_TYPE> & Arow = A[row];
    nonzero_entries[row] = 0;
    //row_plm_cache[row].clear();
    unsigned i = M_table.lookup_product(pi->currMonomial(), t);
    Arow.resize(num_cols - i, F0);
    Arow[0] = pi->currCoeff().value();
    l_head[row] = offset[row] = i;
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
    delete cp;
    ++row;
  }
}

void F4_Reduction_Data::initialize_many(const list<Critical_Pair_Dynamic *> & P) {
  num_cols = M_builder.size();
  num_rows = P.size();
  cout << "Initializing " << num_rows << " x " << num_cols << " basic matrix\n";
  dirty.resize(num_rows, true);
  strategies.resize(num_rows);
  nonzero_entries.resize(num_rows);
  R_built.resize(num_cols);
  num_readers.assign(num_cols, 0);
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
  head.resize(P.size());
  l_head.resize(P.size());
  offset.resize(P.size());
  pref_head.resize(P.size());
  pp_weights.resize(M.size());
  compatible_pps.resize(P.size());
  potential_ideals.resize(P.size());
  unsigned cores = std::thread::hardware_concurrency();
  unsigned num_threads = (cores < num_rows) ? cores : num_rows;
  list<Critical_Pair_Dynamic *> * thread_rows
    = new list<Critical_Pair_Dynamic *>[num_threads];
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
      Monomial * t = new Monomial(pi->currMonomial());
      t->set_monomial_ordering(curr_ord);
      (*t) *= u;
      //cout << "adding " << pi->currMonomial() << " * " << u << " = " << *t << endl;
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
  for (unsigned j = offset[i] + head[i]; j < M.size(); ++j) {
    if (as_poly) {
      if (A[i][j - offset[i]] != 0) {
        cout << " + " << A[i][j - offset[i]] << " " << *M[j];
      }
    } else {
      cout << A[i][j - offset[i]] << ", ";
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

void F4_Reduction_Data::build_reducer(unsigned mi, const Monomial & u) {
  NVAR_TYPE n = Rx.number_of_variables();
  const Prime_Field & F = Rx.ground_field();
  UCOEF_TYPE mod = F.modulus();
  const auto g = R[mi];
  Polynomial_Iterator * gi = g->new_iterator();
  auto & r = R_built[mi];
  vector<COEF_TYPE> r_buf(num_cols - mi, 0);
  while (not gi->fellOff()) {
    const Monomial & t = gi->currMonomial();
    auto j = M_table.lookup_product(u, t);
    r_buf[j - mi] = gi->currCoeff().value();
    gi->moveRight();
  }
  // interreduce the buffer before finalizing the reducer
  for (unsigned ri = r_buf.size() - 1; ri > 0; --ri) { // don't try to reduce by self
    red_lock[mi + ri].lock();
    if (r_buf[ri] != 0 and R[mi + ri] != nullptr) {
      auto & s = R_built[mi + ri]; // reducer for this column
      if (s.size() == 0) {
        Monomial v(n);
        const Monomial & t = *M[mi + ri];
        for (NVAR_TYPE k = 0; k < n; ++k)
          v.set_exponent(k, t[k] - (R[mi + ri]->leading_monomial())[k]);
        build_reducer(mi + ri, v);
      }
      COEF_TYPE a(r_buf[ri]);
      for (auto si : s) {
        r_buf[si.first - mi] += mod - (a * si.second % mod);
        r_buf[si.first - mi] %= mod;
      }
    }
    red_lock[mi + ri].unlock();
  }
  // finalize the reducer
  for (unsigned ri = 0; ri < r_buf.size(); ++ri)
    if (r_buf[ri] != 0)
      r.emplace_back(ri + mi, r_buf[ri]);
  delete gi;
}

void F4_Reduction_Data::reduce_my_rows(
  const vector<int> & my_rows
) {
  NVAR_TYPE n = Rx.number_of_variables();
  const Prime_Field & F = Rx.ground_field();
  Monomial u(n);
  UCOEF_TYPE mod = F.modulus();
  unsigned new_nonzero_entries;
  #define UNREDUCED_REDUCERS true
  if (UNREDUCED_REDUCERS) {
    for (unsigned mi = 0; mi < num_cols; ++mi) {
      // is this monomial reducible?
      if (R[mi] != nullptr) {
        for (int k : my_rows) {
          vector<COEF_TYPE> & Ak = A[k];
          auto off_k = offset[k];
          // do we need to reduce this poly?
          if (off_k <= mi and Ak[mi - off_k] != 0) {
            unsigned i = mi - off_k;
            red_mutex[mi].lock();
            // get reducer for this monomial
            const Abstract_Polynomial * g = R[mi];
            vector<pair<unsigned, COEF_TYPE> > & r = R_built[mi];
            if (r.size() == 0) { // need to create reducer
              Polynomial_Iterator * gi = g->new_iterator();
              //if (k == 0) cout << "reducing " << *M[mi] << " by " << *g << endl;
              // determine multiplier
              const Monomial & t = *M[mi];
              const Monomial & v = gi->currMonomial();
              u.make_product_or_quotient(t, v, false);
              size_t j = mi;
              // loop through g's terms
              while (not gi->fellOff()) {
                const Monomial & t = gi->currMonomial();
                j = M_table.lookup_product(u, t);
                r.emplace_back(j, gi->currCoeff().value());
                gi->moveRight();
              }
              delete gi;
            }
            red_mutex[mi].unlock();
            // prepare for reduction
            new_nonzero_entries = nonzero_entries[k];
            auto sk = strategies[k];
            sk->pre_reduction_tasks(u, *g);
            // determine reduction coefficient
            COEF_TYPE a(mod - ((Ak[i]*F.inverse(r[0].second)) % mod));
            // loop through g's terms
            unsigned & hk = head[k];
            for (auto ri : r) {
              unsigned j = ri.first;
              unsigned l = j - offset[k];
              auto & Akl = Ak[l];
              bool was_zero = Akl == 0;
              Akl += a*ri.second;
              if (was_zero) {
                ++new_nonzero_entries;
                if (hk > l) hk = l;
                was_zero = false;
              }
              if (Akl & (OVERFLOW_MASK != 0)) {
                Akl %= mod;
                if (not was_zero and Akl == 0) {
                  --new_nonzero_entries;
                }
              }
              // advance
            }
            while (hk < Ak.size() and Ak[hk] == 0) {
              ++hk;
              if (hk < Ak.size() and Ak[hk] != 0) {
                Ak[hk] %= mod;
                if (Ak[hk] == 0) {
                  --new_nonzero_entries;
                }
              }
            }
            nonzero_entries[k] = new_nonzero_entries;
            if (nonzero_entries[k] == 0) hk = M.size();
            l_head[k] = hk + offset[k];
            //if (k == 0) { cout << '\t'; print_row(k); }
          }
        }
      }
    }
  } else {
    for (unsigned mi = num_cols - 1; mi < num_cols; --mi) {
      // is this monomial reducible?
      if (R[mi] != nullptr) {
        for (int k : my_rows) {
            vector<COEF_TYPE> & Ak = A[k];
            auto off_k = offset[k];
            // do we need to reduce this poly?
            if (off_k <= mi and Ak[mi - off_k] != 0) {
              unsigned i = mi - off_k;
              // get reducer for this monomial
              const Abstract_Polynomial * g = R[mi];
              // determine multiplier
              Monomial & t = *M[mi];
              for (NVAR_TYPE l = 0; l < n; ++l)
                u.set_exponent(l, t[l] - (g->leading_monomial())[l]);
              red_lock[mi].lock();
              vector<pair<unsigned, COEF_TYPE> > & r = R_built[mi];
              if (r.size() == 0) { // need to create reducer
                build_reducer(mi, u);
              }
              red_lock[mi].unlock();
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
              //if (k == 0) { cout << '\t'; print_row(k); }
            }
        }
      }
    }
  }
  for (auto i: my_rows) {
    auto & Ai = A[i];
    auto & hi = head[i];
    for (auto k = hi; k < Ai.size() and nonzero_entries[i] != 0; ++k) {
      if (Ai[k] != 0) {
        Ai[k] %= mod;
        if (Ai[k] == 0) {
          --nonzero_entries[i];
          if (nonzero_entries[i] == 0)
            hi = M.size();
          else if (k == hi) {
            while (hi < Ai.size() and Ai[hi] == 0) ++hi;
          }
        }
      }
    }
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
  unsigned cores = std::thread::hardware_concurrency() * 3 / 2;
  unsigned num_threads = (cores < num_rows) ? cores : num_rows;
  vector<int> * thread_rows = new vector<int>[num_threads];
  // loop through num_rows
  for (unsigned k = 0; k < num_rows; ++k)
    thread_rows[k % num_threads].push_back(k);
  thread * workers = new thread[num_threads];
  for (unsigned c = 0; c < num_threads; ++c) {
    workers[c] = thread(
        &F4_Reduction_Data::reduce_my_rows, this, std::cref(thread_rows[c])
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

void F4_Reduction_Data::reduce_my_new_rows(
    unsigned i,
    unsigned lhead_i,
    const set<unsigned> & to_reduce
) {
  const Prime_Field & F = Rx.ground_field();
  auto mod = F.modulus();
  const auto & Ai = A[i];
  for (auto j : to_reduce) {
    auto ci = head[i];
    auto & Aj = A[j];
    auto cj = head[i] + offset[i] - offset[j]; // pos in A[j] of A[i]'s head
    // adjust head
    if (offset[j] + head[j] > offset[i] + head[i]) {
      // first make sure row is large enough; if not, resize & copy correctly
      unsigned old_size = A[j].size();
      auto M_size = M.size();
      if (old_size < M_size - (offset[i] + head[i])) {
        auto new_size = M_size - (offset[i] + head[i]);
        A[j].resize(new_size);
        unsigned k = 1;
        for (/* */; k <= old_size; ++k)
          A[j][new_size - k] = A[j][old_size - k];
        for (/* */; k <= new_size; ++k)
          A[j][new_size - k] = 0;
        offset[j] = offset[j] + old_size - new_size;
        cj = 0;
      }
      head[j] = cj;
    }
    auto a = Aj[lhead_i - offset[j]];
    unsigned ops = 0;
    a *= mod - F.inverse(Ai[lhead_i - offset[i]]);
    auto Aj_size = Aj.size();
    while (cj < Aj_size and ops < nonzero_entries[i]) {
      if (Ai[ci] != 0) {
        bool was_zero = (Aj[cj] == 0);
        Aj[cj] += a*Ai[ci]; Aj[cj] %= mod;
        if (was_zero) {
          ++nonzero_entries[j];
        } else if (Aj[cj] == 0) {
          --nonzero_entries[j];
        }
        ++ops;
      }
      ++cj; ++ci;
    }
    unsigned & hj = head[j];
    while (hj < Aj.size() and Aj[hj] == 0) ++hj;
    l_head[j] = hj + offset[j];
    // if (j == 0) { cout << "reduced row " << j << " by new: "; print_row(j); }
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
  unsigned cores = std::thread::hardware_concurrency() * 2;
  unsigned num_threads = (cores < num_rows) ? cores : num_rows;
  set<unsigned> * thread_rows = new set<unsigned>[num_threads];
  thread * workers = new thread[num_threads];
  COEF_TYPE Ai0 = lhead_i; // abs pos of A[i]'s head
  for (unsigned j = 0; j < num_threads; ++j)
    thread_rows[j].clear();
  unsigned k = 0;
  for (unsigned j = 0; j < num_rows; ++j) {
    if (j != i and nonzero_entries[j] != 0
        and Ai0 >= offset[j]
        and (A[j][Ai0 - offset[j]] != 0)
    ) {
      dirty[j] = true;
      thread_rows[k].insert(j);
      ++k; k %= num_threads;
    }
  }
  for (unsigned c = 0; c < num_threads; ++c)
    workers[c] = thread(
        &F4_Reduction_Data::reduce_my_new_rows, this,
        i, lhead_i, std::cref(thread_rows[c])
    );
  for (unsigned c = 0; c < num_threads; ++c)
    workers[c].join();
  delete [] workers;
  delete [] thread_rows;
  //cout << "post reduction:\n";
  //for (auto b: dirty) cout << b << ' '; cout << endl;
  //print_matrix(false);
}

Constant_Polynomial * F4_Reduction_Data::finalize(unsigned i) {
  Constant_Polynomial * result;
  const Prime_Field & F = Rx.ground_field();
  UCOEF_TYPE mod = F.modulus();
  NVAR_TYPE n = M[0]->num_vars();
  vector<COEF_TYPE> & Ai = A[i];
  Monomial * M_final = static_cast<Monomial *>(
      malloc(sizeof(Monomial)*nonzero_entries[i])
  );
  Prime_Field_Element * A_final = static_cast<Prime_Field_Element *>(
      malloc(sizeof(Prime_Field_Element)*nonzero_entries[i])
  );
  COEF_TYPE scale = F.inverse(Ai[l_head[i] - offset[i]]);
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
  result = new Constant_Polynomial(
      nonzero_entries[i],
      Rx,
      M_final, A_final,
      mord
  );
  result->set_strategy(strategies[i]);
  strategies[i] = nullptr;
  for (k = 0; k < nonzero_entries[i]; ++k) M_final[k].deinitialize();
  free(M_final);
  free(A_final);
  return result;
}

void F4_Reduction_Data::monomials_in_row(unsigned i, list<int> & result) const {
  unsigned processed = 0;
  unsigned k = head[i];
  auto Ai = A[i];
  //cout << "head: " << k << "; offset: " << offset[i] << "; maximum: " << M.size() << "; nonzero: " << nonzero_entries[i] << endl;
  //cout << "inserted ";
  while (processed < nonzero_entries[i]) {
    while (Ai[k] == 0) ++k;
    result.push_back(k + offset[i]);
    //cout << k + offset[i] << ", ";
    ++k; ++processed;
  }
  //cout << endl;
}

void F4_Reduction_Data::simplify_identical_rows(set<unsigned> & in_use) {
  const Prime_Field & F = Rx.ground_field();
  UCOEF_TYPE mod = F.modulus();
  set<unsigned> removed;
  // identify redundants
  for (auto i : in_use) {
    if (nonzero_entries[i] != 0) {
      unsigned i0 = offset[i] + head[i];
      for (unsigned j = i + 1; j < number_of_rows(); ++j) {
        if (offset[j] + head[j] == i0 and nonzero_entries[j] == nonzero_entries[i]) {
          bool still_equal = true;
          auto & Ai = A[i], & Aj = A[j];
          auto a = Aj[i0 - offset[j]] * F.inverse(Ai[i0 - offset[i]]) % mod;
          for (unsigned k = i0; still_equal and k < M.size(); ++k) {
            unsigned ik = k - offset[i];
            unsigned jk = k - offset[j];
            still_equal = (a*Ai[ik] % mod ) == Aj[jk];
          }
          if (still_equal) {
            nonzero_entries[j] = 0;
            head[j] = M.size();
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
  const LP_Solver * skel,
  bool & stop,
  vector<bool> & completed
)
{

  int currentLPP_index = F4.head[my_row] + F4.offset[my_row];

  // get the exponent vector of the current LPP, insert it
  const Monomial & currentLPP(*F4.M[currentLPP_index]);
  NVAR_TYPE n = currentLPP.num_vars();
  list<int> initial_candidates;
  initial_candidates.push_back(currentLPP_index);

  // compare other monomials with LPP
  list<int> allPPs;

  if (not stop) {

    F4.monomials_in_row(my_row, allPPs);
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

    if (not stop) {
  
      list<int> & result = F4.compatible_pps[my_row];
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
    //cout << "pushed back " << t << endl;
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
    const list<Abstract_Polynomial *> G,
    const list<Critical_Pair_Dynamic *> & P,
    WGrevlex * curr_ord,
    LP_Solver * & skel,
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
    compatible_pp(first_row, *this, skel, found_single, completed);
  
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
        if (compatible_pps[i].size() > 0 and (not dirty[i])) completed[i] = true;
        else {
        compatible_pps[i].clear();
        //compatible_pp(i, *this, skel, found_single, completed);
        waiters.push_back( std::move(
            async(
                compatible_pp, i, std::ref(*this), skel,
                std::ref(found_single), std::ref(completed)
            )
        ) );
      }
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

  // find current lm and use for reduction
  vector<COEF_TYPE> & Ai = A[winning_row];
  //unsigned j = M_table[winning_lm];
  unsigned j = M_table[potential_ideals[winning_row].front().get_pp()];
  l_head[winning_row] = j;
  static double new_reduction_time = 0;
  time_t start_time = time(nullptr);
  reduce_by_new(winning_row, j, unprocessed);
  time_t end_time = time(nullptr);
  new_reduction_time += difftime(end_time, start_time);
  cout << "spent " << new_reduction_time << " seconds in reducing by new polys\n";
  U.push_back(*M[l_head[winning_row]]);
  cout << "selected " << *M[l_head[winning_row]] << " from row " << winning_row << endl;
  unprocessed.erase(winning_row);
  for (unsigned i = 0; i < number_of_rows(); ++i)
    if (unprocessed.count(i) > 0)
      if (nonzero_entries[i] == 0)
        unprocessed.erase(i);
  delete newideal;
  return winning_row;
}

extern Grading_Order_Data_Allocator<EXP_TYPE> * moda;
extern Grading_Order_Data_Allocator<Monomial> * monoda;

list<Constant_Polynomial *> f4_control(
    const list<Abstract_Polynomial *> &F,
    const bool static_algorithm,
    const unsigned max_refinements,
    const Analysis style
) {
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
    Constant_Polynomial * f = new Constant_Polynomial(*fo);
    f->set_strategy(new Poly_Sugar_Data(f));
    f->strategy()->at_generation_tasks();
    auto * p = new Critical_Pair_Dynamic(f, StrategyFlags::SUGAR_STRATEGY, curr_ord);
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
    cout << "\tdegree: " << operating_degree << endl;
    while (Pnew.empty() and not P.empty()) {
      for (auto pi = P.begin(); pi != P.end(); /* */) { 
        p = *pi;
        if (p->lcm().total_degree() <= operating_degree) {
          report_front_pair(p, StrategyFlags::SUGAR_STRATEGY);
          Pnew.push_back(p);
          auto qi = pi;
          ++qi;
          P.erase(pi);
          pi = qi;
        } else
          ++pi;
      }
      ++operating_degree;
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
          unsigned completed_row = s.select_dynamic_single(
              unprocessed, T, G, P, curr_ord, skel, style
          );
          time_t end_time = time(nullptr);
          dynamic_time += difftime(end_time, start_time);
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
        }
        if (comparisons == max_comparisons) cout << "refinements halted early\n";
      }
      // process remaining pairs statically
      while (unprocessed.size() != 0) {
        unsigned winning_row = *unprocessed.begin();
        if (s.number_of_nonzero_entries(winning_row) != 0) {
          unsigned winning_lm = s.head_monomial_index(winning_row);
          for (auto i: unprocessed) {
            auto nz = s.number_of_nonzero_entries(i);
            if ( (nz != 0)
                 and ( (s.head_monomial_index(i) > winning_lm)
                       or (s.head_monomial_index(i) == winning_lm
                           and nz < s.number_of_nonzero_entries(winning_row)) )
            ) {
              winning_row = i;
              winning_lm = s.head_monomial_index(i);
            }
          }
          s.reduce_by_new(winning_row, winning_lm, unprocessed);
          all_completed_rows.insert(winning_row);
        }
        unprocessed.erase(winning_row);
      }
      for (auto completed_row : all_completed_rows) {
        if (s.number_of_nonzero_entries(completed_row) != 0) {
          Constant_Polynomial * r = s.finalize(completed_row);
          if (ordering_changed) {
            r->set_monomial_ordering(curr_ord);
            for (auto p : P)
              p->change_ordering(curr_ord);
            for (auto & t : T)
              t.set_monomial_ordering(curr_ord);
          }
          cout << "\tadded " << r->leading_monomial() << " from row " << completed_row << endl;
          //r->printlncout();
          //cout << "SANITY CHECK: ordering changed? " << ordering_changed << "; " << curr_ord << endl;
          //for (auto t : T) cout << t << " "; cout << endl;
          very_verbose = false;
          if (very_verbose) { cout << "\tadded "; r->println(); }
          start_time = time(nullptr);
          gm_update_dynamic(P, G, r, StrategyFlags::SUGAR_STRATEGY, curr_ord);
          end_time = time(nullptr);
          gm_time += difftime(end_time, start_time);
          //for (auto g : G) cout << g->leading_monomial() << " "; cout << endl;
          //for (auto p : P) cout << p->how_ordered() << ' '; cout << endl;
          cout << "continuing\n";
        }
      }
    }
    /*list<Constant_Polynomial *> B;
    cout << "basis of degree " << mindeg << endl;
    for (auto g : G) {
      static_cast<Constant_Polynomial *>(g)->set_monomial_ordering(curr_ord);
      B.push_back(static_cast<Constant_Polynomial *>(g));
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
  list<Constant_Polynomial *> B;
  unsigned int num_mons = 0;
  for (Abstract_Polynomial * g : G) {
    g->set_monomial_ordering(curr_ord);
    B.push_back(new Constant_Polynomial(*g));
    num_mons += g->length();
    delete g;
  }
  cout << num_mons << " monomials in basis (possibly counting multiple times)\n";
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
  cout << dynamic_time << " seconds were spent in dynamic overhead\n";
  cout << creation_time << " seconds were spent creating the matrices\n";
  cout << gm_time << " seconds were spent analyzing critical pairs\n";
  return B;
}

#endif