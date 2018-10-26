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
#include "algorithm_buchberger_dynamic.hpp"

#include <thread>
using std::thread;
#include <mutex>
using std::mutex;
#include <bitset>
using std::bitset;

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
using Dynamic_Engine::compatible_pp;
using Dynamic_Engine::verify_and_modify_if_necessary;

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
      if (p->second() != nullptr)
        strategies.back()->pre_reduction_tasks(p->second_multiplier().log(), *(p->second()));
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
    while (not found and g != G.end()) {
      if (mi->first->divisible_by((*g)->leading_monomial())) {
        found = true;
        //cout << "for " << **mi << " selected " << (*g)->leading_monomial() << endl;
        //*ri = *g;
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

void F4_Reduction_Data::initialize_some_rows(
    const list<Critical_Pair_Dynamic *> & P, unsigned row
) {
  const unsigned num_cols = M_builder.size();
  //const Prime_Field & F = P.front()->first()->ground_field();
  const COEF_TYPE F0 = 0;
  for (auto cp : P) {
    auto p = cp->first();
    const Monomial & t = cp->first_multiplier();
    auto pi = p->new_iterator();
    vector<COEF_TYPE> & Arow = A[row];
    nonzero_entries[row] = 0;
    //unsigned i = 0;
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
    delete cp;
    ++row;
  }
}

void F4_Reduction_Data::initialize_many(const list<Critical_Pair_Dynamic *> & P) {
  num_cols = M_builder.size();
  num_rows = P.size();
  cout << "Initializing " << num_rows << " x " << num_cols << " basic matrix\n";
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
  //unsigned row = 0;
  A.resize(P.size());
  head.resize(P.size());
  l_head.resize(P.size());
  offset.resize(P.size());
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
  //
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
  //static unsigned long monomials_processed = 0;
  Polynomial_Iterator * pi = g->new_iterator();
  if (not new_row) pi->moveRight();
  //monomials_processed += g->length();
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

void F4_Reduction_Data::print_row(unsigned i) {
  for (unsigned j = offset[i] + head[i]; j < M.size(); ++j) {
    if (A[i][j - offset[i]] != 0) {
      cout << " + " << A[i][j - offset[i]] << " " << *M[j];
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

double reduction_time = 0;
//mutex red_mutex;

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
          if (k != -1) {
            vector<COEF_TYPE> & Ak = A[k];
            auto off_k = offset[k];
            // do we need to reduce this poly?
            if (off_k <= mi and Ak[mi - off_k] != 0) {
              unsigned i = mi - off_k;
              // get reducer for this monomial
              const Abstract_Polynomial * g = R[mi];
              Polynomial_Iterator * gi = g->new_iterator();
              //if (k == 0) cout << "reducing " << *M[mi] << " by " << *g << endl;
              // determine multiplier
              const Monomial & t = *M[mi];
              const Monomial & v = gi->currMonomial();
              u.make_product_or_quotient(t, v, false);
              red_mutex[mi].lock();
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
                //begin modifications 16 Oct
                Akl += a*ri.second; // Akl %= mod;
                if (was_zero) {
                  ++new_nonzero_entries;
                  if (hk > l) hk = l;
                  was_zero = false;
                }
                if (Akl & OVERFLOW_MASK != 0) {
                  Akl %= mod;
                  if (not was_zero and Akl == 0)
                    --new_nonzero_entries;
                }
                // advance
              }
              while (hk < Ak.size() and Ak[hk] == 0) {
                ++hk;
                if (hk < Ak.size() and Ak[hk] != 0) {
                  Ak[hk] %= mod;
                  if (Ak[hk] == 0) --new_nonzero_entries;
                }
              }
              // end modifications 16 Oct (but see below)
              nonzero_entries[k] = new_nonzero_entries;
              if (nonzero_entries[k] == 0) hk = M.size();
              //if (k == 0) { cout << '\t'; print_row(k); }
            }
          }
        }
      }
    }
  } else {
    for (unsigned mi = num_cols - 1; mi < num_cols; --mi) {
      // is this monomial reducible?
      if (R[mi] != nullptr) {
        for (int k : my_rows) {
          if (k != -1) {
            vector<COEF_TYPE> & Ak = A[k];
            auto off_k = offset[k];
            // do we need to reduce this poly?
            if (off_k <= mi and Ak[mi - off_k] != 0) {
              unsigned i = mi - off_k;
              // get reducer for this monomial
              const Abstract_Polynomial * g = R[mi];
              //if (k == 0) cout << "reducing " << *M[mi] << " by " << *g << endl;
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
  }
  // begin modifications 16 Oct
  //mutex debugging;
  for (auto i: my_rows) {
    auto & Ai = A[i];
    /*debugging.lock();
    cout << "first check of " << i << ": " << head[i] << ", "
         << nonzero_entries[i] << endl;
    check_consistency(i);
    debugging.unlock();*/
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
    /*debugging.lock();
    cout << "second check of " << i << ": " << head[i] << ", "
         << nonzero_entries[i] << endl;
    check_consistency(i);
    debugging.unlock();*/
  }
  // end modifications 16 Oct
}

void F4_Reduction_Data::reduce_by_old() {
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
  //for (unsigned k = 0; k < num_threads; ++k) thread_rows[k].push_back(-1);
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
      if (old_size < M.size() - (offset[i] + head[i]) + 1) {
        A[j].resize(M.size());
        unsigned k = 1;
        for (/* */; k <= old_size; ++k)
          A[j][M.size() - k] = A[j][old_size - k];
        for (/* */; k <= M.size(); ++k)
          A[j][M.size() - k] = 0;
        offset[j] = 0; cj = offset[i] + head[i];
      }
      head[j] = cj;
    }
    auto a = Aj[lhead_i - offset[j]];
    unsigned ops = 0;
    a *= mod - F.inverse(Ai[lhead_i - offset[i]]);
    while (cj < Aj.size() and ops < nonzero_entries[i]) {
      if (Ai[ci] != 0) {
        bool was_zero = (Aj[cj] == 0);
        Aj[cj] += a*Ai[ci]; Aj[cj] %= mod;
        if (was_zero)
          ++nonzero_entries[j];
        else if (Aj[cj] == 0)
          --nonzero_entries[j];
        ++ops;
      }
      ++cj; ++ci;
    }
    unsigned & hj = head[j];
    while (hj < Aj.size() and Aj[hj] == 0) ++hj;
    // if (j == 0) { cout << "reduced row " << j << " by new: "; print_row(j); }
  }
}

void F4_Reduction_Data::reduce_by_new(
    unsigned i, unsigned lhead_i, const set<unsigned> & unprocessed
) {
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
        and Ai0 > offset[j]
        and (A[j][Ai0 - offset[j]] != 0)
    ) {
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
}

vector<Constant_Polynomial *> F4_Reduction_Data::finalize() {
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
      Prime_Field_Element * A_final = static_cast<Prime_Field_Element *>(
          malloc(sizeof(Prime_Field_Element)*nonzero_entries[i])
      );
      COEF_TYPE scale = F.inverse(Ai[head[i]]);
      unsigned k = 0;
      for (unsigned j = head[i]; k < nonzero_entries[i]; ++j) {
        if (Ai[j] != 0) {
          A_final[k].assign(Ai[j]*scale % mod, &F);
          M_final[k].common_initialization();
          M_final[k].initialize_exponents(n);
          //M_final[k].set_monomial_ordering(mord);
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
      for (k = 0; k < nonzero_entries[i]; ++k)
        M_final[k].deinitialize();
      free(M_final);
      free(A_final);
    }
  }
  return result;
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
  COEF_TYPE scale = F.inverse(Ai[l_head[i]]);
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

void F4_Reduction_Data::monomials_in_row(unsigned i, set<Monomial> & result) const {
  unsigned processed = 0;
  unsigned k = head[i];
  auto Ai = A[i];
  while (processed < nonzero_entries[i]) {
    while (Ai[k] == 0) ++k;
    result.insert(*M[k + offset[i]]);
    ++k; ++processed;
  }
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

bool F4_Reduction_Data::verify_and_modify_processed_rows(LP_Solver * skel) {
  bool consistent = true; // innocent until proven guilty
  Ray w = ray_sum(skel->get_rays()); // our tentative ordering
  NVAR_TYPE n = w.get_dimension();
  CONSTR_TYPE *coefficients = new CONSTR_TYPE [n]; // used for coefficients for new constraints
  // cout << "Have ray " << w << endl;
  // loop through all polynomials; verify leading power product is unchanged
  for (auto i : dynamic_processed) {
    // create a ray for the current LPP's exponents
    const Monomial & t = *M[l_head[i]];
    DEG_TYPE td = t.weighted_degree(w.weights());
    // loop through the polynomial's remaining monomials
    //for (poly titer = (*piter)->next; consistent and titer != nullptr; titer = titer->next)
    unsigned k = head[i];
    unsigned mons_considered = 1; // accounting for l_head[i]
    while (mons_considered < nonzero_entries[i]) {
      // don't compare with LPP; that would be Bad (TM)
      if (k + offset[i] != l_head[i])
      {
        ++mons_considered;
        Monomial & u = *M[k + offset[i]];
        // compare weights between a and b; if this fails,
        // recompute the skeleton with a new constraint
        if (td <= u.weighted_degree(w.weights())) {
          if (coefficients == nullptr) // ensure we have space
            coefficients = new CONSTR_TYPE[n];
          for (NVAR_TYPE i = 0; i < n; ++i)
            coefficients[i] = t[i] - u[i];
          Constraint new_constraint(n, coefficients);
          LP_Solver * newskel;
          if (dynamic_cast<Skeleton *>(skel) != nullptr)
            newskel = new Skeleton(*static_cast<Skeleton *>(skel));
          else if (dynamic_cast<GLPK_Solver *>(skel) != nullptr) {
            newskel = new GLPK_Solver(*static_cast<GLPK_Solver *>(skel));
          }
          else if (dynamic_cast<PPL_Solver *>(skel) != nullptr)
            newskel = new PPL_Solver(*static_cast<PPL_Solver *>(skel));
          else
            newskel = nullptr; // this should never happen, of course
          consistent = newskel->solve(new_constraint);
          //cout << "Have ray " << w << endl;
          // if we're consistent, we need to recompute the ordering
          if (consistent) { // if consistent
            w = ray_sum(newskel->get_rays());
            //cout << "Have ray " << w << endl;
            *skel = *newskel;
          } else delete newskel;
        } // if LPP changed
      } // if PP != LPP
      ++k;
    } // loop through PPs
  } // loop through processed
  // cleanup
  delete [] coefficients;
  // finally done
  return consistent;
}

bool F4_Reduction_Data::row_would_change_ordering(unsigned i) {
  auto possibleIdeals = I[i];
  auto currSkel = row_skel[i];
  Ray w = ray_sum(currSkel->get_rays());
  Skeleton * src_skel    = dynamic_cast<Skeleton *>    (currSkel);
  GLPK_Solver * src_GLPK = dynamic_cast<GLPK_Solver *> (currSkel);
  PPL_Solver * src_PPL   = dynamic_cast<PPL_Solver *>  (currSkel);
  // main loop: look for a leading monomial that is compatible with skeleton
  // and other chosen monomials
  bool searching = true;
  const PP_With_Ideal * winner = & (*(possibleIdeals.begin()));
  while (searching) {
    LP_Solver * newSkeleton;
    if (src_skel != nullptr)
      newSkeleton = new Skeleton(*src_skel);
    else if (src_GLPK != nullptr)
      newSkeleton = new GLPK_Solver(*src_GLPK);
    else if (src_PPL != nullptr)
      newSkeleton = new PPL_Solver(*src_PPL);
    vector<Constraint> newvecs;
    set<Monomial> monomials_to_compare;
    for (auto & ideal : I[i]) monomials_to_compare.insert(ideal.get_pp());
    constraints_for_new_pp(*winner, monomials_to_compare, newvecs);
    if (newSkeleton->solve(newvecs)) {
      //cout << "consistent\n";
      if (
          verify_and_modify_if_necessary(newSkeleton, G)
          and verify_and_modify_processed_rows(newSkeleton)
      ) {
        searching = false;
        delete row_skel[i];
        if (src_skel != nullptr)
          row_skel[i] = newSkeleton;
        else if (src_GLPK != nullptr)
          row_skel[i] = newSkeleton;
        else if (src_PPL != nullptr)
          row_skel[i] = newSkeleton;
      }
    } else {
      // this monomial is not, in fact, compatible
      delete newSkeleton;
      possibleIdeals.erase(possibleIdeals.begin());
      winner = & (*possibleIdeals.begin());
    }
  }
  // check whether winner would change ordering
  Ray new_weight = ray_sum(currSkel->get_rays());
  //cout << "Have ray " << w << endl;
  new_weight.simplify_ray();
  bool ordering_changed = false;
  for (
        int i = 0;
        not ordering_changed and i < (int )(new_weight.get_dimension());
        ++i
  ) {
    // w is the old ordering
    ordering_changed = ordering_changed or (new_weight[i] != w[i]);
  }
  cout << "row " << i << "would change ordering? " << ordering_changed << endl;
  return ordering_changed;
}

WGrevlex * F4_Reduction_Data::reduce_and_select_order(
  const list<Monomial> & T, const list<Critical_Pair_Dynamic *> & P,
  LP_Solver * skel
) {
  WGrevlex * result = nullptr;
  //bool ordering_changed = false;
  // set up the logical heads, the ideals, ...
  I.resize(number_of_rows());
  row_skel.resize(number_of_rows(), nullptr);
  Ray w(mord->number_of_weights(), mord->order_weights());
  Dense_Univariate_Integer_Polynomial * hn = hilbert_numerator_bigatti(T);
  for (unsigned i = 0; i < number_of_rows(); ++i) {
    if (dynamic_cast<Skeleton *>(skel) != nullptr)
      row_skel[i] = new Skeleton(*static_cast<Skeleton *>(skel));
    else if (dynamic_cast<GLPK_Solver *>(skel) != nullptr)
      row_skel[i] = new GLPK_Solver(*static_cast<GLPK_Solver *>(skel));
    else if (dynamic_cast<PPL_Solver *>(skel) != nullptr)
      row_skel[i] = new PPL_Solver(*static_cast<PPL_Solver *>(skel));
    l_head[i] = head[i];
    if (nonzero_entries[i] != 0) {
      set<Monomial> all_pps, potential_pps;
      list<Monomial> boundary_pps;
      monomials_in_row(i, all_pps);
      compatible_pp(*M[offset[i] + head[i]], all_pps, potential_pps, boundary_pps, row_skel[i]);
      for (auto t : potential_pps) {
        I[i].emplace_front(t, T, w, P, hn);
      }
      I[i].sort(heuristic_judges_smaller);
      dynamic_unprocessed.insert(i);
    }
  }
  list<Monomial> U(T);
  // main loop
  vector<bool> row_changes_ordering(number_of_rows());
  while (dynamic_unprocessed.size() != 0) {
    for (unsigned i = 0; i < number_of_rows(); ++i)
      row_changes_ordering[i] = false;
    simplify_identical_rows(dynamic_unprocessed);
    set<Monomial> considered_monomials;
    unsigned winning_row = *dynamic_unprocessed.begin();
    for (auto i : dynamic_unprocessed) {
      if (
          considered_monomials.count(I[i].front().get_pp()) == 0
          and row_would_change_ordering(i)
      ) {
        row_changes_ordering[i] = true;
        if (heuristic_judges_smaller(I[i].front(), I[winning_row].front())) {
          winning_row = i;
        }
      }
      considered_monomials.insert(I[i].front().get_pp());
    }
    dynamic_processed.insert(winning_row);
    dynamic_unprocessed.erase(winning_row);
    if (row_changes_ordering[winning_row]) {
      //ordering_changed = true;
      // update skeleton
      skel->copy(row_skel[winning_row]);
      // update ordering
      if (result != nullptr) delete result;
      Ray r(ray_sum(skel->get_rays()));
      r.simplify_ray();
      result = new WGrevlex(r);
      // update row skeletons and PP_With_Ideal's
      U.push_back(I[winning_row].front().get_pp());
      for (auto i : dynamic_unprocessed) {
        delete row_skel[i];
        if (dynamic_cast<Skeleton *>(skel) != nullptr)
          row_skel[i] = new Skeleton(*static_cast<Skeleton *>(skel));
        else if (dynamic_cast<GLPK_Solver *>(skel) != nullptr)
          row_skel[i] = new GLPK_Solver(*static_cast<GLPK_Solver *>(skel));
        else if (dynamic_cast<PPL_Solver *>(skel) != nullptr)
          row_skel[i] = new PPL_Solver(*static_cast<PPL_Solver *>(skel));
        list<PP_With_Ideal> J;
        for (auto PI : I[i])
          J.emplace_back(PI.get_pp(), U, r, P, hn);
        I[i].clear();
        I[i].insert(I[i].begin(), J.begin(), J.end());
      }
      // reconsider compatible monomials?
    }
    set<unsigned> singletons;
    for (auto i : dynamic_unprocessed)
      if (I[i].size() == 1) singletons.insert(i);
    for (auto i : singletons) {
      dynamic_unprocessed.erase(i);
      dynamic_processed.insert(i);
    }
  }
  return result;
}

WGrevlex * F4_Reduction_Data::select_dynamic(
    list<Monomial> & U,
    const list<Abstract_Polynomial *> G,
    const list<Critical_Pair_Dynamic *> & P,
    WGrevlex * curr_ord,
    LP_Solver * & skel
) {
  bool ordering_changed = false;
  WGrevlex * result = curr_ord;
  PP_With_Ideal * newideal = nullptr;
  LP_Solver * winning_skel = nullptr;
  set<unsigned> processed, unprocessed;
  vector<set<Monomial> > potential_pps(number_of_rows());
  for (unsigned i = 0; i < number_of_rows(); ++i) {
    if (nonzero_entries[i] != 0) {
      set<Monomial> all_pps;
      list<Monomial> boundary_pps;
      monomials_in_row(i, all_pps);
      compatible_pp(*M[offset[i] + head[i]], all_pps, potential_pps[i], boundary_pps, skel);
      unprocessed.insert(i);
    }
  }
  Dense_Univariate_Integer_Polynomial *hn = nullptr;
  // select most advantageous unprocessed poly, reduce others
  while (unprocessed.size() != 0) {
    simplify_identical_rows(unprocessed);
    cout << unprocessed.size() << " rows to consider\n";
    unsigned winning_row = *(unprocessed.begin());
    Monomial winning_lm { *M[offset[winning_row] + head[winning_row]] };
    for (unsigned i : unprocessed) {
      // copy skeleton
      if (nonzero_entries[i] != 0 and potential_pps[i].size() > 1) {
        LP_Solver * new_lp;
        list<Monomial> T {U};
        if (dynamic_cast<Skeleton *>(skel) != nullptr)
          new_lp = new Skeleton(*dynamic_cast<Skeleton *>(skel));
        else if (dynamic_cast<GLPK_Solver *>(skel) != nullptr)
          new_lp = new GLPK_Solver(*dynamic_cast<GLPK_Solver *>(skel));
        else if (dynamic_cast<PPL_Solver *>(skel) != nullptr)
          new_lp = new PPL_Solver(*dynamic_cast<PPL_Solver *>(skel));
        // see if we can obtain a better ordering from this polynomial
        vector<Dense_Univariate_Integer_Polynomial*> H(number_of_rows());
        for (auto h : H)
          h = hn == nullptr ?
            nullptr : new Dense_Univariate_Integer_Polynomial(*hn);
        bool new_ordering = false;
        select_monomial(
            potential_pps[i], *potential_pps[i].begin(), U, &H[i], G, P, new_lp,
            new_ordering, heur
        );
        // if we have a new ordering, compare the result to current preference
        // we use PP_With_Ideal to hold relevant dynamic information
        if (not new_ordering) {
          for (auto h : H) delete h;
          delete new_lp;
        } else {
          // first save the index of the leading monomial (need for reduction)
          Monomial row_lm { U.back() };
          U.pop_back();
          // if tempideal beats newideal, reassign newideal
          Ray r(mord->number_of_weights(), mord->order_weights());
          if (newideal == nullptr) { // always true in this case
            ordering_changed = true;
            winning_row = i;
            winning_skel = new_lp;
            delete skel;
            skel = winning_skel;
            newideal = new PP_With_Ideal(row_lm, T, r, P, nullptr);
            newideal->set_hilbert_numerator(H[i]);
          } else {
            PP_With_Ideal * tempideal = new PP_With_Ideal(row_lm, T, r, P, nullptr);
            tempideal->set_hilbert_numerator(H[i]);
          if (not heuristic_judges_smaller(*tempideal, *newideal)) {
              delete new_lp;
              delete tempideal;
            } else { // winner has changed; reduce matrix by it
              ordering_changed = true;
              delete newideal;
              winning_row = i;
              winning_lm = row_lm;
              newideal = tempideal;
              winning_skel = new_lp;
              delete skel;
              skel = winning_skel;
            } // change_winner?
          } // newideal == nullptr?
        } // ordering changed?
      }
    } // loop through rows
    if (ordering_changed) {
      skel = winning_skel;
      result = new WGrevlex(ray_sum(skel->get_rays()));
      cout << "new ordering " << result << endl;
    }
    // find current lm and use for reduction
    unsigned j = offset[winning_row] + head[winning_row];
    DEG_TYPE d = M[j]->weighted_degree(result->order_weights());
    vector<COEF_TYPE> & Ai = A[winning_row];
    for (unsigned k = j; k < M.size(); ++k) {
      if (Ai[k - offset[winning_row]] != 0) {
        DEG_TYPE e = M[k]->weighted_degree(result->order_weights());
        if (e > d) {
          d = e;
          j = k;
        }
      }
    }
    reduce_by_new(winning_row, j, unprocessed);
    unprocessed.erase(winning_row);
    processed.insert(winning_row);
    for (unsigned i = 0; i < number_of_rows(); ++i)
      if (unprocessed.count(i) > 0) {
        if (nonzero_entries[i] == 0)
          unprocessed.erase(i);
        else {
          if (ordering_changed and potential_pps[i].size() > 1) {
            set<Monomial> new_pps;
            list<Monomial> boundary_pps;
            compatible_pp(
              *(potential_pps[i].begin()), potential_pps[i], new_pps, boundary_pps, skel
            );
            if (potential_pps[i].size() != new_pps.size())
              cout << "changed from " << potential_pps[i].size() << " to " << new_pps.size() << endl;
            potential_pps[i] = new_pps;
          }
        }
      }
  }
  delete newideal;
  return result;
}

unsigned F4_Reduction_Data::select_dynamic_single(
    set<unsigned> & unprocessed,
    list<Monomial> & U,
    const list<Abstract_Polynomial *> G,
    const list<Critical_Pair_Dynamic *> & P,
    WGrevlex * curr_ord,
    LP_Solver * & skel
) {
  bool ordering_changed = false;
  PP_With_Ideal * newideal = nullptr;
  LP_Solver * winning_skel = nullptr;
  Dense_Univariate_Integer_Polynomial *hn = nullptr;
  // select most advantageous unprocessed poly, reduce others
  simplify_identical_rows(unprocessed);
  unsigned winning_row = *(unprocessed.begin());
  Monomial winning_lm { *M[offset[winning_row] + head[winning_row]] };
  set<unsigned> remove;
  for (unsigned i : unprocessed) {
    if (nonzero_entries[i] == 0)
      remove.insert(i);
    else {
      set<Monomial> all_pps, potential_pps;
      list<Monomial> boundary_pps;
      monomials_in_row(i, all_pps);
      compatible_pp(*M[offset[i] + head[i]], all_pps, potential_pps, boundary_pps, skel);
      cout << "compatible monomials for row " << i << ": ";
      for (auto t : potential_pps) cout << t << ", "; cout << endl;
      if (potential_pps.size() == 1) {
        // we process rows with only one potential pp first
        // this helps avoid wrong paths
        winning_row = i;
        //winning_lm = *M[offset[winning_row] + head[winning_row]];
        winning_lm = *potential_pps.begin();
        break;
      } else {
        LP_Solver * new_lp;
        list<Monomial> T {U};
        if (dynamic_cast<Skeleton *>(skel) != nullptr) {
          new_lp = new Skeleton(*dynamic_cast<Skeleton *>(skel));
          //cout << "allocated skeleton " << new_lp << endl;
        }
        else if (dynamic_cast<GLPK_Solver *>(skel) != nullptr)
          new_lp = new GLPK_Solver(*dynamic_cast<GLPK_Solver *>(skel));
        else if (dynamic_cast<PPL_Solver *>(skel) != nullptr)
          new_lp = new PPL_Solver(*dynamic_cast<PPL_Solver *>(skel));
        // see if we can obtain a better ordering from this polynomial
        bool new_ordering = false;
        Dense_Univariate_Integer_Polynomial * hn = nullptr;
        select_monomial(
            potential_pps, *potential_pps.begin(), U, & hn, G, P, new_lp,
            new_ordering, heur
        );
        // first save the index of the leading monomial (need for reduction)
        Monomial row_lm { U.back() };
        U.pop_back();
        // if tempideal beats newideal, reassign newideal
        Ray r(ray_sum(new_lp->get_rays()));
        r.simplify_ray();
        if (newideal == nullptr) { // always true in this case
          ordering_changed = new_ordering;
          winning_row = i;
          winning_skel = new_lp;
          //cout << "saving skeleton " << new_lp << endl;
          winning_lm = row_lm;
          newideal = new PP_With_Ideal(row_lm, T, r, P, nullptr);
          newideal->set_hilbert_numerator(hn);
        } else {
          PP_With_Ideal * tempideal = new PP_With_Ideal(row_lm, T, r, P, nullptr);
          tempideal->set_hilbert_numerator(hn);
          if (not heuristic_judges_smaller(*tempideal, *newideal)) {
            delete new_lp;
            //cout << "deleting new skeleton " << new_lp << endl;
            delete tempideal;
          } else { // winner has changed; reduce matrix by it
            ordering_changed = ordering_changed or new_ordering;
            delete newideal;
            winning_row = i;
            winning_lm = row_lm;
            newideal = tempideal;
            delete winning_skel;
            //cout << "deleting old skeleton " << winning_skel << endl;
            winning_skel = new_lp;
          } // change_winner?
        } // newideal == nullptr?
      } // row has more than one potential pp
    } // row has nonzero entries
  } // loop through rows
  if (not ordering_changed and winning_skel == nullptr)
    delete winning_skel;
  else {
    //cout << "deleting " << skel << endl;
    delete skel;
    skel = winning_skel;
  }
  // find current lm and use for reduction
  unsigned j = offset[winning_row] + head[winning_row];
  vector<COEF_TYPE> & Ai = A[winning_row];
  bool searching = true;
  for (unsigned k = j; searching; ++k) {
    if (Ai[k - offset[winning_row]] != 0) {
      searching = (*(M[k])) != winning_lm;
      if (not searching) {
        j = k;
        l_head[winning_row] = k - offset[winning_row];
      }
    }
  }
  static double new_reduction_time = 0;
  time_t start_time = time(nullptr);
  reduce_by_new(winning_row, j, unprocessed);
  time_t end_time = time(nullptr);
  new_reduction_time += difftime(end_time, start_time);
  cout << "spent " << new_reduction_time << " seconds in reducing by new polys\n";
  U.push_back(*M[l_head[winning_row] + offset[winning_row]]);
  cout << "selected " << *M[l_head[winning_row] + offset[winning_row]] << endl;
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

list<Constant_Polynomial *> f4_control(const list<Abstract_Polynomial *> &F) {
  list<Monomial> T;
  Dense_Univariate_Integer_Polynomial * hn = nullptr;
  NVAR_TYPE n = F.front()->number_of_variables();
  //LP_Solver * skel = new Skeleton(n);
  LP_Solver * skel = new PPL_Solver(n);
  //LP_Solver * skel = new GLPK_Solver(n);
  time_t start_f4 = time(nullptr);
  Dynamic_Heuristic heur = Dynamic_Heuristic::ORD_HILBERT_THEN_DEG;
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
  for (Abstract_Polynomial * fo : F)
  {
    Constant_Polynomial * f = new Constant_Polynomial(*fo);
    f->set_strategy(new Poly_Sugar_Data(f));
    f->strategy()->at_generation_tasks();
    P.push_back(new Critical_Pair_Dynamic(f, StrategyFlags::SUGAR_STRATEGY, curr_ord));
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
      }
      else if (p->lcm().total_degree() == mindeg) {
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
      // delete s;
    } else {
      /* OLD WAY
      WGrevlex * new_ord = s.select_dynamic(T, G, P, curr_ord, skel);
      bool ordering_changed = new_ord != curr_ord;
      vector<Constant_Polynomial *> Rvec = s.finalize();
      if (ordering_changed) {
        cout << "new ordering: " << *new_ord << endl;
        for (auto p : P)
          p->change_ordering(new_ord);
        for (auto & t : T)
          t.set_monomial_ordering(new_ord);
        all_orderings_used.push_front(curr_ord);
        curr_ord = new_ord;
      }
      for (auto r : Rvec) {
        if (ordering_changed)
          r->set_monomial_ordering(new_ord);
        T.push_back(r->leading_monomial());
        cout << "\tadded " << r->leading_monomial() << endl;
        very_verbose = false;
        if (very_verbose) { cout << "\tadded "; r->println(); }
        gm_update_dynamic(P, G, r, StrategyFlags::SUGAR_STRATEGY, curr_ord);
      }
      OLD WAY */
      set<unsigned> unprocessed;
      for (unsigned i = 0; i < s.number_of_rows(); ++i)
        if (s.number_of_nonzero_entries(i) != 0)
          unprocessed.insert(i);
      set<unsigned> all_completed_rows;
      bool ordering_changed = false;
      while (unprocessed.size() != 0) {
        time_t start_time = time(nullptr);
        LP_Solver * old_skel = skel;
        unsigned completed_row = s.select_dynamic_single(
            unprocessed, T, G, P, curr_ord, skel
        );
        time_t end_time = time(nullptr);
        dynamic_time += difftime(end_time, start_time);
        all_completed_rows.insert(completed_row);
        ordering_changed |= skel != old_skel;
        if (ordering_changed) {
          Ray w(ray_sum(skel->get_rays()));
          w.simplify_ray();
          WGrevlex * new_ord = new WGrevlex(w);
          cout << "new ordering: " << *new_ord << endl;
          //cout << "skeleton:\n" << *skel << endl;
          all_orderings_used.push_front(curr_ord);
          curr_ord = new_ord;
          /*s.set_ordering(curr_ord);
          for (auto p : P)
            p->change_ordering(curr_ord);
          for (auto & t : T)
            t.set_monomial_ordering(curr_ord);*/
        }
      }
      s.set_ordering(curr_ord);
      for (auto p : P)
        p->change_ordering(curr_ord);
      for (auto & t : T)
        t.set_monomial_ordering(curr_ord);
      for (auto completed_row : all_completed_rows) {
        Constant_Polynomial * r = s.finalize(completed_row);
        if (ordering_changed)
          r->set_monomial_ordering(curr_ord);
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
  //cout << "deleting " << skel << endl;
  cout << number_of_spolys << " s-polynomials computed and reduced\n";
  // cleanup
  cout << G.size() << " polynomials before interreduction\n";
  //check_correctness(G, strategy);
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
