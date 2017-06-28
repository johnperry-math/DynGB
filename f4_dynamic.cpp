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

#include "lp_solver.hpp"
using LP_Solvers::LP_Solver;
#include "skeleton.hpp"
using LP_Solvers::Skeleton;
#include "glpk_solver.hpp"
using LP_Solvers::GLPK_Solver;
#include "ppl_solver.hpp"
using LP_Solvers::PPL_Solver;

using Dynamic_Engine::PP_With_Ideal;
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
    const list<Abstract_Polynomial *> & B
) : G(B), Rx(P.front()->first()->base_ring()) {
  num_cols = 0;
  A.clear();
  head.clear();
  R_build.clear();
  M_build.clear();
  strategies.clear();
  auto mi = M_build.begin();
  auto ri = R_build.begin();
  NVAR_TYPE n = Rx.number_of_variables();
  mord = curr_ord;
  for (auto p : P) {
    Poly_Sugar_Data * new_sugar = new Poly_Sugar_Data(p->first());
    strategies.push_back(new_sugar);
    if (strategies.back() != nullptr) {
      strategies.back()->at_generation_tasks(p->first_multiplier());
      if (p->second() != nullptr)
        strategies.back()->pre_reduction_tasks(p->second_multiplier().log(), *(p->second()));
    }
    add_monomials(curr_ord, mi, ri, p->first(), p->first_multiplier(), true);
    if (p->second() != nullptr) {
      *ri = const_cast<Abstract_Polynomial *>(p->second());
      ++ri; ++mi;
      add_monomials(curr_ord, mi, ri, p->second(), p->second_multiplier());
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
        add_monomials(curr_ord, mi, ri, *g, u);
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
    const list<Critical_Pair_Dynamic *> & P, unsigned row
) {
  const unsigned num_cols = M_build.size();
  const Prime_Field & F = P.front()->first()->ground_field();
  const COEF_TYPE F0 = 0;
  for (auto cp : P) {
    auto p = cp->first();
    const Monomial & t = cp->first_multiplier();
    auto pi = p->new_iterator();
    vector<COEF_TYPE> & Arow = A[row];
    nonzero_entries[row] = 0;
    unsigned i = 0;
    while (not M[i]->like_multiple(pi->currMonomial(), t)) ++i;
    Arow.resize(num_cols - i, F0);
    Arow[0] = pi->currCoeff().value();
    offset[row] = i;
    head[row] = 0;
    nonzero_entries[row] = 1;
    pi->moveRight();
    for (/* */; not pi->fellOff(); ++i) {
      while (not M[i]->like_multiple(pi->currMonomial(), t)) ++i;
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
  num_cols = M_build.size();
  num_rows = P.size();
  strategies.resize(num_rows);
  nonzero_entries.resize(num_rows);
  R_built.resize(num_cols);
  num_readers.assign(num_cols, 0);
  M.clear();
  for (auto mi = M_build.begin(); mi != M_build.end(); ++mi)
    M.push_back(*mi);
  unsigned row = 0;
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
  //reduce_my_rows(thread_rows[0], buffer[0], indices[0]);
  delete [] workers;
  delete [] thread_rows;
  //
}

void F4_Reduction_Data::add_monomials(
    const WGrevlex * curr_ord,
    list<Monomial *>::iterator & t1,
    list<Abstract_Polynomial *>::iterator & r1,
    const Abstract_Polynomial *g,
    const Monomial & u,
    bool new_row
) {
  const_cast<Abstract_Polynomial *>(g)->set_monomial_ordering(curr_ord);
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
  //delete [] mord->order_weights();
  //delete mord;
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
  EXP_TYPE * u = new EXP_TYPE[n]; // exponents of multiplier
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
          for (NVAR_TYPE l = 0; l < n; ++l)
            u[l] = (*M[mi])[l] - (gi->currMonomial())[l];
          red_mutex.lock();
          vector<pair<unsigned, COEF_TYPE> > & r = R_built[mi];
          if (r.size() == 0) { // need to create reducer
            unsigned j = mi;
            // loop through g's terms
            while (not gi->fellOff()) {
              const Monomial & t = gi->currMonomial();
              while (not (M[j]->like_multiple(u, t))) ++j;
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
  delete [] u;
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
    unsigned lhead_i,
    const Prime_Field & F,
    const set<unsigned> & to_reduce
) {
  auto mod = F.modulus();
  const auto & Ai = A[i];
  for (auto j : to_reduce) {
    auto ci = lhead_i - offset[i];
    auto & Aj = A[j];
    auto cj = lhead_i - offset[j]; // pos in A[j] of A[i]'s head
    auto a = Aj[cj];
    unsigned ops = 0;
    a *= mod - F.inverse(Ai[ci]);
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
  }
}

void F4_Reduction_Data::reduce_by_new(unsigned i, unsigned lhead_i) {
  auto F = Rx.ground_field();
  unsigned cores = std::thread::hardware_concurrency() * 2;
  unsigned num_threads = (cores < num_rows) ? cores : num_rows;
  set<unsigned> * thread_rows = new set<unsigned>[num_threads];
  thread * workers = new thread[num_threads];
  //for (unsigned i = 0; i < num_rows; ++i) {
  //  if (nonzero_entries[i] > 0) {
      COEF_TYPE Ai0 = lhead_i; // abs pos of A[i]'s head
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
            i, lhead_i, std::cref(F), std::cref(thread_rows[c])
        );
      for (unsigned c = 0; c < num_threads; ++c)
        workers[c].join();
  //  }
  //}
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
      Monomial * M_final = (Monomial *)malloc(sizeof(Monomial)*nonzero_entries[i]);
      Prime_Field_Element * A_final
          = (Prime_Field_Element *)malloc(sizeof(Prime_Field_Element)*nonzero_entries[i]);
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

WGrevlex * F4_Reduction_Data::select_dynamic(
    list<Monomial> & U,
    const list<Abstract_Polynomial *> G,
    const list<Critical_Pair_Dynamic *> & P,
    WGrevlex * curr_ord,
    LP_Solver * & skel,
    Dynamic_Heuristic heur
) {
  bool ordering_changed = false;
  WGrevlex * result = curr_ord;
  PP_With_Ideal * newideal = nullptr;
  LP_Solver * winning_skel = nullptr;
  set<unsigned> unprocessed;
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
            potential_pps[i], *M[offset[i] + head[i]], U, &H[i], G, P, new_lp,
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
            bool change_winner;
            switch(heur) {
              case Dynamic_Heuristic::ORD_HILBERT_THEN_DEG:
                change_winner = less_by_hilbert_then_degree(*tempideal, *newideal); break;
              case Dynamic_Heuristic::ORD_HILBERT_THEN_LEX:
                change_winner = less_by_hilbert(*tempideal, *newideal); break;
              case Dynamic_Heuristic::DEG_THEN_ORD_HILBERT:
                change_winner = less_by_degree_then_hilbert(*tempideal, *newideal); break;
              case Dynamic_Heuristic::GRAD_HILB_THEN_DEG:
                change_winner = less_by_grad_hilbert_then_degree(*tempideal, *newideal); break;
              case Dynamic_Heuristic::SMOOTHEST_DEGREES:
                change_winner = less_by_smoothest_degrees(*tempideal, *newideal); break;
              case Dynamic_Heuristic::LARGEST_MAX_COMPONENT:
                change_winner = less_by_largest_max_component(*tempideal, *newideal); break;
              case Dynamic_Heuristic::MIN_CRIT_PAIRS:
                change_winner = less_by_num_crit_pairs(*tempideal, *newideal); break;
              case Dynamic_Heuristic::BETTI_HILBERT_DEG:
                change_winner = less_by_betti(*tempideal, *newideal); break;
              case Dynamic_Heuristic::GRAD_BETTI_HILBERT_DEG:
                change_winner = less_by_grad_betti(*tempideal, *newideal); break;
              default:
                change_winner = less_by_hilbert(*tempideal, *newideal); break;
            }
            if (not change_winner) {
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
    // find current lm and use for reduction
    bool found_mon = false;
    unsigned j = offset[winning_row] + head[winning_row];
    for (unsigned k = offset[winning_row]; not found_mon; ++k) {
      if (*M[k] == winning_lm) {
        j = k;
        found_mon = true;
      }
    }
    reduce_by_new(winning_row, j);
    unprocessed.erase(winning_row);
    for (unsigned i = 0; i < number_of_rows(); ++i)
      if (unprocessed.count(i) > 0 and nonzero_entries[i] == 0)
        unprocessed.erase(i);
  }
  if (ordering_changed) {
    skel = winning_skel;
    result = new WGrevlex(ray_sum(skel->get_rays()));
    cout << "new ordering " << result << endl;
  }
  delete newideal;
  return result;
}

list<Constant_Polynomial *> f4_control(const list<Abstract_Polynomial *> &F) {
  list<Monomial> T;
  Dense_Univariate_Integer_Polynomial * hn = nullptr;
  NVAR_TYPE n = F.front()->number_of_variables();
  LP_Solver * skel = new Skeleton(n);
  time_t start_f4 = time(nullptr);
  Dynamic_Heuristic heur = Dynamic_Heuristic::ORD_HILBERT_THEN_DEG;
  cout << "computation started at " << asctime(localtime(&start_f4)) << endl;
  unsigned number_of_spolys = 0;
  double total_time = 0;
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
    F4_Reduction_Data s(curr_ord, Pnew, G);
    time_t end_time = time(nullptr);
    total_time += difftime(end_time, start_time);
    number_of_spolys += Pnew.size();
    Pnew.clear();
    if (not s.is_zero())
      s.reduce_by_old();
    if (s.is_zero()) {
      cout << "\tmatrix reduced to zero\n";
      // delete s;
    } else {
      WGrevlex * new_ord = s.select_dynamic(T, G, P, curr_ord, skel, heur);
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
    }
  }
  delete skel;
  cout << number_of_spolys << " s-polynomials computed and reduced\n";
  // cleanup
  cout << G.size() << " polynomials before interreduction\n";
  //check_correctness(G, strategy);
  G = reduce_basis(G);
  cout << G.size() << " polynomials after interreduction\n";
  for (auto f : Ftemp) delete f;
  list<Constant_Polynomial *> B;
  for (Abstract_Polynomial * g : G) {
    g->set_monomial_ordering(curr_ord);
    B.push_back(new Constant_Polynomial(*g));
    delete g;
  }
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
  cout << "parallel f4 spent " << total_time << " seconds in timed section\n";
  check_correctness(B);
  return B;
}

#endif