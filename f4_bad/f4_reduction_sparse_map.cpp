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
  A.clear();
  I.clear();
  R_build.clear();
  M_build.clear();
  strategies.clear();
  auto mi = M_build.begin();
  auto ri = R_build.begin();
  unsigned cores = std::thread::hardware_concurrency() * 2;
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
  num_cols = M_build.size();
  const Prime_Field & F = P.front()->first()->ground_field();
  const COEF_TYPE F0 = 0;
  num_rows = P.size();
  strategies.resize(num_rows);
  nonzero_entries.resize(num_rows);
  M.clear();
  for (auto mi = M_build.begin(); mi != M_build.end(); ++mi)
    M.push_back(*mi);
  for (unsigned i = 0; i < strategies.size(); ++i)
    strategies[i] = nullptr;
  unsigned row = 0;
  A.resize(P.size());
  I.resize(P.size());
  for (auto cp : P) {
    auto p = cp->first();
    auto t = cp->first_multiplier();
    auto pi = p->new_iterator();
    unordered_map<unsigned, COEF_TYPE> & Arow = A[row];
    unsigned j = 0;
    for (unsigned i = 0; not pi->fellOff(); ++i) {
      while (not M[i]->like_multiple(pi->currMonomial(), t))
        ++i;
      Arow.emplace(i, pi->currCoeff().value());
      pi->moveRight();
      ++nonzero_entries[row];
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
    cout << "A[" << i << "]: { ";
    for (auto a : A[i]) {
      cout << '(' << a.first << ',' << a.second << " ) ";
    }
    cout << "}\n";
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
unsigned num_copies = 0;
unsigned copy_length = 0;
std::mutex my_mutex;

void F4_Reduction_Data::reduce_my_rows(
  const set<unsigned> & my_rows
) {
  vector<unsigned> indices(num_cols);
  NVAR_TYPE n = Rx.number_of_variables();
  const Prime_Field & F = Rx.ground_field();
  EXP_TYPE * u = new EXP_TYPE[n]; // exponents of multiplier
  UCOEF_TYPE mod = F.modulus();
  set<unsigned> unconsidered_keys;
  unsigned new_nonzero_entries;
  for (unsigned k : my_rows) {
    // loop through A[k]'s terms
    unsigned buffer_last = nonzero_entries[k];
    unordered_map<unsigned, COEF_TYPE> & Ak = A[k];
    for (auto a : Ak) unconsidered_keys.insert(a.first);
    while (unconsidered_keys.size() != 0) {
      // do we need to reduce this term, and can we?
      auto mk = min_element(unconsidered_keys.begin(), unconsidered_keys.end());
      unsigned mi = *mk;
      unconsidered_keys.erase(mk);
      COEF_TYPE a = Ak[mi];
      if (R[mi] != nullptr) {
        // get reducer for this monomial
        const Abstract_Polynomial * g = R[mi];
        Polynomial_Iterator * gi = g->new_iterator();
        // set up append
        new_nonzero_entries = nonzero_entries[k];
        // determine multiplier
        for (NVAR_TYPE l = 0; l < n; ++l)
          u[l] = (*M[mi])[l] - (gi->currMonomial())[l];
        // determine reduction coefficient
        COEF_TYPE b(mod - ((a*F.inverse(gi->currCoeff().value())) % mod));
        // loop through g's terms
        unsigned j = mi;
        while (not gi->fellOff()) {
          const Monomial & t = gi->currMonomial();
          while (not M[j]->like_multiple(u, t)) ++j;
          auto aj = Ak.find(j);
          if (aj == Ak.end()) {
            Ak[j] = b*gi->currCoeff().value() % mod;
            gi->moveRight();
            ++new_nonzero_entries;
            unconsidered_keys.insert(j);
          } else { // j == I[k][l]
            aj->second += b*gi->currCoeff().value(); aj->second %= mod;
            // if the current coefficient is zero, do not advance
            if (aj->second == 0) {
              Ak.erase(aj);
              --new_nonzero_entries;
              unconsidered_keys.erase(j);
            }
           // advance
            gi->moveRight();
          }
        }
        if (strategies[k] != nullptr)
          strategies[k]->pre_reduction_tasks(u, *g);
        delete gi;
        nonzero_entries[k] = new_nonzero_entries;
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
  delete [] workers;
  delete [] thread_rows;
}

bool check(unordered_map<unsigned, COEF_TYPE> M) {
  bool result = true;
  for (auto mi = M.begin(); result and mi != M.end(); ++mi)
    result = (mi->second != 0);
  return result;
}

void F4_Reduction_Data::reduce_by_new() {
  const Prime_Field & F = Rx.ground_field();
  UCOEF_TYPE mod = F.modulus();
  // find leading indices
  li.resize(num_rows);
  for (unsigned i = 0; i < num_rows; ++i) {
    if (nonzero_entries[i] == 0)
      li[i] = num_cols;
    else {
      unsigned ti = num_cols;
      for (auto a : A[i])
        if (a.first < ti)
          ti = a.first;
      li[i] = ti;
    }
  }
  // use new rows to reduce other new rows
  for (unsigned i = 0; i < num_rows; ++i) {
    // use a nonzero row
    if (nonzero_entries[i] != 0) {
      COEF_TYPE a = F.inverse(A[i][li[i]]);
      // look for a row that can be reduced by this one
      for (unsigned j = 0; j < num_rows; ++j) {
        if (i != j) { // obviously not the same row
          auto aji = A[j].find(li[i]);
          if (aji != A[j].end()) { // found one
            auto safety = A[j];
            COEF_TYPE b(mod - ((a*aji->second) % mod));
            bool found = false;
            // add A[i]'s terms to A[j]
            for (auto ai : A[i]) {
              // can we find ai in A[j]?
              found = true;
              auto aj = A[j].find(ai.first);
              if (aj == A[j].end()) { // no; create new term
                A[j][ai.first] = (b*ai.second) % mod;
                ++nonzero_entries[j];
              }
              else { // yes; modify coefficient & check for 0
                aj->second += (b*ai.second); aj->second %= mod;
                if (aj->second == 0) {
                  --nonzero_entries[j];
                  // find new leading index
                  if (li[j] == aj->first) {
                    unsigned ti = num_cols;
                    for (auto ak : A[j])
                      if (ak.first < ti and ak.first != aj->first)
                        ti = ak.first;
                    li[j] = ti;
                  }
                  A[j].erase(aj);
                }
              }
            }
          }
        }
      }
    }
  }
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
      Monomial * M_final = (Monomial *)malloc(sizeof(Monomial)*nonzero_entries[i]);
      Prime_Field_Element * A_final
          = (Prime_Field_Element *)malloc(sizeof(Prime_Field_Element)*nonzero_entries[i]);
      unsigned ti = A[i][li[i]];
      COEF_TYPE scale = F.inverse(ti);
      unsigned j = 0;
      for (auto ai : A[i]) {
        A_final[j].assign(((ai.second)*scale) % mod, &F);
        M_final[j].common_initialization();
        M_final[j].initialize_exponents(n);
        M_final[j].set_monomial_ordering(mord);
        M_final[j] = *M[ai.first];
        ++j;
      }
      result.push_back(new Constant_Polynomial(
          nonzero_entries[i],
          Rx,
          M_final, A_final,
          mord
      ));
      result.back()->sort_by_order();
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
      time_t start_time = time(nullptr);
      s.reduce_by_old();
      time_t end_time = time(nullptr);
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
  cout <<  num_copies << " copies took " << total_time << " seconds\n";
  cout << "avg len " << double(copy_length) / double(num_copies) << endl;
  return B;
}

#endif