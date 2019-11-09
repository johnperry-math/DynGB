#ifndef __ALGORITHM_BUCHBERGER_DYNAMIC_CPP_
#define __ALGORITHM_BUCHBERGER_DYNAMIC_CPP_

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

#include "system_constants.hpp"

#include "algorithm_buchberger_basic.hpp"
#include "algorithm_buchberger_dynamic.hpp"

#include "skeleton.hpp"
#include "dynamic_engine.hpp"
using Dynamic_Engine::compatible_pp; using Dynamic_Engine::hilbert_cmp;
#include "particular_orderings.hpp"
#include "reduction_support.hpp"

// instantiation of templated functions

template void sort_pairs_by_strategy<Critical_Pair_Dynamic>(
    list<Critical_Pair_Dynamic *> &
);

template void report_critical_pairs<Critical_Pair_Dynamic>(
    const list<Critical_Pair_Dynamic *>, bool
);

template bool no_triplet<Critical_Pair_Dynamic>(
    const Critical_Pair_Dynamic *, const list<Critical_Pair_Dynamic *>
);

void reduce_over_basis_dynamic(
    Mutable_Polynomial **sp,
    const list<Abstract_Polynomial *>G
) {
  Abstract_Polynomial * g; // used to loop through G
  Mutable_Polynomial * s = *sp; // s-poly
  Mutable_Polynomial * r = s->zero_polynomial(); // remainder / residue
  bool verbose = false;
  bool very_verbose = false;
  // continue reducing until s is zero
  while (!s->is_zero()) {
    if (verbose) cout << "reducing " << s->leading_monomial() << "? "; 
    if (very_verbose) s->println();
    if ((g = find_reducer<list<Abstract_Polynomial *>>(s, G))) {
      if (verbose) cout << "\tyes! by " << g->leading_monomial() << endl;
      if (very_verbose) { cout << "\tyes! by "; g->println(); }
      // DYNAMIC
      if (g->monomial_ordering() != s->monomial_ordering())
        g->set_monomial_ordering(s->monomial_ordering());
      Monomial t(s->leading_monomial());
      t /= g->leading_monomial();
      Prime_Field_Element a = s->leading_coefficient();
      Prime_Field_Element b = g->leading_coefficient();
      a /= b;
      // s = s - atg
      if (s->strategy() != nullptr)
        s->strategy()->pre_reduction_tasks(t, *g);
      s->add_polynomial_multiple(a, t, *g, true);
      if (very_verbose) {
        cout << "\tresults in ";
        if (s->is_zero())
          cout << 0 << endl;
        else
          s->println();
      }
      if (verbose and not s->is_zero())
        cout << "\tresults in " << s->leading_monomial() << endl;
    } else {
      if (verbose or very_verbose) cout << "no!\n";
      // no reducer; move leading term to residue and continue
      Abstract_Polynomial * t = s->detach_head();
      // t should already be smaller than what's already in r, so no danger
      r->add_last(t->leading_coefficient(), t->leading_monomial());
      delete t;
    }
  }
  // move the strategy from s to r
  r->set_strategy(s->strategy());
  s->set_strategy(nullptr);
  delete s;
  *sp = r;
}

void gm_update_dynamic(
    list<Critical_Pair_Dynamic *> & P,
    list<Abstract_Polynomial *> & G,
    Abstract_Polynomial * r,
    StrategyFlags strategy,
    ORDERING_TYPE * ordering
) {
  //cout << "----------------------\n";
  list<Critical_Pair_Dynamic *> C;
  // critical pairs with new polynomial
  for (Abstract_Polynomial * g : G)
    C.push_back(new Critical_Pair_Dynamic(g, r, strategy, ordering));
  // apply Buchberger's lcm criterion to new pairs
  list<Critical_Pair_Dynamic *> D;
  while (C.size() != 0) {
    Critical_Pair_Dynamic * p = C.front();
    C.pop_front();
    if ((p->first()->leading_monomial().is_coprime(
            p->second()->leading_monomial()))
        or (no_triplet(p, C) and no_triplet(p, D))
        )
      D.push_back(p);
    else {
      //cout << "triplet prunes " << *p << endl;
      delete p;
    }
  }
  // apply Buchberger's gcd criterion
  list<Critical_Pair_Dynamic *> E;
  while (D.size() != 0) {
    Critical_Pair_Dynamic * p = D.front();
    D.pop_front();
    if (!(p->first()->leading_monomial().is_coprime(
            p->second()->leading_monomial())))
      E.push_back(p);
    else {
      //cout << "gcd prunes " << *p << endl;
      delete p;
    }
  }
  // apply Buchberger's lcm criterion to old pairs
  list<Critical_Pair_Dynamic *> Q;
  while (P.size() != 0) {
    Critical_Pair_Dynamic * p = P.front();
    P.pop_front();
    if (!(r->leading_monomial() | p->lcm())
          or lcm_alike(p->first()->leading_monomial(), r->leading_monomial(), p)
          or lcm_alike(p->second()->leading_monomial(), r->leading_monomial(), p)
       )
      Q.push_back(p);
    else {
      //cout << "triplet prunes " << *p << endl;
      delete p;
    }
  }
  P = Q;
  // add new pairs to old pairs
  for (Critical_Pair_Dynamic * e : E)
    P.push_back(e);
  /*cout << "All pairs:\n";
  for (auto pi = P.begin(); pi != P.end(); ++pi)
    cout << '\t' << **pi << endl;
  cout << "----------------------\n";*/
  // add new poly to basis
  G.push_back(r);
}

bool none_others_divide(const list<Monomial> & U, const Monomial & t) {
  bool result = true;
  for (const Monomial & u : U) {
    if (not u.is_like(t) and u | t) {
      result = false;
      break;
    }
  }
  return result;
}

list<list<Monomial>> recursive_select(list<list<Monomial>> & Ms) {
  list<list<Monomial>> result;
  if (Ms.size() == 0) {
    list<Monomial> M;
    result.push_back(M);
  } else {
    auto M = Ms.front();
    Ms.pop_front();
    auto lower = recursive_select(Ms);
    for (auto t: M) {
      if (none_others_divide(M, t)) {
        for (auto S : lower) {
          auto tmp(S);
          tmp.push_front(t);
          result.push_back(tmp);
        }
      }
    }
  }
  return result;
}

void initial_analysis(
    const list<Abstract_Polynomial *> & F,
    Monomial_Ordering ** mord,
    LP_Solver * solver
) {
  NVAR_TYPE n = (*F.begin())->leading_monomial().num_vars();
  WT_TYPE * weights = new WT_TYPE[n];
  // arrange possible monomial chains
  list<list<Monomial>> mons;
  list<list<Monomial>> saved_mons;
  for (auto f : F) {
    set<Monomial> T;
    set<Monomial> U;
    list<Monomial> V;
    Polynomial_Iterator * fi = f->new_iterator();
    while (not fi->fellOff()) {
      U.insert((fi->currMonomial()));
      fi->moveRight();
    }
    delete fi;
    compatible_pp(f->leading_monomial(), U, T, V, solver);
    list<Monomial> W;
    for (auto t : T) W.push_back(t);
    mons.push_back(W);
    saved_mons.push_back(W);
  }
  list<list<Monomial>> mon_chains = recursive_select(mons);
  // now check each chain for solvability & desirability
  Dense_Univariate_Integer_Polynomial * win_hn = nullptr;
  Dense_Univariate_Rational_Polynomial * win_hp = nullptr;
  vector<Constraint> all_constraints;
  for (auto M : mon_chains) {
    list<Monomial> T;
    list<Abstract_Polynomial *> G;
    list<Critical_Pair_Dynamic *> P;
    Dense_Univariate_Integer_Polynomial * T_hn = hilbert_numerator_bigatti(M);
    Dense_Univariate_Integer_Polynomial * T_hn2
        = hilbert_second_numerator(n, T_hn);
    Dense_Univariate_Rational_Polynomial * T_hp
        = hilbert_polynomial(n, ideal_dimension(n, T_hn, T_hn2), M, T_hn, T_hn2);
    delete T_hn2;
    if (not (win_hn == nullptr or hilbert_cmp(*win_hn, *win_hp, *T_hn, *T_hp) > 0))
    {
      delete T_hn;
      delete T_hp;
    } else {
      // temporarily solve using GLPK
      GLPK_Solver * tmp_lp = new GLPK_Solver(n);
      set<Constraint> tmp_cnstr;
      bool consistent = true;
      CONSTR_TYPE w[n];
      auto Mi = M.begin();
      for (auto N : saved_mons) {
        auto t = *Mi;
        auto fi = N.begin();
        while (consistent and fi != N.end()) {
          if (not t.is_like(*fi)) {
            for (NVAR_TYPE i = 0; i < n; ++i)
              w[i] = t.degree(i) - fi->degree(i);
            Constraint c(n, w);
            consistent = tmp_lp->solve(c);
            if (consistent) tmp_cnstr.emplace(c);
          }
          ++fi;
        }
        if (not consistent)
          break;
        ++Mi;
      }
      delete tmp_lp;
      if (consistent) {
        all_constraints.clear();
        for (auto c : tmp_cnstr)
          all_constraints.push_back(c);
        if (win_hn != nullptr) {
          delete win_hn;
          delete win_hp;
        }
        win_hn = T_hn;
        win_hp = T_hp;
      }
    }
  }
  solver->solve(all_constraints);
  cout << "prefer ordering with " << solver->get_rays().size() << " rays: ";
  Ray interior_ray = ray_sum(solver->get_rays());
  cout << interior_ray << endl;
  auto wt_ptr = interior_ray.weights();
  for (NVAR_TYPE i = 0; i < n; ++i) weights[i] = wt_ptr[i];
  *mord = new ORDERING_TYPE(n, weights);
  if (win_hn != nullptr) delete win_hn;
  if (win_hp != nullptr) delete win_hp;
}

list<Abstract_Polynomial *> buchberger_dynamic(
    const list<Abstract_Polynomial *> &F,
    SPolyCreationFlags method,
    StrategyFlags strategy,
    WT_TYPE * strategy_weights,
    Dynamic_Heuristic heuristic,
    DynamicSolver solver_type,
    bool analyze_inputs
) {
  unsigned number_of_spolys = 0;
  list<Abstract_Polynomial *> G;
  list<Critical_Pair_Dynamic *> P;
  // set up list of LPPs for selecting monomials
  list<Monomial> currentLPPs;
  // create initial ordering
  NVAR_TYPE n = (*(F.begin()))->base_ring().number_of_variables();
  WT_TYPE * weights = new WT_TYPE [n];
  for (NVAR_TYPE i = 0; i < n; ++i) weights[i] = 1;
  ORDERING_TYPE * curr_ord = new ORDERING_TYPE(n, weights);
  list<ORDERING_TYPE *> all_orderings_used; // so we can free them at the end
  all_orderings_used.push_front(curr_ord);
  const WT_TYPE * wt_ptr;
  // create skeleton
  LP_Solver * solver;
  switch (solver_type) {
    case SKELETON_SOLVER: solver = new Skeleton(n); break;
    case GLPK_SOLVER: solver = new GLPK_Solver(n); break;
    case PPL_SOLVER: solver = new PPL_Solver(n); break;
    default: solver = new Skeleton(n); break;
  }
  Monomial_Ordering * dummy = curr_ord;
  if (analyze_inputs)
    initial_analysis(F, &dummy, solver);
  //G.push_back(*(F.begin()));
  // set up initial ordering
  Abstract_Polynomial * g = new Constant_Polynomial(**(F.begin()));
  cout << "Working with " << g->leading_monomial() << endl;
  switch (strategy) {
  case StrategyFlags::SUGAR_STRATEGY:
    g->set_strategy(new Poly_Sugar_Data(g));
    break;
  case StrategyFlags::WSUGAR_STRATEGY:
    g->set_strategy(new Poly_WSugar_Data(g, strategy_weights));
    break;
  default: break; // includes NORMAL_STRATEGY 
  }
  bool new_world_order = false;
  Dense_Univariate_Integer_Polynomial * hNum = nullptr;
  select_monomial(g, currentLPPs, &hNum, G, P, solver, new_world_order, heuristic);
  G.push_back(g);
  // check ordering & convert if necessary
  if (new_world_order) {
    weights = new WT_TYPE [n];
    cout << "new ordering from " << solver->get_rays().size() << " rays: ";
    Ray interior_ray = ray_sum(solver->get_rays());
    cout << interior_ray << endl;
    wt_ptr = interior_ray.weights();
    for (NVAR_TYPE i = 0; i < n; ++i) weights[i] = wt_ptr[i];
    curr_ord = new ORDERING_TYPE(n, weights);
    all_orderings_used.push_front(curr_ord);
    g->set_monomial_ordering(curr_ord);
  }
  /*cout << "Concluded with " << g->leading_monomial() << endl;
        for (auto r: solver->get_rays()) cout << '\t' << r << endl;
        cout << endl;*/
  // set up critical_pairs
  for (Abstract_Polynomial * f : F)
    if (f != *(F.begin())) {
      Critical_Pair_Dynamic * p = new Critical_Pair_Dynamic(f, strategy, curr_ord);
      P.push_back(p);
    }
  // main loop
  bool verbose = false;
  bool very_verbose = false;
  while (!P.empty()) {
    //cout << "current basis: ";
    //for (auto g : G) cout << g->leading_monomial() << ',';
    //cout << endl;
    cout << "current ordering: " << *curr_ord << endl;
    sort_pairs_by_strategy<Critical_Pair_Dynamic>(P);
    report_critical_pairs(P);
    Critical_Pair_Dynamic * p = P.front();
    P.pop_front();
    report_front_pair(p, strategy);
    cout << "\tweighted degree: " << p->lcm().weighted_degree(weights) << endl;
    /*cout << "\tothers:\n\t\t";
    for (Critical_Pair_Dynamic * pp : P) cout << pp->lcm() << ',' << static_cast<const Pair_Sugar_Data *>(pp->pair_key())->pair_sugar() << ',' << pp->lcm().weighted_degree(weights) << "; ";
    cout << endl;*/
    // make s-poly
    // DYNAMIC
    Mutable_Polynomial * s = p->s_polynomial(method, strategy);
    ++number_of_spolys;
    delete p;
    // cout << "Reducing s-poly "; s->println();
    if (!s->is_zero())
      reduce_over_basis_dynamic(&s, G);
    if (s->is_zero()) {
      cout << "\treduced to zero\n";
      delete s;
    } else {
      Abstract_Polynomial * r;
      r = new Constant_Polynomial(*s);
      // move strategy from s to r
      r->set_strategy(s->strategy());
      s->set_strategy(nullptr);
      delete s;
      cout << "\tadded " << r->leading_monomial() << endl;
      very_verbose = false;
      if (very_verbose) { cout << "\tadded "; r->println(); }
      very_verbose = false;
      // Currently performing "late conversion" of polynomials to new orderings
      // so new_world_order should be useless for the time being
      bool new_world_order = false;
      select_monomial(r, currentLPPs, &hNum, G, P, solver, new_world_order, heuristic);
      if (new_world_order) {
        cout << "new ordering from " << solver->get_rays().size() << " rays: ";
        Ray interior_ray = ray_sum(solver->get_rays());
        cout << interior_ray << endl;
        //for (auto r: solver->get_rays()) cout << '\t' << r << endl;
        //cout << endl;
        curr_ord = new ORDERING_TYPE(interior_ray);
        //cout << skel << endl;
        all_orderings_used.push_front(curr_ord);
        r->set_monomial_ordering(curr_ord);
        cout << "ended with leading monomial " << r->leading_monomial() << endl;
        //cout << "changing " << P.size() << " pairs:\n\t";
        for (Critical_Pair_Dynamic * p : P) {
          //cout << p->lcm() << " ";
          p->change_ordering(curr_ord);
        }
        cout << endl;
      }
      gm_update_dynamic(P, G, r, strategy, curr_ord);
    }
  }
  cout << number_of_spolys << " s-polynomials computed and reduced\n";
  // cleanup
  cout << G.size() << " polynomials before interreduction\n";
  //check_correctness(G, strategy);
  G = reduce_basis(G);
  cout << G.size() << " polynomials after interreduction\n";
  //set<Constant_Polynomial *, smaller_lm> B;
  // sort all polynomials to new ordering
  list<Abstract_Polynomial *> B;
  unsigned long num_mons = 0;
  unsigned long max_mons = 0;
  for (Abstract_Polynomial * g : G)
  {
    unsigned long glen = g->length();
    num_mons += glen;
    if (glen > max_mons) max_mons = glen;
    Constant_Polynomial * b = new Constant_Polynomial(*g);
    b->set_monomial_ordering(curr_ord);
    B.push_back(b);
    delete g;
  }
  cout << "tot # monomials: " << num_mons << endl;
  cout << "max # monomials: " << max_mons << endl;
  cout << "avg # monomials: " << num_mons / B.size() << endl;
  // eliminate old orderings
  all_orderings_used.pop_front();
  // first one should be curr_ord; we do not want to delete that!
  while (all_orderings_used.size() != 0) {
    ORDERING_TYPE * bye_bye_ordering = all_orderings_used.front();
    all_orderings_used.pop_front();
    delete [] bye_bye_ordering->order_weights();
    delete bye_bye_ordering;
  }
  if (hNum != nullptr) delete hNum;
  delete solver;
  return B;
}

#endif