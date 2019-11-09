#ifndef __ALGORITHM_BUCHBERGER_BASIC_CPP_
#define __ALGORITHM_BUCHBERGER_BASIC_CPP_

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

#include "algorithm_buchberger_basic.hpp"

#include "reduction_support.hpp"

template bool no_triplet<Critical_Pair_Basic>(
    const Critical_Pair_Basic *, const list<Critical_Pair_Basic *>
);

template void report_critical_pairs<Critical_Pair_Basic>(
    const list<Critical_Pair_Basic *>, bool
);

template void sort_pairs_by_strategy<Critical_Pair_Basic>(
    list<Critical_Pair_Basic *> &
);

bool lcm_alike(const Monomial & t, const Monomial & u,
               const Critical_Pair_Basic * p)
{
  bool result = true;
  for (NVAR_TYPE i = 0; result and i < t.num_vars(); ++i)
    result = result and p->lcm_degree(i) == ((t.degree(i) >= u.degree(i)) ?
                                                  t.degree(i) : u.degree(i));
  return result;
}

list<Abstract_Polynomial *> reduce_basis(list<Abstract_Polynomial *>G) {
  list <Abstract_Polynomial *> result;
  set <Abstract_Polynomial *> eliminate; // save polynomials for deletion
  // we identify polynomials that need deletion, but wait to delete
  for (Abstract_Polynomial * f : G)
  {
    bool eliminated = false;
    // look for divisors of f in G
    for (Abstract_Polynomial * g : G) {
      if (g != f and g->leading_monomial() | f->leading_monomial()) {
        eliminated = true;
        eliminate.insert(f);
        break;
      }
    }
    // look for divisors of f in result
    if (not eliminated)
      for (Abstract_Polynomial * g : result) {
        if (g != f and g->leading_monomial() | f->leading_monomial()) {
          eliminated = true;
          eliminate.insert(f);
          break;
        }
      }
    if (not eliminated) {
      result.push_back(f);
    }
  }
  // eliminate the eliminated
  for (Abstract_Polynomial * g : eliminate)
    delete g;
  return result;
}

void check_correctness(
    list<Abstract_Polynomial *>G, StrategyFlags strategy, EXP_TYPE max_degree
) {
  cout << "not-so-quick check for correctness\n";
  bool verbose = true;
  for (auto fi = G.begin(); fi != G.end(); ++fi)
    for (auto gi = next(fi); gi != G.end(); ++gi)
    {
      if (max_degree == 0 or (*fi)->leading_monomial().lcm((*gi)->leading_monomial()).total_degree() <= max_degree) {
        Critical_Pair_Basic * p = new Critical_Pair_Basic(*fi, *gi, strategy);
        if (verbose)
          cout << "checking " << p->first()->leading_monomial() << " , "
               << p->second()->leading_monomial() << " : ";
        Mutable_Polynomial * s = p->s_polynomial(
            SPolyCreationFlags::LINKED_LST, strategy
        );
        reduce_over_basis<list<Abstract_Polynomial *>>(&s, G);
        if (s->is_zero()) {
          if (verbose) cout << "checks out\n";
        } else {
          cout << "\tfailure with " << p->first()->leading_monomial() << ','
               << p->second()->leading_monomial() << ':' << s->leading_monomial()
               << endl;
          if (verbose) cout << '\t' << *s << endl << endl;
        }
        delete s;
      }
    }
}

void gm_update(
    list<Critical_Pair_Basic *> & P,
    list<Abstract_Polynomial *> & G,
    Abstract_Polynomial * r,
    StrategyFlags strategy
) {
  //cout << "----------------------\n";
  list<Critical_Pair_Basic *> C;
  // critical pairs with new polynomial
  for (Abstract_Polynomial * g : G)
    C.push_back(new Critical_Pair_Basic(g, r, strategy));
  // apply Buchberger's lcm criterion to new pairs
  list<Critical_Pair_Basic *> D;
  while (C.size() != 0) {
    Critical_Pair_Basic * p = C.front();
    C.pop_front();
    if ((p->first()->leading_monomial().is_coprime(
            p->second()->leading_monomial()))
        or (no_triplet(p, C) and no_triplet(p, D))
        )
      D.push_back(p);
    else {
      delete p;
      //cout << "triplet prunes " << *p << endl;
    }
  }
  // apply Buchberger's gcd criterion
  list<Critical_Pair_Basic *> E;
  while (D.size() != 0) {
    Critical_Pair_Basic * p = D.front();
    D.pop_front();
    if (!(p->first()->leading_monomial().is_coprime(
            p->second()->leading_monomial())))
      E.push_back(p);
    else {
      delete p;
      // cout << "gcd prunes " << *p << endl;
    }
  }
  // apply Buchberger's lcm criterion to old pairs
  list<Critical_Pair_Basic *> Q;
  while (P.size() != 0) {
    Critical_Pair_Basic * p = P.front();
    P.pop_front();
    bool crit1 = !(r->leading_monomial() | p->lcm());
    bool crit2 = lcm_alike(p->first()->leading_monomial(), r->leading_monomial(), p);
    bool crit3 = p->second() != nullptr and
          lcm_alike(p->second()->leading_monomial(), r->leading_monomial(), p);
    if ( crit1 or crit2 or crit3)
      Q.push_back(p);
    else {
      if (p->s_polynomial() != nullptr) delete p->s_polynomial();
      delete p;
      //cout << "triplet prunes " << *p << endl;
    }
  }
  P = Q;
  // add new pairs to old pairs
  for (Critical_Pair_Basic * e : E)
    P.push_back(e);
  /*cout << "All pairs:\n";
  for (auto pi = P.begin(); pi != P.end(); ++pi)
    cout << '\t' << **pi << endl;
  cout << "----------------------\n";*/
  // add new poly to basis
  G.push_back(r);
}

void report_basis(
    list<Abstract_Polynomial *> G,
    bool verbose,
    bool very_verbose
) {
  cout << G.size() << " polys in basis\n";
  if (verbose) {
    for (Abstract_Polynomial * g : G) cout << g->leading_monomial() << '\t';
    cout << endl;
  }
  if (very_verbose) {
    for (Abstract_Polynomial * g : G) g->println();
  }
}

void report_front_pair(Critical_Pair_Basic *p, StrategyFlags strategy) {
  cout << "processing pair: " << *p << endl;
  if (
      strategy == StrategyFlags::SUGAR_STRATEGY or
      strategy == StrategyFlags::WSUGAR_STRATEGY
  ) {
    cout << "\tsugar: "
         << (static_cast<const Pair_Sugar_Data *>((p->pair_key()))->pair_sugar())
         << endl;
  }
}

list<Abstract_Polynomial *> buchberger(
    const list<Abstract_Polynomial *> &F,
    SPolyCreationFlags method,
    StrategyFlags strategy,
    WT_TYPE * strategy_weights
) {
  unsigned number_of_spolys = 0;
  list<Abstract_Polynomial *> G;
  list<Critical_Pair_Basic *> P;
  // set up basis with generators
  for (Abstract_Polynomial * fo : F)
  {
    Constant_Polynomial * f = new Constant_Polynomial(*fo);
    switch(strategy) {
      case StrategyFlags::NORMAL_STRATEGY: break; // don't need polynomial data
      case StrategyFlags::SUGAR_STRATEGY:
        f->set_strategy(new Poly_Sugar_Data(f));
        break;
      case StrategyFlags::WSUGAR_STRATEGY:
        f->set_strategy(new Poly_WSugar_Data(f, strategy_weights));
        break;
      default: break; // assume normal strategy
    }
    if (f->strategy() != nullptr) { f->strategy()->at_generation_tasks(); }
    P.push_back(new Critical_Pair_Basic(f, strategy));
  }
  // main loop
  bool verbose = false;
  bool very_verbose = false;
  while (!P.empty()) {
    sort_pairs_by_strategy(P);
    report_critical_pairs(P);
    Critical_Pair_Basic * p = P.front();
    report_front_pair(p, strategy);
    P.pop_front();
    // make s-poly
    Mutable_Polynomial * s = p->s_polynomial(method, strategy);
    ++number_of_spolys;
    if (p->is_generator())
      delete p->first();
    delete p;
    // cout << "Reducing s-poly "; s->println();
    if (!s->is_zero())
      reduce_over_basis<list<Abstract_Polynomial *>>(&s, G);
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
      if (very_verbose) { cout << "\tadded "; r->println(); }
      gm_update(P, G, r, strategy);
    }
  }
  cout << number_of_spolys << " s-polynomials computed and reduced\n";
  // cleanup
  cout << G.size() << " polynomials before interreduction\n";
  //check_correctness(G, strategy);
  G = reduce_basis(G);
  cout << G.size() << " polynomials after interreduction\n";
  //set<Constant_Polynomial *, smaller_lm> B;
  list<Abstract_Polynomial *> B;
  unsigned long num_mons = 0;
  unsigned long max_mons = 0;
  for (Abstract_Polynomial * g : G) {
    unsigned long glen = g->length();
    num_mons += glen;
    if (glen > max_mons) max_mons = glen;
    B.push_back(new Constant_Polynomial(*g));
  }
  cout << "tot # monomials: " << num_mons;
  cout << "max # monomials: " << max_mons;
  cout << "avg # monomials: " << num_mons / B.size() << endl;
  return B;
}

#endif
