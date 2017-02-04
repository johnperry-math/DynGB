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

F4_Reduction_Data::F4_Reduction_Data(
    Critical_Pair_Basic & p,
    list<Abstract_Polynomial *> & B
) : G(B) {
  head = len = 0;
  A = nullptr;
  switch(p.first()->strategy()->type()) {
    case NORMAL_STRATEGY: break; /* nothing to do */
    case SUGAR_STRATEGY :
      strategy = new Poly_Sugar_Data(p.first());
      break;
    case WSUGAR_STRATEGY: {
      Poly_WSugar_Data * sd
          = static_cast<Poly_WSugar_Data *>(p.first()->strategy());
      const WT_TYPE * w;
      if (sd == nullptr) w = nullptr;
      else
        w = (static_cast<const Weighted_Ordering *>(
              p.first()->monomial_ordering())
            )->order_weights();
      strategy = new Poly_WSugar_Data(p.first(), w);
      break;
    }
    default: break; /* default to normal strategy */
  }
  if (strategy != nullptr) {
    strategy->at_generation_tasks(p.first_multiplier());
    if (p.second() != nullptr)
      strategy->pre_reduction_tasks(p.second_multiplier(), *(p.second()));
  }
  mord = p.first()->monomial_ordering();
  M_build = new list<Monomial *>();
  R_build = new list<Abstract_Polynomial *>();
  // add monomials of first
  list<Abstract_Polynomial *>::iterator ri = R_build->begin();  
  add_monomials(nullptr, ri, p.first(), &p.first_multiplier());
  list<Monomial *>::iterator mi = M_build->begin();
  if (p.second() != nullptr) {
    *ri = const_cast<Abstract_Polynomial *>(p.second());
    add_monomials(&mi, ri, p.second());
    ++ri; ++mi;
  }
  // for each monomial, find an appropriate reducer
  for ( /* already initialized */ ; mi != M_build->end(); ++mi, ++ri) {
    list<Abstract_Polynomial *>::iterator g = G.begin();
    bool found = ((*ri) != nullptr);
    while (not found and g != G.end()) {
      if ((**mi).divisible_by((*g)->leading_monomial())) {
        found = true;
        *ri = *g;
        add_monomials(&mi, ri, *g);
        g = G.end();
      }
      ++g;
    }
  }
  initialize(const_cast<Abstract_Polynomial *>(p.first()), p.first_multiplier());
}

void F4_Reduction_Data::initialize(Abstract_Polynomial * p, const Monomial & t) {
  len = M_build->size();
  Prime_Field & F = p->ground_field();
  head = 0;
  M = new Monomial *[len] { nullptr };
  A = (Prime_Field_Element *)malloc(sizeof(Prime_Field_Element)*len);
  for (unsigned i = 0; i < len; ++i) {
    A[i].assign(0, &F);
  }
  Polynomial_Iterator * pi = p->new_iterator();
  list<Monomial *>::iterator mi = M_build->begin();
  nonzero_entries = 0;
  for (unsigned i = 0; i < len; ++i) {
    M[i] = *mi;
    if ((not pi->fellOff()) and (**mi).like_multiple(pi->currMonomial(), t)) {
      A[i] = pi->currCoeff();
      pi->moveRight();
      ++nonzero_entries;
    }
    ++mi;
  }
  R.resize(R_build->size());
  unsigned i = 0;
  for (Abstract_Polynomial * r : *R_build)
    R[i++] = r;
  delete pi;
  delete M_build;
  delete R_build;
  Rx = &(p->base_ring());
}

void F4_Reduction_Data::add_monomials(
    list<Monomial *>::iterator * ti,
    list<Abstract_Polynomial *>::iterator & rp,
    const Abstract_Polynomial *g,
    const Monomial * u
) {
  NVAR_TYPE n = g->leading_monomial().num_vars();
  if (ti == nullptr) { // g is a generator
    Polynomial_Iterator * pi = g->new_iterator();
    while (not (pi->fellOff())) {
      Monomial * t = new Monomial(pi->currMonomial());
      (*t) *= (*u);
      M_build->push_back(t);
      R_build->push_back(nullptr);
      pi->moveRight();
    }
    delete pi;
    rp = R_build->begin();
  } else { // g is not a generator; advance and compare while inserting
    list<Abstract_Polynomial *>::iterator ri(rp);
    list<Monomial *>::iterator tp = *ti;
    Monomial * v = new Monomial(**tp);
    (*v) /= g->leading_monomial();
    ++tp; ++ri;
    Polynomial_Iterator * pi = g->new_iterator(); pi->moveRight();
    while (not (pi->fellOff())) {
      Monomial * t = new Monomial(pi->currMonomial());
      (*t) *= (*v);
      while (tp != M_build->end() and *t < **tp) {
        ++tp; ++ri;
      }
      if (tp == M_build->end()) {
        M_build->push_back(t);
        R_build->push_back(nullptr);
      } else if (*t == **tp)
        delete t;
      else {
        M_build->insert(tp, t);
        R_build->insert(ri, nullptr);
      }
      pi->moveRight();
    }
    delete pi;
  }
}

F4_Reduction_Data::~F4_Reduction_Data() {
  for (unsigned i = 0; i < len; ++i)
    delete M[i];
  delete [] M;
  free(A);
  if (strategy != nullptr) delete strategy;
}

void F4_Reduction_Data::list_reducers() {
  for (unsigned i = 0; i < len; ++i) {
    cout << *(M[i]) << " to be reduced by ";
    if (R[i] == nullptr)
      cout << "none\n";
    else
      cout << R[i]->leading_monomial() << endl;
  }
}

bool F4_Reduction_Data::is_zero() {
  return nonzero_entries == 0;
}

Constant_Polynomial * F4_Reduction_Data::finalize() {
  Constant_Polynomial * result;
  if (nonzero_entries == 0)
    result = nullptr;
  else {
    NVAR_TYPE n = M[0]->num_vars();
    Monomial * M_final = (Monomial *)malloc(sizeof(Monomial)*nonzero_entries);
    Prime_Field_Element * A_final
        = (Prime_Field_Element *)malloc(sizeof(Prime_Field_Element)*nonzero_entries);
    unsigned i = 0;
    Prime_Field_Element scale(A[head].inverse(), A[head].field());
    for (unsigned j = head; j < len; ++j) {
      if (not A[j].is_zero()) {
        A_final[i] = A[j]*scale;
        M_final[i].common_initialization();
        M_final[i].initialize_exponents(n);
        M_final[i].set_monomial_ordering(mord);
        M_final[i] = *M[j];
        ++i;
      }
    }
    result = new Constant_Polynomial(
        nonzero_entries,
        *Rx,
        M_final, A_final,
        mord
    );
    free(M_final);
    free(A_final);
  }
  return result;
}

void F4_Reduction_Data::reduce() {
  NVAR_TYPE n = M[0]->num_vars();
  EXP_TYPE * u = new EXP_TYPE[n]; // exponents of multiplier
  // loop through our terms
  for (unsigned i = head; i < len; ++i) {
    //u /= u;
    // do we need to reduce this term, and can we?
    if ((not A[i].is_zero()) and R[i] != nullptr) {
      // get reducer for this monomial
      const Abstract_Polynomial * g = R[i];
      Polynomial_Iterator * gi = g->new_iterator();
      // determine multiplier
      for (NVAR_TYPE k = 0; k < n; ++k) {
        u[k] = (*M[i])[k] - (gi->currMonomial())[k];
      }
      // determine reduction coefficient
      Prime_Field_Element a(A[i]*gi->currCoeff().inverse());
      // loop through g's terms
      // by construction, monomials of u*g[i] should already appear in M,
      // so the line marked *** SHOULD NOT need a guard
      unsigned j = i;
      while (not gi->fellOff()) {
        const Monomial & t = gi->currMonomial();
        while (not M[j]->like_multiple(u, t)) ++j;
        bool was_zero = A[j].is_zero();
        A[j] -= a*gi->currCoeff();
        if (was_zero and not A[j].is_zero())
          ++nonzero_entries;
        else if (not was_zero and A[j].is_zero())
          --nonzero_entries;
        ++j;
        gi->moveRight();
      }
      if (strategy != nullptr)
        strategy->pre_reduction_tasks(u, *g);
      advance_head();
      delete gi;
    }
  }
  delete [] u;
}

list<Constant_Polynomial *> f4_control(
    const list<Abstract_Polynomial *> &F,
    int method,
    unsigned strategy,
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
      case NORMAL_STRATEGY: break; // don't need polynomial data
      case SUGAR_STRATEGY:
        f->set_strategy(new Poly_Sugar_Data(f));
        break;
      case WSUGAR_STRATEGY:
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
    //Mutable_Polynomial * s = p->s_polynomial(method, strategy);
    F4_Reduction_Data s(*p, G);
    //s.list_reducers();
    ++number_of_spolys;
    // cout << "Reducing s-poly "; s->println();
    if (not s.is_zero())
      s.reduce();
    if (s.is_zero()) {
      cout << "\treduced to zero\n";
      // delete s;
    } else {
      Abstract_Polynomial * r = s.finalize();
      // move strategy from s to r
      r->set_strategy(s.get_strategy());
      s.clear_strategy();
      //delete s;
      cout << "\tadded " << r->leading_monomial() << endl;
      if (very_verbose) { cout << "\tadded "; r->println(); }
      gm_update(P, G, r, strategy);
    }
    if (p->is_generator())
      delete p->first();
    delete p;
  }
  cout << number_of_spolys << " s-polynomials computed and reduced\n";
  // cleanup
  cout << G.size() << " polynomials before interreduction\n";
  //check_correctness(G, strategy);
  G = reduce_basis(G);
  cout << G.size() << " polynomials after interreduction\n";
  list<Constant_Polynomial *> B;
  for (Abstract_Polynomial * g : G) {
    B.push_back(new Constant_Polynomial(*g));
    //if (F.find(g) == F.end()) delete g;
  }
  return B;
}

#endif