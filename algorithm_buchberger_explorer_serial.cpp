#ifndef __ALGORITHM_BUCHBERGER_EXPLORER_CPP_
#define __ALGORITHM_BUCHBERGER_EXPLORER_CPP_

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

#include "hilbert_functions.hpp"

#include "algorithm_buchberger_basic.hpp"

#include "algorithm_buchberger_explorer.hpp"

#include "reduction_support.hpp"

/**
  @brief compares the Hilbert function at position @c i (&ldquo;older&rdquo;)
    with the Hilbert function at position @c j (&ldquo;newer&rdquo;).
    Returns @c true iff j&rsquo;s Hilbert function is measurably better than
    i&rsquo;s.
*/
bool newer_HF_smaller(Monomial & t, unsigned i, Monomial & u, unsigned j,
    Dense_Univariate_Rational_Polynomial ** HP,
    Dense_Univariate_Integer_Polynomial ** HFN
) {
  // this adapts LessByHilbert in dynamic_engine.cpp
  // in this case, we want the opposite result,
  // so the return statement at the end is negated
  bool result;
  NVAR_TYPE n = t.num_vars();
  // first check the coefficients of the Hilbert polynomial
  Dense_Univariate_Rational_Polynomial HPdiff(*(HP[i]));
  HPdiff -= *(HP[j]);
  if (not HPdiff.is_zero())
    result = (HPdiff.numerator(HPdiff.degree()) < 0);
  else // use Hilbert series
  {
    Dense_Univariate_Integer_Polynomial * h1 = HFN[i];
    Dense_Univariate_Integer_Polynomial * h2 = HFN[j];
    DEG_TYPE i = 0;
    for ( /* already initialized */ ;
          i <= h1->degree() and i <= h2->degree() and (*h1)[i] == (*h2)[i];
          i++)
    { /* taken care of in loop */ }
    if (i > h1->degree())
    {
      if (i > h2->degree())
      { // the numerators are equal; break tie via lex
        int i = 0;
        while (i < n and t[i] == u[i]) ++i;
        if (i == n) result = false;
        else result = (t.degree(i) > u.degree(i));
      }
      else
        result = true;
    }
    else
    {
      if (i > h2->degree()) result = false;
      else result = (*h1)[i] < (*h2)[i];
    }
  }
  //cout << "\tfirst less than second? " << result << endl;
  return not result;
}

typedef Critical_Pair_Basic * Critical_Pair_Basic_Ptr;

list<Constant_Polynomial *> buchberger_explorer(
    const list<Abstract_Polynomial *> &F,
    SPolyCreationFlags method,
    StrategyFlags strategy,
    WT_TYPE * strategy_weights,
    const unsigned number_to_advance
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
  Critical_Pair_Basic ** Pp = new Critical_Pair_Basic_Ptr[number_to_advance];
  Dense_Univariate_Rational_Polynomial ** HP
      = new Dense_Univariate_Rational_Polynomial *[number_to_advance];
  Dense_Univariate_Integer_Polynomial ** HFN
      = new Dense_Univariate_Integer_Polynomial *[number_to_advance];
  Dense_Univariate_Integer_Polynomial ** HSN
      = new Dense_Univariate_Integer_Polynomial *[number_to_advance];
  list<Monomial> T;
  Mutable_Polynomial * s;
  // create, reduce s-polynomials
  while (!P.empty()) {
    sort_pairs_by_strategy(P);
    report_critical_pairs(P);
    // first select number_to_advance pairs for reduction, store in Pp, reduce
    unsigned i = 0;
    for (/* already initialized */; !P.empty() and i < number_to_advance; ++i) {
      Pp[i] = P.front();
      report_front_pair(Pp[i], strategy);
      P.pop_front();
      // make s-poly
      if ((s = Pp[i]->s_polynomial()) == nullptr)
        s = Pp[i]->s_polynomial(method, strategy);
      // cout << "Reducing s-poly "; s->println();
      if (!s->is_zero())
        reduce_over_basis(&s, G);
      Pp[i]->set_spoly(s);
    }
    // create, compare Hilbert functions
    while (i < number_to_advance) { Pp[i] = nullptr; ++i; }
    unsigned winning_index = number_to_advance + 1;
    for (i = 0; i < number_to_advance and Pp[i] != nullptr; ++i) {
      if (Pp[i] != nullptr) {
        s = Pp[i]->s_polynomial();
        if (s->is_zero()) {
          cout << "\treduced to zero\n";
          HP[i] = nullptr;
          HFN[i] = HSN[i]= nullptr;
        } else {
          T.push_back(s->leading_monomial());
          unsigned n = T.front().num_vars();
          HFN[i] = hilbert_numerator_bigatti(T);
          HSN[i] = hilbert_second_numerator(n, HFN[i]);
          unsigned d = ideal_dimension(n, HFN[i], HSN[i]);
          HP[i] = hilbert_polynomial(n, d, T, HFN[i], HSN[i]);
          T.pop_back();
          if (winning_index > number_to_advance)
            winning_index = i;
          else {
            if (newer_HF_smaller(
                Pp[winning_index]->s_polynomial()->leading_monomial(),
                winning_index,
                Pp[i]->s_polynomial()->leading_monomial(),
                i,
                HP, HFN))
              winning_index = i;
          }
        }
      }
    }
    // move result of winning reduction to basis; return others to P
    cout << "sorting results\n";
    for (i = 0; i < number_to_advance and Pp[i] != nullptr; ++i) {
      delete HP[i]; delete HFN[i]; delete HSN[i];
      if (i != winning_index or Pp[i]->s_polynomial()->is_zero()) {
        // not the winning pair: move to front of list of pairs
        // (for further reduction)
        if (not Pp[i]->s_polynomial()->is_zero())
          P.push_front(Pp[i]);
        else {
          ++number_of_spolys;
          cout << "\treduced to zero\n";
          s = Pp[i]->s_polynomial();
          //cout << "deleting " << s << endl;
          delete s;
          delete Pp[i];
        }
      } else { // winning pair
        ++number_of_spolys;
        cout << "selected "
             << Pp[i]->first()->leading_monomial() << ','
             << ((Pp[i]->is_generator()) ? 0 : Pp[i]->second()->leading_monomial())
             << endl;
        if (Pp[i]->is_generator()) delete Pp[i]->first();
        s = Pp[i]->s_polynomial();
        Abstract_Polynomial * r = new Constant_Polynomial(*s);
        delete Pp[i];
        // move strategy from s to r
        r->set_strategy(s->strategy());
        s->set_strategy(nullptr);
        //cout << "deleting " << s << endl;
        delete s;
        cout << "\tadded " << r->leading_monomial() << endl;
        if (very_verbose) { cout << "\tadded "; r->println(); }
        gm_update(P, G, r, strategy);
      }
    }
  }
  delete [] Pp;
  delete [] HP; delete [] HFN; delete [] HSN;
  cout << number_of_spolys << " s-polynomials computed and reduced\n";
  // cleanup
  cout << G.size() << " polynomials before interreduction\n";
  //check_correctness(G, strategy);
  G = reduce_basis(G);
  cout << G.size() << " polynomials after interreduction\n";
  //set<Constant_Polynomial *, smaller_lm> B;
  list<Constant_Polynomial *> B;
  for (Abstract_Polynomial * g : G) {
    B.push_back(new Constant_Polynomial(*g));
    //if (F.find(g) == F.end()) delete g;
  }
  return B;
}

#endif