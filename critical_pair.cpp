#ifndef __CRITICAL_PAIR_CPP_
#define __CRITICAL_PAIR_CPP_

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

#include "critical_pair.hpp"

#include "polynomial.hpp"
#include "strategies.hpp"
#include "sugar_strategy.hpp"
#include "weighted_sugar_strategy.hpp"
#include "particular_orderings.hpp"

ostream & operator<<(ostream & os, const Critical_Pair_Basic & p) {
  os << "(";
  p.p->leading_monomial().print(false, os, p.p->base_ring().name_list());
  os << ',';
  if (p.q == nullptr) cout << 0;
  else p.q->leading_monomial().print(false, os, p.q->base_ring().name_list());
  os << ';';
  p.tpq.print(false, os, p.p->base_ring().name_list());
  os << ')';
  return os;
}

ostream & operator<<(ostream & os, const Critical_Pair_Dynamic & p) {
  os << static_cast<const Critical_Pair_Basic &>(p);
  os << p.ordering;
  return os;
}

Critical_Pair_Basic::Critical_Pair_Basic(
    Abstract_Polynomial * f,
    StrategyFlags strategy
) : tpq(f->leading_monomial()), p(f), q(nullptr), s(nullptr),
  tp(tpq.num_vars()), tq(tpq.num_vars())
{
  tpq.set_monomial_ordering(f->monomial_ordering());
  tp.set_monomial_ordering(f->monomial_ordering());
  tp.set_monomial_ordering(f->monomial_ordering());
  switch(strategy) {
  case StrategyFlags::NORMAL_STRATEGY:
    key = new Normal_Strategy(*this);
    break;
  case StrategyFlags::SUGAR_STRATEGY :
    key = new Pair_Sugar_Data(*this);
    break;
  case StrategyFlags::WSUGAR_STRATEGY:
    key = new Pair_WSugar_Strategy(*this);
    break;
  default: key = new Normal_Strategy(*this); break;
  }
}

Critical_Pair_Basic::Critical_Pair_Basic(
    Abstract_Polynomial * f,
    Abstract_Polynomial * g,
    StrategyFlags strategy
) : tpq(f->leading_monomial().lcm(g->leading_monomial())),
    p((f->leading_monomial().larger_than(g->leading_monomial())) ? f : g ),
    q((f->leading_monomial().larger_than(g->leading_monomial())) ? g : f ),
    s(nullptr),
    tp(tpq),
    tq(tpq)
{
  tpq.set_monomial_ordering(f->monomial_ordering());
  tp.set_monomial_ordering(f->monomial_ordering());
  tp.set_monomial_ordering(f->monomial_ordering());
  tp /= p->leading_monomial();
  tq /= q->leading_monomial();
  switch(strategy) {
  case StrategyFlags::NORMAL_STRATEGY:
    key = new Normal_Strategy(*this);
    break;
  case StrategyFlags::SUGAR_STRATEGY :
    key = new Pair_Sugar_Data(*this);
    break;
  case StrategyFlags::WSUGAR_STRATEGY:
    key = new Pair_WSugar_Strategy(*this);
    break;
  default: key = new Normal_Strategy(*this); break;
  }
}

Critical_Pair_Basic::Critical_Pair_Basic(const Critical_Pair_Basic & other) :
    tpq(other.tpq), tp(other.tp), tq(other.tq),
    p(other.p), q(other.q), s(other.s),
    key(other.key)
{ /* done */ }

Mutable_Polynomial * Critical_Pair_Basic::s_polynomial(
    SPolyCreationFlags method, StrategyFlags strategy
) {
  if (s != nullptr) return s;
  if (key != nullptr) key->pre_spolynomial_tasks();
  // create multiple of first() according to the indicated method
  switch (method) {
  case SPolyCreationFlags::LINKED_LST:
    s = new Polynomial_Linked_List(*p);
    break;
  case SPolyCreationFlags::GEOBUCKETS:
    s = new Polynomial_Geobucket(*p);
    break;
  case SPolyCreationFlags::DOUBLE_BUF:
    s = new Double_Buffered_Polynomial(*p);
    break;
  default:
    s = new Polynomial_Geobucket(*p);
    break;
  }
  //cout << "s-poly: " << *s << endl;
  switch(strategy) {
  case StrategyFlags::NORMAL_STRATEGY:
    break; /* nothing to do */
  case StrategyFlags::SUGAR_STRATEGY :
    s->set_strategy(new Poly_Sugar_Data(s)); break;
  case StrategyFlags::WSUGAR_STRATEGY: {
    Poly_WSugar_Data * sd = static_cast<Poly_WSugar_Data *>(p->strategy());
    const WT_TYPE * w;
    if (sd == nullptr) w = nullptr;
    else
      w = (static_cast<const Weighted_Ordering *>(p->monomial_ordering()))
            ->order_weights();
    s->set_strategy(new Poly_WSugar_Data(s, w));
    break;
  }
  default: break; /* default to normal strategy */
  }
  if (s->strategy() != nullptr) { s->strategy()->at_generation_tasks(tp); }
  s->multiply_by_monomial(tp);
  //cout << "after multiply by monomial " << tp << " : " << *s << endl;
  if (q != nullptr) {
    Prime_Field_Element a = s->leading_coefficient();
    if (not q->leading_coefficient().is_one())
      s->multiply_by_scalar(q->leading_coefficient());
    if (s->strategy() != nullptr) { s->strategy()->pre_reduction_tasks(tq, *q); }
    s->add_polynomial_multiple(a, tq, *q, true);
  }
  //cout << "after adding other poly: " << *s << endl;
  return s;
}

Mutable_Polynomial * Critical_Pair_Dynamic::s_polynomial(
  SPolyCreationFlags method, StrategyFlags strategy
) {
  if (s != nullptr) return s;
  if (p->monomial_ordering() != ordering)
    p->set_monomial_ordering(ordering);
  if (q != nullptr and q->monomial_ordering() != ordering)
    q->set_monomial_ordering(ordering);
  switch(strategy) {
    case StrategyFlags::NORMAL_STRATEGY: break; /* nothing to do */
    case StrategyFlags::SUGAR_STRATEGY: break; /* nothing to do? */
    case StrategyFlags::WSUGAR_STRATEGY: {
      if (p->strategy() == nullptr)
        p->set_strategy(new Poly_WSugar_Data(p, ordering->order_weights()));
      else {
        Poly_WSugar_Data * sd = static_cast<Poly_WSugar_Data *>(p->strategy());
        sd->change_weights(ordering->order_weights());
      }
      if (q != nullptr and q->strategy() == nullptr)
        q->set_strategy(new Poly_WSugar_Data(q, ordering->order_weights()));
      else if (q != nullptr) {
        Poly_WSugar_Data * sd = static_cast<Poly_WSugar_Data *>(q->strategy());
        sd->change_weights(ordering->order_weights());
      }
    }
    default: break; /* default to normal strategy */
  }
  return Critical_Pair_Basic::s_polynomial(method, strategy);
};

#endif