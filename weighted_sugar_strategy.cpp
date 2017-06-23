#ifndef __SUGAR_STRATEGY_CPP_
#define __SUGAR_STRATEGY_CPP_

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

#include "polynomial.hpp"
#include "strategies.hpp"
#include "critical_pair.hpp"
#include "weighted_sugar_strategy.hpp"
#include "particular_orderings.hpp"

Poly_WSugar_Data::Poly_WSugar_Data(
    const Abstract_Polynomial * poly, const WT_TYPE *w
) {
  p = poly;
  weights = w;
  sugar = p->weighted_degree(w);
}

bool Poly_WSugar_Data::equivalent(const Poly_Strategy_Data & sd) const {
  const Poly_WSugar_Data * psd = static_cast<const Poly_WSugar_Data *>(&sd);
  return sugar == psd->sugar;
}

bool Poly_WSugar_Data::first_larger(const Poly_Strategy_Data & sd) const {
  const Poly_WSugar_Data * psd = static_cast<const Poly_WSugar_Data *>(&sd);
  return sugar > psd->sugar;
}

void Poly_WSugar_Data::at_generation_tasks() {
  sugar = 0;
  for (
      Polynomial_Iterator * pi = p->new_iterator();
      not pi->fellOff(); pi->moveRight()
  ) {
    if (pi->currMonomial().total_degree() > sugar)
      sugar = pi->currMonomial().weighted_degree(weights);
  }
}

void Poly_WSugar_Data::at_generation_tasks(const Monomial & t) {
  sugar = 0;
  Polynomial_Iterator * pi = p->new_iterator();
  for ( /* already initialized */ ; not pi->fellOff(); pi->moveRight()) {
    if (
        pi->currMonomial().weighted_degree(weights)
          + t.weighted_degree(weights) > sugar
    ) {
      sugar = pi->currMonomial().weighted_degree(weights)
                + t.weighted_degree(weights);
    }
  }
  delete pi;
}

void Poly_WSugar_Data::pre_reduction_tasks(
  const EXP_TYPE * u, const Abstract_Polynomial & g
) {
  Poly_WSugar_Data * sd = static_cast<Poly_WSugar_Data *>(g.strategy());
  DEG_TYPE d = sd->sugar;
  const Monomial & t = g.leading_monomial();
  for (NVAR_TYPE i = 0; i < t.num_vars(); ++i) d += weights[i]*u[i];
  if (sugar < d) sugar = d;
}

DEG_TYPE Poly_WSugar_Data::poly_sugar() const { return sugar; }

const WT_TYPE * Poly_WSugar_Data::weighting() const { return weights; }

ostream & operator <<(ostream &os, const Poly_WSugar_Data &sd) {
  os << "Weighted sugar strategy recording sugar " << sd.sugar;
  return os;
}

Pair_WSugar_Strategy::Pair_WSugar_Strategy(Critical_Pair_Basic & cpb)
{
  cp = &cpb;
  if (cpb.first()->strategy() == nullptr) {
    const CachedWGrevlex_Ordering * mord
        = dynamic_cast<const CachedWGrevlex_Ordering *>(cpb.first()->monomial_ordering());
    const WT_TYPE * w = (mord == nullptr) ? nullptr : mord->order_weights();
    sugar = cpb.first()->weighted_degree(w);
  }
  else {
    Poly_WSugar_Data * wdata
        = dynamic_cast<Poly_WSugar_Data *>(cpb.first()->strategy());
    const WT_TYPE * w = (wdata == nullptr) ? nullptr : wdata->weighting();
    sugar = wdata->poly_sugar() + cpb.first_multiplier().weighted_degree(w);
  }
  if (cpb.second() != nullptr) {
    const Abstract_Polynomial * q = cpb.second();
    if (false) cout << "second: " << *q << endl;
    Poly_WSugar_Data * wdata
        = dynamic_cast<Poly_WSugar_Data *>(q->strategy());
    const WT_TYPE * w = (wdata == nullptr) ? nullptr : wdata->weighting();
    DEG_TYPE second_sugar = wdata->poly_sugar()
          + cpb.second_multiplier().weighted_degree(w);
    if (second_sugar > sugar) sugar = second_sugar;
  }
}

bool Pair_WSugar_Strategy::equivalent(const Pair_Strategy_Data & sd) const {
  const Pair_WSugar_Strategy * nsd = static_cast<const Pair_WSugar_Strategy *>(&sd);
  bool result = (sugar == nsd->sugar)
    and (cp->lcm() == nsd->cp->lcm())
    and (cp->first_multiplier() == nsd->cp->first_multiplier());
  if (result) {
    if (cp->second() == nullptr) result = (nsd->cp->second() == nullptr);
    else {
      if (nsd->cp->second() != nullptr) result = false;
      else 
      result = ((cp->second() == nsd->cp->second())
          and (cp->second_multiplier() == nsd->cp->second_multiplier()));
    }
  }
  return result;
}

bool Pair_WSugar_Strategy::first_larger(const Pair_Strategy_Data & sd) const {
  const Pair_WSugar_Strategy * nsd = static_cast<const Pair_WSugar_Strategy *>(&sd);
  bool result = false;
  if (sugar > nsd->sugar)
    result = true;
  else if (sugar == nsd->sugar and cp->lcm() > nsd->cp->lcm())
    result = true;
  else if (sugar == nsd->sugar and cp->lcm().is_like(nsd->cp->lcm())) {
    if (cp->first_multiplier() > nsd->cp->first_multiplier())
      result = true;
    else if (
        cp->first_multiplier().is_like(nsd->cp->first_multiplier())
    ) {
      if (cp->second() == nullptr) result = false;
      else if (nsd->cp->second() == nullptr)
        result = true;
      else
        result = (cp->second_multiplier() > nsd->cp->second_multiplier());
    }
  }
  return result;
}

DEG_TYPE Pair_WSugar_Strategy::pair_sugar() { return sugar; }

#endif