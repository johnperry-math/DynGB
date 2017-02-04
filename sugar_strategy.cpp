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
* Foobar is distributed in the hope that it will be useful,                   *
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
#include "sugar_strategy.hpp"

Poly_Sugar_Data::Poly_Sugar_Data(const Abstract_Polynomial * poly){
  p = poly;
  sugar = p->standard_degree();
}

bool Poly_Sugar_Data::equivalent(const Poly_Strategy_Data & sd) const {
  const Poly_Sugar_Data * psd = static_cast<const Poly_Sugar_Data *>(&sd);
  return sugar == psd->sugar;
}

bool Poly_Sugar_Data::first_larger(const Poly_Strategy_Data & sd) const {
  const Poly_Sugar_Data * psd = static_cast<const Poly_Sugar_Data *>(&sd);
  return sugar > psd->sugar;
}

void Poly_Sugar_Data::at_generation_tasks() {
  sugar = 0;
  Polynomial_Iterator * pi = p->new_iterator();
  for ( /* already initialized */; not pi->fellOff(); pi->moveRight() ) {
    if (pi->currMonomial().total_degree() > sugar)
      sugar = pi->currMonomial().total_degree();
  }
  delete pi;
}

void Poly_Sugar_Data::at_generation_tasks(const Monomial & t) {
  sugar = 0;
  Polynomial_Iterator * pi = p->new_iterator();
  for ( /* already initialized */ ; not pi->fellOff(); pi->moveRight()) {
    if (pi->currMonomial().total_degree() + t.total_degree() > sugar)
      sugar = pi->currMonomial().total_degree() + t.total_degree();
  }
  delete pi;
}

void Poly_Sugar_Data::pre_reduction_tasks(
  const EXP_TYPE * u, const Abstract_Polynomial & g
) {
  const Monomial & t = g.leading_monomial();
  Poly_Sugar_Data * sd = static_cast<Poly_Sugar_Data *>(g.strategy());
  DEG_TYPE d = sd->sugar;
  for (NVAR_TYPE i = 0; i < t.num_vars(); ++i) d += u[i];
  if (sugar < d) sugar = d;
}

DEG_TYPE Poly_Sugar_Data::poly_sugar() const { return sugar; }

ostream & operator <<(ostream &os, const Poly_Sugar_Data &sd) {
  os << "Sugar strategy recording sugar " << sd.sugar;
  return os;
}

Pair_Sugar_Data::Pair_Sugar_Data(Critical_Pair_Basic & cpb)
{
  cp = &cpb;
  if (cpb.first()->strategy() == nullptr)
    sugar = cpb.first()->standard_degree();
  else
    sugar = ((Poly_Sugar_Data *)(cpb.first()->strategy()))->poly_sugar()
              + cpb.first_multiplier().total_degree();
  if (cpb.second() != nullptr) {
    DEG_TYPE second_sugar =
        ((Poly_Sugar_Data *)(cpb.second()->strategy()))->poly_sugar()
      + cpb.second_multiplier().total_degree();
    sugar = (sugar >= second_sugar) ? sugar : second_sugar;
  }
}

bool Pair_Sugar_Data::equivalent(const Pair_Strategy_Data & sd) const {
  const Pair_Sugar_Data * nsd = static_cast<const Pair_Sugar_Data *>(&sd);
  bool result;
  result = (sugar == nsd->sugar)
    and (cp->lcm() == nsd->cp->lcm())
    and (cp->first()->leading_monomial() == nsd->cp->first()->leading_monomial());
  if (result) {
    if (cp->second() == nullptr) result = (nsd->cp->second() == nullptr);
    else {
      if (nsd->cp->second() != nullptr) { result = false; }
      else
        result = ((cp->second() == nsd->cp->second())
             and (cp->second()->leading_monomial()
                    == nsd->cp->second()->leading_monomial()));
    }
  }
  return result;
}

bool Pair_Sugar_Data::first_larger(const Pair_Strategy_Data & sd) const {
  const Pair_Sugar_Data * nsd = static_cast<const Pair_Sugar_Data *>(&sd);
  bool result = false;
  if (sugar > nsd->sugar)
    result = true;
  else if (sugar == nsd->sugar and cp->lcm() > nsd->cp->lcm())
    result = true;
  else if (sugar == nsd->sugar and cp->lcm().is_like(nsd->cp->lcm())) {
    if (cp->first()->leading_monomial() > nsd->cp->first()->leading_monomial())
      result = true;
    else if (
        cp->first()->leading_monomial().is_like(
          nsd->cp->first()->leading_monomial()
        )
    ) {
      if (cp->second() == nullptr) result = false;
      else if (nsd->cp->second() == nullptr) result = true;
      else if (
          cp->second()->leading_monomial() > nsd->cp->second()->leading_monomial()
      ) {
        result = true;
      }
    }
  }
  return result;
}

#endif