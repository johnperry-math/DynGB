#ifndef __NORMAL_STRATEGY_CPP_
#define __NORMAL_STRATEGY_CPP_

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
#include "normal_strategy.hpp"

Normal_Strategy::Normal_Strategy(Critical_Pair_Basic & cpb)
{
  cp = &cpb;
}

bool Normal_Strategy::equivalent(const Pair_Strategy_Data & sd) const {
  const Normal_Strategy * nsd = static_cast<const Normal_Strategy *>(&sd);
  if (cp->second() == nullptr or nsd->cp->second() == nullptr)
    return false;
  else {
    return (cp->lcm() == nsd->cp->lcm())
      and (cp->first()->leading_monomial() == nsd->cp->first()->leading_monomial())
      and (cp->second()->leading_monomial()
              == nsd->cp->second()->leading_monomial());
  } 
}

bool Normal_Strategy::first_larger(const Pair_Strategy_Data & sd) const {
  const Normal_Strategy * nsd = static_cast<const Normal_Strategy *>(&sd);
  bool result = false;
  if (cp->lcm() > nsd->cp->lcm())
    result = true;
  else if (cp->lcm().is_like(nsd->cp->lcm())) {
    if (cp->first()->leading_monomial() > nsd->cp->first()->leading_monomial())
      result = true;
    else if (
        cp->first()->leading_monomial().is_like(
          nsd->cp->first()->leading_monomial()
        )
    ) {
      if (cp->second() != nullptr and nsd->cp->second() != nullptr)
        result = (cp->second()->leading_monomial()
                    > nsd->cp->second()->leading_monomial());
    }
  }
  return result;
}

ostream & operator <<(ostream &os, const Normal_Strategy &ns) {
  os << "Normal Strategy recording lcm " << ns.cp->lcm();
  return os;
}

#endif