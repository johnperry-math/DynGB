#ifndef __NORMAL_STRATEGY_HPP_
#define __NORMAL_STRATEGY_HPP_

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

#include "strategies.hpp"

class Critical_Pair_Basic;

/**
  @ingroup strategygroup
  @class Normal_Strategy
  @brief ordering critical pairs using the normal strategy
  @author John Perry
  @date 2016
*/
class Normal_Strategy : public Pair_Strategy_Data {
public:
  /** @name Construction */
  ///@{
  /**
    @brief all the information we need is in @c cpb already so no additional
      processing is necessary
  */
  explicit Normal_Strategy(Critical_Pair_Basic & cpb);
  ///@}
  /** @name Comparison */
  ///@{
  /** @brief implementation of equivalent() */
  virtual bool equivalent(const Pair_Strategy_Data & sd) const;
  /** @brief implementation of first_larger() */
  virtual bool first_larger(const Pair_Strategy_Data & sd) const;
  ///@}
  /** @name Basic properties */
  ///@{
  virtual StrategyFlags type() { return StrategyFlags::NORMAL_STRATEGY; }
  ///@}
  /** @name I/O */
  ///@{
  friend ostream & operator <<(ostream &, const Normal_Strategy &);
  ///@}
protected:
  /** @brief the critical pair to which this @c Normal_Strategy belongs */
  Critical_Pair_Basic * cp;
};

#endif