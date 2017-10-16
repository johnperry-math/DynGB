#ifndef __SUGAR_STRATEGY_
#define __SUGAR_STRATEGY_

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
#include "polynomial.hpp"

class Critical_Pair_Basic;

/**
  @ingroup strategygroup
  @class Poly_Sugar_Data
  @author John Perry
  @date 2016
  @brief polynomial-related data for a sugar strategy
*/
class Poly_Sugar_Data : public Poly_Strategy_Data {
public:
  /**
    @brief records @c poly as the reference for @c this
    @param poly the polynomial whose sugar data @p this should be
  */
  Poly_Sugar_Data(const Abstract_Polynomial * poly);
  /**
    @brief returns @c true iff the sugars are equal
    @param sd strategy data containing sugar
    @return @p true if and only if the sugar of @p this is equivalent to the
      sugar of @p sd
  */
  virtual bool equivalent(const Poly_Strategy_Data & sd) const;
  /**
    @brief returns @c true iff @c this sugar is larger than @c sd &rsquo;s
    @param sd strategy data containing sugar
    @return @p true if and only if @p this sugar is larger than @p sd &rsquo;s
  */
  virtual bool first_larger(const Poly_Strategy_Data & sd) const;
  /**
    @brief sets the sugar to the largest total degree of a monomial in
      the assigned previously polynomial
  */
  virtual void at_generation_tasks();
  /**
    @brief sets the sugar to the largest total degree of a monomial in
      the product of the monomial and the previously assigned polynomial
    @param t a Monomial to multiply to @p this
  */
  virtual void at_generation_tasks(const Monomial & t);
  /**
    @brief re-evaluates the sugar, if need be
    @param u exponent vector of a Monomial to multiply to @p g,
        then subtract from @p this
    @param g a polynomial
  */
  virtual void pre_reduction_tasks(
    const EXP_TYPE * u, const Abstract_Polynomial & g
  );
  /** @name Basic properties */
  ///@{
  /** @brief the sugar itself */
  DEG_TYPE poly_sugar() const;
  /** @brief the strategy type */
  virtual StrategyFlags type() { return StrategyFlags::SUGAR_STRATEGY; }
  ///@}
  /** @name Modification */
  ///@{
  /** @brief for those times when a different sugar is appropriate */ 
  void force_sugar(DEG_TYPE new_sugar) { sugar = new_sugar; }
  ///@}
  /** @name I/O */
  ///@{
  friend ostream & operator <<(ostream &, const Poly_Sugar_Data &);
  ///@}
protected:
  /** @brief the polynomial&rsquo;s sugar */
  unsigned long long sugar;
};

/**
  @ingroup strategygroup
  @class Pair_Sugar_Data
  @brief ordering critical pairs using the sugar strategy
  @author John Perry
  @date 2016
*/
class Pair_Sugar_Data : public Pair_Strategy_Data {
public:
  /** @name Construction */
  ///@{
  /**
    @brief all the information we need is in @c cpb already so no additional
      processing is necessary
  */
  Pair_Sugar_Data(Critical_Pair_Basic & cpb);
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
  DEG_TYPE pair_sugar() const { return sugar; }
  ///@}
protected:
  /** @brief the critical pair to which this @c Normal_Strategy belongs */
  Critical_Pair_Basic * cp;
  /** @brief the pair*rsquo;s sugar */
  DEG_TYPE sugar;
};

#endif