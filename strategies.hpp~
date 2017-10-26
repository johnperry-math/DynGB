#ifndef __STRATEGIES_HPP_
#define __STRATEGIES_HPP_

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

#include "monomial.hpp"

class Abstract_Polynomial;

/**
  @enum StrategyFlags
  @author John Perry
  @date 2016
  @brief flag indicating which strategy to use for computation
  @ingroup strategygroup
*/
enum class StrategyFlags {
  NORMAL_STRATEGY = 1,
  SUGAR_STRATEGY,
  WSUGAR_STRATEGY
};

/**
  @class Poly_Strategy_Data
  @author John Perry
  @date 2016
  @brief polynomial-related strategy data
  @ingroup strategygroup
  @details Some strategies will not need this class, but others will,
    such as sugar (to store a polynomial&rsquo;s sugar).
*/
class Poly_Strategy_Data {
public:
  /** @name Destruction */
  ///@{
  virtual ~Poly_Strategy_Data() { };
  ///@}
  /** @name Basic properties */
  ///@{
  virtual StrategyFlags type() = 0;
  ///@}
  /** @name Comparison */
  ///@{
  /**
    @brief should return @c true iff strategy considers @c this and other equivalent
  */
  virtual bool equivalent(const Poly_Strategy_Data &) const = 0;
  /** @brief alias for equivalent() */
  bool operator==(const Poly_Strategy_Data &other) const;
  /**
    @brief should return @c true iff strategy considers @c this larger than other
  */
  virtual bool first_larger(const Poly_Strategy_Data &) const = 0;
  /** @brief alias for first_larger() */
  bool operator >(const Poly_Strategy_Data &other) const;
  /** @brief is @c this larger than or equivalent to other? */
  bool operator >=(const Poly_Strategy_Data &other) const;
  /** @brief is @c this smaller than other? */
  bool operator <(const Poly_Strategy_Data &other) const;
  /** @brief is @c this smaller than or equivalent to other? */
  bool operator <=(const Poly_Strategy_Data &other) const;
  ///@}
  /** @name Computation */
  ///@{
  /**
    @brief hook called while first generating polynomial
    @details The default is to do nothing, which is good for the normal strategy.
      Other strategies, however, may perform some record-keeping; the sugar
      strategy, for instance, will compute the polynomial&rsquo;s sugar.
  */
  virtual void at_generation_tasks() { }
  /**
    @brief hook called while first generating a polynomial \e multiple
    @param t a Monomial that we are about to multiply to create an s-polynomial
    @details The default is to do nothing, which is good for the normal strategy.
      Other strategies, however, may perform some record-keeping; the sugar
      strategy, for instance, will compute the polynomial&rsquo;s sugar.
  */
  virtual void at_generation_tasks(const Monomial & t) { }
  /**
    @brief hook called while finding a reducer
    @param r polynomial that is to be reduced
    @param g polynomial that reduces @c r
    @return @p true if and only if it is valid to reduce @p r by @p g
    @details The default is to do nothing, which is good for the normal and sugar
      strategies. Other strategies, however, may impose constraints on the
      reduction; involutive and signature-based reduction, for instance, both
      forbid certain reductions.
  */
  virtual bool valid_reduction(
      const Abstract_Polynomial & r, const Abstract_Polynomial & g
  ) const;
  /**
  */
  void pre_reduction_tasks(const Monomial & u, const Abstract_Polynomial & g) {
    pre_reduction_tasks(u.log(), g);
  }
  /**
    @brief hook called immediately before performing reduction
    @param g polynomial that reduces @c r (where @c r is @c this)
    @param u exponents of monomial satisfying @f$u\textrm{lm}(g)=\textrm{lm}(r)@f$
    @details The default is to do nothing, which is good for the normal strategy.
      Other strategies, however, may require some processing before reduction;
      the sugar strategy is an example, as reduction may increase
      a polynomial&rsquo;s sugar
  */
  virtual void pre_reduction_tasks(
      const EXP_TYPE * u, const Abstract_Polynomial & g
  ) { }
  ///@}
  /** @name I/O */
  ///@{
  /** @brief print strategy-related data in the polynomial */
  friend ostream & operator <<(ostream &, const Poly_Strategy_Data &);
  ///@}
protected:
  /** @brief the polynomial to which this strategy applies */
  const Abstract_Polynomial * p;
};

/**
  @class Pair_Strategy_Data
  @author John Perry
  @date 2016
  @brief Structure for sorting critical pairs.
  @ingroup strategygroup
*/
class Pair_Strategy_Data {
public:
  /** @name Destruction */
  ///@{
  virtual ~Pair_Strategy_Data();
  ///@}
  /** @name Comparison */
  ///@{
  /** @brief should return @c true iff @c this and other are equivalent */
  virtual bool equivalent(const Pair_Strategy_Data &) const = 0;
  /** @brief alias for equivalent() */
  bool operator ==(const Pair_Strategy_Data & sd) const;
  /** @brief should return @c true iff @c this is larger than other */
  virtual bool first_larger(const Pair_Strategy_Data &) const = 0;
  /** @brief alias for first_larger() */
  bool operator >(const Pair_Strategy_Data & sd) const;
  /** @brief is @c this larger than or equivalent to other? */
  bool operator >=(const Pair_Strategy_Data & sd) const;
  /** @brief is @c this smaller than other? */
  bool operator <(const Pair_Strategy_Data & sd) const;
  /** @brief is @c this smaller than or equivalent to other? */
  bool operator <=(const Pair_Strategy_Data & sd) const;
  ///@}
  /** @name Computation */
  ///@{
  /**
    @brief hook called immediately before computing a new s-polynomiald
    @details The default is to do nothing, which is good for the normal strategy.
      Other strategies, however, may require some processing before reduction;
      the sugar strategy is an example, as an s-polynomial needs to record
      a new polynomial&rsquo;s sugar
  */
  virtual void pre_spolynomial_tasks() const;
  ///@}
};

#endif