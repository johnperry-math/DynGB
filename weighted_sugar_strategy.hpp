#ifndef __WEIGHTED_SUGAR_STRATEGY_HPP_
#define __WEIGHTED_SUGAR_STRATEGY_HPP_

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

#include "system_constants.hpp"

#include "strategies.hpp"
#include "polynomial.hpp"

class Critical_Pair_Basic;

/**
  @ingroup strategygroup
  @class Poly_WSugar_Data
  @author John Perry
  @date 2016
  @brief polynomial-related data for a weighted sugar strategy
  @details This strategy is based on a weight vector,
    rather than simply adding exponents. While the weight vector \e can be
    that of a weighted ordering, one is also free to use a different weight
    altogether. Simply supply the relevant data.
*/
class Poly_WSugar_Data : public Poly_Strategy_Data {
public:
  /**
    @brief records \c poly as the reference for \c this
    @details This expects the weights to be as long as the number of monomials.
      It copies the list, so you can subsequently modify it without subsequently
      affecting correctness.
    @warning The weights are not copied. Please do not modify them, as behavior
      is hard to predict in this case. It is predictable, but it&rsquo;s likely
      highly undesirable. Better to create a new Poly_WSugar_Data instead.
    @param poly polynomial whose sugar @c this is
    @param w weight vector to compute weighted sugar
  */
  Poly_WSugar_Data(const Abstract_Polynomial * poly, const WT_TYPE * w);
  /**
    @brief returns \c true iff the sugars are equal
    @param sd sugar data to compare with @c this
    @return true if and only if @c this has sugar comparable to @p sd
  */
  virtual bool equivalent(const Poly_Strategy_Data & sd) const;
  /**
    @brief returns \c true iff \c this sugar is larger than \c sd &rsquo;s
    @param sd sugar data to compare with @c this
    @return @p true if and only if @c this has larger sugar than @p sd
  */
  virtual bool first_larger(const Poly_Strategy_Data & sd) const;
  /**
    @brief sets the sugar to the largest weighted degree of a monomial in
      the assigned previously polynomial
  */
  virtual void at_generation_tasks();
  /**
    @brief sets the sugar to the largest weighted degree of a monomial in
      the product of the monomial and the previously assigned polynomial
    @param t a Monomial to multiply to @c this
  */
  virtual void at_generation_tasks(const Monomial & t);
  /**
    @brief re-evaluates the sugar, if need be
    @param u monomial to multiply to @p g, then reduce @p this by the product
    @param g a polynomial whose leading monomial divides @p this
  */
  virtual void pre_reduction_tasks(
    const EXP_TYPE * u, const Abstract_Polynomial & g
  );
  /** @name Basic properties */
  ///@{
  /** @brief type of strategy */
  StrategyFlags type() { return WSUGAR_STRATEGY; }
  /** @brief the sugar itself */
  DEG_TYPE poly_sugar() const;
  /** @brief the weights used to compute the sugar */
  const WT_TYPE * weighting() const;
  ///@}
  /** @name Modification */
  ///@{
  /** @brief changes the weights used to compute the sugar to \c w */
  void change_weights(const WT_TYPE * w) { weights = w; }
  ///@}
  /** @name I/O */
  ///@{
  friend ostream & operator <<(ostream &, const Poly_WSugar_Data &);
  ///@}
protected:
  /** @brief the polynomial&rsquo;s sugar */
  DEG_TYPE sugar;
  /** @brief pointer to the weights */
  const WT_TYPE * weights;
};

/**
  @ingroup strategygroup
  @class Pair_WSugar_Strategy
  @brief ordering critical pairs using the weighted sugar strategy
  @author John Perry
  @date 2016
*/
class Pair_WSugar_Strategy : public Pair_Strategy_Data {
public:
  /** @name Construction */
  ///@{
  /**
    @brief creates a pair whose weighted sugar is the maximum of that
      of the products of the polynomials in the critical pair
  */
  Pair_WSugar_Strategy(Critical_Pair_Basic &);
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
  /** @brief the weighted sugar */
  DEG_TYPE pair_sugar();
  ///@}
  /**
    @name Modification
    @brief Useful for when the ordering changes.
  */
  ///@{
  virtual void adjust_sugar(DEG_TYPE new_sugar) { sugar = new_sugar; }
  ///@}
protected:
  /** @brief the critical pair to which this \c Normal_Strategy belongs */
  Critical_Pair_Basic * cp;
  /** @brief the pair*rsquo;s sugar */
  DEG_TYPE sugar;
};

#endif