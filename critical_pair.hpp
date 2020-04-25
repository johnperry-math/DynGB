#ifndef __CRITICAL_PAIR_HPP_
#define __CRITICAL_PAIR_HPP_

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
#include "polynomial_geobucket.hpp"
#include "polynomial_linked_list.hpp"
#include "polynomial_double_buffered.hpp"

#include "strategies.hpp"
#include "sugar_strategy.hpp"
#include "normal_strategy.hpp"
#include "weighted_sugar_strategy.hpp"

#include "particular_orderings.hpp"

// forward declaration
class Pair_Sugar_Data;
class Poly_Sugar_Data;

/**
  @defgroup strategygroup Strategies
  @brief strategies for Gr&ouml;bner basis computation
*/

/**
  \enum SPolyCreationFlags
  @author John Perry
  @date 2015
  @brief flag indicating which structure to use for an s-polynomial
  @ingroup GBComputation
  @details The following options are currently available:
    - LINKED_LST polynomials are created and manipulated in linked list format
    - GEOBUCKETS polynomials are created and manipulated in geobucket format
      @cite YanGeobuckets
    - DOUBLE_BUF polynomials are created and manipulated in a novel,
      doubled-buffered format
*/
enum class SPolyCreationFlags {
  MIN_SPCREATE_FLAG = 0,
  LINKED_LST,
  GEOBUCKETS,
  DOUBLE_BUF,
  MAX_SPCREATE_FLAG
};

/**
  @class Critical_Pair_Basic
  @author John Perry
  @date 2015
  @brief Controls the creation of s-polynomials.
  @ingroup GBComputation
  
  @details This class encapsulates the information necessary to create
  an @f$S@f$-polynomial, including the actual generation of
  an @f$S@f$-polynomial, though it does not retain that information.
*/
class Critical_Pair_Basic {
public:
  /** @name Construction */
  ///@{
  /** @brief create critical pair (f,0) for initial polynomial */
  Critical_Pair_Basic(Abstract_Polynomial * f, StrategyFlags strategy);
  /** @brief create critical pair (f,g) for two polynomials */
  Critical_Pair_Basic(
      Abstract_Polynomial * f, Abstract_Polynomial * g, StrategyFlags strategy
  );
  /** @brief copy constructor */
  explicit Critical_Pair_Basic(const Critical_Pair_Basic &);
  ///@}
  /** @name Destruction */
  ///@{
  virtual ~Critical_Pair_Basic() { if (key != nullptr) delete key; }
  ///@}
  /** @name Basic properties */
  ///@{
  /** @brief whether this pair is from a generator */
  inline bool is_generator() const { return q == nullptr; }
  /** @brief first polynomial in pair */
  inline const Abstract_Polynomial * first() const { return p; }
  /** @brief second polynomial in pair */
  inline const Abstract_Polynomial * second() const { return q; }
  /** @brief lcm of leading monomials of polynomials */
  inline const Monomial & lcm() const { return tpq; }
  /** @brief degree of ith variable in lcm */
  inline unsigned lcm_degree(unsigned i) const { return tpq.degree(i); }
  /** @brief monomial needed to multiply first polynomial to lcm() */
  inline const Monomial & first_multiplier() const { return tp; }
  /** @brief monomial needed to multiply second polynomial to lcm() */
  inline const Monomial & second_multiplier() const { return tq; }
  /** @brief strategy used for comparison of pairs */
  inline const Pair_Strategy_Data * pair_key() const { return key; }
  /**
    @brief to use only if s-polynomial is already computed by another method
    @warning If you have not already created an s-polynomial using one of
        SPolyCreationFlags, then this returns @c nullptr and is useless.
    @return the s-polynomial of this pair
  */
  virtual Mutable_Polynomial * s_polynomial() { return s; }
  ///@}
  /** @name Computation */
  ///@{
  /** @brief creates the s-polynomial for first() and second() */
  virtual Mutable_Polynomial * s_polynomial(
      SPolyCreationFlags method, StrategyFlags strategy
  ) ;
  ///@}
  /** @name Modification */
  ///@{
  /** @brief sets the s-polynomial to @c p, for explorer methods*/
  virtual void set_spoly(Mutable_Polynomial * p) { s = p; }
  /** @brief in case you want to swap the polynomials, for whatever reason */
  void swap() {
    auto temp = p;
    p = q;
    q = temp;
    auto temp_t = tp;
    tp = tq;
    tq = temp_t;
  }
  ///@}
  /** @name I/O */
  ///@{
  friend ostream & operator<<(ostream &, const Critical_Pair_Basic &);
  ///@}
protected:
  /** @brief first polynomial in the critical pair */
  Abstract_Polynomial * p;
  /** @brief second polynomial in the critical pair */
  Abstract_Polynomial * q;
  /** @brief S-polynomial */
  Mutable_Polynomial * s;
  /** @brief lcm of leading monomials of @f$p@f$ and @f$q@f$ */
  Monomial tpq;
  /** @brief monomial multiple of @f$p@f$ that produces @f$S@f$-polynomial */
  Monomial tp;
  /** @brief monomial multiple of @f$q@f$ that produces @f$S@f$-polynomial */
  Monomial tq;
  /** @brief strategy used to sort critical pairs */
  Pair_Strategy_Data * key = nullptr;
};

/**
  @class Critical_Pair_Dynamic
  @author John Perry
  @date 2016
  @brief Controls the creation of s-polynomials, specialized for the dynamic
      algorithm.
  @ingroup GBComputation
  
  @details This class encapsulates the information necessary to create
  an @f$S@f$-polynomial, including the actual generation of
  an @f$S@f$-polynomial, though it does not retain that information.
  The main difference with @c Critical_Pair_Basic is the enabling of late
  re-sorting.
*/
class Critical_Pair_Dynamic : public Critical_Pair_Basic {
public:
  /** @name Construction */
  ///@{
  /**
    @brief create critical pair (f,0) for initial polynomial,
    with given ordering
  */
  Critical_Pair_Dynamic(
      Abstract_Polynomial * f, StrategyFlags strategy,
      Weighted_Ordering * how_to_order
  ) : Critical_Pair_Basic(f, strategy) {
    ordering = how_to_order;
  };
  /** @brief copy constructor */
  Critical_Pair_Dynamic(const Critical_Pair_Dynamic & other) :
      Critical_Pair_Basic(other), ordering(other.ordering)
  { /* done */ }
  /**
    @brief create critical pair (f,g) for two polynomials, with given ordering
  */
  Critical_Pair_Dynamic(
      Abstract_Polynomial * f, Abstract_Polynomial * g, StrategyFlags strategy,
      Weighted_Ordering * how_to_order
  ) : Critical_Pair_Basic(f, g, strategy) {
    ordering = how_to_order;
    tpq.set_monomial_ordering(ordering);
    tp.set_monomial_ordering(ordering);
    tq.set_monomial_ordering(ordering);
    delete key;
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
  ///@}
  // No need for a new destructor (I hope)
  /** @name Computation */
  ///@{
  /**
    @brief creates the s-polynomial for first() and second()
    @details If either first() or second() was sorted using an ordering
        different from the one assigned the pair, it is sorted anew.
    @param method a flag for what structure to use while reducing the
        s-polynomial; see SPolyCreationFlags
    @param strategy a flag for which strategy to use in reduction; see
        StrategyFlags
    @return the generated s-polynomial
  */
  virtual Mutable_Polynomial * s_polynomial(
      SPolyCreationFlags method, StrategyFlags strategy
  );
  ///@}
  /** @name Basic properties */
  ///@{
  /** @brief the ordering associated with this pair */
  Weighted_Ordering * how_ordered() const { return ordering; }
  ///@}
  /** @name Modification */
  ///@{
  /**
    @brief change the ordering associated with this pair
    @details This is necessary at least for critical pairs where one polynomial
        is a generator that we have not yet computed.
        We delay re-sorting the polynomial until the s-polynomial&rsquo;s
        computation.
    @param new_order the new ordering for this critical pair
  */
  void change_ordering(Weighted_Ordering * new_order) {
    ordering = new_order;
    // the leading term may change when an input polynomial is involved,
    // so we need to reconsider the leading monomial & "lcm" here
    if (is_generator()) {
      p->set_monomial_ordering(ordering);
      tpq = p->leading_monomial();
    }
    tpq.set_monomial_ordering(ordering);
    tp.set_monomial_ordering(ordering);
    tq.set_monomial_ordering(ordering);
    Pair_WSugar_Strategy * strat = dynamic_cast<Pair_WSugar_Strategy *>(key);
    if (strat != nullptr) {
      DEG_TYPE new_sugar = first()->weighted_degree(new_order->order_weights());
      new_sugar += first_multiplier().weighted_degree(new_order->order_weights());
      if (second() != nullptr) {
        DEG_TYPE second_sugar = first()->weighted_degree(new_order->order_weights());
        second_sugar += second_multiplier().weighted_degree(new_order->order_weights());
        if (second_sugar > new_sugar) new_sugar = second_sugar;
      }
      strat->adjust_sugar(new_sugar);
    }
  }
  ///@}
  /** @name I/O */
  ///@{
  /** @brief outputs the pair, and the ordering */
  friend ostream & operator<<(ostream &, const Critical_Pair_Dynamic &);
  ///@}
protected:
  /**
    @brief the @c Monomial_Ordering assigned to this pair &mdash; might disagree
      with that of polynomials; they must be updated immediately if this is a
      generator; otherwise, update can wait until the critical pair is computed.
  */
  Weighted_Ordering * ordering;
};

#endif