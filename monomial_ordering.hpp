#ifndef __MONOMIAL_ORDERING_HPP_
#define __MONOMIAL_ORDERING_HPP_

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

#include <cstdlib>

#include "system_constants.hpp"

class Monomial;

/**
  @defgroup orderinggroup Monomial Orderings
  @brief classes defining monomial orderings
*/

/**
  @ingroup orderinggroup
  @class Monomial_Order_Data
  @author John Perry
  @date 2015
  @brief data for a monomial ordering: optional, but stored in @c Monomial
*/
class Monomial_Order_Data {
public:
  /** @name Construction */
  ///@{
  /** @brief default clone returns @c nullptr */
  virtual Monomial_Order_Data * clone() { return nullptr; }
  ///@}
  /** @name Basic data */
  ///@{
  /** @brief default value is useless; orderings that supply gradings should redefine */
  virtual DEG_TYPE grading(NVAR_TYPE) const { return 0; }
  ///@}
  /** @name Destruction */
  ///@{
  /** @brief does nothing but guarantee polymorphism (stupid, stupid C++) */
  virtual ~Monomial_Order_Data() {}
  ///@}
};

/**
  @ingroup orderinggroup
  @class Monomial_Ordering
  @author John Perry
  @date 2015
  @brief interface to a monomial ordering
  @warning Avoid changing the monomial ordering data in the comparison functions,
  as the user may wish to compare
  two monomials according to a different ordering than the current value.
  This is the reason those functions are marked to leave the Monomials constant.
  The expected behavior is that first_larger() does whatever the ordering
  wants, so if you need monomial data check first to decide whether it exists,
  and applies to this ordering!
  Use set_data() if you want to change the monomials first.

  @warning There is no strict need to check whether the monomials
  are associated to the same ordering, but if the ordering uses
  Monomial_Order_Data one would be foolhardy not to check first.
*/
class Monomial_Ordering {
public:
  /** @name Destruction */
  ///@{
  /** @brief needs virtual destructor for polymorphic @c delete */
  virtual ~Monomial_Ordering();
  ///@}
  /** @name Utility */
  ///@{
  /**
    @brief sets monomial ordering&rsquo;s data; default is to do nothing
    @details Child classes that override this function are strongly recommended
      to use set_ordering_degree() of the Monomial class to set a primary
      degree. For weighted/graded degree orderings, this typically improves
      performance nontrivially.
  */
  virtual void set_data(Monomial &) const;
  ///@}
  /** @name Comparison */
  ///{@
  /**
    @param t a Monomial to compare to @f$ u @f$
    @param u a Monomial to compare to @f$ t @f$
    @return 0 if the Monomials are like; negative if smaller, positive
      if larger -- for efficiency, you probably want to redefine this
  */
  virtual int cmp(const Monomial & t, const Monomial & u) const = 0;/*{
    int result = 0;
    if (first_larger(t, u)) result = 1;
    else if (first_larger(u, t)) result = -1;
    else result = 0;
    return result;
  }*/
  /**
    @return @c true iff the first Monomial is larger than the second
  */
  virtual bool first_larger(const Monomial &, const Monomial &) const = 0;
  /**
    @return @c true iff the first Monomial is larger or equal to the second
    @see first_larger()
    @warning Do not override unless you know what you&rsquo;re doing
  */
  bool first_larger_or_equal(const Monomial &, const Monomial &) const;
  /**
    @return @c true iff the first Monomial is smaller than the second
    @warning Do not override unless you know what you&rsquo;re doing
  */
  virtual bool first_smaller(const Monomial &, const Monomial &) const = 0;
  /**
    @return @c true iff the first Monomial is smaller or equal
      to the second
    @see first_larger()
    @warning Do not override unless you know what you&rsquo;re doing
  */
  bool first_smaller_or_equal(const Monomial &, const Monomial &) const;
  /**
    @brief returns @c true iff the first Monomial is larger than the specified
      multiple of the second
  */
  virtual bool first_larger_than_multiple(
      const Monomial &, const Monomial &, const Monomial &
  ) const = 0;
  /**
    @return @c true iff the first Monomial is larger or equal to 
      the specified multiple of the second
    @see first_larger()
    @warning Do not override unless you know what you&rsquo;re doing
  */
  bool first_larger_or_equal_than_multiple(
      const Monomial &, const Monomial &, const Monomial &
  ) const;
  /**
    @return @c true iff the first Monomial is smaller than the specified
      multiple of the second
    @warning Do not override unless you know what you&rsquo;re doing
  */
  bool first_smaller_than_multiple(
      const Monomial &, const Monomial &, const Monomial &
  ) const;
  /**
    @return @c true iff the first Monomial is smaller or equal to
      the specified multiple of the second
    @see first_larger()
    @warning Do not override unless you know what you&rsquo;re doing
  */
  bool first_smaller_or_equal_than_multiple(
      const Monomial &, const Monomial &, const Monomial &
  ) const;
  ///@}
};

/**
  @ingroup orderinggroup
  @class Weighted_Ordering
  @author John Perry
  @date 2016
  @brief interface to a weighted monomial ordering
  @see Monomial_Ordering
  @details This class adds all of one method to Monomial_Ordering.
*/
class Weighted_Ordering : public Monomial_Ordering {
public:
  /** @name Basic properties */
  ///@{
  /** @brief returns the weights used by this orderings */
  virtual const WT_TYPE * order_weights() const = 0;
  /** @brief returns the number of weights (same as number of indeterminates) */
  virtual NVAR_TYPE number_of_weights() const = 0;
  ///@}
};


#endif