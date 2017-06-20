#ifndef __POLYNOMIAL_RING_H_
#define __POLYNOMIAL_RING_H_

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

/**
  @class Polynomial_Ring
  @author John Perry
  @date 2015
  @brief Encapsulates information about a polynomial ring for easy access:
    ground field, number of indeterminates, &hellip;
  @ingroup polygroup
*/

#include <string>

using std::string;

#include "system_constants.hpp"

#include "fields.hpp"

#include "indeterminate.hpp"

class Indeterminate;

class Polynomial_Ring {
public:
  /** @name Construction */
  ///@{
  /**
    @brief Initialize the ring for the given field and number of indeterminates
    @details The array of names is copied,
      so you can reuse or discard the memory as you prefer.
    @param num_vars number of indeterminates in this Polynomial_Ring
    @param field ground field of this Polynomial_Ring
    @param var_names an optional array of @p num_vars names to assign to the
      indeterminates
  */
  Polynomial_Ring(
      NVAR_TYPE num_vars, Prime_Field & field, string * var_names = nullptr
  );
  ///@}
  /** @name Destruction */
  ///@{
  /** @brief Deletes the names. */
  ~Polynomial_Ring();
  ///@}
  /** @name Modification */
  ///@{
  /**
    @brief sets the names of the indeterminates, if you do not want the default
    @details Returns @c True if and only if successful.
      The names are not copied, so please do not discard them until you are done
      with this ring, or else have reassigned the names.
    @warning If @c length is not the same as the number of indeterminates,
      things can go very, very badly.
    @param new_names names for the indeterminates
    @param length number of names in @p new_names
    @return whether all supplied names were set
  */
  bool set_names(string * new_names, NVAR_TYPE length);
  ///@}
  /** @name Basic properties */
  ///@{
  /** @brief number of indeterminates (variables) in the ring */
  virtual NVAR_TYPE number_of_variables() const;
  /** @brief an array of the ring&rsquo;s indeterminates; use only to biuld polynomials */
  virtual Indeterminate * indeterminates();
  /** @brief ground field */
  virtual Prime_Field & ground_field() const;
  /** @brief name of the @f$i@f$th indeterminate */
  virtual const string name(NVAR_TYPE i) const;
  /** @brief names of all the variabiles */
  virtual const string * name_list() const { return names; }
  ///@}
protected:
  /** @brief the ring's ground field */
  Prime_Field & F;
  /**
    @brief the number of indeterminates for every monomial in the ring
    @warning The system does not verify that the number of indeterminates
      in a monomial agrees with the number specified in this ring.
      This information is easily available through the number_of_variables()
      command.
  */
  NVAR_TYPE n;
  /** @brief optional names for the variables */
  string * names;
};

#endif
