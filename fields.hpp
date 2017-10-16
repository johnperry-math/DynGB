#ifndef __FIELDS_H_
#define __FIELDS_H_

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

#include <vector>
using std::vector;
#include <iostream>
using std::ostream;

#include "system_constants.hpp"

/**
  @defgroup FieldsGroup Fields
  @brief classes for performing field arithmetic
*/

// forward declaration
class Prime_Field_Element;

/**
  @class Prime_Field
  @author John Perry
  @date 2015
  @brief Information necessary for a field modulo a prime.
  @ingroup FieldsGroup
  
  @details This class encapsulates the information necessary for a field modulo a prime:
  both the modulus @f$m@f$ and a list of inverses of non-zero elements.
  The constructors do not verify that @f$m@f$ is prime,
  but expect misbehavior if not.

  @warning The constructor does not check whether the supplied modulus is prime.
  If @f$m@f$ is not prime, behavior is undefined, and probably spectacularly bad.
*/
class Prime_Field
{
public:
  /** @name Construction */
  ///@{
  /**
    @brief Please read detailed information.
    @param modulus All computation in this field is done modulo this number;
      it should be prime, but we do not check this. See this warning.
    @param show_modulus indicates whether the modulus is shown when printing an
      element of this field. You probably don't want this.

    @warning This constructor does not check whether @p modulus is prime.
    If @f$m@f$ is not prime, behavior is undefined,
    and probably spectacularly bad.
  */
  Prime_Field(UCOEF_TYPE modulus, bool show_modulus = false);
  /**
    @brief Copy constructor: copies previously-computed inverses.
    @param F the source to copy

    Allocates new data, so you can discard <c>F</c> later if you want.
  */
  explicit Prime_Field(const Prime_Field & F);
  ///@}
  /** @name Destruction */
  ///@{
  ~Prime_Field();
  ///@}
  /**
    @name Basic properties
    The following functions give information about the prime field,
    but do not modify it.
  */
  ///@{
  /** @brief Returns the field's modulus. */
  inline unsigned modulus() const { return m; }
  /**
    @brief Returns the inverse of @f$a@f$, modulo @f$m@f$.
    @param a the element of this field whose inverse we desire
    @return @f$ a^{-1} @f$
    @details Looks up the inverse in a table.
  */
  COEF_TYPE inverse(COEF_TYPE a) const { return Fi[a]; };
  /** @brief &ldquo;unity&rdquo; is the multiplicative identity. */
  Prime_Field_Element unity() const;
  /** @brief &ldquo;zero&rdquo; is the additive identity. */
  Prime_Field_Element zero() const;
  ///@}
  /** \name I/O */
  ///@{
  /** @brief determines whether to print the modulus */
  void set_print_modulus(bool b);
  /** @brief indicates whether to print the modulus */
  bool get_print_modulus() const;
  ///@}
  /** @name Modification */
  Prime_Field & operator = (const Prime_Field & other) {
    m = other.m;
    print_modulus = other.print_modulus;
    Fi = other.Fi;
    return *this;
  }
  ///@{
  ///@}
protected:
  /** @brief characteristic/modulus of the field */
  UCOEF_TYPE m;
  /**
    @brief for @f$i\neq0@f$, @f$Fi_i@f$ is multiplicative inverse
      of @f$i@f$, mod @f$m@f$
  */
  vector<COEF_TYPE> Fi;
  /**
    @brief determines whether a coefficient's modulus is printed
  */
  bool print_modulus;
};

/**
  @class Prime_Field_Element
  @author John Perry
  @date 2015
  @brief Element of a field of prime characteristic.
  @ingroup FieldsGroup

  @details This class encapsulates an element of a field of prime characteristic
  (Prime_Field).

  @warning Do not delete the prime field that <c>this</c> lives in while
  <c>this</c> is still active.
  Behavior is unpredictable in this circumstance,
  as the field is necessary for some record keeping, e.g.,
  to look up multiplicative inverses.
*/
class Prime_Field_Element
{
public:
  /** @name Construction */
  ///@{
  /**
    @brief Constructs a prime field element in the specified field,
    with the value 0.
    @param field this element&rsquo;s parent
    @details A pointer to <c>field</c> is attached to <c>this</c>
    in order to find inverses during arithmetic.
    @warning Do not delete this field while <c>this</c> is still active.
    Behavior is unpredictable in this circumstance.
  */
  explicit Prime_Field_Element(const Prime_Field *field);
  /**
    @brief Constructs a prime field element in the specified field,
    with the specified value.
    @param value the element&rsquo;s value
    @param field the element&rsquo;s base field
    @details A pointer to <c>field</c> is attached to <c>this</c>
    in order to find inverses during arithmetic.
    @warning Do not delete this field while <c>this</c> is still active.
    Behavior is unpredictable in this circumstance.
  */
  Prime_Field_Element(COEF_TYPE value, const Prime_Field *field);
  /** @brief constructs a prime field element from another */
  Prime_Field_Element(const Prime_Field_Element & other) {
    a = other.a; F = other.F; m = other.m;
  }
  ///@}
  /** @name Basic properties
    The following functions give information about the monomial,
    but do not modify it.
  */
  ///@{
  /** @brief The value of the element. This always satisfies @f$0\leq a\leq m@f$. */
  COEF_TYPE value() const;
  /** @brief The field's modulus. */
  unsigned modulus() const;
  /** @brief The field this element lies in. */
  const Prime_Field * field() const;
  /** @brief <c>true</c> iff <c>this</c> and <c>b</c> have the same modulus. */
  bool like(const Prime_Field_Element & b) const;
  /** @brief Returns the multiplicative inverse of <c>this</c>. */
  COEF_TYPE inverse() const;
  /** @brief Is <c>this</c> the additive identity? */
  inline bool is_zero() const { return a == 0; };
  /** @brief Is <c>this</c> the multiplicative identity? */
  bool is_one() const;
  ///@}
  /** @name Modification */
  ///@{
  /** @brief assignment to other */
  Prime_Field_Element & operator = (const Prime_Field_Element & other) {
    a = other.a; F = other.F; m = other.m;
    return *this;
  }
  /** @brief for initializing arrays of Prime_Field_Element */
  void assign(COEF_TYPE val, const Prime_Field * K) {
    a = val; F = K; m = F->modulus();
  }
  ///@}
  /**
    @name Computation
    Computes something, and may modify <c>this</c>.
  */
  ///@{
  /** @brief &ldquo;Negates&rdquo; <c>this</c>. */
  void negate();
  /**
    @brief increases the value of <c>this</c>.
    @param other element to add
  */
  void operator +=(const Prime_Field_Element &other);
  /**
    @brief decreases the value of <c>this</c>.
    @param other element to subtract
  */
  void operator -=(const Prime_Field_Element &other);
  /**
    @brief multiply the value of <c>this</c> by @p other
    @param other element to multiply
  */
  void operator *=(const Prime_Field_Element &other);
  /**
    @brief Changes the value of <c>this</c>.
    @param b scalar to multiply
    @details This function is useful for avoiding the construction
    of a prime field element when you know what you want to multiply.
  */
  void operator *=(const COEF_TYPE b);
  /**
    @brief Changes the value of <c>this</c>.
    @param other element to multiply
    @details Multiplies by multiplicative inverse of <c>other</c>. */
  void operator /=(const Prime_Field_Element &other);
  ///@}
  /** @name Friend functions for computation
    Will not modify <c>this</c>.
  */
  ///@{
  /** @brief Does not modify <c>this</c>. */
  friend Prime_Field_Element operator+(const Prime_Field_Element &,
                                        const Prime_Field_Element &);
  /** @brief Does not modify <c>this</c>. */
  friend Prime_Field_Element operator-(const Prime_Field_Element &,
                                        const Prime_Field_Element &);
  /** @brief Does not modify <c>this</c>. */
  friend Prime_Field_Element operator*(const Prime_Field_Element &,
                                        const Prime_Field_Element &);
  /** @brief Does not modify <c>this</c>. */
  friend Prime_Field_Element operator+(const Prime_Field_Element &, const int);
  /** @brief Does not modify <c>this</c>. */
  friend Prime_Field_Element operator-(const Prime_Field_Element &, const int);
  /** @brief Does not modify <c>this</c>. */
  friend Prime_Field_Element operator*(const Prime_Field_Element &, const int);
  /**
    @brief Does not modify <c>this</c>. If you want to modify <c>this</c>,
    see negate().
  */
  friend Prime_Field_Element operator-(const Prime_Field_Element &);
  ///@}
  /** @name I/O
    @brief Will not modify <c>this</c>.
  */
  ///@{
  friend ostream & operator << (ostream &, const Prime_Field_Element &);
  ///@}
protected:
  /** @brief the field this element lives in; used to find inverses */
  const Prime_Field *F;
  /**
    @brief the number&rsquo;s value; descendants should ensure @f$0\leq a<m@f$
  */
  COEF_TYPE a;
  /**
    @brief the number&rsquo; modulus, stored here to avoid the expense of
      accessing it in @f$F@f$.
  */
  unsigned m;
};

#endif