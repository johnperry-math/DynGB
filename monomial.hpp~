#ifndef __MONOMIAL_H_
#define __MONOMIAL_H_

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

#include <string>
#include <cstdlib>
#include <iostream>
#include <initializer_list>

#include "system_constants.hpp"

#include "goda.hpp"
#include "monomial_ordering.hpp"

extern Monomial_Ordering * generic_grevlex_ptr;

using std::ostream;
using std::cout;
using std::initializer_list;
using std::endl;
using std::string;

#define EPACK false /** whether to pack exponents */
#define PACKSIZE uint8_t /** size of packed exponent -- do not change unless you are willing to change MASKALL in <c>.cpp</c> file (and maybe more) */
#define NPACK sizeof(unsigned long long)/sizeof(PACKSIZE)    /** number of exponents to pack */

/**
  @class Monomial
  @author John Perry
  @date 2015
  @brief Implementation of monomials.
  @ingroup polygroup
  
  @details A monomial is a product of powers, for instance @f$x_1^2x_3@f$.
  (We do not include coefficients in our definition of monomials.)
  This class encapsulates basic monomial functionality necessary
  for computing Gr&ouml;bner bases.

  @warning When initializing an array of @c Monomial , be sure to use the
  initialization functions to format and allocate space.
  
  @warning Any monomial requires that you specify a number of variables.
  At the present time, the limit on this is 64 (for the sake of masking;
  though an earlier verison of the code did not have this limit).
  Most operations require that the number of variables be the same,
  but the code does not check for this.
  For what I am designing this for, this is acceptable,
  but if you intend to work with monomials
  that contain a different number of variables,
  check to make sure that any monomial arithmetic involves monomials
  with the same number of variables. 
*/
class Monomial {
public:
  /** @name Construction */
  ///@{
  /** @brief things all @c Monomial initializers must do */
  inline void common_initialization(const Monomial_Ordering * ord = nullptr) {
    exponents = nullptr; ordering_data = nullptr; ordering = ord;
  }
  /** @brief clears ordering data */
  inline void clear_ordering_data() { ordering_data = nullptr; }
  /**
    @brief allocates memory for exponents
    This is useful when you want to allocate an array of monomials.
  */
  void initialize_exponents(NVAR_TYPE number_of_vars);
  /**
    @brief deallocates memory for exponents
    This is useful when you want to deallocate an array of monomials.
  */
  void deinitialize();
  /** @brief change @f$i@f$th exponent to @f$e@f$ */
  void set_exponent(NVAR_TYPE i, DEG_TYPE e);
  /** @brief The first constructor is equivalent to instantiating 1. */
  Monomial(
      NVAR_TYPE number_of_vars,
      const Monomial_Ordering * order = generic_grevlex_ptr
  );
  /** @brief Copy constructor, allocates new space for exponents. */
  Monomial(const Monomial &other);
  /**
    @brief Constructor from initializer list.
    @param powers the Monomial&rsquo;s exponents
    @param order a Monomial_Ordering in the appropriate number of variables
    @details The list should contains the powers of the exponents, in order.
  */
  Monomial(
      initializer_list<EXP_TYPE> powers,
      const Monomial_Ordering * order = generic_grevlex_ptr
  );
  /**
    @brief Constructor from array of exponents.
    @param size number of variables in the Monomial
    @param powers array of exponents &mdash; use at least @p n unless you like
        to crash
    @param order a Monomial_Ordering for @p n variables
    @details Copies exponents so you can reuse yours.
  */
  Monomial(
      NVAR_TYPE size, const EXP_TYPE *powers,
      const Monomial_Ordering * order = generic_grevlex_ptr
  );
  ///@}
  /** @name Destruction */
  ///@{
  /** @brief Destructor deallocates exponents. */
  ~Monomial();
  ///@}
  /** @name Basic properties
    The following functions give information about the monomial,
    but do not modify it.
  */
  ///@{
  /** @brief number of variables */
  inline NVAR_TYPE num_vars() const { return n; }
  /** @brief all exponents 0? */
  bool is_one() const;
  /** @brief Degree of @f$i@f$th variable. */
  DEG_TYPE degree(NVAR_TYPE i) const;
  /** @brief Degree of @f$i@f$th variable. Synonymous with degree().*/
  DEG_TYPE operator[](NVAR_TYPE i) const;
  /**
    @brief Sum of exponents of the first @p m variables.
    @details If @c m is zero (default), computes for all variables.
    @param m compute up to @p m indeterminates; the default value of 0 computes
        all
    @return total standard degree of this Monomial
  */
  DEG_TYPE total_degree(NVAR_TYPE m=0) const;
  /**
    @return weighted sum of first @p m exponents, using given @p weights
    @param weights the weights for each indeterminate: the first weight will
        be multiplied by the first exponent, the second weight by the second
        exponent, and so forth
    @param m compute up to @p m indeterminates; the default value of 0 computes
        all
    @details If @p weights is @c nullptr then returns total_degree(m).
        If @c m is zero (default), computes for all indeterminates.
  */
  DEG_TYPE weighted_degree(const WT_TYPE *weights, NVAR_TYPE m=0) const;
  /** @brief Direct access to the exponents, for whatever reason. */
  const EXP_TYPE * log() const { return exponents; }
  ///@}
  /**
    @name Comparison
    Compares this monomial to another.
  */
  ///@{
  /** @brief the Monomial_Ordering associated with this Monomial */
  inline const Monomial_Ordering * monomial_ordering() const { return ordering; }
  /** @brief the Monomial_Order_Data associated with this Monomial */
  inline Monomial_Order_Data * monomial_ordering_data() const {
    return ordering_data;
  }
  /** @brief sets the Monomial_Ordering associated with this Monomial */
  void set_monomial_ordering(const Monomial_Ordering * mord);
  /** @brief sets the Monomial_Order_Data associated with this Monomial */
  void set_ordering_data(Monomial_Order_Data * mordat);
  /** @brief sets the ordering degree for this Monomial */
  inline void set_ordering_degree(DEG_TYPE d) { ord_degree = d; }
  /** @brief returns the ordering degree for this Monomial */
  inline DEG_TYPE ordering_degree() const { return ord_degree; }
  /** @brief equal/alike? */
  bool operator ==(const Monomial &other) const;
  /** @brief unequal/unlike? */
  bool operator !=(const Monomial &other) const;
  /**
    @brief <c>true</c> iff <c>this</c> has no common factor with <c>other</c>.
  */
  bool is_coprime(const Monomial &other) const;
  /** @brief Returns 0 if like, negative if smaller */
  int cmp(const Monomial & u) const {
    return monomial_ordering()->cmp(*this, u);
  }
  /** @brief Have same variables, same powers? Synonymous with operator==(). */
  bool is_like(const Monomial &other) const;
  /** @brief Have same variables, same powers? Comparison of exponent vectors. */
  bool is_like(const EXP_TYPE * e) const {
    bool result = true;
    for (unsigned i = 0; result and i < n; ++i) result = exponents[i] == e[i];
    return result;
  }
  /** @brief is <c>this</c> like @f$uv@f$? */
  bool like_multiple(const Monomial &u, const Monomial &v) const;
  /**
    @brief is <c>this</c> like @f$uv@f$?
      (where @f$u@f$ is a monomial whose exponents are those of <c>e</c>
  */
  bool like_multiple(EXP_TYPE * e, const Monomial & v) const;
  /**
    @brief compares monomial with @f$u@f$ according to monomial ordering
  */
  bool larger_than(const Monomial &) const;
  /**
    @brief Compares monomial with @f$u@f$ according to monomial ordering.
  */
  bool operator >(const Monomial & u) const;
  /**
    @brief Compares monomial with @f$u@f$ according to monomial ordering.
  */
  bool operator >=(const Monomial & u) const;
  /**
    @brief Compares monomial with @f$u@f$ according to monomial ordering.
  */
  bool operator <(const Monomial & u) const;
  /**
    @brief Compares monomial with @f$u@f$ according to monomial ordering.
  */
  bool operator <=(const Monomial & u) const;
  /**
    @brief Compares monomial with @f$uv@f$ according to monomial ordering.
  */
  bool larger_than_multiple(
      const Monomial & u, const Monomial &v
  ) const;
  /** @brief Divisible by <c>other</c>? */
  bool divisible_by(const Monomial &other) const;
  /** @brief operator for divisibility */
  bool operator |(const Monomial &other) const;
  ///@}
  /**
    @name Computation
    Computes something, and may modify <c>this</c>.
  */
  ///@{
  /** @brief assignment */
  Monomial & operator =(const Monomial &other);
  /** @brief Multiply <c>this</c> by <c>other</c>. */
  Monomial & operator *=(const Monomial &other);
  /** @brief Return result of <c>this</c> by <c>other</c>. */
  Monomial operator *(const Monomial &other) const;
  /**
    @brief divide @c this by @p other
    @param other Monomial to divide @c this by
    @warning Does not check if <c>this</c> is divisible by other.
      See divisible_by(), which is the tool to use.
    @return @c true if and only if @c this and @p other have the same number
      of variables
  */
  bool operator /=(const Monomial &other);
  /** @brief Least common multiple: largest exponents. */
  Monomial lcm(const Monomial & u) const;
  /** @brief Greatest common divisor: smallest exponents. */
  Monomial gcd(const Monomial & u) const;
  /**
    @brief colon operator: exponents needed to make @f$u@f$ divisible by @c this
  */
  Monomial colon(const Monomial & u) const;
  ///@}
  /** @name I/O */
  ///@{
  /**
    @brief prints the monomial with the given names; adds
      a newline if the boolean is true
  */
  void print(bool=true, ostream & = cout, const string * names = nullptr) const;
  /**
    @brief equivalent to @c print() with default values
  */
  void printlncout() const { print(); }
  /**
    @brief essentially <c>u.print(false, ostream)</c>
  */
  friend ostream & operator << (ostream &, const Monomial & u);
  ///@}
  /** @name Memory management */
  ///@{
  /** @brief requests memory form Monomial's Grading_Order_Data_Allocator */
  void * operator new(size_t);
  /** @brief returns data to Monomial's Grading_Order_Data_Allocator */
  void operator delete(void *);
  ///@}
protected:
  /** @brief number of variables */
  NVAR_TYPE n;
  /** @brief has size n */
  EXP_TYPE * exponents;
  /**
    @brief degree associated with monomial ordering (used for faster comparisons)
  */
  DEG_TYPE ord_degree;
  /** @brief Monomial_Ordering associated with this polynomial */
  const Monomial_Ordering * ordering;
  /** @brief optional data for a monomial ordering */
  Monomial_Order_Data * ordering_data;
  /** @brief divisibility mask (up to 64 variables)  */
  uint64_t mask;
  #if EPACK
  /** @brief equality mask, a packing of exponents */
  uint64_t emask;
  #endif
};

#endif