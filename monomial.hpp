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
#include <bitset>
using std::bitset;
#include <cstdlib>
#include <iostream>
#include <initializer_list>

#include "system_constants.hpp"

#include "goda.hpp"
#include "indeterminate.hpp"
#include "monomial_ordering.hpp"

class Indeterminate;

extern Monomial_Ordering * generic_grevlex_ptr;

using std::ostream;
using std::cout;
using std::initializer_list;
using std::endl;
using std::string;

#define EPACK false /** whether to pack exponents */
#define PACKSIZE uint8_t /** size of packed exponent -- do not change unless you are willing to change MASKALL in <c>.cpp</c> file (and maybe more) */
#define NPACK sizeof(unsigned long long)/sizeof(PACKSIZE)    /** number of exponents to pack */

/** @brief used to locate a variable's exponent in the Monomial data structure */
struct Exponent_Location_Struct {
  /**
    @brief location of the desired exponent, if already set
    (check against @c already_set !)'
    @warning If @c already_set is @c false, this data is @b invalid !
  */
  NVAR_TYPE loc = 0;
  /**
    @brief @c true if and only if the degree of this exponent has been set
    @warning Ignore the value of @c loc if @c already_set is @c false!
  */
  bool already_set;
};

typedef struct Exponent_Location_Struct Exponent_Location;

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
  /**
    @brief "product" constructor when monomial's expoennts already initialized
    @param product @c true if and only if we want a product
    @param t @c Monomial to multiply
    @param u @c Monomial to multiply
  */
  void make_product_or_quotient(
      const Monomial & t, const Monomial & u, bool product = true
  );
  /**
    @brief "Product" constructor: constructs a product or quotient
    @param product @c true if and only if we want a product
    @param t @c Monomial to multiply
    @param u @c Monomial to multiply
  */
  Monomial(const Monomial & t, const Monomial & u, bool product = true);
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
  /**
    @brief access to the exponents, for whatever reason.
    @return an exponent vector @c e where <c>e[i]</c> is
      the degree of the monomial in the <c>i</c>th variable
    @details This is a new array.
      Clients should delete it once they are done using it.
  */
  EXP_TYPE * log() const {
    EXP_TYPE * result = new EXP_TYPE[n] { 0 };
    for (NVAR_TYPE i = 0; i < last; i += 2)
      result[exponents[i]] = exponents[i + 1];
    return result;
  }
  /**
    @brief access to the exponents, packed; use in conjunction with @c packed_size()
    @return array of exponents in the form @f$(i_1,e_1,\ldots,i_k,e_k)@f$
      such that @f$ x_{i_j}^e_{i_j} @f$ divides @c this
    @warning In the current implementation, this provides direct access to the
      underlying data structure. <b>Do not modify or delete.</b>
  */
  const EXP_TYPE * packed_log() const {
    return exponents;
  }
  /** @brief indicates the size of the result of @c packed_log() */
  NVAR_TYPE packed_size() const { return last; }
  /** @brief exponent mask */
  const bitset<MASK_SIZE> mask() const {
    bitset<MASK_SIZE> result;
    for (NVAR_TYPE i = 0; i < last; i += 2)
      result[exponents[i]] = 1;
    return result;
  }
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
  bool like_multiple(const EXP_TYPE * e, const Monomial & v) const;
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
  /** @brief remove zero exponents */
  void compress();
  /** @brief assignment */
  Monomial & operator =(const Monomial &other);
  /** @brief Multiply <c>this</c> by <c>other</c>. */
  Monomial & operator *=(const Monomial &other);
  /** @brief Return result of <c>this</c> by <c>other</c>. */
  Monomial operator *(const Monomial &other) const;
  /** @brief Return result of <c>this</c> by <c>other</c>. */
  Monomial operator *(const Indeterminate &other) const;
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
  /**
    @brief locate the indicated exponent
    @return information on the exponent's location; @see Exponent_Location_Struct
  */
  Exponent_Location find_exponent(NVAR_TYPE) const;
  /** @brief make the last variable invalid */
  void make_last_invalid() { last = 0; }
  /**
    @brief check if exponents are valid
    @return @c true if and only if this monomial has valid data
    @details Currently this is a test to see if the polynomial has any exponents
    defined at all. It will in fact give the opposite result of @c is_one().
  */
  bool valid_exponents() const { return last > 0; }
  /** @brief make space for a new exponent, starting at given position */
  void shift_exponents_right(NVAR_TYPE);
  /** @brief number of variables */
  NVAR_TYPE n;
  /** @brief has size n */
  EXP_TYPE * exponents;
  /** @brief last valid exponent */
  NVAR_TYPE last = 0;
  /**
    @brief degree associated with monomial ordering (used for faster comparisons)
  */
  DEG_TYPE ord_degree;
  /** @brief Monomial_Ordering associated with this polynomial */
  const Monomial_Ordering * ordering;
  /** @brief optional data for a monomial ordering */
  Monomial_Order_Data * ordering_data;
};

#endif