#ifndef __PARTICULAR_ORDERINGS_HPP_
#define __PARTICULAR_ORDERINGS_HPP_

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

#include <exception>

#include "system_constants.hpp"

#include "monomial.hpp"

#include "lp_solver.hpp"
using LP_Solvers::Ray;

/**
  @ingroup orderinggroup
  @class Generic_Grevlex
  @author John Perry
  @date 2015
  @brief generic grevlex ordering, works with any number of variables

  @details The difference between Generic_Grevlex and Grevlex_Ordering is that the
  former doesn&rsquo;t track any Monomial_Order_Data, while the latter relies
  on it.
  The latter should, as a result, be more time efficient, while the former is
  more space efficient. Generic_Grevlex will also order monomials which by
  design have a different number of variables, though that should not in
  general be something one encounters, and it&rsquo;s rather dangerous with
  first_larger_than_multiple().
*/
class Generic_Grevlex : public Monomial_Ordering {
public:
  /** @name Comparison */
  ///@{
  /** @brief returns @c true iff the first Monomial is larger than the second */
  virtual bool first_larger(const Monomial &, const Monomial &) const override;
  /** @brief returns @c true iff the first Monomial is smaller than the second */
  virtual bool first_smaller(const Monomial &, const Monomial &) const override;
  /**
    @brief returns @c true iff the first Monomial is larger than the product
      of the second and the third
  */
  virtual bool first_larger_than_multiple(
      const Monomial &, const Monomial &, const Monomial &
  ) const override;
  /**
    @param t a Monomial, to compare to @f$ u @f$
    @param u a Monomial, to compare to @f$ t @f$
    @return 1 if @f$ t > u @f$ , -1 if @f$ t < u @f$ , and 0 otherwise
  **/
  virtual int cmp(const Monomial & t, const Monomial & u) const override {
    int result = 0;
    DEG_TYPE a = t.ordering_degree();
    DEG_TYPE b = u.ordering_degree();
    if (a > b) result = 1;
    else if (a < b) result = -1;
    else if (first_larger(t, u)) result = 1;
    else if (first_larger(u, t)) result = -1;
    else result = 0;
    return result;
  }
  /**
    @brief sets the Monomial&rsquo;s @c monomial_ordering_data
  */
  ///@}
  /** @name Utility */
  ///@{
  virtual void set_data(Monomial & t) const override;
  ///@}
};

extern Monomial_Ordering * generic_grevlex_ptr;

/**
  @ingroup orderinggroup
  @class Grevlex_Order_Data
  @author John Perry
  @date 2015
  @brief data for the grevlex monomial ordering

  @details The data involves an array of @f$n@f$ <c>DEG_TYPE</c>,
  where the first entry is the sum of the first @f$n@f$ variables,
  the second entry is the sum of all but the last variable, etc.
*/
class Grevlex_Order_Data : public Monomial_Order_Data {
public:
  /** @name Construction */
  ///@{
  /**
    @brief creates an array of partial weights of @c t
    @param t a Monomial whose parial weights @c this caches
  */
  Grevlex_Order_Data(const Monomial & t);
  /** @brief copy constructor */
  Grevlex_Order_Data(const Grevlex_Order_Data &);
  /** @brief clone &ldquo;constructor&rdquo; */
  virtual Grevlex_Order_Data * clone() override;
  ///@}
  /** @name Destruction */
  ///@{
  /** @brief deletes the array creates by the constructor */
  ///@{
  ~Grevlex_Order_Data();
  ///@}
  /** @name Basic properties */
  ///@{
  /** @brief returns the sum of the first @f$i@f$ variables&rsquo; exponents */
  DEG_TYPE operator [] (NVAR_TYPE i) const;
  virtual DEG_TYPE grading(NVAR_TYPE i) const override { return gradings[i]; }
  ///@}
  /** @name Computation */
  ///@{
  /**
    @brief assigns gradings to a pre-allocated array
    @warning This does not create the array if it does not exist already!
  */
  void assign_gradings(const Monomial &);
  ///@}
protected:
  /** @brief list of partial sums of exponents */
  DEG_TYPE *gradings;
  /** @brief length of @c gradings */
  const NVAR_TYPE number_of_gradings;
};

/**
  @ingroup orderinggroup
  @class Grevlex_Ordering
  @author John Perry
  @date 2015
  @brief the grevlex ordering for a specified number of variables

  @details The grevlex ordering first compares the sums of the exponents,
  then the sums of all but the last exponent,
  then the sums of all but the last two exponents, and so forth,
  until either the sums differ or it runs out of variables.
*/
class Grevlex_Ordering : public Monomial_Ordering {
public:
  /** @name Construction */
  ///@{
  /**
    @brief creates a grevlex ordering specific to the specified number of variables
  */
  Grevlex_Ordering(NVAR_TYPE number_of_variables);
  ///@}
  /** @name Comparison */
  ///@{
  /**
    @brief returns @c true iff @f$t>u@f$ by sums of successively fewer exponents
  */
  virtual bool first_larger(const Monomial & t, const Monomial & u) const override;
  /**
    @brief returns @c true iff @f$t< u@f$ by sums of successively fewer exponents
  */
  virtual bool first_smaller(const Monomial & t, const Monomial & u) const override;
  /**
    @param t a Monomial, to compare to @f$ uv @f$
    @param u a Monomial, to multiply to @f$ v @f$
    @param v a Monomial, to multiply to @f$ u @f$
    @return @c true iff @f$t>uv@f$
  */
  virtual bool first_larger_than_multiple(
      const Monomial & t, const Monomial & u, const Monomial & v
  ) const override;
  /**
    @return the weighted sum of the first i exponents
    @warning Be sure that @c t has the correct ordering!
    @param t a Monomial whose partial degree we want
    @param i index of the indeterminate to which we compute the degree
  */
  DEG_TYPE partial_degree(const Monomial & t, NVAR_TYPE i) const;
  /**
    @param t a Monomial, to compare to @f$ u @f$
    @param u a Monomial, to compare to @f$ t @f$
    @return 1 if @f$ t > u @f$ , -1 if @f$ t < u @f$ , and 0 otherwise
  **/
  virtual int cmp(const Monomial & t, const Monomial & u) const override {
    int result = 0;
    if (first_larger(t, u)) result = 1;
    else if (first_larger(u, t)) result = -1;
    else result = 0;
    return result;
  }
  ///@}
  /** @name Utility */
  ///@{
  /**
    @brief computes the sum of the first i exponents
  */
  DEG_TYPE compute_ith_weight(const Monomial & t, NVAR_TYPE i) const;
  /**
    @brief sets the Monomial&rsquo;s @c monomial_ordering_data
  */
  virtual void set_data(Monomial & t) const override;
  ///@}
protected:
  /** @brief the number of variables, which should remain constant */
  const NVAR_TYPE n;
};

/**
  @ingroup orderinggroup
  @class WGrevlex
  @author John Perry
  @date 2016
  @brief the grevlex ordering for a specified number of variables

  @details The grevlex ordering first compares the sums of the exponents,
  then the sums of all but the last exponent,
  then the sums of all but the last two exponents, and so forth,
  until either the sums differ or it runs out of variables.
*/
class WGrevlex : public Weighted_Ordering {
public:
  /** @name Construction */
  ///@{
  /**
    @brief Creates a grevlex ordering specific to the specified number of
        variables, with the given weights. The final parameter indicates
        whether to apply the weights when breaking ties.
  */
  WGrevlex(NVAR_TYPE, WT_TYPE *, bool=true);
  /**
    @brief Creates a weighted grevlex ordering using the specified ray.
      The final parameter indicates whether to apply the weights
      when breaking ties.
  */
  WGrevlex(Ray r, bool=true);
  ///@}
  /** @name Basic properties */
  ///@{
  /** @brief this weighted ordering&rsquo;s weights */
  virtual const WT_TYPE * order_weights() const override { return weights; }
  /** @brief returns the number of weights (same as number of indeterminates) */
  virtual NVAR_TYPE number_of_weights() const override { return n; }
  ///@}
  /** @name Comparison */
  ///@{
  /**
    @brief returns @c true iff @f$t>u@f$ by weighted sums of
        successively fewer exponents
  */
  virtual bool first_larger(const Monomial & t, const Monomial & u) const override;
  /**
    @brief returns @c true iff @f$t< u@f$ by weighted sums of
        successively fewer exponents
  */
  virtual bool first_smaller(const Monomial & t, const Monomial & u) const override;
  /**
    @brief returns @c true iff @f$t>uv@f$ by weighted sums of
        successively fewer exponents
  */
  virtual bool first_larger_than_multiple(
      const Monomial & t, const Monomial & u, const Monomial & v
  ) const override;
  /**
    @return the weighted sum of the first i exponents
    @warning Be sure that @c t has the correct ordering!
    @param t a Monomial whose partial degree we want
    @param i index of the indeterminate to which we compute the degree
  */
  DEG_TYPE partial_degree(const Monomial & t, NVAR_TYPE i) const;
  /**
    @param t a Monomial, to compare to @f$ u @f$
    @param u a Monomial, to compare to @f$ t @f$
    @return 1 if @f$ t > u @f$ , -1 if @f$ t < u @f$ , and 0 otherwise
  **/
  virtual int cmp(const Monomial & t, const Monomial & u) const override {
    int result = 0;
    DEG_TYPE a = t.ordering_degree();
    DEG_TYPE b = u.ordering_degree();
    if (a > b) result = 1;
    else if (a < b) result = -1;
    else if (first_larger(t, u)) result = 1;
    else if (first_larger(u, t)) result = -1;
    else result = 0;
    return result;
  }
  ///@}
  /** @name Utility */
  ///@{
  /**
    @brief computes the weighted sum of the first i exponents
  */
  DEG_TYPE compute_ith_weight(const Monomial & t, NVAR_TYPE i) const;
  /**
    @brief sets the Monomial&rsquo;s @c monomial_ordering_data
  */
  virtual void set_data(Monomial & t) const override;
  ///@}
  /** @name I/O */
  ///@{
  friend ostream & operator<<(ostream & os, const WGrevlex &word) {
    os << "WGOrd( ";
    for (unsigned i = 0; i < word.n; ++i) {
      os << word.weights[i] << ' ';
    }
    os << ")";
    return os;
  }
  ///@}
protected:
  /** @brief the number of variables, which should remain constant */
  const NVAR_TYPE n;
  /** @brief the weights for this ordering */
  WT_TYPE * weights;
  /** @brief whether to apply the weights to all the variables */
  bool thorough_weighting;
};

/**
  @ingroup orderinggroup
  @class Lex_Ordering
  @author John Perry
  @date 2015
  @brief the lex ordering for a specified number of variables

  @details The lex ordering is a dictionary ordering. The first variable is
    considered largest, and monomials with a larger degree in the first
    variable will be considered larger than monomials with a smaller degree in
    the first variable, regardless of the overall degree in all variables.
*/
class Lex_Ordering : public Monomial_Ordering {
public:
  /** @name Construction */
  ///@{
  /**
    @brief creates a lex ordering specific to @f$n@f$ variables
    @param number_of_variables number of variables this ordering should check
  */
  Lex_Ordering(NVAR_TYPE number_of_variables);
  ///@}
  /** @name Comparison */
  ///@{
  /**
    @return @c true iff @f$t>u@f$
    @param t a Monomial, to compare to @f$ u @f$
    @param u a Monomial, to compare to @f$ t @f$
  */
  virtual bool first_larger(const Monomial & t, const Monomial & u) const override;
  /**
    @return @c true iff @f$t< u@f$
    @param t a Monomial, to compare to @f$ u @f$
    @param u a Monomial, to compare to @f$ t @f$
  */
  virtual bool first_smaller(const Monomial & t, const Monomial & u) const override;
  /**
    @param t a Monomial, to compare to @f$ uv @f$
    @param u a Monomial, to multiply to @f$ v @f$
    @param v a Monomial, to multiply to @f$ u @f$
    @return @c true iff @f$t>uv@f$
  */
  virtual bool first_larger_than_multiple(
      const Monomial & t, const Monomial & u, const Monomial & v
  ) const override;
  /**
    @param t a Monomial, to compare to @f$ u @f$
    @param u a Monomial, to compare to @f$ t @f$
    @return 1 if @f$ t>u @f$ , -1 if @f$ t < u @f$ , and 0 otherwise
  */
  virtual int cmp(const Monomial & t, const Monomial & u) const override {
    int result = 0;
    if (first_larger(t, u)) result = 1;
    else if (first_larger(u, t)) result = -1;
    else result = 0;
    return result;
  }
  ///@}
protected:
  /** @brief the number of variables, which should remain constant */
  const NVAR_TYPE n;
};

/**
  @ingroup orderinggroup
  @class WGrevlex_Order_Data
  @author John Perry
  @date 2015
  @brief data for the weighted grevlex monomial ordering

  @details The data involves an array of @f$n@f$ <c>DEG_TYPE</c>,
  where the first entry is a weighted sum of the first @f$n@f$ variables,
  the second entry is the \e ordinary sum of all but the last variable, etc.
*/
class WGrevlex_Order_Data : public Monomial_Order_Data {
public:
  /** @name Construction */
  ///@{
  /**
    @brief creates an array of partial weights of @c t
    @warning Assign the correct ordering to @c t first!
    @param t a Monomial whose weights @c this will cache
  */
  WGrevlex_Order_Data(Monomial & t);
  /** @brief copy constructor */
  WGrevlex_Order_Data(const WGrevlex_Order_Data &);
  /** @brief clone constructor */
  virtual WGrevlex_Order_Data * clone() override;
  ///@}
  /** @name Destruction */
  ///@{
  /** @brief deletes the array of partial weights */
  ~WGrevlex_Order_Data();
  ///@}
  /** @name Basic properties */
  ///@{
  /** @name returns the weighted sum of the first @f$i@f$ variables */
  DEG_TYPE operator [] (NVAR_TYPE i) const;
  inline virtual DEG_TYPE grading(NVAR_TYPE i) const override { return gradings[i]; }
  ///@}
  /** @name Computation */
  ///@{
  /**
    @brief assigns gradings to a pre-allocated array
    @warning This does not create the array if it does not exist already!
    @param t a Monomial whose gradings we want
  */
  void assign_gradings(Monomial &t);
  ///@}
  /** @name Memory management */
  ///@{
  /** @brief requests memory form WGrevlex_Ordering's Grading_Order_Data_Allocator */
  void * operator new(size_t);
  /** @brief returns data to WGrevlex_Ordering's Grading_Order_Data_Allocator */
  void operator delete(void *);
  ///@}
protected:
  /** @brief array of partial weighted sums of exponents */
  DEG_TYPE *gradings;
  /** @brief length of @c gradings */
  const NVAR_TYPE number_of_gradings;
};

/**
  @ingroup orderinggroup
  @class CachedWGrevlex_Ordering
  @author John Perry
  @date 2015
  @brief the weighted grevlex ordering for a specified number of variables,
    with cached weights for each monomial

  @details The weighted grevlex ordering first compares weighted sums
  of all exponents, then the sums of all but the last exponent,
  then the sums of all but the last two exponents, and so forth,
  until either the sums differ or it runs out of variables.
  Whether the weights are applied to more than one sum can be decided
  at construction.
*/
class CachedWGrevlex_Ordering : public Weighted_Ordering {
public:
  /** @name Construction */
  ///@{
  /**
    @brief creates a weighted grevlex ordering specific to @f$n@f$ variables,
      using the weights specified by @f$w@f$
    @param thorough determines whether the weights are applied to subsequent
      orderings; if @c false, the second, third, etc. sums are ordinary sums of
      the exponents, rather than weighted sums
    @param number_of_variables the number of variables this will compare;
      @c weights should have this length!
    @param w the weights to multiply to each exponent
  */
  CachedWGrevlex_Ordering(
      NVAR_TYPE number_of_variables, WT_TYPE * w, bool thorough=true
  );
  ///@}
  /** @name Basic properties */
  ///@{
  /** @brief the weights that define this ordering */
  virtual const WT_TYPE * order_weights() const override;
  /** @brief returns the number of weights (same as number of indeterminates) */
  virtual NVAR_TYPE number_of_weights() const override { return n; }
  ///@}
  /** @name Comparison */
  ///@{
  /** @brief resturns 0 if they are alike; positive if first larger; negative otherwise */
  virtual int cmp(const Monomial &, const Monomial &) const override;
  /**
    @brief returns @c true iff @f$t>u@f$ by sums of successively fewer exponents
  */
  virtual bool first_larger(const Monomial & t, const Monomial & u) const override;
  /**
    @brief returns @c true iff @f$t< u@f$ by sums of successively fewer exponents
  */
  virtual bool first_smaller(const Monomial & t, const Monomial & u) const override;
  /**
    @brief returns @c true iff @f$t>u@f$ by sums of successively fewer exponents
  */
  virtual bool first_larger_than_multiple(
      const Monomial & t, const Monomial & u, const Monomial & v
  ) const override;
  ///@}
  /** @name Utility */
  ///@{
  /**
    @return the weighted sum of the first i exponents
    @param t a Monomial whose degree we want
    @param i index to an exponent
    @warning Be sure that @c t has the correct ordering!
  */
  inline DEG_TYPE partial_degree(const Monomial & t, NVAR_TYPE i) const {
    return t.monomial_ordering_data()->grading(i);
  }
  /**
    @return the sum of the first i exponents
    @param t a Monomial whose degree we want
    @param i index to an exponent
  */
  DEG_TYPE compute_ith_weight(const Monomial & t, NVAR_TYPE i) const;
  /**
    @brief sets the Monomial&rsquo;s @c monomial_ordering_data
  */
  virtual void set_data(Monomial & t) const override;
  ///@}
  /** @name I/O */
  ///@{
  friend ostream & operator<<(ostream & os, const CachedWGrevlex_Ordering &word) {
    os << "WGOrd( ";
    for (unsigned i = 0; i < word.n; ++i) {
      os << word.weights[i] << ' ';
    }
    os << ")";
    return os;
  }
  ///@}
protected:
  /** @brief the number of variables, which should remain constant */
  const NVAR_TYPE n;
  /** @brief the weights that define this ordering */
  const WT_TYPE * weights;
  /** @brief whether to apply the weights to more than the first sum */
  const bool fully_apply;
};

/**
  @ingroup orderinggroup
  @class Nonsingular_Matrix_Ordering_Exception
  @author John Perry
  @date 2016
  @brief exceptions for Matrix_Ordering
*/
class Nonsingular_Matrix_Ordering_Exception : public std::exception {
  virtual const char * what() const throw() override;
};

/**
  @ingroup orderinggroup
  @class Matrix_Ordering
  @author John Perry
  @date 2016
  @brief orderings defined by a nonsingular matrix
*/
class Matrix_Ordering : Monomial_Ordering {
public:
  /**
    @param rows number of rows desired in the matrix
    @param cols number of columns desired in the matrix
    @param data should contain @f$ rows \times cols @f$ elements
    @brief checks that @c data defines a nonsingular matrix, and sets things up
  */
  Matrix_Ordering(NVAR_TYPE rows, NVAR_TYPE cols, const WT_TYPE **data);
  /**
    @param t a Monomial, to compare to @f$ u @f$
    @param u a Monomial, to compare to @f$ t @f$
    @return @c true iff first Monomial is larger than second
  */
  virtual bool first_larger(const Monomial & t, const Monomial & u) const override;
  /**
    @param t a Monomial, to compare to @f$ u @f$
    @param u a Monomial, to compare to @f$ t @f$
    @return @c true iff first Monomial is smaller than second
  */
  virtual bool first_smaller(const Monomial & t, const Monomial & u) const override;
  /**
    @param t a Monomial, to compare to @f$ uv @f$
    @param u a Monomial, to multiply to @f$ v @f$
    @param v a Monomial, to multiply to @f$ u @f$
    @return @c true iff first Monomial is larger than product of second
      and third
  */
  virtual bool first_larger_than_multiple(
      const Monomial & t, const Monomial & u, const Monomial & v
  ) const override;
  /**
    @param t a Monomial, to compare to @f$ u @f$
    @param u a Monomial, to compare to @f$ t @f$
    @return 1 if @f$ t>u @f$ , -1 if @f$ t < u @f$ , and 0 otherwise
  */
  virtual int cmp(const Monomial & t, const Monomial & u) const override {
    int result = 0;
    if (first_larger(t, u)) result = 1;
    else if (first_larger(u, t)) result = -1;
    else result = 0;
    return result;
  }
protected:
  /** @brief the number of rows */
  const NVAR_TYPE m;
  /** @brief the number of columns */
  const NVAR_TYPE n;
  /** @brief the matrix that defines this ordering */
  const WT_TYPE ** W;
};

#endif