#ifndef __LP_SOLVER_HPP_
#define __LP_SOLVER_HPP_

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

#include <set>
using std::set;
#include <vector>
using std::vector;
#include <iostream>
using std::ostream; using std::cout; using std::endl;

#include "system_constants.hpp"

#include "monomial.hpp"

/**
  @defgroup CLSSolvers Constrained Linear System Solvers
  @brief classes that solve constrained linear systems
  @details Classes in this group solve constrained linear systems;
  that is, they solve systems of the form
  @f$\left\{\sum_{j=1}^na_{ij}x_j\leq b_i\right\}_{i=1}^m@f$ (where the inequality
  may or may not be strict).
*/

/**
  @namespace LP_Solvers
  @ingroup CLSSolvers 
  @brief classes that solve constrained linear systems
  @details Classes in this group solve constrained linear systems;
  that is, they solve systems of the form
  @f$\left\{\sum_{j=1}^na_{ij}x_j\leq b_i\right\}_{i=1}^m@f$ (where the inequality
  may or may not be strict).
*/
namespace LP_Solvers {

/**
  @brief memory manager for ray entries
  @ingroup memorygroup
  @details Automatically initialized, but clients need to call the destructor
    when finished.
*/
extern Grading_Order_Data_Allocator<DEG_TYPE> * doda;

/**
  @brief a constraint @f$ c_1 x_1 + \ldots + c_n x_n \geq 0 @f$
  @author John Perry
  \version 1.0
  @date October 2014
  @copyright The University of Southern Mississippi
  @details This class encapsulates a simple constraint for a skeleton; that is,
  an inequality of the form @f$ c_1 x_1 + \ldots + c_n x_n \geq 0 @f$.
  Constraints can be ordered lexicographically using the less-than operator,
  allowing for their inclusion in ordered collections, such as sets.
  @ingroup CLSSolvers
*/
class Constraint
{

public:

  /** @name Construction */
  ///@{

  /**
    @brief Initialize constraint to the given coefficients.

    The resulting constraint is @f$ c_1x_1 + \cdots + c_nx_n \geq 0, @f$
    where @f$ c_i @f$ is the coefficient of @f$ x_i @f$.
    @param num_variables length of coeffs
    @param coeffs copies this array of coefficients
    @pre the size of the array needs to be at least as long as the dimension!
  */
  Constraint(NVAR_TYPE, CONSTR_TYPE []);

  /**
    @brief Initialize constraint to the given coefficients.

    The resulting constraint is @f$ c_1x_1 + \cdots + c_nx_n \geq 0, @f$
    where @f$ c_i @f$ is the coefficient of @f$ x_i @f$.
    @param coeffs copies thiis vector of coefficients
    @post @c nvars will have the value equal to `coeffs.size()`
  */
  Constraint(vector<CONSTR_TYPE> &);

  /**
    @brief Copies the coefficients of the other constraint,
    including the allocation of new memory.
  */
  Constraint(const Constraint &);

  ///@}

  /** @name Destruction */
  ///@{

  /**
    @brief Deletes memory allocated by the constructor.

    Currently, that means it deletes an array created by the constructors.
  */
  ~Constraint();

  ///@}

  /** @name Basic properties */
  ///@{

  /**
    @brief Returns the number of variables in the constraint.
  */
  inline NVAR_TYPE get_number_of_variables() const { return nvars; }

  /**
    @brief Returns the coefficient indicated. Numbering starts at 0.
  */
  inline CONSTR_TYPE operator[](NVAR_TYPE index) const { return coefficients[index]; };

  /** @brief Returns the coefficients that determine this constraints. */
  inline const CONSTR_TYPE * coeffs() const { return coefficients; }

  ///@}

  /** @name Comparisons */
  ///@{

  friend bool operator<(const Constraint & a, const Constraint & b);

  friend bool operator==(const Constraint & a, const Constraint & b);

  friend bool operator!=(const Constraint & a, const Constraint & b);

  ///@}

  /** @name I/O */
  ///@{

  friend ostream & operator<<(ostream & os, const Constraint & c);

  ///@}

private:

  /** number of variables/coefficients */
  NVAR_TYPE nvars;
  /** coefficients */
  CONSTR_TYPE * coefficients;

};

/**
  @brief Lexicographic comparison of constraints.
  @warning This is unsafe when number of variables is not the same.
    It does not check, since the assumption is that you know what you're doing.
  @param a first constraint
  @param b second constraint
  @return @c true if and only if @p a and @p b have different entries
*/
bool operator<(const Constraint & a, const Constraint & b);

/** @brief check for constraint equality
  @warning This is unsafe when number of variables is not the same.
    It does not check, since the assumption is that you know what you're doing.
  @param a first constraint
  @param b second constraint
  @return @c true if and only if @p a and @p b have the same entries
*/
bool operator==(const Constraint & a, const Constraint & b);

/** @brief check for constraint inequality
  @warning This is unsafe when number of variables is not the same.
    It does not check, since the assumption is that you know what you're doing.
  @param a first constraint
  @param b second constraint
  @return @c true if and only if @p a and @p b have different entries
*/
bool operator!=(const Constraint & a, const Constraint & b);

/** @brief print a representation of the constraint to the stream

  @details Output is of the form @f$ c_1 x_1 + \ldots + c_n x_n @f$ ,
    where @f$ c_i @f$ is the coefficient of @f$ x_i @f$.
  
  @param os output stream
  @param c constraint to print
  @return the output stream
*/
ostream & operator<<(ostream & os, const Constraint & c);

/**
  @brief a ray defined by nonnegative coordinates @f$(a_1,\ldots,a_n)@f$
  @author John Perry
  \version 1.0
  @date October 2014
  @copyright The University of Southern Mississippi
  @ingroup CLSSolvers
  @details This class encapsulates a ray, one major part of the definition of a skeleton.
  Rays can be initialized to a particular set of coefficients, or to a particular
  axis (which is then translated into the corresponding coefficients).

  A special feature is that a rays can track the constraints known to be
  active at the ray, allowing for more efficient computation in the double
  description method. Adding known constraints can be done with or without
  checking whether the constraint actually is active, so this should be done
  with care.
*/
class Ray
{

public:

  /** @name Construction */
  ///@{

  /**
    @brief Creates a ray with the given number of variables, all set to 0.

    The optional second argument specifies a direction, and sets that coordinate
    to 1. In this case, there is no need to set the ray's known active constraints,
    as this is known and populated automatically.
    @pre The dimension should be greater than zero. While the direction need not
      be specified&hellip; (see postcondition)
    @post &hellip;the result when the direction is zero is a zero ray.
      If the direction is @f$ i @f$, then the result is the @f$i@f$th canonical vector.
  */
  Ray(NVAR_TYPE, long = -1);

  /**
    @brief Creates a ray with the given number of variables,
    with coordinates set to the value of the array.
    @pre the size of the array needs to be at least as long
      as the number of variables!
  */
  Ray(NVAR_TYPE, const RAYENT_TYPE []);

  /**
    @brief Creates a ray with the given number of variables,
    with coordinates set to the value of the array.
    @pre the size of the array needs to be at least as long
      as the number of variables!
  */
  Ray(NVAR_TYPE, const EXP_TYPE []);

  /**
    @brief Creates a ray whose coordinates are given by the vector.
    @post The dimension of this ray will equal the number of entries in the vector,
      and the values of their entries will be equal.
  */
  Ray(const vector<RAYENT_TYPE> &);

  /**
    @brief Copies the coordinates of the other ray.

    Allocates new memory, and copies the active constraints.
  */
  Ray(const Ray &);

  ///@}

  /** @name Destruction */
  ///@{

  /**
    @brief Deletes memory allocated by the constructor.

    Currently, that means it deletes `coords`.
  */
  ~Ray();

  ///@}

  /** @name Basic properies */
  ///@{
  
  /** @brief Returns the dimension of this ray. */
  inline NVAR_TYPE get_dimension() const { return dim; };

  /** @brief Returns the entry indicated. Numbering starts at 0. */
  inline RAYENT_TYPE operator[](NVAR_TYPE index) const { return coords[index]; };

  /** @brief Returns the weights. */
  inline const RAYENT_TYPE * weights() const { return coords; };

  /** @brief Synonym for []. I have no idea why I added this. */
  inline const RAYENT_TYPE coordinate(NVAR_TYPE index) { return coords[index]; };

  /**
    @brief Returns `true` if and only if the hyperplane is active at this ray.
    @param hyperplane a constraint; we would like to know whether @p this lies
        on it
    @return true if and only if @p this lies on @p constraint
    @details Practically speaking, if the hyperplane is defined by the vector
    @f$ \mathbf c @f$ and the ray is defined by @f$ \mathbf r @f$ ,
    this function returns true if and only if @f$ c\cdot r = 0 @f$.
  */
  inline bool is_active_at(const Constraint &hyperplane) const
  {
    return 0 == obtain_dot_product(hyperplane);
  };

  /**
    @brief Returns `true` if and only if this ray is above the hyperplane.
    @param hyperplane a constraint; we would like to know whether @p this is
        above it
    @return true if and only if @p this is above @p constraint
    @details Practically speaking, if the hyperplane is defined by the vector
    @f$ \mathbf c @f$ and the ray is defined by @f$ \mathbf r @f$ ,
    this function returns true if and only if @f$ c\cdot r > 0 @f$.
  */
  inline bool is_above(Constraint &hyperplane)
  {
    return 0 < obtain_dot_product(hyperplane);
  };

  /**
    @brief Returns `true` if and only if this ray is below the hyperplane.
    @param hyperplane a constraint; we would like to know whether @p this is
        below it
    @return true if and only if @p this is below @p constraint
    @details Practically speaking, if the hyperplane is defined by the vector
    @f$ \mathbf c @f$ and the ray is defined by @f$ \mathbf r @f$ ,
    this function returns true if and only if @f$ c\cdot r < 0 @f$.
  */
  inline bool is_below(Constraint &hyperplane)
  {
    return 0 > obtain_dot_product(hyperplane);
  };

  ///@}

  /** @name Comparison */
  ///@{

  friend bool operator==(const Ray & a, const Ray & b);

  friend bool operator!=(const Ray & a, const Ray & b);

  friend bool operator < (const Ray & a, const Ray & b);

  ///@}

  /** @name Computation */
  ///@{

  /**
    @brief Convenience function to compute dot product between ray
        and the given constraint.
    @return the dot product of @p this and @p constraint
    @warning This is unsafe when dimension is not the same.
      It does not check, since the assumption is that you know what you're doing.
  */
  DOTPROD_TYPE obtain_dot_product(const Constraint &) const;

  ///@}

  /** @name Modification */
  ///@{

  /**
    @brief Simplifies the ray by dividing its components by the least common
    denominator.
  */
  void simplify_ray();

  /**
    @brief Assignment operator; assigns the value of `other` to `this`.
    @return @p this
    @warning This is unsafe when dimension is not the same.
      It does not check, since the assumption is that you know what you're doing.
  */
  Ray & operator=(const Ray &);

  /**
    @brief Swap two rays of equal dimension by swapping their data,
      avoiding memory reallocation.
    @warning This is unsafe when dimension is not the same.
      It does not check, since the assumption is that you know what you're doing.
  */
  void swap(Ray &);

  ///@}

  /** @name I/O */
  ///@{

  friend ostream & operator<<(ostream & os, const Ray & r);

  ///@}

private:

  NVAR_TYPE dim; /**< number of entries in @c coords */

  RAYENT_TYPE * coords; /**< coordinates of the ray */

};

/**
  @brief Output is of the form @f$(r_1, \ldots, r_n)@f$.
  @param os output stream
  @param r ray to print
  @return the output stream
*/
ostream & operator<<(ostream & os, const Ray & r);

/**
  @brief indicates whether the two rays are equal
  @warning This is unsafe when number of variables is not the same.
    It does not check, since the assumption is that you know what you're doing.
  @param a first ray
  @param b second ray
  @return @c true if and only if the rays&rsquo; entries have the same values
  @details Notice that the rays can point in the same direction,
    but still be considered unequal.
*/
bool operator==(const Ray & a, const Ray & b);

/**
  @brief Indicates whether the two rays are unequal.
  @warning This is unsafe when number of variables is not the same.
    It does not check, since the assumption is that you know what you're doing.
  @param a first ray
  @param b second ray
  @return @c true if and only if the rays&rsquo; entries have different values
  @details Notice that the rays can point in the same direction,
    but still be considered unequal.
*/
bool operator!=(const Ray & a, const Ray & b);

/**
  @brief Lexicographic comparison of rays.
  @warning This is unsafe when dimension is not the same.
    It does not check, since the assumption is that you know what you're doing.
  @param a first ray
  @param b second ray
  @return @c true if and only if @p a is lexicographically smaller than @p b
*/
bool operator < (const Ray & a, const Ray & b);

/**
  \ingroup CLSSolvers
  @brief Multiply every coordinate in the given ray by the given scalar.
  @return a copy of the ray, scaled by the requesed amount
*/
Ray operator*(const RAYENT_TYPE, const Ray &);

/**
  @ingroup CLSSolvers
  @brief compute the dot product of the specified rays, one of which is a vector
  @return the dot product
  @warning This is unsafe when dimension is not the same.
    It does not check, since the assumption is that you know what you're doing.
*/
RAYENT_TYPE operator*(const Ray &, const vector<long> &);

/**
  @ingroup CLSSolvers
  @brief compute the dot product of the specified rays, one of which is a vector
  @return the dot product
  @warning This is unsafe when dimension is not the same.
    It does not check, since the assumption is that you know what you're doing.
*/
RAYENT_TYPE operator*( const vector<long> &, const Ray &);

/**
  \ingroup CLSSolvers
  @brief Add the two rays.
  @return the sum of the two rays
  @warning This is unsafe when dimension is not the same.
    It does not check, since the assumption is that you know what you're doing.
*/
Ray operator+(const Ray &, const Ray &);

/**
  \ingroup CLSSolvers
  @brief Subtract the two rays.
  @return the difference of the two rays
  @warning This is unsafe when dimension is not the same.
    It does not check, since the assumption is that you know what you're doing.
*/
Ray operator-(const Ray &, const Ray &);

/**
  \ingroup CLSSolvers
  @brief Add all the rays in a set.
  @return a ray that is the sum of all rays in the given set
  @warning This is unsafe when dimension is not the same.
    It does not check, since the assumption is that you know what you're doing.
*/
Ray ray_sum(const set<Ray> &);

/**
  \ingroup CLSSolvers
  @brief Compute the dot product on the rays.
  @return the dot product of the rays
  @warning This is unsafe when dimension is not the same.
    It does not check, since the assumption is that you know what you're doing.
*/
RAYENT_TYPE operator*(const Ray &, const Ray &);

/**
  \ingroup CLSSolvers
  @brief Compute the dot product between the ray and the constraint.
  @param r a ray
  @param c a constraint
  @return the dot product of the ray and the constraint
  @warning This is unsafe when dimension is not the same.
    It does not check, since the assumption is that you know what you're doing.
*/
inline DOTPROD_TYPE operator*(const Ray &r, const Constraint &c)
{ return r.obtain_dot_product(c); }

/**
  \ingroup CLSSolvers
  @brief Compute the dot product between the ray and the constraint.
  @param c a Constraint
  @param r a ray
  @return the dot product of the ray and the constraint
  @warning This is unsafe when dimension is not the same.
    It does not check, since the assumption is that you know what you're doing.
*/
inline DOTPROD_TYPE operator*(Constraint &c, Ray &r)
{ return r.obtain_dot_product(c); }

/**
  @brief exact or approximate polyhedral cone solution,
      with methods allowing definition and refinement
  @author John Perry
  \version 1.0
  @date January 2017
  @copyright The University of Southern Mississippi
  @ingroup CLSSolvers
  @details  This class encapsulates the skeleton of a polyhedral cone,
  defined by a sequence of inequalities of the form
  @f$ c_1 x_1 + \cdots c_n x_n \geq 0 @f$.

  @warning Some classes may provide only an <i>approximate</i> cone; see, for
    example, GLPK_Solver. In addition, Clients must ensure two things.
      -# The rays must
         have the same number @f$ m @f$ of dimensions, constraints must have
         the same number @f$ n @f$ of variables, and @f$ m=n @f$. Violating any
         of these three conditions will lead to undesirable behavior.
      -# When refining the cone, it is essential to check that
         the return value of solve() is @c true; for if it is not,
         then the cone is no longer be consistent.
         Please read the relevant documentation.
*/

class LP_Solver {
public:
  /** @name Construction */
  ///@{
  /**
    @brief performs a deep copy, similar to a copy constructor
    @return @c true iff copying was successful
    @warning Do not mix-and-match solvers. At the present time, a PPL_Solver
        is not equipped to copy a GLPK_Solver, or vice versa.
        (This doesn't even make sense between exact and approximate solvers.)
  */
  virtual bool copy(const LP_Solver *) = 0;
  ///@}
  /** @name Destruction */
  ///@{
  /** @brief the default destructor does nothing (this is an abstract class) */
  virtual ~LP_Solver() { }
  ///@}
  /** @name Modification */
  ///@{
  /**
    @brief Adds the indicated constraint (singular!) and re-computes the solution.

    \return @c true if and only if the new constraint is consistent with the
      current constraints

    @warning Checking the return value is crucial!
      If the function returns @c false, you have an inconsistent system!
      While the <i>present</i> cone will remain consistent,
      <b>the function will not roll back previous changes you have made</b>,
      so if you want to iterate again,
      your best bet is to copy the skeleton, and try that copy.
      Accept the new constraints only if that copy succeeds,
      in which case, you might as well discard the original, and keep the copy.
  */
  virtual bool solve(const Constraint &) = 0;
  /**
    @brief Adds the indicated constraints (plural!) and re-computes the solution.

    \return @c true if and only if the new constraints are consistent with the
      current constraints

    @warning Checking the return value is crucial!
      If the function returns @c false, you have an inconsistent system!
      While the <i>present</i> cone will remain consistent,
      <b>the function will not roll back previous changes you have made</b>,
      so if you want to iterate again,
      your best bet is to copy the skeleton, and try that copy.
      Accept the new constraints only if that copy succeeds,
      in which case, you might as well discard the original, and keep the copy.
  */
  virtual bool solve(const vector<Constraint> &) = 0;

  ///@}

  /** @name Basic properies */
  ///@{
  /** @brief Returns the dimension of the underlying vector space. */
  virtual NVAR_TYPE get_dimension() const = 0;
  /** @brief Returns the number of rays defining the skeleton. */
  virtual unsigned long get_number_of_rays() { return rays.size(); }
  /**
    @brief Returns rays that define a skeleton.
    @return a set of Ray that define the given skeleton
    @details When using an approximate solver such as GLPK_Solver,
      this will give only an approximate skeleton.
  */
  virtual const set<Ray> & get_rays() const;
  /**
    @brief returns the number of constraints used by the skeleton
    @return number of constraints
  */
  virtual unsigned long get_number_of_constraints() const = 0;
  ///@}

  /** @name Computation */
  ///@{
  /** @brief tests for consistency of a constraint generated by two monomials. */
  virtual inline bool makes_consistent_constraint(
      const Monomial & t, const Monomial & u, bool show_data = false
  ) const {
    bool inconsistent = true;
    for (
          auto riter = rays.begin();
          inconsistent and riter != rays.end();
          ++riter
    ) {
      int d = 0;
      for (int i = 0; i < riter->get_dimension(); ++i)
        d += (*riter)[i] * (t.degree(i) - u.degree(i));
      if (d > 0)
        inconsistent = false;
      if (show_data) {
        cout << d << ' ';
        if (!inconsistent) cout << *riter << endl;
      }
    }
    if (show_data) cout << endl;
    return not inconsistent;
  }
  ///@}
protected:
  set<Ray> rays; /**< the skeleton (may be approximate, depending on solver) */
};

}

#endif