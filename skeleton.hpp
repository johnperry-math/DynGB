#ifndef SKELETON_H
#define SKELETON_H

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
  @brief implementation of classes for double description method
  @author John Perry
  \version 1.0
  @date October 2014
  @copyright The University of Southern Mississippi
  @ingroup CLSSolvers
*/

#include <set>
#include <vector>
#include <iostream>
using std::cout; using std::endl;

#include "lp_solver.hpp"
#include "polynomial.hpp"

#include "system_constants.hpp"

namespace LP_Solvers {

/**
  @brief an edge @f$(r_1,r_2)@f$ connecting the two rays @f$ r_1 @f$ and @f$ r_2 @f$
  @author John Perry
  \version 1.0
  @date October 2014
  @copyright The University of Southern Mississippi
  @ingroup CLSSolvers
  @details  This class encapsulates an edge, the other major part of a skeleton.
  Edges describe how the rays of the skeleton are connected.
  Edges are ordered, so that the smaller ray always comes first.

  @warning An edge's rays should have the same dimension.
    To start with, it doesn't make mathematical sense to &ldquo;join&rdquo;
    two rays of different dimension.
    Moreover, comparison of edges requires comparison of rays,which requires
    that the rays have the same dimension.
    (But you wouldn't be dumb enough to do this in the first place.)
*/
class Edge {

public:

  /** @name Construction */
  ///@{

  /** @brief Creates a new edge that joins the two rays. */
  Edge(const Ray &, const Ray &);

  /** @brief Copies the rays in @c other to two new rays. */
  Edge(const Edge &);

  ///@}

  /** @name Destruction */
  ///@{

  ~Edge() {} /**< Does nothing beyond what the compiler would do. */

  ///@}

  /** @name Basic properties */
  ///@{

  /** @brief Returns the first ray listed in this edge. */
  inline Ray get_first_ray() const { return first; };

  /** @brief Returns the second ray listed in this edge. */
  inline Ray get_second_ray() const { return second; };

  ///@}

  /** @name Comparison */
  ///@{
  
  friend bool operator==(const Edge & e1, const Edge & e2);

  friend bool operator<(const Edge & e1, const Edge & e2);

  ///@}

  /** @name I/O */
  ///@{

  friend ostream & operator<<(ostream &ostr, const Edge &e);

  ///@}

  /** @name Modification */
  ///@{

  /** @brief Assignment operator */
  Edge & operator=(const Edge &);

  ///@}

private:

  Ray first, second; /**< the rays defining this edge */

};

/**
  @brief Equal if and only if.
  @param e1 an edge
  @param e2 an edge
  @return @c true if and only if the edges have identical entries
*/
bool operator==(const Edge &e1, const Edge &e2);

/**
  @brief Compares two edges lexicographically.

  @details If the first ray in @c this edge is smaller,
  then @c this edge is smaller.
  Otherwise, if the first rays are equal, and the second ray in @c this edge
  is smaller, then @c this edge is smaller.

  @return @c true if and only if the first edge is lexicographically smaller

  @param e1 an edge
  @param e2 an edge
*/
bool operator<(const Edge & e1, const Edge & e2);

///@}

/**
  @brief Output has the form @f$ \{ \mathbf{r}_1, \mathbf{r}_2 \} @f$
  where @f$ \mathbf{r}_1 @f$ is the first ray in this edge, etc.
  @param ostr output stream to write to
  @param e edge to write
  @return the output stream
*/
ostream & operator<<(ostream &ostr, const Edge &e);

/**
  \ingroup CLSSolvers
  @brief counts the number of constraints active in both sets
  @return the number of constraints common to both sets.
*/
int number_of_common_constraints(bool *, bool *, unsigned);


/**
  \ingroup CLSSolvers
  @brief indicates which constraints are active for both sets
  @return the intersection between the given sets of constraints.
*/
vector<bool> intersections_of_active_constraints(bool *, bool *, unsigned);

/**
  \ingroup CLSSolvers
  @brief indicates which constraints are active for both sets
  @param a first set
  @param b second set
  @param result set where @c true occurs only if it does in both @p a and @p b
  @param m number of entries in @p a and @p b
  @return the intersection between the given sets of constraints.
*/
void intersections_of_active_constraints(
    bool * a, bool * b, bool * result, unsigned m
);


/**
  \ingroup CLSSolvers
  @brief determines whether the first set of active constraints is a subset
    of the second
  @return @c true if and only if the first set is a subset of the second.
*/
bool is_first_subset_of_second(bool *, bool *, unsigned);

/**
  \ingroup CLSSolvers
  @brief computes the union of the specified edge sets
  @return union of the specified edge sets
*/
set<Edge> union_of_edge_sets(const set<Edge> &, const set<Edge> &);

/**
  @brief skeleton of a polyhedral cone, with methods allowing definition and refinement
  @author John Perry
  \version 1.1
  @date October 2014
  @copyright The University of Southern Mississippi
  @ingroup CLSSolvers
  @details  This class implements the Double Description Method,
  an iterative algorithm for computing the skeleton of a cone.
  This particular version uses Zolotykh's GraphAdj criterion
  @cite Zolotych_DoubleDescription.
  The iterative nature means that the cone can be updated with new constraints,
  passed to that algorithm, and the skeleton will be automatically recomputed.
*/
class Skeleton : public LP_Solver {

public:

  /** @name Construction */
  ///@{

  /** @brief Initialization common to all constructors */
  void common_initialization(NVAR_TYPE);

  /**
    @brief Constructs a basic skeleton in the given number of dimensions,
    initialized to the axes, or (equivalently) to the set of constraints
    @f$ x_i \geq 0 @f$.

    The rays are informed of their active constraints.
    @pre the argument should be at least two
    @post the skeleton of the positive orthant
  */
  explicit Skeleton(NVAR_TYPE);

  /**
    @brief Constructs a skeleton described by the given system of constraints.

    Practically speaking, it first generates a basic skeleton,
    then iterates on the given constraints.
    @pre `u.size() == v.size()` for all `u`, `v` in the vector
    @post unless the system supplied was inconsistent, a valid skeleton of the
        corresponding polyhedral cone
    @warning Your program will almost certainly fail
    if you do not respect the precondition.
  */
  Skeleton(NVAR_TYPE, const vector<Constraint> &);

  /** @brief Performs a deep copy of `other`. */
  explicit Skeleton(const Skeleton &);

  virtual bool copy(const LP_Solver *) override;

  ///@}

  /** @name Destruction */
  ///@{

  /**
    @brief Currently does nothing the compiler wouldn't do.
  */
  virtual ~Skeleton();

  ///@}

  /** @name Basic properties */
  ///@{

  inline NVAR_TYPE get_dimension() const override { return dim; };

  /** @brief Returns the number of edges defining the skeleton. */
  inline unsigned long get_number_of_edges() { return edges.size(); };

  /** @brief Returns the edges that define the skeleton. */
  inline set<Edge> get_edges() { return edges; };

  /** @brief Returns the number of constraints defining the skeleton. */
  inline unsigned long get_number_of_constraints() const override {
    return constraints.size();
  };

  /** @brief Returns the constraints that define the skeleton. */
  inline const vector<Constraint> & get_constraints() {
    return constraints;
  };

  /** @brief Returns the indicated constraint. Numbering starts at 0. */
  inline const Constraint & get_constraint(int index) {
    return constraints[index];
  };

  ///@}

  /** @name Computation */
  ///@{

  /** @brief returns the set of constraints in the skeleton active at @c u */
  //inline vector<bool> which_constraints_active_at(const Ray & u) const {
  inline void which_constraints_active_at(const Ray & u, bool * result) const {
    for (unsigned i = 0; i < constraints.size(); ++i) {
      if (u.is_active_at(constraints[i]))
        result[i] = true;
      else
        result[i] = false;
    }
  }

  /** @brief tests for consistency of a potentially new constraint. */
  inline bool is_consistent(const Constraint & c) const {
    bool inconsistent = true;
    for (
          auto riter = rays.begin();
          inconsistent and riter != rays.end();
          ++riter
    ) {
      if (((*riter) * c) > 0)
        inconsistent = false;
    }
    return not inconsistent;
  }

  ///@}

  /** @name Modification */
  ///@{

  /**
    @brief Adds the indicated constraints (plural!) and re-computes the skeleton.

    \return `true` if and only if the new constraints are consistent with the
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
  virtual bool solve(const vector<Constraint> &) override;

  /**
    @brief Adds the indicated constraint (singular!) and re-computes the skeleton.

    \return `true` if and only if the new constraint is consistent with the
      current constraints

    @warning Checking the return value is crucial!
      If the function returns `false`, you have an inconsistent system!
      While the <i>present</i> cone will remain consistent,
      <b>the function will not roll back previous changes you have made</b>,
      so if you want to iterate again,
      your best bet is to copy the skeleton, and try that copy.
      Accept the new constraints only if that copy succeeds,
      in which case, you might as well discard the original, and keep the copy.
  */
  virtual bool solve(const Constraint &) override;

  /**
    @brief Re-computes the edges in the skeleton
    using Zolotych's `GraphAdj` algorithm and returns the result.
  */
  set<Edge> adjacencies_by_graphs(const set<Ray> &);

  /**
    @brief Assignment operator; empties current set & copies from other.
  */
  Skeleton & operator=(const Skeleton &);

  ///@}

  /** @name I/O */
  ///@{

  friend ostream & operator<<(ostream & os, const Skeleton & s);

  ///@}


private:

  int dim; /**< dimension of skeleton */

  set<Edge> edges; /**< edges defining skeleton */

  vector<Constraint> constraints; /**< constraints defining skeleton */

};

/**
  @brief prints out the constraints, then the rays, then the edges of @p s.
  @param os output stream to print to
  @param s skeleton to print
  @return the output stream
*/
ostream & operator<<(ostream & os, const Skeleton & s);

}

#endif
