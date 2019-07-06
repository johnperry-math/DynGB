#ifndef __DYNAMIC_ENGINE_H__
#define __DYNAMIC_ENGINE_H__

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

#include <cstddef>

#include "monomial.hpp"
#include "monomial_ideal.hpp"
#include "hilbert_functions.hpp"
#include "polynomial.hpp"
#include "critical_pair.hpp"
#include "lp_solver.hpp"
using LP_Solvers::LP_Solver;
using LP_Solvers::Ray; using LP_Solvers::Constraint;
#include "skeleton.hpp"
using LP_Solvers::Skeleton;
#include "glpk_solver.hpp"
using LP_Solvers::GLPK_Solver;
#include "ppl_solver.hpp"
using LP_Solvers::PPL_Solver;

#include <map>
using std::map;
#include <vector>
using std::vector;
#include <iostream>
using std::cout; using std::endl;

/**
  @namespace Dynamic_Engine
  @brief Namespace for functions directly related to dynamic computation
      of a Gr&ouml;bner basis.
*/
namespace Dynamic_Engine {

/**
  @ingroup GBComputation
  @author John Perry
  @date 2014-2017
  @brief values for the heuristic, as in <c>select_monomial()</c>
  @details Pass this value to <c>select_monomial()</c>. The methods are:
    -# `ORD_HILBERT_THEN_DEG` considers the standard Hilbert polynomial,
       breaking ties by the standard Hilbert series,
       breaking remaining ties by the (weighted) degree.
    -# `DEG_THEN_ORD_HILBERT` considers the degree of the leading power products,
       breaking ties by the standard Hilbert polynomial,
       breaking remaining ties by the standard Hilbert series.
    -# `GRAD_HILB_THEN_DEG` considers the standard Hilbert polynomial,
       breaking ties by the graded Hilbert series,
       breaking remaining ties by the (weighted) degree.
       (The graded Hilbert polynomial is actuall a quasi-polynomial,
       and we do not yet have a criterion to apply to this.)
    -# `DEG_THEN_GRAD_HILB` considers the degree of the leading power products,
       breaking ties by the standard Hilbert polynomial,
       breaking remaining ties by the graded Hilbert series.
    -# `MIN_CRIT_PAIRS` tries to minimize the number of new critical pairs.
    -# `BETTI_HILBERT_DEG` tries to obtain a set of incremental Betti numbers
       with the largest values, comparing from the bottom and stopping when
       not equal
    -# `BIG_BETTI_HILBERT_DEG` tries to obtain a set of incremental Betti
       with the largest values, comparing from the top and stopping when not
       equal
    -# `GRAD_BETTI_HILBERT_DEG` tries to obtain a set of incremental Betti
       numbers with the largest values, comparing from the bottom and stopping
       when not equal
*/
enum class Dynamic_Heuristic
{
  MIN_HEURISTIC = 0,
  ORD_HILBERT_THEN_LEX, ORD_HILBERT_THEN_DEG, GRAD_HILB_THEN_DEG,
  MIN_CRIT_PAIRS, GRAD_MIN_CRIT_PAIRS,
  BETTI_HILBERT_DEG, BIG_BETTI_HILBERT_DEG, GRAD_BETTI_HILBERT_DEG,
  SMOOTHEST_DEGREES, LARGEST_MAX_COMPONENT,
  DEG_THEN_GRAD_HILB, DEG_THEN_ORD_HILBERT,
  EVIL_RANDOM,
  MAX_HEURISTIC
};

/**
  @ingroup GBComputation
  @author John Perry
  @date 2017
  @brief compares two Hilbert functions, using polynomial and numerator
  @return -1 if the first Hilbert function is smaller than the second,
    0 if they are indistinguishable, and
    1 if the second is smaller than the first. For details, see
    @cite CaboaraDynAlg.
  @param hn1 unreduced Hilbert numerator for first function
  @param hn2 unreduced Hilbert numerator for second function
  @param hp1 Hilbert polynomial for first function
  @param hp2 Hilbert polynomial for second function
  @warning None of these should be @c nullptr.
*/
int hilbert_cmp(
  const Dense_Univariate_Integer_Polynomial  & hn1,
  const Dense_Univariate_Rational_Polynomial & hp1,
  const Dense_Univariate_Integer_Polynomial  & hn2,
  const Dense_Univariate_Rational_Polynomial & hp2
);

/**
  @ingroup GBComputation
  @author John Perry
  @date 2014
  @brief Used to associate a potential leading power product
  with the resulting monomial ideal if it were chosen
  as the actual leading power product.
*/
class PP_With_Ideal {
  public:
    /** @name Construction */
    ///@{
    /**
      @brief Construct a monomial/ideal pair.
      @param u proposed new Monomial for ideal
      @param F current ideal
      @param w current monomial ordering
      @param P current list of critical pairs
      @param h unreduced Hilbert numerator of @p F (does not verify correctness)
    */
    PP_With_Ideal(
        const Monomial & u, const list<Monomial> & F, const Ray & w,
        const list<Critical_Pair_Dynamic *> & P,
        const Dense_Univariate_Integer_Polynomial * h = nullptr
    ) : t(u), ordering(w), num_new_pairs(-1), min_deg(-1), max_deg(-1),
        I(u.num_vars(), F, h), pairs(P)
    {
      I.add_generator(u);
    };
    /** @brief copy constructor */
    PP_With_Ideal(const PP_With_Ideal & PI)
    : t(PI.t), ordering(PI.ordering), num_new_pairs(PI.num_new_pairs),
      min_deg(PI.min_deg), max_deg(PI.max_deg), I(PI.I), pairs(PI.pairs)
    { }
    ///@}
    /** @name Destruction */
    ///@{
    /** @brief does nothing the default wouldn't do */
    ~PP_With_Ideal() { }
    ///@}
    /** @name Basic properties */
    ///@{
    /** @brief the leading monomial being added to the ideal */
    inline const Monomial & get_pp() const { return t; };
    /** @brief the old ideal of leading monomials */
    inline const Monomial_Ideal & get_ideal() const { return I; };
    /** @brief the current monomial ordering */
    inline const Ray get_ordering() const { return ordering; };
    /**
      @brief the incremental Betti numbers obtained by adding the monomial to
      the ideal
    */
    inline const map<DEG_TYPE, unsigned long> & get_inc_betti(bool graded = false)
    {
      if (graded)
        return I.inc_betti(ordering.weights());
      else
        return I.inc_betti();
    }
    /**
      @brief the Hilbert numerator obtained by adding the monomial to the ideal
      (numerator is <i>not</i> reduced)
    */
    inline Dense_Univariate_Integer_Polynomial * get_hilbert_numerator(
        bool graded = false
    ) const {
      if (graded)
        return I.hilbert_numerator(ordering.weights());
      else 
        return I.hilbert_numerator();
    }
    /**
      @brief the Hilbert numerator obtained by adding the monomial to the ideal
      (numerator <i>is</i> reduced)
    */
    inline Dense_Univariate_Integer_Polynomial * get_hilbert_reduced_numerator(
        bool graded = false
    ) const {
      if (graded)
        return I.reduced_hilbert_numerator(ordering.weights());
      else
        return I.reduced_hilbert_numerator();
    }
    /**
      @brief the Hilbert polynomial obtained by adding the monomial to the ideal
    */
    inline Dense_Univariate_Rational_Polynomial * get_hilbert_polynomial() const {
      return I.hilbert_poly();
    };
    /**
      @brief estimate of the number of new critical pairs generated by adding
      the monomial to the ideal
    */
    inline int how_many_new_pairs() const { return num_new_pairs; }
    /**
      @brief the degree of the new critical pairs generated by adding the
      monomial to the ideal
    */
    inline int deg_of_new_pairs() const { return min_deg; }
    /**
      @brief computes the difference in degree between the first and last
      monomials of the ideal
    */
    inline int get_difference_in_degree() const {
      if (min_deg < 0)
      {
        const list<Monomial> & G = I.generators();
        min_deg = max_deg = G.front().weighted_degree(ordering.weights());
        for (const Monomial & g : G) {
          int d = g.weighted_degree(ordering.weights());
          if (d < min_deg) min_deg = d;
          if (d > max_deg) max_deg = d;
        }
      }
      return max_deg - min_deg;
    }
    ///@}
    /** @name Modification */
    ///@{
    /** @brief Computes the number of critical pairs the monomial would add */
    void compute_number_new_pairs() const;
    /**
      @brief assigns a value to the hilbert numerator when it&rsquo;s
        already known
    */
    void set_hilbert_numerator(Dense_Univariate_Integer_Polynomial *h) {
      I.set_hilbert_numerator(h);
    }
    ///@}
  protected:
    /** @brief the last monomial added to @f$I@f$ */
    Monomial t;
    /** @brief the ideal of leading terms */
    mutable Monomial_Ideal I;
    /**
      @brief the list of critical pairs of @f$I@f$ at this point in the algorithm
    */
    const list<Critical_Pair_Dynamic *> & pairs;
    /** @brief the current ordering of the Gr&ouml;bner basis computation */
    Ray ordering;
    /** @brief estimate of number of new pairs */
    mutable int num_new_pairs;
    /** @brief minimum weighted degree of monomials in ideal */
    mutable int min_deg;
    /** @brief minimum weighted degree of monomials in ideal */
    mutable int max_deg;
};

/**
  @ingroup GBComputation
  @author John Perry
  @date 2014
  @brief Compute the compatible leading monomials of a polynomial.
  @param currentLPP the current leading power product of the polynomial
  @param allPPs set of all power products of the polynomial
  @param result set of power products of the polynomial compatible with `bndrys`
  @param skel existing skeleton that defines currently-compatible orderings
  @param boundary_mons boundary monomials (no apparent purpose at the moment)
*/
void compatible_pp(
  const Monomial & currentLPP,            // the current LPP
  const set<Monomial> &allPPs,    // the monomials to consider;
                                  // some will be removed
  set<Monomial> &result,          // returned as PPs for Hilbert function
                                  // ("easy" (& efficient?) to extract exps
  list<Monomial> &boundary_mons,  // boundary monomials
  const LP_Solver *skel                 // used for alternate refinement
);

/**
  @ingroup GBComputation
  @author John Perry
  @date 2014
  @brief Verifies that the leading power products of the current basis
  remain compatible with the proposed refinement of ordering.
  @param skel the skeleton corresponding to the choice of @f$w@f$
  @param currentPolys the current basis
  @return @p true if and only if we were able to make this skeleton
    (and thus the corresponding choice of a leading monomial) consistent
    with the existing basis.
*/
bool verify_and_modify_if_necessary(
  LP_Solver *skel,
  const list<Abstract_Polynomial *> &currentPolys
);

/**
  @ingroup GBComputation
  @author John Perry
  @date 2014
  @brief Create constraints for a candidate LPP.
  @param pp_I pair of PP with the ideal it would have.
  @param monomials_for_comparison monomials used to generate constraints with LPP
  @param result the new constraints
*/
void constraints_for_new_pp(
  const PP_With_Ideal & pp_I,
  const set<Monomial> & monomials_for_comparison,
  vector<Constraint> & result
);

/**
  @ingroup GBComputation
  @author John Perry
  @date 2014
  @brief Selects a leading power product for a polynomial.
  @details Applies a particular Dynamic_Heuristic
  (default is `ORD_HILBERT_THEN_DEG`)
  and ensures compatibility with previous choices of LPP for other monomials.
  @param r the polynomial in need of a new choice of LPP
  @param CurrentLPPs the current choices of LPPs for `CurrentPolys`
  @param CurrentPolys the current basis of the ideal
  @param critpairs the current list of critical pairs
  @param currSkel the current skeleton, corresponding to the choices `CurrentLPPs`
  @param ordering_changed whether the monomial selected changes the ordering
  @param method the method to apply; see `Dynamic_Heuristic`
  @param current_hilbert_numerator Hilbert numerator for CurrentLPPs (changes to
      match new monomial, hence the double reference)
*/
void select_monomial(
    Abstract_Polynomial * r,                          // changes
    list<Monomial> &CurrentLPPs,      // changes
    Dense_Univariate_Integer_Polynomial ** current_hilbert_numerator,
    const list<Abstract_Polynomial *> &CurrentPolys,
    const list<Critical_Pair_Dynamic *> & critpairs,
    LP_Solver * currSkel,                // possibly changes
    bool &ordering_changed,
    Dynamic_Heuristic method = Dynamic_Heuristic::ORD_HILBERT_THEN_DEG
);

/**
  @ingroup GBComputation
  @brief same as
    <c>select_monomial(Abstract_Polynomial *, list<Monomial> &, Dense_Univariate_Integer_Polynomial **, const list<Abstract_Polynomial *> &, const list<Critical_Pair_Dynamic *> &,LP_Solver *, bool &, Dynamic_Heuristic)</c>
    but with the monomials already extracted
  @param r the polynomial in need of a new choice of LPP;
      this will change to the monomials that are actually compatible
  @param currentLPP the current leading monomial
  @param CurrentLPPs the current choices of LPPs for `CurrentPolys`
  @param current_hilbert_numerator Hilbert numerator for CurrentLPPs (changes to
      match new monomial, hence the double reference)
  @param CurrentPolys the current basis of the ideal
  @param critpairs the current list of critical pairs
  @param currSkel the current skeleton, corresponding to the choices `CurrentLPPs`
  @param ordering_changed whether the monomial selected changes the ordering
  @param method the method to apply; see `Dynamic_Heuristic`
  @author John Perry
  @date 2017
*/
void select_monomial(
    const set<Monomial> & r,
    const Monomial & currentLPP,
    list<Monomial> &CurrentLPPs,      // changes
    Dense_Univariate_Integer_Polynomial ** current_hilbert_numerator,
    const list<Abstract_Polynomial *> &CurrentPolys,
    const list<Critical_Pair_Dynamic *> & critpairs,
    LP_Solver * currSkel,                // possibly changes
    bool &ordering_changed,
    Dynamic_Heuristic method = Dynamic_Heuristic::ORD_HILBERT_THEN_DEG
);

/**
  @brief compares the desirability of two potential expansions
    of a monomial ideal, selecting one at random
  @param a,b ideals to compare
  @return @c true iff the first ideal is considered smaller
  @details This heuristic is useful to show that on sufficiently complex
    systems the other heuristics actually do something intelligent.
    On some simple systems this heuristic can sometimes beat the others,
    because there are not very many bases possible, and the greedy algorithms
    inevitable fail to find the "best." Once you get out of those, however,
    this heuristic inevitably turns disastrous.
*/
bool less_by_random(const PP_With_Ideal & a, const PP_With_Ideal & b);

/**
  @brief compares the desirability of two potential expansions
    of a monomial ideal via the Hilbert heuristic
  @param a,b ideals to compare
  @return @c true iff the first ideal is considered smaller
  @details The Hilbert heuristic compares first the Hilbert polynomial,
    then the Hilbert numerator. See @cite CaboaraDynAlg.
    In the case of a tie, it prefers the monomial that is smaller
    with respect to the current ordering.
*/
bool less_by_hilbert (PP_With_Ideal &a, PP_With_Ideal &b);

/**
  @brief compares the desirability of two potential expansions
    of a monomial ideal by a &ldquo;smoothest degree&rdquo; heuristic
  @param a,b ideals to compare
  @return @c true iff the first ideal is considered smaller
  @details The smoothest degree heuristic looks for a minimal change
    in weighted degree.
    It aims for polynomials that are &ldquo;close to homogeneous.&rdquo;
*/
bool less_by_smoothest_degrees (PP_With_Ideal &a, PP_With_Ideal &b);

/**
  @brief compares the desirability of two potential expansions
    of a monomial ideal via the largest maximal component heuristic
  @param a,b ideals to compare
  @return @c true iff the first ideal is considered smaller
  @details The largest maximal component heuristic aims to maximize
    the subset of the support that has largest degree.
    This aims for polynomials that are &ldquo;close to regular.&rdquo;
*/
bool less_by_largest_max_component (PP_With_Ideal &a, PP_With_Ideal &b);

/**
  @brief compares the desirability of two potential expansions
    of a monomial ideal via the basic Betti heuristic
  @param a,b ideals to compare
  @return @c true iff the first ideal is considered smaller
  @details The Betti heuristic compares the number of critical pairs.
    In the case of a tie, it prefers the monomial that is smaller
    with respect to the current ordering.
*/
bool less_by_num_crit_pairs (PP_With_Ideal &a, PP_With_Ideal &b);

/**
  @brief compares the desirability of two potential expansions
    of a monomial ideal via the Hilbert+degree heuristic
  @param a,b ideals to compare
  @return @c true iff the first ideal is considered smaller
  @details The Hilbert heuristic compares first the Hilbert polynomial,
    then the Hilbert numerator. See @cite CaboaraDynAlg.
    It tries to break ties by total degree of the monomials.
    In the case of a tie, it prefers the monomial that is smaller
    with respect to the current ordering.
*/
bool less_by_hilbert_then_degree(PP_With_Ideal &a, PP_With_Ideal &b);

/**
  @brief compares the desirability of two potential expansions
    of a monomial ideal via the graded Hilbert+weighted degree heuristic
  @param a,b ideals to compare
  @return @c true iff the first ideal is considered smaller
  @details The Hilbert heuristic compares first the Hilbert polynomial,
    then the graded Hilbert numerator.
    It tries to break ties using the weighted degree of the monomials.
    In the case of a tie, it prefers the monomial that is smaller
    with respect to the current ordering.
*/
bool less_by_grad_hilbert_then_degree(PP_With_Ideal &a, PP_With_Ideal &b);

/**
  @brief compares the desirability of two potential expansions
    of a monomial ideal via the degree+Hilbert heuristic
  @param a,b ideals to compare
  @return @c true iff the first ideal is considered smaller
  @details This heuristic compares first the total degrees of the monomials.
    In case of a tie, it applies the Hilbert heuristic,
    which compares first the Hilbert polynomial,
    then the Hilbert numerator. See @cite CaboaraDynAlg.
    In the case of a tie, it prefers the monomial that is smaller
    with respect to the current ordering.
*/
bool less_by_degree_then_hilbert(PP_With_Ideal &a, PP_With_Ideal &b);

/**
  @brief compares the desirability of two potential expansions
    of a monomial ideal via the weighted degree+Hilbert heuristic
  @param a,b ideals to compare
  @return @c true iff the first ideal is considered smaller
  @details This heuristic first compares the monomials&rsquo; weighted degrees.
    In case of a tie, it applies the Hilbert heuristic,
    which compares first the Hilbert polynomial,
    then the Hilbert numerator. See @cite CaboaraDynAlg.
    In the case of a tie, it prefers the monomial that is smaller
    with respect to the current ordering.
*/
bool less_by_wdegree_then_hilbert(PP_With_Ideal &a, PP_With_Ideal &b);

/**
  @brief compares the desirability of two potential expansions
    of a monomial ideal via the weighted degree+graded Hilbert heuristic
  @param a,b ideals to compare
  @return @c true iff the first ideal is considered smaller
  @details This heuristic first compares the monomials&rsquo; weighted degrees.
    In case of a tie, it applies the Hilbert heuristic,
    which compares first the Hilbert polynomial,
    then the Hilbert numerator. See @cite CaboaraDynAlg.
    In the case of a tie, it prefers the monomial that is smaller
    with respect to the current ordering.
*/
bool less_by_degree_then_grad_hilbert(PP_With_Ideal &a, PP_With_Ideal &b);

/**
  @brief compares the desirability of two potential expansions
    of a monomial ideal via the standard graded Betti heuristic
  @param a,b ideals to compare
  @return @c true iff the first ideal is considered smaller
  @details The graded Betti heuristic compares the Betti numbers
    from the standard grading, breaking ties with the Hilbert+degree heuristic.
    For the Betti part, it proceeds through the list of critical pairs
    at each degree, quitting when it either runs out or finds a difference.
    In this case it prefers the ideal that has more pairs at lower degree.
    In the case of a tie, it prefers the monomial that is smaller
    with respect to the current ordering.
*/
bool less_by_betti (PP_With_Ideal & a, PP_With_Ideal & b);

/**
  @brief compares the desirability of two potential expansions
    of a monomial ideal via the &ldquo;Big Betti&rdquo; heuristic
  @param a,b ideals to compare
  @return @c true iff the first ideal is considered smaller
  @details The Hilbert heuristic compares first the Hilbert polynomial,
    then the Hilbert numerator. See @cite CaboaraDynAlg.
    In the case of a tie, it prefers the monomial that is smaller
    with respect to the current ordering.
*/
bool less_by_big_betti (PP_With_Ideal & a, PP_With_Ideal & b);

/**
  @brief the graded Betti heuristic with a non-standard grading
  @param a,b ideals to compare
  @return @c true iff the first ideal is considered smaller
*/
bool less_by_grad_betti (PP_With_Ideal & a, PP_With_Ideal & b);


}

#endif