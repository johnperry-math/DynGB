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
* Foobar is distributed in the hope that it will be useful,                   *
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
#include "skeleton.hpp"
#include "glpk_solver.hpp"
#include "ppl_solver.hpp"

#include <map>
#include <vector>
#include <iostream>

using std::map; using std::vector; using std::cout; using std::endl;

/**
  @ingroup GBComputation
  @author John Perry
  @date 2014-2017
  @details Pass this value to SelectMonomial. The methods are:
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
enum class DynamicHeuristic
{
  MIN_HEURISTIC = 0,
  ORD_HILBERT_THEN_LEX, ORD_HILBERT_THEN_DEG, GRAD_HILB_THEN_DEG,
  MIN_CRIT_PAIRS, GRAD_MIN_CRIT_PAIRS,
  BETTI_HILBERT_DEG, BIG_BETTI_HILBERT_DEG, GRAD_BETTI_HILBERT_DEG,
  SMOOTHEST_DEGREES, LARGEST_MAX_COMPONENT,
  DEG_THEN_GRAD_HILB, DEG_THEN_ORD_HILBERT,
  MAX_HEURISTIC
};

/**
  @ingroup GBComputation
  @author John Perry
  @date 2017
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
int hilbertCmp(
  const Dense_Univariate_Integer_Polynomial  & hn1,
  const Dense_Univariate_Rational_Polynomial & hp1,
  const Dense_Univariate_Integer_Polynomial  & hn2,
  const Dense_Univariate_Rational_Polynomial & hp2
);

/**
  @ingroup GBComputation
  @author John Perry
  @date 2014
  Used to associate a potential leading power product
  with the resulting monomial ideal if it were chosen
  as the actual leading power product.
*/
class PPWithIdeal {
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
    PPWithIdeal(
        Monomial u, const list<Monomial> & F, ::ray & w,
        const list<Critical_Pair_Dynamic *> & P,
        const Dense_Univariate_Integer_Polynomial * h = nullptr
    ) : t(u), ordering(w), num_new_pairs(-1), min_deg(-1), max_deg(-1),
        I(u.num_vars(), F, h), pairs(P)
    {
      I.add_generator(u);
    };
    /** @brief copy constructor */
    PPWithIdeal(const PPWithIdeal & PI)
    : t(PI.t), ordering(PI.ordering), num_new_pairs(PI.num_new_pairs),
      min_deg(PI.min_deg), max_deg(PI.max_deg), I(PI.I), pairs(PI.pairs)
    { }
    ///@}
    /** @name Destruction */
    ///@{
    /** @brief does nothing the default wouldn't do */
    ~PPWithIdeal() { }
    ///@}
    /** @name Basic properties */
    ///@{
    /** @brief the leading monomial being added to the ideal */
    inline const Monomial & getPP() const { return t; };
    /** @brief the old ideal of leading monomials */
    inline const Monomial_Ideal & getIdeal() const { return I; };
    /** @brief the current monomial ordering */
    inline ::ray getOrdering() { return ordering; };
    /**
      @brief the incremental Betti numbers obtained by adding the monomial to
      the ideal
    */
    inline const map<DEG_TYPE, unsigned long> & getIncBetti(bool graded = false)
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
    inline Dense_Univariate_Integer_Polynomial * getHilbertNumerator(
        bool graded = false
    ) {
      if (graded)
        return I.hilbert_numerator(ordering.weights());
      else 
        return I.hilbert_numerator();
    }
    /**
      @brief the Hilbert numerator obtained by adding the monomial to the ideal
      (numerator <i>is</i> reduced)
    */
    inline Dense_Univariate_Integer_Polynomial * getHilbertReducedNumerator(
        bool graded = false
    ) {
      if (graded)
        return I.reduced_hilbert_numerator(ordering.weights());
      else
        return I.reduced_hilbert_numerator();
    }
    /**
      @brief the Hilbert polynomial obtained by adding the monomial to the ideal
    */
    inline Dense_Univariate_Rational_Polynomial * getHilbertPolynomial() {
      return I.hilbert_poly();
    };
    /**
      @brief estimate of the number of new critical pairs generated by adding
      the monomial to the ideal
    */
    inline int howManyNewPairs() const { return num_new_pairs; }
    /**
      @brief the degree of the new critical pairs generated by adding the
      monomial to the ideal
    */
    inline int degOfNewPairs() const { return min_deg; }
    /**
      @brief computes the difference in degree between the first and last
      monomials of the ideal
    */
    inline int getDifferenceInDegree()
    {
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
    void computeNumberNewPairs();
    ///@}
  protected:
    /** @brief the last monomial added to \f$I\f$ */
    Monomial t;
    /** @brief the ideal of leading terms */
    Monomial_Ideal I;
    /**
      @brief the list of critical pairs of \f$I\f$ at this point in the algorithm
    */
    const list<Critical_Pair_Dynamic *> & pairs;
    /** @brief the current ordering of the Gr&ouml;bner basis computation */
    ::ray ordering;
    /** @brief estimate of number of new pairs */
    int num_new_pairs;
    /** @brief minimum weighted degree of monomials in ideal */
    int min_deg;
    /** @brief minimum weighted degree of monomials in ideal */
    int max_deg;
};

/**
  @ingroup GBComputation
  @author John Perry
  @date 2014
  Compute the compatible leading monomials of a polynomial.
  \param currentLPP the current leading power product of the polynomial
  \param allPPs set of all power products of the polynomial
  \param bndrys boundary (or, &ldquo;corner&rdquo;) vectors of current cone
  \param result set of power products of the polynomial compatible with `bndrys`
  \param skel existing skeleton that defines currently-compatible orderings
  \param boundary_mons boundary monomials (no apparent purpose at the moment)
*/
void compatiblePP(
  Monomial currentLPP,            // the current LPP
  const set<Monomial> &allPPs,    // the monomials to consider;
                              // some will be removed
  const set<::ray> &bndrys,     // known boundary vectors
  set<Monomial> &result,          // returned as PPs for Hilbert function
                              // ("easy" (& efficient?) to extract exps
  set<Monomial> &boundary_mons,   // boundary monomials
  LP_Solver *skel              // used for alternate refinement
);

/**
  @ingroup GBComputation
  @author John Perry
  @date 2014
  Verifies that the leading power products of the current basis
  remain compatible with the proposed refinement of ordering.
  \param skel the skeleton corresponding to the choice of \f$w\f$
  \param currentPolys the current basis
  @return @p true if and only if we were able to make this skeleton
    (and thus the corresponding choice of a leading monomial) consistent
    with the existing basis.
*/
bool verifyAndModifyIfNecessary(
  LP_Solver *skel,
  const list<Abstract_Polynomial *> &currentPolys
);

/**
  @ingroup GBComputation
  @author John Perry
  @date 2014
  Create constraints for a candidate LPP.
  \param pp_I pair of PP with the ideal it would have.
  \param monomialsForComparison monomials used to generate constraints with LPP
  \param result the new constraints
*/
void ConstraintsForNewPP(
  const PPWithIdeal &pp_I,
  const set<Monomial> &monomialsForComparison,
  vector<constraint> &result
);

/**
  @ingroup GBComputation
  @author John Perry
  @date 2014
  @brief Selects a leading power product for a polynomial.
  @details Applies a particular DynamicHeuristic
  (default is `ORD_HILBERT_THEN_DEG`)
  and ensures compatibility with previous choices of LPP for other monomials.
  \param r the polynomial in need of a new choice of LPP
  \param CurrentLPPs the current choices of LPPs for `CurrentPolys`
  \param CurrentPolys the current basis of the ideal
  @param critpairs the current list of critical pairs
  \param currSkel the current skeleton, corresponding to the choices `CurrentLPPs`
  @param ordering_changed whether the monomial selected changes the ordering
  \param method the method to apply; see `DynamicHeuristic`
  @param current_hilbert_numerator Hilbert numerator for CurrentLPPs (changes to
      match new monomial, hence the double reference)
*/
void SelectMonomial(
    Abstract_Polynomial * r,                          // changes
    list<Monomial> &CurrentLPPs,      // changes
    Dense_Univariate_Integer_Polynomial ** current_hilbert_numerator,
    const list<Abstract_Polynomial *> &CurrentPolys,
    const list<Critical_Pair_Dynamic *> & critpairs,
    LP_Solver * currSkel,                // possibly changes
    bool &ordering_changed,
    DynamicHeuristic method = DynamicHeuristic::ORD_HILBERT_THEN_DEG
);

#endif