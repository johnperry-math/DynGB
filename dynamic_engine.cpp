#ifndef __DYNAMIC_ENGINE_C__
#define __DYNAMIC_ENGINE_C__

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

#include <iostream>
#include <list>

using std::cout; using std::endl;
using std::list;

#include "betti.hpp"
#include "monomial.hpp"
#include "polynomial.hpp"
#include "dynamic_engine.hpp"

int hilbertCmp(
  const Dense_Univariate_Integer_Polynomial & hn1,
  const Dense_Univariate_Rational_Polynomial & hp1,
  const Dense_Univariate_Integer_Polynomial & hn2,
  const Dense_Univariate_Rational_Polynomial & hp2
) {
  int result = 0;
  Dense_Univariate_Rational_Polynomial HPdiff(hp1);
  HPdiff -= hp2;
  if (not HPdiff.is_zero()) {
    if (HPdiff.numerator(HPdiff.degree()) < 0) result = -1;
    else result = 1;
  } else { // use Hilbert series
    DEG_TYPE i = 0;
    for (/* */; i <= hn1.degree() and i <= hn2.degree() and hn1[i] == hn2[i]; i++)
    { /* taken care of in loop */ }
    if (i > hn1.degree()) {
      if (i > hn2.degree()) // functions are equal
        result = 0;
      else
        result = 1;
    } else {
      if (i > hn2.degree()) result = 1;
      else {
        if (hn1[i] < hn2[i]) result = -1;
        else result = 1;
      }
    }
  }
  return result;
}

bool LessByHilbert (PPWithIdeal &a, PPWithIdeal &b)
{
  //cout << "Less by Hilbert then Lex\n";
  bool result;
  int n = a.getPP().num_vars();
  // first check the coefficients of the Hilbert polynomial
  int hilcheck = hilbertCmp(
    *a.getHilbertNumerator(), *a.getHilbertPolynomial(),
    *b.getHilbertNumerator(), *b.getHilbertPolynomial()
  );
  if (hilcheck == 0) // the numerators are equal; break tie via current ordering
    result = (a.getPP() < b.getPP());
  else {
    if (hilcheck == -1) result = true;
    else result = false;
  }
  //cout << "\tfirst less than second? " << result << endl;
  return result;
};

bool LessBySmoothestDegrees (PPWithIdeal &a, PPWithIdeal &b)
{
  return a.getDifferenceInDegree() < b.getDifferenceInDegree();
};

bool LessByLargestMaxComponent (PPWithIdeal &a, PPWithIdeal &b)
{
  return a.getDifferenceInDegree() < b.getDifferenceInDegree();
}

void PPWithIdeal::computeNumberNewPairs()
{
  int n = t.num_vars();
  int m = I.size();
  const list<Monomial> T = I.generators();
  const vector<Monomial> I(T.begin(), T.end());
  num_new_pairs = min_deg = 0;
  bool * keepers = new bool [m];
  for (int i = 0; i < m; ++i) keepers[i] = true;
  // first main loop: apply Buchberger's lcm criterion to new pairs
  for (int i = 0; i < m; ++i)
  {
    // if gcd(t,Ii) == 1 then skip i for the time being
    if (not t.is_coprime(I[i]))
    {
      bool has_divisor = false;
      for (int j=0; (not has_divisor) and j < m; ++j)
      {
        if (i != j and keepers[j])
        {
          // if some j satisfies lcm(t,Ij) | lcm(t,Ii) then do not count i
          has_divisor = true;
          for (int k=1; has_divisor and k <= n; ++k)
            if // deg(lcm(t,lm(Ii))) <  deg(lcm(t,lm(Ij))) ?
               (((t[k] > I[i][k]) ? t[k] : I[i][k])
                 < ((t[k] > I[j][k]) ? t[k] : I[j][k]))
              has_divisor = false;
        }
      }
      if (has_divisor) keepers[i] = false;
    }
  }
  // second main loop: apply Buchberger's gcd criterion to new pairs, count survivors
  for (int i = 0; i < m; ++i)
  {
    if (keepers[i] and not t.is_coprime(I[i]))
    {
      int new_deg = 0;
      // determine deg(lcm(t,Si))
      for (int k=1; k <= n; ++k)
        new_deg += (t[k] > I[i][k]) ? t[k] : I[i][k];
      if (min_deg == 0 or min_deg > new_deg)
      {
        min_deg = new_deg; num_new_pairs = 1;
      }
      else if (min_deg == new_deg)
      {
        ++num_new_pairs;
        //cout << '\t'; pWrite(pHead(strat->S[i]));
      }
    }
  }
  delete [] keepers;
  // cout << "we get " << num_new_pairs << " from "; pWrite(t);
  // third main loop: apply Buchberger's lcm criterion to old pairs, UNLESS
  // all three lcm's are equal
  for (Critical_Pair_Dynamic * p : pairs) {
      const Monomial & u = p->lcm();
      int new_deg = 0;
      for (int k=1; k <= n; ++k) new_deg += u[k];
      // for (int k=1; k <= n; ++k) new_deg += pGetExp(u,k) * Rx->wvhdl[0][k];
      // no point continuing if it wouldn't change min_deg
      if (min_deg == 0 or new_deg <= min_deg)
      {
        // see if new poly divides this one
        bool has_divisor = true;
        // see if Buchberger triple has same lcm
        bool all_equal = true;
        const Monomial & p1 = p->first()->leading_monomial();
        const Monomial & p2 = (p->second() == nullptr) ?
            p->first()->leading_monomial()
          : p->second()->leading_monomial();
        for (int k=1; has_divisor and all_equal and k <= n; ++k)
        {
          if (t[k] > u[k]) has_divisor = false;
          else
          {
            // check lcm(t,lm(p1)) == lcm(t,lm(p2)) == lcm(lm(p1),lm(p2)) in xk
            int a = (t[k] > p1[k]) ? t[k] : p1[k];
            int b = (t[k] > p2[k]) ? t[k] : p2[k];
            int c = (p1[k] > p2[k]) ? p1[k] : p2[k];
            all_equal = (a == c) and (b == c);
          }
        }
        if (not has_divisor or all_equal)
        {
          if (min_deg == 0 or min_deg > new_deg)
          {
            min_deg = new_deg; num_new_pairs = 1;
          }
          else // the only reason we'd be here is if min_deg == new_deg
          {
            ++num_new_pairs;
            //cout << '\t'; pWrite(u);
          }
        }
      }
  }
  // cout << " which makes " << num_new_pairs << " pairs total at degree " << min_deg << ".\n";
}

bool LessByNumCritPairs (PPWithIdeal &a, PPWithIdeal &b)
{
  bool result;
  NVAR_TYPE n = a.getIdeal().number_of_variables();
  // first check if the number of critical pairs has been computed
  if (a.howManyNewPairs() < 0) a.computeNumberNewPairs();
  if (b.howManyNewPairs() < 0) b.computeNumberNewPairs();
  /*if (a.degOfNewPairs() < b.degOfNewPairs())
    result = true;
  else if (a.degOfNewPairs() > b.degOfNewPairs())
    result = false;*/
  // at this point, the degrees of the new pairs will be equal
  /*else*/ if (a.howManyNewPairs() > b.howManyNewPairs())
    result = true;
  else if (a.howManyNewPairs() < b.howManyNewPairs())
    result = false;
  else // the numerators are equal; break tie via monomial ordering
    result = a.getPP() < b.getPP();
  //cout << "\tfirst less than second? " << result << endl;
  return result;
};

bool LessByHilbertThenDegree(PPWithIdeal &a, PPWithIdeal &b)
{
  //cout << "Less by Hilbert then Deg\n";
  bool result;
  int n = a.getIdeal().number_of_variables();
  // first check the coefficients of the Hilbert polynomial
  int hilcheck = hilbertCmp(
      *a.getHilbertNumerator(), *a.getHilbertPolynomial(),
      *b.getHilbertNumerator(), *b.getHilbertPolynomial()
  );
  if (hilcheck == -1) result = true;
  else if (hilcheck == 1) result = false;
  else { // the numerators are equal; break tie via monomial degree
    if (a.getPP().total_degree() < b.getPP().total_degree())
      result = true;
    else if (a.getPP().total_degree() > b.getPP().total_degree())
      result = false;
    else
      result = (a.getPP() < b.getPP());
  }
  return result;
};

bool LessByGradHilbertThenDegree(PPWithIdeal &a, PPWithIdeal &b) {
  bool result;
  int n = a.getIdeal().number_of_variables();
  // first check the coefficients of the Hilbert polynomial
  Dense_Univariate_Rational_Polynomial * hp1 = a.getHilbertPolynomial();
  Dense_Univariate_Rational_Polynomial * hp2 = b.getHilbertPolynomial();
  int hilcheck = hilbertCmp(
    *a.getHilbertNumerator(true), *hp1,
    *b.getHilbertNumerator(true), *hp2
  );
  if (hilcheck == -1) result = true;
  else if (hilcheck == 1) result = false;
  else if (hilcheck == 0) { // break tie via monomial degree
    if (a.getPP().weighted_degree(a.getOrdering().weights())
          < b.getPP().weighted_degree(b.getOrdering().weights()))
      result = true;
    else if (a.getPP().weighted_degree(a.getOrdering().weights())
                > b.getPP().weighted_degree(b.getOrdering().weights()))
      result = false;
    else
      result = a.getPP() < b.getPP();
  }
  return result;
};

bool LessByDegreeThenHilbert(PPWithIdeal &a, PPWithIdeal &b)
{
  //cout << "Less by Deg then Hilbert\n";
  bool result;
  int n = a.getPP().num_vars();
  // first check the weighted degree
  if (a.getPP().total_degree() < b.getPP().total_degree())
    result = true;
  else if (a.getPP().total_degree() > b.getPP().total_degree())
    result = false;
  else {
    // now check the coefficients of the Hilbert polynomial
    Dense_Univariate_Rational_Polynomial HPdiff(*(a.getHilbertPolynomial()));
    HPdiff -= *(b.getHilbertPolynomial());
    if (not HPdiff.is_zero())
      result = (HPdiff.numerator(HPdiff.degree()) < 0);
    else // use Hilbert series
    {
      Dense_Univariate_Integer_Polynomial * h1 = a.getHilbertNumerator();
      Dense_Univariate_Integer_Polynomial * h2 = b.getHilbertNumerator();
      DEG_TYPE i = 0;
      for ( /* already initialized */ ;
           i <= h1->degree() and i <= h2->degree() and (*h1)[i] == (*h2)[i];
           i++)
      { /* taken care of in loop */ }
      if (i > h1->degree())
      {
        if (i > h2->degree())
        {
          int i = 0;
          while (i < n and a.getPP()[i] == b.getPP()[i]) ++i;
          if (i == n) result = false;
          else
            result = (a.getPP()[i] < b.getPP()[i]);
        }
        else
          result = true;
      }
      else
      {
        if (i > h2->degree()) result = false;
        else result = (*h1)[i] < (*h2)[i];
      }
    }
  }
  return result;
};

bool LessByWDegreeThenHilbert(PPWithIdeal &a, PPWithIdeal &b)
{
  //cout << "Less by Deg then Hilbert\n";
  bool result;
  int n = a.getPP().num_vars();
  // first check the weighted degree
  if (a.getPP().weighted_degree(a.getOrdering().weights())
        < b.getPP().weighted_degree(b.getOrdering().weights()))
    result = true;
  else if (a.getPP().weighted_degree(a.getOrdering().weights())
        > b.getPP().weighted_degree(b.getOrdering().weights()))
    result = false;
  else {
    // now check the coefficients of the Hilbert polynomial
    Dense_Univariate_Rational_Polynomial HPdiff
        = a.getHilbertPolynomial() - b.getHilbertPolynomial();
    if (not HPdiff.is_zero())
      result = (HPdiff.numerator(HPdiff.degree()) >= 0);
    else // use Hilbert series
    {
      Dense_Univariate_Integer_Polynomial * h1 = a.getHilbertNumerator();
      Dense_Univariate_Integer_Polynomial * h2 = b.getHilbertNumerator();
      DEG_TYPE i = 0;
      for ( /* already initialized */ ;
           i <= h1->degree() and i <= h2->degree() and (*h1)[i] == (*h2)[i];
           i++)
      { /* taken care of in loop */ }
      if (i > h1->degree())
      {
        if (i > h2->degree())
        {
          int i = 0;
          while (i < n and a.getPP()[i] == b.getPP()[i]) ++i;
          if (i == n) result = false;
          else
            result = (a.getPP()[i] < b.getPP()[i]);
        }
        else
          result = true;
      }
      else
      {
        if (i > h2->degree()) result = false;
        else result = (*h1)[i] < (*h2)[i];
      }
    }
  }
  return result;
};

bool LessByDegreeThenGradHilbert(PPWithIdeal &a, PPWithIdeal &b)
{
  bool result;
  int n = a.getPP().num_vars();
  // first check the weighted degree
  if (a.getPP().weighted_degree(a.getOrdering().weights())
        < b.getPP().weighted_degree(b.getOrdering().weights()))
    result = true;
  else if (a.getPP().weighted_degree(a.getOrdering().weights())
        > b.getPP().weighted_degree(b.getOrdering().weights()))
    result = false;
  else {
    // now check the coefficients of the Hilbert polynomial
    Dense_Univariate_Rational_Polynomial HPdiff(*(a.getHilbertPolynomial()));
    HPdiff -= *(b.getHilbertPolynomial());
    if (not HPdiff.is_zero())
      result = (HPdiff.numerator(HPdiff.degree()) < 0);
    else // use Hilbert series
    {
      Dense_Univariate_Integer_Polynomial * h1 = a.getHilbertNumerator(true);
      Dense_Univariate_Integer_Polynomial * h2 = b.getHilbertNumerator(true);
      DEG_TYPE i = 0;
      for ( /* already initialized */ ;
           i <= h1->degree() and i <= h2->degree() and (*h1)[i] == (*h2)[i];
           i++)
      { /* taken care of in loop */ }
      if (i > h1->degree())
      {
        if (i > h2->degree())
        {
          int i = 0;
          while (i < n and a.getPP()[i] == b.getPP()[i]) ++i;
          if (i == n) result = false;
          else result = (a.getPP()[i] < b.getPP()[i]);
        }
        else
          result = true;
      }
      else
      {
        if (i > h2->degree()) result = false;
        else result = (*h1)[i] < (*h2)[i];
      }
    }
  }
  return result;
};

bool LessByBetti (PPWithIdeal & a, PPWithIdeal & b) {
  bool result;
  const map<DEG_TYPE, unsigned long> & Ba = a.getIncBetti();
  const map<DEG_TYPE, unsigned long> & Bb = b.getIncBetti();
  auto Bai = Ba.begin();
  auto Bbi = Bb.begin();
  for (/* */ ;
       Bai != Ba.end() and Bbi != Bb.end()
          and (Bai->first == Bbi->first and Bai->second == Bbi->second);
       ++Bai, ++Bbi
  ) { }
  if (Bai == Ba.end()) {
    if (Bbi == Bb.end())
      result = LessByHilbertThenDegree(a, b);
    else
      result = false;
  } else {
    if (Bbi == Bb.end()) {
      result = true;
    } else {
      if (Bai->first < Bbi->first)
        result = false;
      else if (Bai->first > Bbi->first)
        result = true;
      else if (Bai->second < Bbi->second)
        result = false;
      else
        result = a.getPP() < b.getPP();
    }
  }
  return result;
}

bool LessByBigBetti (PPWithIdeal & a, PPWithIdeal & b) {
  bool result;
  const map<DEG_TYPE, unsigned long> & Ba = a.getIncBetti();
  const map<DEG_TYPE, unsigned long> & Bb = b.getIncBetti();
  auto Bai = Ba.end();
  auto Bbi = Bb.end();
  for (/* */ ;
       Bai != Ba.begin() and Bbi != Bb.begin()
          and (Bai->first == Bbi->first and Bai->second == Bbi->second);
       --Bai, --Bbi
  ) { }
  if (Bai == Ba.begin()) {
    if (Bbi == Bb.begin())
      result = LessByHilbertThenDegree(a, b);
    else
      result = true;
  } else {
    if (Bbi == Bb.begin()) {
      result = false;
    } else {
      if (Bai->first < Bbi->first)
        result = false;
      else if (Bai->first > Bbi->first)
        result = true;
      else if (Bai->second < Bbi->second)
        result = false;
      else
        result = a.getPP() < b.getPP();
    }
  }
  return result;
}

bool LessByGradBetti (PPWithIdeal & a, PPWithIdeal & b) {
  bool result;
  const map<DEG_TYPE, unsigned long> & Ba = a.getIncBetti(true);
  const map<DEG_TYPE, unsigned long> & Bb = b.getIncBetti(true);
  auto Bai = Ba.begin();
  auto Bbi = Bb.begin();
  for (/* */ ;
       Bai != Ba.end() and Bbi != Bb.end()
          and (Bai->first == Bbi->first and Bai->second == Bbi->second);
       ++Bai, ++Bbi
  ) { }
  if (Bai == Ba.end()) {
    if (Bbi == Bb.end())
      result = LessByGradHilbertThenDegree(a, b);
    else
      result = false;
  } else {
    if (Bbi == Bb.end()) {
      result = true;
    } else {
      if (Bai->first < Bbi->first)
        result = false;
      else if (Bai->first > Bbi->first)
        result = true;
      else if (Bai->second < Bbi->second)
        result = false;
      else
        result = a.getPP() < b.getPP();
    }
  }
  return result;
}

void compatiblePP(
  Monomial currentLPP,            // the current LPP
  const set<Monomial> &allPPs,    // the monomials to consider;
                                      // some will be removed
  const set<::ray> &bndrys,     // known boundary vectors
  set<Monomial> &result,          // returned as PPs for Hilbert function
                                      // ("easy" (& efficient?) to extract exps
  set<Monomial> &boundary_mons,   // boundary monomials
  LP_Solver *skel              // used for alternate refinement
)
{
  // get the exponent vector of the current LPP, insert it
  NVAR_TYPE n = currentLPP.num_vars();
  ::ray aray(n, currentLPP.log());
  set<Monomial> initial_candidates;
  initial_candidates.insert(currentLPP);
  // compare other monomials with LPP
  /*for (
        auto b_ptr = allPPs.begin();
        b_ptr != allPPs.end();
        ++b_ptr
      )
  {
    // take the dot product of each monomial's exponents with pp,
    // add the ones that pass muster to initial_candidates
    //cout << '\t'; p_Write(*b_ptr, Rx);
    p_GetExpV(*b_ptr, b, Rx);
    for (unsigned long i = 1; i <= n; ++i) { along[i-1] = b[i]; }
    ::ray bray(n, along);
    bool searching = true;
    for (
         auto w_ptr = bndrys.begin();
         searching and w_ptr != bndrys.end();
         ++w_ptr
        )
    {
      // check b against a with all rays
      //cout << "checking " << bray << " against " << aray << " via " << *w_ptr << ": " << (*w_ptr) * bray << ' ' << (*w_ptr) * aray << endl;
      if ((*w_ptr) * bray > (*w_ptr) * aray)
      {
        // only need one success
        initial_candidates.insert(*b_ptr);
        cout << "succeeded with " << bray << endl;
        searching = false;
      }
    }
  }*/
  for (const Monomial & b : allPPs)
    if (not b.is_like(currentLPP) and skel->makes_consistent_constraint(b, currentLPP))
      initial_candidates.insert(b);
  for (const Monomial & t : initial_candidates)
    result.insert(t);
  // (new) alternate refinement: compare remaining monomials against each other,
  // using skeleton to ensure consistency
  // this approach should be more efficient than the one below, yet equivalent to it
  for (const Monomial & t : initial_candidates)
  {
    // cout << "testing for consistency: "; pWrite(*t);
    bool good_constraints = true;
    /*for (set<Monomial>::const_iterator u = initial_candidates.begin();
         good_constraints and u != initial_candidates.end();
         ++u)
      if (not t.is_like(*u))
        if (not skel->makes_consistent_constraint(t,*u))
          good_constraints = false;*/
    for (const Monomial & u : initial_candidates)
      if (not t.is_like(u))
        if (not skel->makes_consistent_constraint(t, u)) {
          good_constraints = false;
          break;
        }
    if (good_constraints)
    {
      boundary_mons.insert(t);
      // cout << "\tconsistent!\n";
    }
  }
  // second refinement: compare remaining monomials with each other
  /*for (auto p_ptr = initial_candidates.begin(); p_ptr != initial_candidates.end(); ++p_ptr)
  {
    p_GetExpV(*p_ptr, a, Rx);
    for (unsigned long i = 1; i <= n; ++i) { along[i-1] = a[i]; }
    ::ray pray(n, along);
    bool cleared_the_hurdle = false;
    for (auto w_ptr = bndrys.begin();
         not cleared_the_hurdle and w_ptr != bndrys.end(); ++w_ptr)
    {
      bool maybe_this_vector = true;
      ::ray w = *w_ptr;
      for (auto q_ptr = initial_candidates.begin();
           maybe_this_vector and not cleared_the_hurdle and q_ptr != initial_candidates.end(); ++q_ptr)
        if (*p_ptr != *q_ptr)
        {
          p_GetExpV(*q_ptr, b, Rx);
          for (unsigned long i = 1; i <= n; ++i) { along[i-1] = b[i]; }
          ::ray qray(n, along);
          // guard against invalid exclusions
          unsigned long long wt = (*w_ptr) * qray;
          if ((((*w_ptr) * pray) <= wt) and wt != 0) { maybe_this_vector = false; }
        }
      cleared_the_hurdle = maybe_this_vector;
      if (not cleared_the_hurdle) boundary_mons.insert(*p_ptr);
    }
    if (cleared_the_hurdle) result.insert(*p_ptr);
  } */
  //cout << "boundary monomials:\n";
  //for (auto piter = boundary_mons.begin(); piter != boundary_mons.end(); ++piter) pWrite(*piter);
  //cout << "compatible monomials:\n";
  //for (auto piter = result.begin(); piter != result.end(); ++piter) pWrite(*piter);
  //result = initial_candidates;
}

bool verifyAndModifyIfNecessary(
  LP_Solver *skel,
  const list<Abstract_Polynomial *> &currentPolys
)
{
  bool consistent = true; // innocent until proven guilty
  ::ray w = ray_sum(skel->get_rays()); // our tentative ordering
  // cout << "Have ray " << w << endl;
  NVAR_TYPE n = w.get_dimension();
  RAYENT_TYPE *entries = new RAYENT_TYPE [n]; // used for entries for new rays
  CONSTR_TYPE *coefficients = new CONSTR_TYPE [n]; // used for coefficients for new constraints
  // loop through all polynomials; verify leading power product is unchanged
  for (
       auto piter = currentPolys.begin(); piter != currentPolys.end(); ++piter
       // next line commented out for GLPK_Solver (approx skel causes issues...)
       //Abstract_Polynomial * pp : currentPolys
      )
  {
    // create a ray for the current LPP's exponents
    Abstract_Polynomial * pp = *piter;
    const Monomial & t = pp->leading_monomial();
    for (NVAR_TYPE i = 0; i < n; ++i) entries[i] = t[i];
    ::ray a(n, entries);
    // loop through the polynomial's remaining monomials
    //for (poly titer = (*piter)->next; consistent and titer != nullptr; titer = titer->next)
    Polynomial_Iterator * ti;
    for (
          ti = pp->new_iterator();
          consistent and not ti->fellOff();
          ti->moveRight()
    ) {
      // don't compare with LPP; that would be Bad (TM)
      if (pp->leading_monomial() != ti->currMonomial())
      {
        // create a ray for the PP's exponents
        for (NVAR_TYPE i = 0; i < n; ++i) entries[i] = ti->currMonomial()[i];
        ::ray b(n, entries);
        // compare weights between a and b; if this fails,
        // recompute the skeleton with a new constraint
        if (a*w <= b*w)
        {
          if (coefficients == nullptr) // ensure we have space
            coefficients = new CONSTR_TYPE[n];
          for (NVAR_TYPE i = 0; i < n; ++i)
            coefficients[i] = a[i] - b[i];
          constraint new_constraint(n, coefficients);
          LP_Solver * newskel;
          if (dynamic_cast<skeleton *>(skel) != nullptr)
            newskel = new skeleton(*static_cast<skeleton *>(skel));
          else if (dynamic_cast<GLPK_Solver *>(skel) != nullptr) {
            newskel = new GLPK_Solver(*static_cast<GLPK_Solver *>(skel));
            //piter = currentPolys.begin(); break;
          }
          else if (dynamic_cast<PPL_Solver *>(skel) != nullptr)
            newskel = new PPL_Solver(*static_cast<PPL_Solver *>(skel));
          consistent = newskel->solve(new_constraint);
          w = ray_sum(skel->get_rays());
          //cout << "Have ray " << w << endl;
          // if we're consistent, we need to recompute the ordering
          if (consistent and a*w > b*w)
          {
            w = ray_sum(skel->get_rays());
            //cout << "Have ray " << w << endl;
            *skel = *newskel;
          } // if consistent
          else consistent = false;
          delete newskel;
        } // if LPP changed
      } // if PP != LPP
    } // loop through PPs
    delete ti;
  } // loop through polys
  // cleanup
  delete [] entries;
  delete [] coefficients;
  // finally done
  return consistent;
}

void ConstraintsForNewPP(
  const PPWithIdeal &I,
  const set<Monomial> &monomialsForComparison,
  vector<constraint> &result
)
{
  // setup
  NVAR_TYPE n = I.getIdeal().number_of_variables();
  //int * a = new int[n];   // space for exponent vectors
  //int * b = new int[n];
  const EXP_TYPE * a, * b;
  CONSTR_TYPE *c = new CONSTR_TYPE[n];  // space for coefficients of constraint
  a = I.getPP().log();  // exponent vector of candidate
  // loop through exponent vectors of other 
  for (const Monomial & t : monomialsForComparison)
  {
    // insert only different PPs (since I->t should also be in that set)
    if (t != I.getPP())
    {
      b = t.log();
      for (NVAR_TYPE i = 0; i < n; ++i) c[i] = a[i] - b[i];
      result.push_back(constraint(n,c));
    }
  }
  delete [] c;
  //delete [] b;
  //delete [] a;
}

void SelectMonomial(
    Abstract_Polynomial * r,                    // changes
    list<Monomial> & CurrentLPPs,       // changes
    Dense_Univariate_Integer_Polynomial ** current_hilbert_numerator,
    const list<Abstract_Polynomial *> & CurrentPolys,
    const list<Critical_Pair_Dynamic *> & crit_pairs,
    LP_Solver * currSkel,                        // possibly changes
    bool & ordering_changed,
    DynamicHeuristic method
)
{
  //cout << "entering selmon\n";
  //cout << "skeleton before: " << currSkel << endl;
  skeleton * src_skel = dynamic_cast<skeleton *>(currSkel);
  GLPK_Solver * src_GLPK = dynamic_cast<GLPK_Solver *>(currSkel);
  PPL_Solver * src_PPL = dynamic_cast<PPL_Solver *>(currSkel);
  ::ray w = ray_sum(currSkel->get_rays());
  //cout << "Have ray " << w << endl;
  vector<WT_TYPE> ord(w.get_dimension());
  for (NVAR_TYPE i = 0; i < w.get_dimension(); ++i) { ord.push_back(w[i]); }
  const Monomial & currentLPP = r->leading_monomial();
  //cout << "comparing against: "; p_Write(currentLPP, Rx);
  set<Monomial> allPPs, boundaryPPs, compatiblePPs;
  // transform monomials into exponent vectors
  Polynomial_Iterator * ti;
  for (
       ti = r->new_iterator();
       not ti->fellOff();
       ti->moveRight()
  ) {
    allPPs.insert(ti->currMonomial());
  }
  delete ti;
  // loop through all exponent vectors
  cout << allPPs.size() << " possible monomials\n";
  compatiblePP(currentLPP, allPPs, currSkel->get_rays(), compatiblePPs, boundaryPPs, currSkel);
  cout << compatiblePPs.size() << " compatible monomials\n";
  //for (auto piter = compatiblePPs.begin(); piter != compatiblePPs.end(); ++piter)
  //  p_Write(*piter, Rx);
  // list possible future ideals, sort by Hilbert Function
  list<PPWithIdeal> possibleIdealsBasic;
  for (const Monomial & t : compatiblePPs)
  {
    PPWithIdeal newIdeal(t, CurrentLPPs, w, crit_pairs, *current_hilbert_numerator);
    possibleIdealsBasic.push_back(newIdeal);
  }
  //cout << "heuristic: " << method << endl;
  switch(method)
  {
    case DynamicHeuristic::ORD_HILBERT_THEN_LEX:
      possibleIdealsBasic.sort(LessByHilbert);
      break;
    case DynamicHeuristic::ORD_HILBERT_THEN_DEG:
      possibleIdealsBasic.sort(LessByHilbertThenDegree);
      break;
    case DynamicHeuristic::DEG_THEN_ORD_HILBERT:
      possibleIdealsBasic.sort(LessByDegreeThenHilbert);
      break;
    case DynamicHeuristic::GRAD_HILB_THEN_DEG:
      possibleIdealsBasic.sort(LessByGradHilbertThenDegree);
      break;
    case DynamicHeuristic::DEG_THEN_GRAD_HILB:
      possibleIdealsBasic.sort(LessByDegreeThenGradHilbert);
      break;
    case DynamicHeuristic::SMOOTHEST_DEGREES:
      possibleIdealsBasic.sort(LessBySmoothestDegrees);
      break;
    case DynamicHeuristic::LARGEST_MAX_COMPONENT:
      possibleIdealsBasic.sort(LessByLargestMaxComponent);
      break;
    case DynamicHeuristic::MIN_CRIT_PAIRS:
      possibleIdealsBasic.sort(LessByNumCritPairs);
      break;
    case DynamicHeuristic::BETTI_HILBERT_DEG:
      possibleIdealsBasic.sort(LessByBetti);
      break;
    case DynamicHeuristic::GRAD_BETTI_HILBERT_DEG:
      possibleIdealsBasic.sort(LessByGradBetti);
      break;
    default: possibleIdealsBasic.sort(LessByHilbert);
  }
  PPWithIdeal * winner = & possibleIdealsBasic.front();
  bool searching = true;
  if (possibleIdealsBasic.size() != 1)
  {
    // test each combination of LPPs for consistency
    // one of them must work (current LPP, if nothing else -- see previous case) 
    set<Monomial> PPunion;
    for (const Monomial & t : compatiblePPs)
      PPunion.insert(t);
    for (const Monomial & t : boundaryPPs)
      PPunion.insert(t);
    for (PPWithIdeal & I : possibleIdealsBasic) {
      //cout << currSkel << endl;
      LP_Solver * newSkeleton;
      if (src_skel != nullptr)
        newSkeleton = new skeleton(*src_skel);
      else if (src_GLPK != nullptr)
        newSkeleton = new GLPK_Solver(*src_GLPK);
      else if (src_PPL != nullptr)
        newSkeleton = new PPL_Solver(*src_PPL);
      vector<constraint> newvecs;
      //cout << "testing " << I.getPP() << endl;
      ConstraintsForNewPP(I, PPunion, newvecs);
      if (newSkeleton->solve(newvecs))
      {
        //cout << "consistent\n";
        if (verifyAndModifyIfNecessary(newSkeleton, CurrentPolys))
        {
          searching = false;
          if (src_skel != nullptr)
            src_skel -> copy(newSkeleton);
          else if (src_GLPK != nullptr)
            src_GLPK -> copy(newSkeleton);
          else if (src_PPL != nullptr)
            src_PPL -> copy(newSkeleton);
          winner = & I;
          delete newSkeleton;
          break;
        }
      }
      else
      {
        //cout << "inconsistent\n";
        // cout << newSkeleton;
        // this monomial is not, in fact, compatible
        compatiblePPs.erase(I.getPP());
      }
      delete newSkeleton;
    }
  }
  else if (possibleIdealsBasic.size() == 1 and compatiblePPs.size() != 1)
  {
    vector<constraint> newvecs;
    ConstraintsForNewPP(*(possibleIdealsBasic.begin()), compatiblePPs, newvecs);
    currSkel->solve(newvecs);
    verifyAndModifyIfNecessary(currSkel, CurrentPolys);
  }
    
  // set marked lpp, new Hilbert numerator
  CurrentLPPs.push_back(winner->getPP());
  if (*current_hilbert_numerator != nullptr) delete *current_hilbert_numerator;
  *current_hilbert_numerator
      = new Dense_Univariate_Integer_Polynomial(*(winner->getHilbertNumerator()));
  // TODO: delete elements of allPPs (not clear how: elements are in a set)
  ::ray new_weight = ray_sum(currSkel->get_rays());
  //cout << "Have ray " << w << endl;
  new_weight.simplify_ray();
  for (
        int i = 0;
        not ordering_changed and i < (int )(new_weight.get_dimension());
        ++i
  ) {
    // w is the old ordering
    ordering_changed = ordering_changed or (new_weight[i] != w[i]);
  }
  cout << "ordering changed? " << ordering_changed << endl;
  //cout << "finished with "; pWrite(t);
  /*for (unsigned long i = 0; i < CurrentLPPs.size(); ++i) pWrite(CurrentLPPs[i]);
  cout << endl;
  cout << "skeleton after:\n";
  cout << currSkel; */
  //cout << "returning from selmon\n";
}

#endif