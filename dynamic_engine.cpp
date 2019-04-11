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
* DynGB is distributed in the hope that it will be useful,                    *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with DynGB. If not, see <http://www.gnu.org/licenses/>.               *
\*****************************************************************************/

#include <iostream>
using std::cout; using std::endl;
#include <list>
using std::list;
#include <cstdlib>
#include <ctime>

#include "betti.hpp"
#include "monomial.hpp"
#include "polynomial.hpp"
#include "dynamic_engine.hpp"

namespace Dynamic_Engine {

bool less_by_random(const PP_With_Ideal & a, const PP_With_Ideal & b) {
  static WT_TYPE * weights = nullptr;
  const Monomial & t = a.get_pp();
  const Monomial & u = b.get_pp();
  NVAR_TYPE n = t.num_vars();
  if (weights == nullptr) {
    srand(time(nullptr));
    weights = new WT_TYPE[n];
    cout << "evil random: ";
    for (NVAR_TYPE i = 0; i < n; ++i) {
      weights[i] = rand() % 1000;
      cout << weights[i] << ' ';
    }
    cout << endl;
  }
  bool result = false;
  if (t.weighted_degree(weights) < u.weighted_degree(weights))
    result = true;
  return result;
}

int hilbert_cmp(
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

bool less_by_hilbert (PP_With_Ideal &a, PP_With_Ideal &b) {
  //cout << "Less by Hilbert then Lex\n";
  bool result;
  // first check the coefficients of the Hilbert polynomial
  int hilcheck = hilbert_cmp(
    *a.get_hilbert_numerator(), *a.get_hilbert_polynomial(),
    *b.get_hilbert_numerator(), *b.get_hilbert_polynomial()
  );
  if (hilcheck == 0) // the numerators are equal; break tie via current ordering
    result = (a.get_pp() < b.get_pp());
  else {
    if (hilcheck == -1) result = true;
    else result = false;
  }
  //cout << "\tfirst less than second? " << result << endl;
  return result;
};

bool less_by_smoothest_degrees (PP_With_Ideal &a, PP_With_Ideal &b)
{
  return a.get_difference_in_degree() < b.get_difference_in_degree();
};

bool less_by_largest_max_component (PP_With_Ideal &a, PP_With_Ideal &b) {
  return a.get_difference_in_degree() < b.get_difference_in_degree();
}

void PP_With_Ideal::compute_number_new_pairs() const {
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

bool less_by_num_crit_pairs (PP_With_Ideal &a, PP_With_Ideal &b)
{
  bool result;
  // first check if the number of critical pairs has been computed
  if (a.how_many_new_pairs() < 0) a.compute_number_new_pairs();
  if (b.how_many_new_pairs() < 0) b.compute_number_new_pairs();
  /*if (a.degOfNewPairs() < b.degOfNewPairs())
    result = true;
  else if (a.degOfNewPairs() > b.degOfNewPairs())
    result = false;*/
  // at this point, the degrees of the new pairs will be equal
  /*else*/ if (a.how_many_new_pairs() > b.how_many_new_pairs())
    result = true;
  else if (a.how_many_new_pairs() < b.how_many_new_pairs())
    result = false;
  else // the numerators are equal; break tie via monomial ordering
    result = a.get_pp() < b.get_pp();
  //cout << "\tfirst less than second? " << result << endl;
  return result;
};

bool less_by_hilbert_then_degree(PP_With_Ideal &a, PP_With_Ideal &b)
{
  //cout << "Less by Hilbert then Deg\n";
  bool result;
  // first check the coefficients of the Hilbert polynomial
  int hilcheck = hilbert_cmp(
      *a.get_hilbert_numerator(), *a.get_hilbert_polynomial(),
      *b.get_hilbert_numerator(), *b.get_hilbert_polynomial()
  );
  if (hilcheck == -1) result = true;
  else if (hilcheck == 1) result = false;
  else { // the numerators are equal; break tie via monomial degree
    if (a.get_pp().total_degree() < b.get_pp().total_degree())
      result = true;
    else if (a.get_pp().total_degree() > b.get_pp().total_degree())
      result = false;
    else {
      //cout << "breaking a Hilbert & degree tie by ordering: " << a.get_pp() << " < " << b.get_pp() << "? ";
      result = (a.get_pp() < b.get_pp());
      //cout << result << endl;
    }
  }
  return result;
};

bool less_by_grad_hilbert_then_degree(PP_With_Ideal &a, PP_With_Ideal &b) {
  bool result = false;
  // first check the coefficients of the Hilbert polynomial
  Dense_Univariate_Rational_Polynomial * hp1 = a.get_hilbert_polynomial();
  Dense_Univariate_Rational_Polynomial * hp2 = b.get_hilbert_polynomial();
  int hilcheck = hilbert_cmp(
    *a.get_hilbert_numerator(true), *hp1,
    *b.get_hilbert_numerator(true), *hp2
  );
  if (hilcheck == -1) result = true;
  else if (hilcheck == 1) result = false;
  else if (hilcheck == 0) { // break tie via monomial degree
    if (a.get_pp().weighted_degree(a.get_ordering().weights())
          < b.get_pp().weighted_degree(b.get_ordering().weights()))
      result = true;
    else if (a.get_pp().weighted_degree(a.get_ordering().weights())
                > b.get_pp().weighted_degree(b.get_ordering().weights()))
      result = false;
    else
      result = a.get_pp() < b.get_pp();
  }
  return result;
};

bool less_by_degree_then_hilbert(PP_With_Ideal &a, PP_With_Ideal &b)
{
  //cout << "Less by Deg then Hilbert\n";
  bool result;
  int n = a.get_pp().num_vars();
  // first check the weighted degree
  if (a.get_pp().total_degree() < b.get_pp().total_degree())
    result = true;
  else if (a.get_pp().total_degree() > b.get_pp().total_degree())
    result = false;
  else {
    // now check the coefficients of the Hilbert polynomial
    Dense_Univariate_Rational_Polynomial HPdiff(*(a.get_hilbert_polynomial()));
    HPdiff -= *(b.get_hilbert_polynomial());
    if (not HPdiff.is_zero())
      result = (HPdiff.numerator(HPdiff.degree()) < 0);
    else // use Hilbert series
    {
      Dense_Univariate_Integer_Polynomial * h1 = a.get_hilbert_numerator();
      Dense_Univariate_Integer_Polynomial * h2 = b.get_hilbert_numerator();
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
          while (i < n and a.get_pp()[i] == b.get_pp()[i]) ++i;
          if (i == n) result = false;
          else
            result = (a.get_pp()[i] < b.get_pp()[i]);
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

bool less_by_wdegree_then_hilbert(PP_With_Ideal &a, PP_With_Ideal &b)
{
  //cout << "Less by Deg then Hilbert\n";
  bool result;
  int n = a.get_pp().num_vars();
  // first check the weighted degree
  if (a.get_pp().weighted_degree(a.get_ordering().weights())
        < b.get_pp().weighted_degree(b.get_ordering().weights()))
    result = true;
  else if (a.get_pp().weighted_degree(a.get_ordering().weights())
        > b.get_pp().weighted_degree(b.get_ordering().weights()))
    result = false;
  else {
    // now check the coefficients of the Hilbert polynomial
    Dense_Univariate_Rational_Polynomial HPdiff(*a.get_hilbert_polynomial());
    HPdiff -= *b.get_hilbert_polynomial();
    if (not HPdiff.is_zero())
      result = (HPdiff.numerator(HPdiff.degree()) >= 0);
    else // use Hilbert series
    {
      Dense_Univariate_Integer_Polynomial * h1 = a.get_hilbert_numerator();
      Dense_Univariate_Integer_Polynomial * h2 = b.get_hilbert_numerator();
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
          while (i < n and a.get_pp()[i] == b.get_pp()[i]) ++i;
          if (i == n) result = false;
          else
            result = (a.get_pp()[i] < b.get_pp()[i]);
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

bool less_by_degree_then_grad_hilbert(PP_With_Ideal &a, PP_With_Ideal &b)
{
  bool result;
  int n = a.get_pp().num_vars();
  // first check the weighted degree
  if (a.get_pp().weighted_degree(a.get_ordering().weights())
        < b.get_pp().weighted_degree(b.get_ordering().weights()))
    result = true;
  else if (a.get_pp().weighted_degree(a.get_ordering().weights())
        > b.get_pp().weighted_degree(b.get_ordering().weights()))
    result = false;
  else {
    // now check the coefficients of the Hilbert polynomial
    Dense_Univariate_Rational_Polynomial HPdiff(*(a.get_hilbert_polynomial()));
    HPdiff -= *(b.get_hilbert_polynomial());
    if (not HPdiff.is_zero())
      result = (HPdiff.numerator(HPdiff.degree()) < 0);
    else // use Hilbert series
    {
      Dense_Univariate_Integer_Polynomial * h1 = a.get_hilbert_numerator(true);
      Dense_Univariate_Integer_Polynomial * h2 = b.get_hilbert_numerator(true);
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
          while (i < n and a.get_pp()[i] == b.get_pp()[i]) ++i;
          if (i == n) result = false;
          else result = (a.get_pp()[i] < b.get_pp()[i]);
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

bool less_by_betti (PP_With_Ideal & a, PP_With_Ideal & b) {
  bool result;
  const map<DEG_TYPE, unsigned long> & Ba = a.get_inc_betti();
  const map<DEG_TYPE, unsigned long> & Bb = b.get_inc_betti();
  auto Bai = Ba.begin();
  auto Bbi = Bb.begin();
  for (/* */ ;
       Bai != Ba.end() and Bbi != Bb.end()
          and (Bai->first == Bbi->first and Bai->second == Bbi->second);
       ++Bai, ++Bbi
  ) { }
  if (Bai == Ba.end()) {
    if (Bbi == Bb.end())
      result = less_by_hilbert_then_degree(a, b);
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
        result = a.get_pp() < b.get_pp();
    }
  }
  return result;
}

bool less_by_big_betti (PP_With_Ideal & a, PP_With_Ideal & b) {
  bool result;
  const map<DEG_TYPE, unsigned long> & Ba = a.get_inc_betti();
  const map<DEG_TYPE, unsigned long> & Bb = b.get_inc_betti();
  auto Bai = Ba.end();
  auto Bbi = Bb.end();
  for (/* */ ;
       Bai != Ba.begin() and Bbi != Bb.begin()
          and (Bai->first == Bbi->first and Bai->second == Bbi->second);
       --Bai, --Bbi
  ) { }
  if (Bai == Ba.begin()) {
    if (Bbi == Bb.begin())
      result = less_by_hilbert_then_degree(a, b);
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
        result = a.get_pp() < b.get_pp();
    }
  }
  return result;
}

bool less_by_grad_betti (PP_With_Ideal & a, PP_With_Ideal & b) {
  bool result;
  const map<DEG_TYPE, unsigned long> & Ba = a.get_inc_betti(true);
  const map<DEG_TYPE, unsigned long> & Bb = b.get_inc_betti(true);
  auto Bai = Ba.begin();
  auto Bbi = Bb.begin();
  for (/* */ ;
       Bai != Ba.end() and Bbi != Bb.end()
          and (Bai->first == Bbi->first and Bai->second == Bbi->second);
       ++Bai, ++Bbi
  ) { }
  if (Bai == Ba.end()) {
    if (Bbi == Bb.end())
      result = less_by_grad_hilbert_then_degree(a, b);
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
        result = a.get_pp() < b.get_pp();
    }
  }
  return result;
}

void compatible_pp(
  const Monomial & currentLPP,            // the current LPP
  const set<Monomial> & allPPs,   // the monomials to consider; some removed
  set<Monomial> &result,          // returned as PPs for Hilbert function
                                  // ("easy" (& efficient?) to extract exps
  list<Monomial> &boundary_mons,   // boundary monomials
  const LP_Solver *skel                 // used for alternate refinement
)
{
  // known boundary vectors
  const set<Ray> &bndrys = skel->get_rays();
  // get the exponent vector of the current LPP, insert it
  NVAR_TYPE n = currentLPP.num_vars();
  Ray aray(n, currentLPP.log());
  list<Monomial> initial_candidates;
  initial_candidates.push_back(currentLPP);
  // compare other monomials with LPP
  for (const Monomial & b : allPPs)
    if (not b.is_like(currentLPP) and skel->makes_consistent_constraint(b, currentLPP))
      initial_candidates.push_back(b);
  //for (const Monomial & t : initial_candidates)
  //  result.insert(t);
  // (new) alternate refinement: compare remaining monomials against each other,
  // using skeleton to ensure consistency
  // this approach should be more efficient than the one below, yet equivalent to it
  for (const Monomial & t : initial_candidates)
  {
    cout << "testing for consistency: " << t << endl;
    bool good_constraints = true;
    /*for (set<Monomial>::const_iterator u = initial_candidates.begin();
         good_constraints and u != initial_candidates.end();
         ++u)
      if (not t.is_like(*u))
        if (not skel->makes_consistent_constraint(t,*u))
          good_constraints = false;*/
    for (const Monomial & u : initial_candidates) {
      cout << "\tagainst " << u << endl;
      if (not t.is_like(u))
        if (not skel->makes_consistent_constraint(t, u)) {
          good_constraints = false;
          cout << "\tinconsistent\n";
          break;
        }
    }
    if (good_constraints)
    {
      //boundary_mons.push_back(t);
      result.insert(t);
      cout << "\tconsistent!\n";
    }
  }
}

bool verify_and_modify_if_necessary(
  LP_Solver *skel,
  const list<Abstract_Polynomial *> &currentPolys
)
{
  bool consistent = true; // innocent until proven guilty
  Ray w { ray_sum(skel->get_rays()) }; // our tentative ordering
  cout << "Have ray " << w << endl;
  NVAR_TYPE n = w.get_dimension();
  RAYENT_TYPE *entries = new RAYENT_TYPE [n]; // used for entries for new rays
  CONSTR_TYPE *coefficients = new CONSTR_TYPE [n]; // used for coefficients for new constraints
  // loop through all polynomials; verify leading power product is unchanged
  for (
       auto piter = currentPolys.begin(); consistent and piter != currentPolys.end(); ++piter
       // next line commented out for GLPK_Solver (approx skel causes issues...)
       //Abstract_Polynomial * pp : currentPolys
      )
  {
    // create a ray for the current LPP's exponents
    Abstract_Polynomial * pp = *piter;
    const Monomial & t = pp->leading_monomial();
    for (NVAR_TYPE i = 0; i < n; ++i) entries[i] = t[i];
    Ray a(n, entries);
    // loop through the polynomial's remaining monomials
    Polynomial_Iterator * ti;
    for (
          ti = pp->new_iterator();
          consistent and not ti->fellOff();
          ti->moveRight()
    ) {
      // don't compare with LPP; that would be Bad (TM)
      if (t != ti->currMonomial())
      {
        //cout << "\tagainst " << ti->currMonomial() << ": ";
        // create a ray for the PP's exponents
        for (NVAR_TYPE i = 0; i < n; ++i) entries[i] = ti->currMonomial()[i];
        Ray b(n, entries);
        // compare weights between a and b; if this fails,
        // recompute the skeleton with a new constraint
        //cout << a*w << ',' << b*w << endl;
        if (a*w <= b*w)
        {
          //cout << "ray fails to preserve " << t << " > " << ti->currMonomial() << " ; adding constraint and trying to recover\n";
          if (coefficients == nullptr) // ensure we have space
            coefficients = new CONSTR_TYPE[n];
          for (NVAR_TYPE i = 0; i < n; ++i)
            coefficients[i] = a[i] - b[i];
          Constraint new_constraint(n, coefficients);
          LP_Solver * newskel = nullptr;
          if (dynamic_cast<Skeleton *>(skel) != nullptr)
            newskel = new Skeleton(*static_cast<Skeleton *>(skel));
          else if (dynamic_cast<GLPK_Solver *>(skel) != nullptr)
            newskel = new GLPK_Solver(*static_cast<GLPK_Solver *>(skel));
          else if (dynamic_cast<PPL_Solver *>(skel) != nullptr)
            newskel = new PPL_Solver(*static_cast<PPL_Solver *>(skel));
          consistent = newskel->solve(new_constraint);
          auto w_tmp { ray_sum(newskel->get_rays()) };
          //cout << "Have ray " << w << endl;
          // if we're consistent, we need to recompute the ordering
          //cout << "\t\t" << a*w << ',' << b*w << endl;
          if (consistent and a*w_tmp > b*w_tmp)
          {
            skel->copy(newskel);
            piter = currentPolys.begin();
            w = w_tmp;
            cout << " recovered with << " << w_tmp << endl;
            break;
          } // if consistent
          else {
            consistent = false;
            cout << t << " fails again with ordering " << w_tmp << endl;
          }
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
  cout << "is consistent? " << consistent << endl;
  return consistent;
}

void constraints_for_new_pp(
  const PP_With_Ideal &I,
  const set<Monomial> &monomials_for_comparison,
  vector<Constraint> &result
)
{
  // setup
  NVAR_TYPE n = I.get_ideal().number_of_variables();
  //int * a = new int[n];   // space for exponent vectors
  //int * b = new int[n];
  const EXP_TYPE * a, * b;
  CONSTR_TYPE *c = new CONSTR_TYPE[n];  // space for coefficients of constraint
  a = I.get_pp().log();  // exponent vector of candidate
  // loop through exponent vectors of other 
  for (const Monomial & t : monomials_for_comparison)
  {
    // insert only different PPs (since I->t should also be in that set)
    if (t != I.get_pp())
    {
      //cout << "adding constraint " << I.get_pp() << " > " << t << endl;
      b = t.log();
      for (NVAR_TYPE i = 0; i < n; ++i) c[i] = a[i] - b[i];
      result.push_back(Constraint(n,c));
    }
  }
  delete [] c;
  //delete [] b;
  //delete [] a;
}

void select_monomial(
    Abstract_Polynomial * r,                    // changes
    list<Monomial> & CurrentLPPs,       // changes
    Dense_Univariate_Integer_Polynomial ** current_hilbert_numerator,
    const list<Abstract_Polynomial *> & CurrentPolys,
    const list<Critical_Pair_Dynamic *> & crit_pairs,
    LP_Solver * currSkel,                        // possibly changes
    bool & ordering_changed,
    Dynamic_Heuristic method
) {
  // transform monomials into exponent vectors
  set<Monomial> allPPs;
  Polynomial_Iterator * ti;
  for (
       ti = r->new_iterator();
       not ti->fellOff();
       ti->moveRight()
  ) {
    allPPs.insert(ti->currMonomial());
  }
  delete ti;
  select_monomial(
      allPPs, r->leading_monomial(), CurrentLPPs, current_hilbert_numerator,
      CurrentPolys, crit_pairs, currSkel, ordering_changed, method
  );
}

void select_monomial(
    const set<Monomial> & allPPs,
    const Monomial & currentLPP,
    list<Monomial> & CurrentLPPs,       // changes
    Dense_Univariate_Integer_Polynomial ** current_hilbert_numerator,
    const list<Abstract_Polynomial *> & CurrentPolys,
    const list<Critical_Pair_Dynamic *> & crit_pairs,
    LP_Solver * currSkel,                        // possibly changes
    bool & ordering_changed,
    Dynamic_Heuristic method
)
{
  //cout << "entering selmon\n";
  //cout << "skeleton before: " << currSkel << endl;
  Skeleton * src_skel    = dynamic_cast<Skeleton *>    (currSkel);
  GLPK_Solver * src_GLPK = dynamic_cast<GLPK_Solver *> (currSkel);
  PPL_Solver * src_PPL   = dynamic_cast<PPL_Solver *>  (currSkel);
  Ray w = ray_sum(currSkel->get_rays());
  //cout << "Have ray " << w << endl;
  vector<WT_TYPE> ord(w.get_dimension());
  for (NVAR_TYPE i = 0; i < w.get_dimension(); ++i) { ord.push_back(w[i]); }
  //cout << "comparing against: "; p_Write(currentLPP, Rx);
  list<Monomial> boundaryPPs;
  set<Monomial> compatible_pps;
  // loop through all exponent vectors
  cout << allPPs.size() << " possible monomials\n";
  compatible_pp(currentLPP, allPPs, compatible_pps, boundaryPPs, currSkel);
  cout << compatible_pps.size() << " compatible monomials\n";
  // list possible future ideals, sort by Hilbert Function
  list<PP_With_Ideal> possibleIdealsBasic;
  for (const Monomial & t : compatible_pps)
  {
    PP_With_Ideal newIdeal(t, CurrentLPPs, w, crit_pairs, *current_hilbert_numerator);
    possibleIdealsBasic.push_back(newIdeal);
    //cout << "pushed back " << t << endl;
  }
  //cout << "heuristic: " << method << endl;
  switch(method)
  {
    case Dynamic_Heuristic::ORD_HILBERT_THEN_LEX:
      possibleIdealsBasic.sort(less_by_hilbert);
      break;
    case Dynamic_Heuristic::ORD_HILBERT_THEN_DEG:
      possibleIdealsBasic.sort(less_by_hilbert_then_degree);
      break;
    case Dynamic_Heuristic::DEG_THEN_ORD_HILBERT:
      possibleIdealsBasic.sort(less_by_degree_then_hilbert);
      break;
    case Dynamic_Heuristic::GRAD_HILB_THEN_DEG:
      possibleIdealsBasic.sort(less_by_grad_hilbert_then_degree);
      break;
    case Dynamic_Heuristic::DEG_THEN_GRAD_HILB:
      possibleIdealsBasic.sort(less_by_degree_then_grad_hilbert);
      break;
    case Dynamic_Heuristic::SMOOTHEST_DEGREES:
      possibleIdealsBasic.sort(less_by_smoothest_degrees);
      break;
    case Dynamic_Heuristic::LARGEST_MAX_COMPONENT:
      possibleIdealsBasic.sort(less_by_largest_max_component);
      break;
    case Dynamic_Heuristic::MIN_CRIT_PAIRS:
      possibleIdealsBasic.sort(less_by_num_crit_pairs);
      break;
    case Dynamic_Heuristic::BETTI_HILBERT_DEG:
      possibleIdealsBasic.sort(less_by_betti);
      break;
    case Dynamic_Heuristic::GRAD_BETTI_HILBERT_DEG:
      possibleIdealsBasic.sort(less_by_grad_betti);
      break;
    case Dynamic_Heuristic::EVIL_RANDOM:
      possibleIdealsBasic.sort(less_by_random);
      break;
    default: possibleIdealsBasic.sort(less_by_hilbert);
  }
  PP_With_Ideal * winner = & possibleIdealsBasic.front();
  if (possibleIdealsBasic.size() != 1)
  {
    // test each combination of LPPs for consistency
    // one of them must work (current LPP, if nothing else -- see previous case) 
    set<Monomial> PPunion;
    for (const Monomial & t : compatible_pps)
      PPunion.insert(t);
    for (const Monomial & t : boundaryPPs)
      PPunion.insert(t);
    for (PP_With_Ideal & I : possibleIdealsBasic) {
      //cout << currSkel << endl;
      LP_Solver * newSkeleton;
      if (src_skel != nullptr)
        newSkeleton = new Skeleton(*src_skel);
      else if (src_GLPK != nullptr)
        newSkeleton = new GLPK_Solver(*src_GLPK);
      else if (src_PPL != nullptr)
        newSkeleton = new PPL_Solver(*src_PPL);
      vector<Constraint> newvecs;
      //cout << "testing " << I.get_pp() << endl;
      /*cout << '\t' << *I.get_hilbert_polynomial() << endl;
      cout << '\t' << *I.get_hilbert_numerator() << endl;*/
      constraints_for_new_pp(I, PPunion, newvecs);
      if (newSkeleton->solve(newvecs))
      {
        //cout << I.get_pp() << " is consistent\n";
        if (verify_and_modify_if_necessary(newSkeleton, CurrentPolys))
        {
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
        //cout << I.get_pp() << "is inconsistent\n";
        // cout << newSkeleton;
        // this monomial is not, in fact, compatible
        compatible_pps.erase(I.get_pp());
      }
      delete newSkeleton;
    }
  }
  else if (possibleIdealsBasic.size() == 1 and compatible_pps.size() != 1)
  {
    vector<Constraint> newvecs;
    constraints_for_new_pp(*(possibleIdealsBasic.begin()), compatible_pps, newvecs);
    currSkel->solve(newvecs);
    verify_and_modify_if_necessary(currSkel, CurrentPolys);
  }
    
  // set marked lpp, new Hilbert numerator
  CurrentLPPs.push_back(winner->get_pp());
  if (*current_hilbert_numerator != nullptr) delete *current_hilbert_numerator;
  *current_hilbert_numerator
      = new Dense_Univariate_Integer_Polynomial(*(winner->get_hilbert_numerator()));
  // TODO: delete elements of allPPs (not clear how: elements are in a set)
  Ray new_weight = ray_sum(currSkel->get_rays());
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

}

#endif