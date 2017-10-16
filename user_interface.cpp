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

#include <list>
#include <string>
#include <iostream>

using std::list;
using std::string;
using std::cin; using std::cout; using std::endl;

#include "system_constants.hpp"

#include "monomial.hpp"
#include "dynamic_engine.hpp"
#include "polynomial_array.hpp"

#include "algorithm_buchberger_dynamic.hpp"

/**
  @ingroup utils
  @author John Perry
  @date 2017
  @brief reads ideal definition from @c cin, then reads options for computation,
    writes Gr&ouml;bner basis to @c cout
  @return 0 if success, nonzero otherwise
  @details The definition of the ideal should take place as follows, <b>in this
    order</b>. Bad Things Will Happen<sup>TM</sup> if you do not heed.
    The good news is that the program prompts the user.
      -# a prime number, which defines the base field
      -# the number of indeterminates
      -# whether to specify the indeterminates&rsquo; names
      -# if so, the indeterminates&rsquo; names
      -# the number of polynomials
      -# the definition of the polynomials (each terminated by a new line)
      -# static or dynamic computation
      -# if dynamic,
        -# which solver (LP_Solver)
        -# which heuristic (Dynamic_Heuristic)
        -# whether to perform a global analysis of the input
    The code has some other options, but they are mere curiosities for the time
    being.
*/
void user_interface() {
  UCOEF_TYPE p;
  cout << "prime number for base field: ";
  cin >> p;
  NVAR_TYPE n;
  cout << "number of indeterminates: ";
  cin >> n;
  char specify_names;
  cout << "do you want to specify the indeterminates' names? (y/n) ";
  cin >> specify_names;
  string * names = new string[n];
  Prime_Field F(p);
  Polynomial_Ring * P;
  if (not specify_names) {
    P = new Polynomial_Ring(n, F);
    cout << "indeterminate names are: ";
    for (NVAR_TYPE i = 0; i < n; ++i) {
      cout << P->name(i);
      names[i] = P->name(i);
    }
  } else {
    for (NVAR_TYPE i = 0; i < n; ++i)
      cin >> names[i];
    P = new Polynomial_Ring(n, F, names);
  }
  unsigned m;
  cout << "number of polynomials: ";
  cin >> m;
  cin.ignore();
  list<Abstract_Polynomial *> I;
  for (unsigned np = 0; np < m; ++np) {
    cout << "please enter polynomial #" << np << ": (use spaces to separate factors) ";
    string inpoly;
    getline(cin, inpoly);
    list<Monomial> M; // monomials
    list<Prime_Field_Element> A; // coefficients
    unsigned nt = 0; // number of terms
    unsigned i = 0;
    while (i < inpoly.size()) {
      bool reading_term = true;
      COEF_TYPE a = 1;
      EXP_TYPE exp[n];
      for (NVAR_TYPE i = 0; i < n; ++i) exp[i] = 0;
      while (inpoly[i] == ' ') ++i;
      while (reading_term and i < inpoly.size()) {
        unsigned j = i;
        if (
            inpoly[i] == '+' or inpoly[i] == '-' or
            (inpoly[i] >= '0' and inpoly[i] <= '9')
        ) {
          bool positive = true;
          if (inpoly[i] == '-') positive = false;
          if (inpoly[i] == '+' or inpoly[i] == '-') ++i;
          while (inpoly[i] == ' ') ++i;
          if (inpoly[i] >= '0' and inpoly[i] <= '9') {
            j = i;
            while (inpoly[j] >= '0' and inpoly[j] <= '9') ++j;
            a *= stol(inpoly.substr(i, j - i));
            if (not positive) a *= -1;
            i = j;
          }
        }
        else if (inpoly[i] == '*')
          ++i;
        else if (inpoly[i] == '\0')
          ++i;
        else { // should be characters
          unsigned j = i;
          while ((inpoly[j] >= 'a' and inpoly[j] <= 'z') or
                 (inpoly[j] >= 'A' and inpoly[j] <= 'Z') or
                 (inpoly[j] >= '0' and inpoly[j] <= '9'))
            ++j;
          string var_name = inpoly.substr(i, j - i);
          bool found = false;
          unsigned k = 0;
          for (/* */; not found and k < n; ++k)
            if (var_name.compare(names[k]) == 0) {
              found = true;
              --k;
            }
          i = j;
          while (inpoly[i] == ' ') ++i;
          if (inpoly[i] != '^')
            exp[k] += 1;
          else {
            ++i;
            while (inpoly[i] == ' ') ++i;
            unsigned j = i;
            while (inpoly[j] >= '0' and inpoly[j] <= '9') ++j;
            exp[k] += stol(inpoly.substr(i, j - i));
            i = j;
          }
        }
        while (inpoly[i] == ' ') ++i;
        if (inpoly[i] == '+' or inpoly[i] == '-') reading_term = false;
      }
      ++nt;
      Monomial t(n, exp);
      M.push_back(t);
      A.emplace_back(a, &F);
    }
    Abstract_Polynomial * p = new Constant_Polynomial(*P, M, A);
    p->sort_by_order();
    cout << "read " << *p << endl;
    I.push_back(p);
  }
  string computation;
  while (
      computation.compare("s") and computation.compare("d")
      and computation.compare("static") and computation.compare("dynamic")
  ) {
    cout << "static (s) or dynamic (d) computation? ";
    getline(cin, computation);
  }
  list<Constant_Polynomial *> B;
  bool dynamic = not (computation.compare("d") and computation.compare("dynamic"));
  if (not dynamic) {
    B = buchberger(I);
  } else {
    DynamicSolver solver;
    string solver_choice;
    while (
        solver_choice.compare("skel") and solver_choice.compare("glpk")
        and solver_choice.compare("ppl") and solver_choice.compare("skeleton")
    ) {
      cout << "which solver? ([skel]eton, glpk, ppl) ";
      getline(cin, solver_choice);
    }
    if (not solver_choice.compare("skel") or not solver_choice.compare("skeleton"))
      solver = SKELETON_SOLVER;
    else if (not solver_choice.compare("glpk"))
      solver = GLPK_SOLVER;
    else if (not solver_choice.compare("ppl"))
      solver = PPL_SOLVER;
    string heur_choice;
    Dynamic_Heuristic heuristic;
    while (
      heur_choice.compare("h") and heur_choice.compare("gh")
      and heur_choice.compare("b") and heur_choice.compare("bb")
      and heur_choice.compare("gb") and heur_choice.compare("c")
      and heur_choice.compare("d") and heur_choice.compare("r")
    ) {
      cout << "available heuristics are\n";
      cout << "\t h = hilbert with standard grading\n";
      cout << "\tgh = hilbert with order-based grading\n";
      cout << "\t b = betti with standard grading\n";
      cout << "\tbb = \"big\" betti with standard grading\n";
      cout << "\tgb = betti with order-based grading\n";
      cout << "\t c = minimal number of critical pairs\n";
      cout << "\t d = minimal degree, ties broken by hilbert\n";
      cout << "\t r = random choice\n";
      cout << "which heuristic would you like? ";
      getline(cin, heur_choice);
    }
    if      (not heur_choice.compare( "h")) heuristic = Dynamic_Heuristic::ORD_HILBERT_THEN_DEG;
    else if (not heur_choice.compare("gh")) heuristic = Dynamic_Heuristic::GRAD_HILB_THEN_DEG;
    else if (not heur_choice.compare( "b")) heuristic = Dynamic_Heuristic::BETTI_HILBERT_DEG;
    else if (not heur_choice.compare("bb")) heuristic = Dynamic_Heuristic::BIG_BETTI_HILBERT_DEG;
    else if (not heur_choice.compare("gb")) heuristic = Dynamic_Heuristic::GRAD_BETTI_HILBERT_DEG;
    else if (not heur_choice.compare( "c")) heuristic = Dynamic_Heuristic::MIN_CRIT_PAIRS;
    else if (not heur_choice.compare( "d")) heuristic = Dynamic_Heuristic::DEG_THEN_ORD_HILBERT;
    else if (not heur_choice.compare( "r")) heuristic = Dynamic_Heuristic::EVIL_RANDOM;
    string whether_analysis;
    bool analyze_first = false;
    while (whether_analysis.compare("y") and whether_analysis.compare("n")) {
      cout << "perform global analysis at the outset? (y or n) ";
      cin >> whether_analysis;
    }
    if (not whether_analysis.compare("y")) analyze_first = true;
    B = buchberger_dynamic(
        I, SPolyCreationFlags::GEOBUCKETS, StrategyFlags::SUGAR_STRATEGY,
        nullptr, heuristic, solver, analyze_first
    );
  }
  //check_correctness(B, StrategyFlags::NORMAL_STRATEGY);
  cout << "have basis with " << B.size() << " elements:\n";
  for (auto b : B) {
    cout << b->leading_monomial() << ", ";
    delete b;
  }
  cout << endl;
  for (auto p : I) delete p;
  delete P;
  delete [] names;
}

int main() {
  user_interface();
}